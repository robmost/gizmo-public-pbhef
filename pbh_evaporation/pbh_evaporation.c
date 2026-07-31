#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file pbh_evaporation.c
 *  \brief DM density calculation around particles
 *
 *    This file contains a loop modeled on the standard gas density computation which
 *    determines the DARK MATTER density around a given set of particles and adjusts the smoothing length for this
*     calculation.
 *    The dark matter density is used to set the energy injection due to primordial black hole (PBH) evaporation.
 *
 * This file was written by Robert Mostoghiu Paun, for GIZMO, based on Florian List's dark matter annihilation feedback routine.
 *
 * Method 1 - receiver-based approach (activate using PBHEF)
  - PBH (traced by N-Body DM particles) density is calculated at each gas particle using smoothing length HsmlPBH
  - from PBH density, PBH evaporation rate at gas particle is calculated
  - energy injection at each gas particle
 * Method 2 - donor-based approach (PBHEF=2)
  - the loop runs the other way round: HsmlPBH is a gas search radius and DensityPBH the gas density at the DM particle
  - the evaporation power follows from the DM particle mass alone, and is shared over its gas neighbors
 *
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially (condensed, different criteria for kernel lengths, optimizations,
 * rewritten parallelism, new physics included, new variable/memory conventions added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#ifdef PBHEF

struct kernel_density /*! defines a number of useful variables we will use below */
{
  double dp[3],dv[3],r, wk, dwk, hinv, hinv3, hinv4, mj_wk, mj_dwk_r;
};


/*! elements active in the density loop below: gas for the receiver-based method, DM for the donor-based one */
static int pbhef_density_isactive(int n)
{

    if(P[n].TimeBin < 0) {return 0;}
    if(P[n].Mass <= 0) {return 0;}
#if (PBHEF == 1)
    if(P[n].Type != 0){return 0;}  /* only gas particles */
#elif (PBHEF == 2)
	if(P[n].Type != 1){return 0;}  /* only DM particles */
#endif
    return 1;
}


#if (PBHEF == 1)
#define PBH_NGB_BITMASK 2  /* gas searches for DM neighbors */
#else
#define PBH_NGB_BITMASK 1  /* DM searches for gas neighbors */
#endif


/*! largest kernel length allowed for the DM neighbor search. All.MaxHsml defaults to an effectively
    infinite value, which lets a particle in a DM-poor region walk an enormous part of the tree, so we
    also keep the search within a small multiple of the particle's own kernel, as disp_density() does */
static double pbhef_return_maxhsml(int i)
{
    double maxsoft = All.MaxHsml;
#if (PBHEF == 1)
    if(PPP[i].Hsml > 0) {maxsoft = DMIN(maxsoft, 10.*PPP[i].Hsml);}
#else
    double soft = ForceSoftening_KernelRadius(i); /* a DM particle has no gas kernel of its own */
    if(soft > 0) {maxsoft = DMIN(maxsoft, 100.*soft);}
#endif
    return maxsoft;
}


#define CORE_FUNCTION_NAME pbhef_density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME pbhef_density_particle2in    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME pbhef_density_out2particle  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(pbhef_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/*! this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat HsmlPBH;
  int NodeList[NODELISTLENGTH];
}
 *DATAIN_NAME, *DATAGET_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static void pbhef_density_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    in->HsmlPBH = P[i].HsmlPBH;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    for(k=0;k<3;k++) {if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}}
}

/*! this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble NgbPBH;
    MyLongDouble RhoPBH;
    MyLongDouble DhsmlNgbPBH;
    MyLongDouble Particle_DivVelPBH;
}
 *DATARESULT_NAME, *DATAOUT_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static void pbhef_density_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(P[i].NumNgbPBH, out->NgbPBH, mode);
	ASSIGN_ADD(P[i].DensityPBH, out->RhoPBH, mode);
	ASSIGN_ADD(P[i].DhsmlNgbFactorPBH, out->DhsmlNgbPBH, mode);
	ASSIGN_ADD(P[i].Particle_DivVelPBH, out->Particle_DivVelPBH, mode);
}


/*! This function represents the core of the initial dm kernel-identification and volume computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine should in general contain no writes to shared memory. for optimization reasons, a couple of such writes have been included here in the sub-code for some sink routines -- those need to be handled with special care, both for thread safety and because of iteration. in general writes to shared memory in density.c are strongly discouraged -- */
int pbhef_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, k, n, startnode, numngb_inbox, listindex = 0; double r2, h2, u, mass_j;
    struct kernel_density kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0) {pbhef_density_particle2in(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    h2 = local.HsmlPBH * local.HsmlPBH; kernel_hinv(local.HsmlPBH, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlPBH, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, PBH_NGB_BITMASK);
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Mass <= 0) continue;
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                if(r2 < h2) /* this loop is only considering particles inside local.HsmlPBH, i.e. seen-by-main */
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);

                    out.NgbPBH += kernel.wk;
                    out.RhoPBH+= kernel.mj_wk;
                    out.DhsmlNgbPBH += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);

                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        if(P[j].Type == 0) {for(k=0;k<3;k++) {kernel.dv[k] = local.Vel[k] - SphP[j].VelPred[k];}} // VelPred is gas-only
                        else {for(k=0;k<3;k++) {kernel.dv[k] = local.Vel[k] - P[j].Vel[k];}}

                        NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,kernel.dv,1); /* wrap velocities for shearing boxes if needed */

                        out.Particle_DivVelPBH -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */

                    } // kernel.r > 0
                } // if(r2 < h2)
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }// while(startnode)
    if(mode == 0) {pbhef_density_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


/*! This function computes the local DM neighbor kernel for each active DM or hydro element, the number of neighbours in the current kernel radius, and the divergence
 * and rotation of the velocity field.  This is used then to compute the effective volume of the element in MFM/MFV/SPH-type methods, which is then used to
 * update volumetric quantities like density and pressure. The routine iterates to attempt to find a target kernel size set adaptively -- see code user guide for details
 */
void pbhef_density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_PBHEFDMDENSMISC] += measure_time(); double t00_truestart = my_second(); MyFloat *LeftPBH, *RightPBH; double fac, fac_lim, desnumngb, desnumngbdev; long long ntot;
    int i, k, npleft, iter=0, redo_particle, particle_set_to_minhsml_flag = 0, particle_set_to_maxhsml_flag = 0;
    LeftPBH = (MyFloat *) mymalloc("LeftPBH", NumPart * sizeof(MyFloat));
    RightPBH = (MyFloat *) mymalloc("RightPBH", NumPart * sizeof(MyFloat));

    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if(pbhef_density_isactive(i)) {

            P[i].NumNgbPBH = 0;
            LeftPBH[i] = RightPBH[i] = 0;

            double maxsoft = pbhef_return_maxhsml(i); /* before the first pass, need to ensure the particles do not exceed the maximum Hsml allowed */
            if((P[i].HsmlPBH <= 0) || !isfinite(P[i].HsmlPBH) || (P[i].HsmlPBH > 0.99*maxsoft)) {P[i].HsmlPBH = 0.99*maxsoft;} /* don't set to exactly maxsoft because our looping below won't treat this correctly */

        }} /* done with intial zero-out loop */

    /* allocate buffers to arrange communication */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        double tstart = my_second(), tend;

        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {

            if(pbhef_density_isactive(i))  // This makes sure that for method 1 (2), only gas particles (DM particles) are treated
            {
                if(P[i].NumNgbPBH > 0)
                {
                    P[i].DhsmlNgbFactorPBH *= P[i].HsmlPBH / (NUMDIMS * P[i].NumNgbPBH);
                    P[i].Particle_DivVelPBH /= P[i].NumNgbPBH;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    P[i].NumNgbPBH *= NORM_COEFF * pow(P[i].HsmlPBH,NUMDIMS);
                } else {
                    P[i].NumNgbPBH = P[i].DhsmlNgbFactorPBH = P[i].Particle_DivVelPBH = 0;
                }

                // inverse of fluid volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(P[i].DhsmlNgbFactorPBH > -0.9) {P[i].DhsmlNgbFactorPBH = 1 / (1 + P[i].DhsmlNgbFactorPBH);} else {P[i].DhsmlNgbFactorPBH = 1;} /* note: this would be -1 if only a single particle at zero lag is found */
                P[i].Particle_DivVelPBH *= P[i].DhsmlNgbFactorPBH;

                double minsoft = All.MinHsml;
                double maxsoft = pbhef_return_maxhsml(i);
                desnumngb = All.DesNumNgb; // Leave it constant
                /* the DM density only sets a heating rate that is linear in it, and a kernel of DesNumNgb neighbours
                   already carries ~1/sqrt(DesNumNgb) shot noise, so converging to All.MaxNumNgbDeviation (tuned for the
                   hydro volume partition, and ~0.16% for the usual settings) buys no accuracy and costs many iterations */
                desnumngbdev = DMAX(All.MaxNumNgbDeviation, 0.1*desnumngb);

                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}

                redo_particle = 0; // set redo_particle = 0, check if it needs to be set to 1 in the following

                /* check if we are in the 'normal' range between the max/min allowed values */
                if((P[i].NumNgbPBH < (desnumngb - desnumngbdev) && P[i].HsmlPBH < 0.999*maxsoft) ||
                   (P[i].NumNgbPBH > (desnumngb + desnumngbdev) && P[i].HsmlPBH > 1.001*minsoft))
                    {redo_particle = 1;}

                /* check maximum kernel size allowed */
                particle_set_to_maxhsml_flag = 0;
                if((P[i].HsmlPBH >= 0.999*maxsoft) && (P[i].NumNgbPBH < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P[i].HsmlPBH == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P[i].HsmlPBH = maxsoft;
                        particle_set_to_maxhsml_flag = 1;
                    }
                }

                /* check minimum kernel size allowed */
                particle_set_to_minhsml_flag = 0;
                if((P[i].HsmlPBH <= 1.001*minsoft) && (P[i].NumNgbPBH > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P[i].HsmlPBH == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P[i].HsmlPBH = minsoft;
                        particle_set_to_minhsml_flag = 1;
                    }
                }

                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("PBHEF loop parameters:\n i=%d task=%d ID=%llu iter=%d Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)",
                               i, ThisTask, (unsigned long long) P[i].ID, iter, P[i].Type, P[i].HsmlPBH, P[i].DhsmlNgbFactorPBH, LeftPBH[i], RightPBH[i],
                               (float) P[i].NumNgbPBH, RightPBH[i] - LeftPBH[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
						PRINT_WARNING("SLOW CONVERGENCE IN PBHEF LOOP!");
                    }

                    /* need to redo this particle */
                    npleft++;

                    if(LeftPBH[i] > 0 && RightPBH[i] > 0)
                        if((RightPBH[i] - LeftPBH[i]) < 1.0e-3 * LeftPBH[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            redo_particle = 0;
                            continue;
                        }

                    if((particle_set_to_maxhsml_flag==0)&&(particle_set_to_minhsml_flag==0))
                    {
                        if(P[i].NumNgbPBH < (desnumngb - desnumngbdev)) {LeftPBH[i] = DMAX(P[i].HsmlPBH, LeftPBH[i]);}
                        else
                        {
                            if(RightPBH[i] != 0) {if(P[i].HsmlPBH < RightPBH[i]) {RightPBH[i] = P[i].HsmlPBH;}} else {RightPBH[i] = P[i].HsmlPBH;}
                        }

                        // right/left define upper/lower bounds from previous iterations
                        if(RightPBH[i] > 0 && LeftPBH[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(RightPBH[i]/LeftPBH[i]);}
                            if(P[i].NumNgbPBH > 1)
                            {
                                double jumpvar = P[i].DhsmlNgbFactorPBH * log( desnumngb / P[i].NumNgbPBH ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                P[i].HsmlPBH *= exp(jumpvar);
                            } else {
                                P[i].HsmlPBH *= 2.0;
                            }
                            if((P[i].HsmlPBH<RightPBH[i])&&(P[i].HsmlPBH>LeftPBH[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(P[i].HsmlPBH > RightPBH[i] / hfac) {P[i].HsmlPBH = RightPBH[i] / hfac;}
                                    if(P[i].HsmlPBH < LeftPBH[i] * hfac) {P[i].HsmlPBH = LeftPBH[i] * hfac;}
                                }
                            } else {
                                if(P[i].HsmlPBH>RightPBH[i]) P[i].HsmlPBH=RightPBH[i];
                                if(P[i].HsmlPBH<LeftPBH[i]) P[i].HsmlPBH=LeftPBH[i];
                                P[i].HsmlPBH = pow(P[i].HsmlPBH * LeftPBH[i] * RightPBH[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(RightPBH[i] == 0 && LeftPBH[i] == 0)
                            {
                                char buf[1000]; sprintf(buf, "RightPBH[i] == 0 && LeftPBH[i] == 0 && P[i].HsmlPBH=%g\n", P[i].HsmlPBH); terminate(buf);
                            }

                            if(RightPBH[i] == 0 && LeftPBH[i] > 0)
                            {
                                if (P[i].NumNgbPBH > 1)
                                    {fac_lim = log( desnumngb / P[i].NumNgbPBH ) / NUMDIMS;} // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if((P[i].NumNgbPBH < 2*desnumngb)&&(P[i].NumNgbPBH > 0.1*desnumngb))
                                {
                                    double slope = P[i].DhsmlNgbFactorPBH;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P[i].DhsmlNgbFactorPBH==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac < fac_lim+0.231)
                                    {
                                        P[i].HsmlPBH *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        P[i].HsmlPBH *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    {P[i].HsmlPBH *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }

                            if(RightPBH[i] > 0 && LeftPBH[i] == 0)
                            {
                                if(P[i].NumNgbPBH > 1)
                                    {fac_lim = log( desnumngb / P[i].NumNgbPBH ) / NUMDIMS;} // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if(fac_lim < -1.535) {fac_lim = -1.535;} // decreasing N_ngb by factor ~100

                                if((P[i].NumNgbPBH < 2*desnumngb)&&(P[i].NumNgbPBH > 0.1*desnumngb))
                                {
                                    double slope = P[i].DhsmlNgbFactorPBH;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P[i].DhsmlNgbFactorPBH==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac > fac_lim-0.231)
                                    {
                                        P[i].HsmlPBH *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        {P[i].HsmlPBH *= exp(fac_lim-0.231);} // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    {P[i].HsmlPBH *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minhsml_flag], RightPBH && LeftPBH > 0
                    } // closes redo_particle, neither maxHsml or minHsml
                    /* resets for max/min values */
                    if(P[i].HsmlPBH < minsoft) {P[i].HsmlPBH = minsoft;}
                    if(particle_set_to_minhsml_flag==1) {P[i].HsmlPBH = minsoft;}
                    if(P[i].HsmlPBH > maxsoft) {P[i].HsmlPBH = maxsoft;}
                    if(particle_set_to_maxhsml_flag==1) {P[i].HsmlPBH = maxsoft;}
                } // redo particle
                else
                {
                    P[i].TimeBin = -P[i].TimeBin - 1; // Mark as inactive
                    redo_particle = 0;
                }
            } //  if(pbhef_density_isactive(i)), active particle
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i]), end DM loop

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 10) {PRINT_STATUS("PBHEF ngb iteration %d: need to repeat for %d%09d particles", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("PBHEF failed to converge in neighbour iteration in pbhef_density()\n"); fflush(stdout); endrun(1156);}
        }
    }
    while(ntot > 0);

    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(RightPBH); myfree(LeftPBH);

    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }


    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_PBHEFDMDENSCOMPUTE] += timecomp; CPU_Step[CPU_PBHEFDMDENSWAIT] += timewait;
    CPU_Step[CPU_PBHEFDMDENSCOMM] += timecomm; CPU_Step[CPU_PBHEFDMDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

/*! This function computes the alpha coefficient for PBH evaporation, according to the Mosbech et al. (2022) analytical fit.
    alpha_eff is defined there (eq. 9) as (G^2/hbar c^4) * M0^3 / (3 tau), with tau the numerically computed lifetime, so it
    is an effective coefficient averaged over the life of the black hole and is a function of the -initial- mass alone. That
    is what makes the constant-alpha solution M(t) = (M0^3 - 3 alpha hbar c^4 t / G^2)^(1/3) give the correct lifetime, and
    it is why this must not be re-evaluated on the current mass as the black hole evaporates. */
double pbhef_calculate_alpha(double m_pbh_initial_grams)
{
    const double alpha_highmass = 2.011e-4; /* the large-mass limit of the fit */
    double alpha = alpha_highmass;
    if (m_pbh_initial_grams < 1.0e18)
    {
        const double c1 = -0.3015;
        const double c2 = 0.3113;
        const double p_exponent = -0.0008;
        alpha = c1 + c2 * pow(m_pbh_initial_grams, p_exponent);
    }
    /* the low-mass branch is a small difference of two numbers near 0.3: it already falls below the large-mass limit at
       about 1.0e17 g and turns negative above about 2.3e17 g, well before the 1e18 g switch quoted in the paper. Since
       alpha_eff = M0^3/(3 tau) is positive by construction, take the larger of the two branches, which puts the switch at
       the point where they cross and leaves the fit continuous and always positive. */
    return DMAX(alpha, alpha_highmass);
}


#ifndef PBHEF_NO_MASS_LOSS
/*! Integrand for the elapsed cosmic time, dt = da / (a H(a)) */
static double pbhef_time_integ(double a, void *param)
{
    return 1 / (hubble_function(a) * a);
}
#endif


/*! Rate of change of the cube of the PBH mass, which is a constant.
 *  From dM/dt = -(hbar c^4/G^2) alpha / M^2 we get d(M^3)/dt = -3 (hbar c^4/G^2) alpha, and alpha is held fixed at its
 *  initial-mass value (see pbhef_calculate_alpha above), so M^3 is exactly linear in cosmic time:
 *      M(t)^3 = M0^3 - 3 (hbar c^4/G^2) alpha t
 *  All.PBH_EvaporationConstant holds hbar c^6/G^2, so divide by c^2 to get hbar c^4/G^2. */
static double pbhef_mass3_decay_rate(void)
{
    return 3.0 * (All.PBH_EvaporationConstant / (C_LIGHT_CODE * C_LIGHT_CODE)) * All.PBH_Alpha;
}


/*! Tabulate the cosmic time elapsed since All.TimeBegin as a function of the scale factor. That is the only quantity that
 *  needs tabulating: given it, the PBH mass follows in closed form (see pbhef_mass3_decay_rate above), so there is no need
 *  to integrate the mass loss itself. The grid is uniform in log(a) and the integration is done with the same GSL routine
 *  and the same conventions as init_drift_table(), so the lookup below mirrors get_drift_factor(). */
void pbhef_init_mass_evolution(void)
{
#ifndef PBHEF_NO_MASS_LOSS
    if(All.ComovingIntegrationOn)
    {
#define PBH_WORKSIZE 100000
        int i; double result, abserr;
        double logTimeBegin = log(All.TimeBegin), logTimeMax = log(All.TimeMax);
        gsl_function F; gsl_integration_workspace *workspace;
        workspace = gsl_integration_workspace_alloc(PBH_WORKSIZE);
        F.function = &pbhef_time_integ;
        for(i = 0; i < PBH_TABLE_SIZE; i++)
        {
            gsl_integration_qag(&F, exp(logTimeBegin),
                                exp(logTimeBegin + ((logTimeMax - logTimeBegin) / PBH_TABLE_SIZE) * (i + 1)), 0,
                                1.0e-8, PBH_WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
            PBH_Table_Time[i] = result;
        }
        gsl_integration_workspace_free(workspace);
#undef PBH_WORKSIZE
    }

    if(ThisTask == 0)
    {
        double final_mass; pbhef_get_current_mass(All.TimeMax, &final_mass);
        printf("PBH mass loss: initial mass %g, mass at the end of the run (a=%g) %g [code units]\n", All.PBH_InitialMass, All.TimeMax, final_mass);
        if(final_mass <= 0) {printf("PBH mass loss: the black holes evaporate completely before the end of the run, after which there is no heating.\n");}
        printf("\n");
    }
#else
    if(ThisTask == 0) {printf("PBH mass loss not enabled: the PBH mass is fixed at the initial value.\n\n");}
#endif
}

/*! Constant part of the heating rate. The rate per unit gas mass is
 *  f0 * (rho_dm/rho_gas) * C * alpha / (m0 * m(t)^2), and everything except the local density ratio is the
 *  same for every cell on a given step, so this is evaluated once per step rather than once per cell.
 *  The power of m(t) is -2 rather than -3 because the number of black holes is conserved as they lose
 *  mass, so the mass fraction in black holes, f(t), itself scales as m(t)/m0. The donor-based method uses
 *  the same quantity: times a DM particle's mass it gives the power of the black holes it carries.
 */
double pbhef_heating_prefactor(void)
{
    if(!pbhef_is_active()) {return 0;}
    double current_pbh_mass; pbhef_get_current_mass(All.Time, &current_pbh_mass);
    return All.PBH_MassFraction * All.PBH_EvaporationConstant * All.PBH_Alpha / (All.PBH_InitialMass * current_pbh_mass * current_pbh_mass);
}

/*! the prefactor for the current time, remembered between calls: the donor loops want it per particle, and
    evaluating it means a lookup of the black hole mass. Same static caching get_drift_factor() uses. Safe
    because the callers are either unthreaded or ask for it once before their threads start */
double pbhef_prefactor_cached(void)
{
    static double time_last = -1, prefactor_last = 0;
    if(All.Time != time_last) {time_last = All.Time; prefactor_last = pbhef_heating_prefactor();}
    return prefactor_last;
}



#if (PBHEF == 1)
/*! Receiver-based energy injection for a single gas cell: adds the PBH evaporation heating to its internal
 *  energy rate. Called from hydro_final_operations_and_cleanup(), which must have converted
 *  DtInternalEnergy to specific energy units first, and which zeroes it again for cells that are frozen or
 *  decoupled further down the same loop.
 */
void pbhef_receiver_inject(int i, double heating_prefactor)
{
    if(SphP[i].Density <= 0) {return;}
    double dm_dens_over_gas_dens = P[i].DensityPBH / SphP[i].Density; // dimensionless ratio of DM density to gas density
    double heat_source = heating_prefactor * dm_dens_over_gas_dens;

    SphP[i].DtInternalEnergy += heat_source; // add internal energy created by PBH evaporation
    SphP[i].PBHEF_Dtu += heat_source;

#ifdef PBHEF_DEBUG
    if(P[i].ID == All.PBH_EnergyID) {
        double current_pbh_mass; pbhef_get_current_mass(All.Time, &current_pbh_mass);
        printf(" ..PBHEF (after injection): i=%d, Type=%d, ID=%llu, atime=%g, InternalEnergy=%g, rhoDM=%g, rho=%g,\n"
               "                            dm_dens_over_gas_dens=%g, f=%g, C=%g, alpha=%g,\n"
               "                            PBHm0=%g, PBHm(t)=%g, heat_source=%g, PBHEF_Dtu=%g, DtInternalEnergy=%g\n",
          i, P[i].Type, (unsigned long long) P[i].ID, All.cf_atime, SphP[i].InternalEnergy, P[i].DensityPBH*All.cf_a3inv, SphP[i].Density*All.cf_a3inv,
          dm_dens_over_gas_dens, All.PBH_MassFraction, All.PBH_EvaporationConstant, All.PBH_Alpha,
          All.PBH_InitialMass, current_pbh_mass, heat_source, SphP[i].PBHEF_Dtu, SphP[i].DtInternalEnergy);
        fflush(stdout);
    }
#endif
}
#endif


/*! Returns 0 if the PBH evaporation heating rate is identically zero for this step, so that callers can
 *  skip the DM neighbour search and the energy injection altogether. This happens if the module is enabled
 *  but given a zero PBH mass fraction, or once the black holes have fully evaporated. Note that the DM
 *  density is then left at its last computed value rather than being updated.
 */
int pbhef_is_active(void)
{
    if(All.PBH_MassFraction <= 0) {return 0;}
    if(All.PBH_Alpha <= 0) {return 0;}
    double current_pbh_mass; pbhef_get_current_mass(All.Time, &current_pbh_mass);
    if(current_pbh_mass <= 0) {return 0;}
    return 1;
}


/*! Elapsed cosmic time since All.TimeBegin, at scale factor a. In a non-cosmological run All.Time already is the time,
 *  so this is just the difference. In a cosmological run it comes from the table built above, using the same log(a)
 *  indexing as get_drift_factor(). */
static double pbhef_elapsed_time(double a)
{
#ifndef PBHEF_NO_MASS_LOSS
    if(!All.ComovingIntegrationOn) {return DMAX(a - All.TimeBegin, 0);}

    double logTimeBegin = log(All.TimeBegin), logTimeMax = log(All.TimeMax);
    if(a <= All.TimeBegin || logTimeMax <= logTimeBegin) {return 0;}

    double u = (log(a) - logTimeBegin) / (logTimeMax - logTimeBegin) * PBH_TABLE_SIZE;
    int i1 = (int) u; if(i1 >= PBH_TABLE_SIZE) {i1 = PBH_TABLE_SIZE - 1;}
    if(i1 <= 1) {return u * PBH_Table_Time[0];}
    return PBH_Table_Time[i1 - 1] + (PBH_Table_Time[i1] - PBH_Table_Time[i1 - 1]) * (u - i1);
#else
    return 0;
#endif
}


/*! Current PBH mass. M^3 is exactly linear in cosmic time for a fixed alpha, so this is evaluated in closed form rather
 *  than interpolated: only the time itself comes from a table. A black hole that has evaporated returns a mass of zero. */
void pbhef_get_current_mass(double a, double *mass_out)
{
#ifndef PBHEF_NO_MASS_LOSS
    double mass3 = All.PBH_InitialMass * All.PBH_InitialMass * All.PBH_InitialMass - pbhef_mass3_decay_rate() * pbhef_elapsed_time(a);
    if(mass3 <= 0) {*mass_out = 0;} else {*mass_out = cbrt(mass3);}
#else
    *mass_out = All.PBH_InitialMass;
#endif
}


#if (PBHEF == 2)

/*! Donor-based energy injection, after List et al. 2019. A DM particle carries f0*M_i/m0 black holes, each
    radiating C*alpha/m(t)^2, so its power is the heating prefactor times its own mass, with no dependence
    on the DM density around it. That power is shared over its gas neighbors with the mass-weighted kernel
    weights of eq. 5, normalised by their sum, which is the gas density pbhef_density() already put in
    DensityPBH; per unit receiver mass the share is therefore W(r,h)/DensityPBH and the shares sum to one.
    The share is stored as a rate, in the slot for the donor's own timebin, because donor and receivers are
    usually not active on the same step. Slots are cleared for the bins active on this step, which is safe
    since every donor in an active bin is about to deposit again. */

static double PBH_RateIntended;   /*!< rate the active donors mean to inject, this task */
static double PBH_RateCoupled;    /*!< rate that actually landed on gas, this task */
static long PBH_NumUncapped;      /*!< active receivers found on a longer step than their donor, this task */

/*! total power of the black holes carried by DM particle i */
double pbhef_donor_power(int i)
{
    return pbhef_prefactor_cached() * P[i].Mass;
}

/*! DM particles that are active and have gas to give their energy to */
static int pbhef_donor_isactive(int n)
{
    if(!pbhef_density_isactive(n)) {return 0;}
    if(P[n].HsmlPBH <= 0 || P[n].DensityPBH <= 0) {return 0;} /* no gas inside the kernel, so no receivers */
    if(pbhef_donor_power(n) <= 0) {return 0;}
    return 1;
}

/*! gas cells that may take the energy: the same ones hydro_final_operations_and_cleanup() does not zero.
    A skipped neighbor's share is not handed to the others, so it shows up in the budget report below */
static int pbhef_donor_can_receive(int j)
{
    if(P[j].Type != 0 || P[j].Mass <= 0) {return 0;}
#ifdef BOX_BND_PARTICLES
    if(P[j].ID == 0) {return 0;} /* frozen boundary particle */
#endif
#ifdef GALSF_SUBGRID_WINDS
    if(SphP[j].DelayTime > 0) {return 0;} /* briefly decoupled wind particle */
#endif
    return 1;
}


#define CORE_FUNCTION_NAME pbhef_donor_evaluate /* name of the 'core' function doing the actual inter-neighbor operations */
#define INPUTFUNCTION_NAME pbhef_donor_particle2in    /* name of the function which loads the element data needed */
#define OUTPUTFUNCTION_NAME pbhef_donor_out2particle  /* name of the function which takes the data returned from other processors */
#define CONDITIONFUNCTION_FOR_EVALUATION if(pbhef_donor_isactive(i)) /* elements which are 'active' for this loop */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below */

/*! variables sent -from- the donating DM particle */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3];
    MyFloat HsmlPBH;
    MyDouble Power;       /* total power of the black holes this DM particle carries */
    MyDouble GasDensity;  /* normalisation: the mass-weighted kernel sum over its gas neighbors */
    int Timebin;          /* which slot the receivers should store the rate in */
    int NodeList[NODELISTLENGTH];
}
 *DATAIN_NAME, *DATAGET_NAME;

static void pbhef_donor_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->HsmlPBH = P[i].HsmlPBH;
    in->GasDensity = P[i].DensityPBH;
    in->Power = pbhef_donor_power(i);
    in->Timebin = P[i].TimeBin;
}

/*! variables sent -back to- the donating DM particle */
static struct OUTPUT_STRUCT_NAME
{
    MyDouble RateCoupled; /* rate that actually landed on the gas, for the budget report */
    int MaxTimebin;       /* largest bin this donor may take, from its receivers */
}
 *DATARESULT_NAME, *DATAOUT_NAME;

/*! nothing needs to go back onto the DM particle itself: just accumulate the budget */
static void pbhef_donor_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    #pragma omp atomic
    PBH_RateCoupled += out->RateCoupled;
#ifdef PBHEF_LIMIT_DM_TIMESTEP
    if(out->MaxTimebin < P[i].PBHEF_MaxTimebin) {P[i].PBHEF_MaxTimebin = out->MaxTimebin;}
#endif
}

/*! gives each gas neighbor of one DM particle its share of the evaporation rate. NOTE unlike the density
    loops here this writes to -other- elements, so those writes are atomic as in mechanical_fb.c */
int pbhef_donor_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, startnode, numngb_inbox, listindex = 0, n_uncapped = 0; double dp[3], r2, h2, r, u, wk, dwk, hinv, hinv3, hinv4;
    struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    out.MaxTimebin = TIMEBINS;
    if(mode == 0) {pbhef_donor_particle2in(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Power <= 0 || local.GasDensity <= 0 || local.HsmlPBH <= 0)
    {
        if(mode == 0) {pbhef_donor_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
        return 0;
    }
    h2 = local.HsmlPBH * local.HsmlPBH; kernel_hinv(local.HsmlPBH, &hinv, &hinv3, &hinv4);
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode; /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlPBH, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 1); // gas only
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(!pbhef_donor_can_receive(j)) {continue;}
                dp[0] = local.Pos[0] - P[j].Pos[0];
                dp[1] = local.Pos[1] - P[j].Pos[1];
                dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(dp[0],dp[1],dp[2],1);
                r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
                if(r2 < h2) /* same condition as the density loop, so the weights sum to GasDensity */
                {
                    r = sqrt(r2); u = r * hinv;
                    kernel_main(u, hinv3, hinv4, &wk, &dwk, -1); /* -1: only wk is needed */
                    if(wk <= 0) {continue;}

                    double dtu_j = local.Power * wk / local.GasDensity; /* specific energy rate for this cell */

                    /* A receiver on a longer step than its donor would integrate this rate for too long,
                       so scale it down by the ratio of the two steps. The cap below is meant to prevent
                       that, but it only takes effect from the next step, so the case is reachable. */
                    if(local.Timebin < P[j].TimeBin)
                    {
                        dtu_j *= (double)(((integertime) 1) << local.Timebin) / (double)(((integertime) 1) << P[j].TimeBin);
                        if(TimeBinActive[P[j].TimeBin]) {n_uncapped++;} /* cap should have applied by now */
                    }

                    #pragma omp atomic
                    SphP[j].PBHEF_DtuBin[local.Timebin] += (float)dtu_j;

                    /* keep the receiver on a step no longer than the donor's */
                    if(P[j].PBHEF_MaxTimebin > local.Timebin)
                    {
                        #pragma omp critical
                        {if(P[j].PBHEF_MaxTimebin > local.Timebin) {P[j].PBHEF_MaxTimebin = local.Timebin;}}
                    }
#ifdef PBHEF_LIMIT_DM_TIMESTEP /* and, optionally, the donor within a few bins of its receivers */
                    {int cap = P[j].TimeBin + PBHEF_LIMIT_DM_TIMESTEP; if(cap < out.MaxTimebin) {out.MaxTimebin = cap;}}
#endif

                    out.RateCoupled += dtu_j * P[j].Mass;
                } // if(r2 < h2)
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    } // while(startnode)
    if(n_uncapped) {
        #pragma omp atomic
        PBH_NumUncapped += n_uncapped;
    }
    if(mode == 0) {pbhef_donor_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


/*! top-level driver: must run after pbhef_density(), which sets HsmlPBH and DensityPBH */
void pbhef_donor_feedback(void)
{
    int i, b;
    /* clear the slots whose donors will deposit again. must cover every gas cell, not just the active
       ones, since an active donor also feeds sleeping receivers */
    for(i = 0; i < N_gas; i++) {for(b = 0; b < TIMEBINS; b++) {if(TimeBinActive[b]) {SphP[i].PBHEF_DtuBin[b] = 0;}}}

    double prefactor = pbhef_prefactor_cached(); /* also warms the cache before the threaded loop below */
    PBH_RateCoupled = 0; PBH_RateIntended = 0; PBH_NumUncapped = 0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(pbhef_density_isactive(i)) {PBH_RateIntended += prefactor * P[i].Mass;}
    }

    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */

    /* total the slots into the rate the rest of the code reads; only occupied bins can hold anything.
       The total is then coupled the way the receiver mode couples its own, by adding to DtInternalEnergy
       so that kicks.c can fold it into the implicit cooling solve. Only cells that are about to
       integrate their rate may take it: hydro_force() has just rebuilt DtInternalEnergy for those, while
       an inactive cell still carries the value it was given when it last woke and would be paid twice.
       Cells the hydro loop excludes are already absent here, since pbhef_donor_can_receive() rejects
       them before anything is deposited. */
    for(i = 0; i < N_gas; i++)
    {
        double dtu = 0;
        for(b = 0; b < TIMEBINS; b++) {if(TimeBinCount[b]) {dtu += SphP[i].PBHEF_DtuBin[b];}}
        SphP[i].PBHEF_Dtu = dtu;
        if(dtu > 0 && TimeBinActive[P[i].TimeBin]) {SphP[i].DtInternalEnergy += dtu;}
    }

    /* the shares sum to one, so these should agree to round-off. this is not conservative by
       construction though: it relies on the slot bookkeeping and the capping, so it is worth watching */
    double budget_local[2] = {PBH_RateIntended, PBH_RateCoupled}, budget[2];
    MPI_Allreduce(&budget_local[0], &budget[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(budget[0] > 0) {PRINT_STATUS(" ..PBHEF donor budget: intended=%g coupled=%g (fraction=%.6f)", budget[0], budget[1], budget[1]/budget[0]);}

    /* an active receiver on a longer step than its donor means the cap has not held. its rate is scaled
       down so no energy is lost, but it is a sign the coupling is not doing its job */
    long nbad_local = PBH_NumUncapped, nbad; MPI_Allreduce(&nbad_local, &nbad, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if(nbad) {PRINT_WARNING("PBHEF donor: %ld active gas cells were on a longer timestep than a donor feeding them", nbad);}

    CPU_Step[CPU_PBHEFDONOR] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif /* PBHEF == 2 */


#endif /* PBHEF */
