#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
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
 * Method 1 - receiver-based approach (activate using PBH_EVAPORATION_FEEDBACK)
  - PBH (traced by N-Body DM particles) density is calculated at each gas particle using smoothing length HsmlDM
  - from PBH density, PBH evaporation rate at gas particle is calculated
  - energy injection at each gas particle
 * Method 2 - donor-based approach, also uses the functionality in this file in order to determine the PBH density around DM particles.
 *
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially (condensed, different criteria for kernel lengths, optimizations,
 * rewritten parallelism, new physics included, new variable/memory conventions added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#if defined(PBH_EVAPORATION_FEEDBACK) || defined(PBH_EVAPORATION_FEEDBACK_DM)

struct kernel_density /*! defines a number of useful variables we will use below */
{
  double dp[3],dv[3],r, wk, dwk, hinv, hinv3, hinv4, mj_wk, mj_dwk_r;
};


/*! routine to determine if a given element is actually going to be active in the density subroutines below to calculate HsmlDM and rhoDM */
int dm_density_isactive(int n)
{

    if(P[n].TimeBin < 0) {return 0;}
    if(P[n].Mass <= 0) {return 0;}
#if defined(PBH_EVAPORATION_FEEDBACK)
    if(P[n].Type != 0){return 0;}  /* only gas particles */
#elif defined(PBH_EVAPORATION_FEEDBACK_DM)
	if(P[n].Type != 1){return 0;}  /* only DM particles */
#endif
    return 1;
}


#define CORE_FUNCTION_NAME dm_density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME dmkerneldensity_particle2in    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME dmkerneldensity_out2particle  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(dm_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'dm_density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/*! this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat HsmlDM;
  int NodeList[NODELISTLENGTH];
}
 *DATAIN_NAME, *DATAGET_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void dmkerneldensity_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    in->HsmlDM = P[i].HsmlDM;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    for(k=0;k<3;k++) {if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}}
}

/*! this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble NgbDM;
    MyLongDouble RhoDM;
    MyLongDouble DhsmlNgbDM;
    MyLongDouble Particle_DivVelDM;
}
 *DATARESULT_NAME, *DATAOUT_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void dmkerneldensity_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(P[i].NumNgbDM, out->NgbDM, mode);
	ASSIGN_ADD(P[i].DensityDM, out->RhoDM, mode);
	ASSIGN_ADD(P[i].DhsmlNgbFactorDM, out->DhsmlNgbDM, mode);
	ASSIGN_ADD(P[i].Particle_DivVelDM, out->Particle_DivVelDM, mode);
}


/*! This function represents the core of the initial dm kernel-identification and volume computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine should in general contain no writes to shared memory. for optimization reasons, a couple of such writes have been included here in the sub-code for some sink routines -- those need to be handled with special care, both for thread safety and because of iteration. in general writes to shared memory in density.c are strongly discouraged -- */
int dm_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, startnode, numngb_inbox, listindex = 0; double r2, h2, u, mass_j;
    struct kernel_density kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0) {dmkerneldensity_particle2in(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    h2 = local.HsmlDM * local.HsmlDM; kernel_hinv(local.HsmlDM, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlDM, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 2); // search for DM particles only
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
                if(r2 < h2) /* this loop is only considering particles inside local.HsmlDM, i.e. seen-by-main */
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);

                    out.NgbDM += kernel.wk;
                    out.RhoDM+= kernel.mj_wk;
                    out.DhsmlNgbDM += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);

                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        kernel.dv[0] = local.Vel[0] - P[j].Vel[0];  // we use Vel here to estimate divergence since VelPred is not calculated for DM particles
                        kernel.dv[1] = local.Vel[1] - P[j].Vel[1];
                        kernel.dv[2] = local.Vel[2] - P[j].Vel[2];

                        NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,kernel.dv,1); /* wrap velocities for shearing boxes if needed */

                        out.Particle_DivVelDM -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */

                    } // kernel.r > 0
                } // if(r2 < h2)
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }// while(startnode)
    if(mode == 0) {dmkerneldensity_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


/*! This function computes the local DM neighbor kernel for each active DM or hydro element, the number of neighbours in the current kernel radius, and the divergence
 * and rotation of the velocity field.  This is used then to compute the effective volume of the element in MFM/MFV/SPH-type methods, which is then used to
 * update volumetric quantities like density and pressure. The routine iterates to attempt to find a target kernel size set adaptively -- see code user guide for details
 */
void dm_density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_PBHEFDMDENSMISC] += measure_time(); double t00_truestart = my_second(); MyFloat *LeftDM, *RightDM; double fac, fac_lim, desnumngb, desnumngbdev; long long ntot;
    int i, k, npleft, iter=0, redo_particle, particle_set_to_minhsmlDM_flag = 0, particle_set_to_maxhsmlDM_flag = 0;
    LeftDM = (MyFloat *) mymalloc("LeftDM", NumPart * sizeof(MyFloat));
    RightDM = (MyFloat *) mymalloc("RightDM", NumPart * sizeof(MyFloat));

    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if(dm_density_isactive(i)) {

            P[i].NumNgbDM = 0;
            LeftDM[i] = RightDM[i] = 0;

            double maxsoft = All.MaxHsml; /* before the first pass, need to ensure the particles do not exceed the maximum Hsml allowed */
            if((P[i].HsmlDM < 0) || !isfinite(P[i].HsmlDM) || (P[i].HsmlDM > 0.99*maxsoft)) {P[i].HsmlDM = 0.99*maxsoft;} /* don't set to exactly maxsoft because our looping below won't treat this correctly */

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

            if(dm_density_isactive(i))  // This makes sure that for method 1 (2), only gas particles (DM particles) are treated
            {
                if(P[i].NumNgbDM > 0)
                {
                    P[i].DhsmlNgbFactorDM *= P[i].HsmlDM / (NUMDIMS * P[i].NumNgbDM);
                    P[i].Particle_DivVelDM /= P[i].NumNgbDM;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    P[i].NumNgbDM *= NORM_COEFF * pow(P[i].HsmlDM,NUMDIMS);
                } else {
                    P[i].NumNgbDM = P[i].DhsmlNgbFactorDM = P[i].Particle_DivVelDM = 0;
                }

                // inverse of fluid volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(P[i].DhsmlNgbFactorDM > -0.9) {P[i].DhsmlNgbFactorDM = 1 / (1 + P[i].DhsmlNgbFactorDM);} else {P[i].DhsmlNgbFactorDM = 1;} /* note: this would be -1 if only a single particle at zero lag is found */
                P[i].Particle_DivVelDM *= P[i].DhsmlNgbFactorDM;

                double minsoft = All.MinHsml;
                double maxsoft = All.MaxHsml;
                desnumngb = All.DesNumNgb; // Leave it constant
				desnumngbdev = All.MaxNumNgbDeviation;

                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}

                redo_particle = 0; // set redo_particle = 0, check if it needs to be set to 1 in the following

                /* check if we are in the 'normal' range between the max/min allowed values */
                if((P[i].NumNgbDM < (desnumngb - desnumngbdev) && P[i].HsmlDM < 0.999*maxsoft) ||
                   (P[i].NumNgbDM > (desnumngb + desnumngbdev) && P[i].HsmlDM > 1.001*minsoft))
                    {redo_particle = 1;}

                /* check maximum kernel size allowed */
                particle_set_to_maxhsmlDM_flag = 0;
                if((P[i].HsmlDM >= 0.999*maxsoft) && (P[i].NumNgbDM < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P[i].HsmlDM == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsmlDM_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P[i].HsmlDM = maxsoft;
                        particle_set_to_maxhsmlDM_flag = 1;
                    }
                }

                /* check minimum kernel size allowed */
                particle_set_to_minhsmlDM_flag = 0;
                if((P[i].HsmlDM <= 1.001*minsoft) && (P[i].NumNgbDM > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P[i].HsmlDM == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsmlDM_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P[i].HsmlDM = minsoft;
                        particle_set_to_minhsmlDM_flag = 1;
                    }
                }

                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("PBHEF loop parameters:\n i=%d task=%d ID=%llu iter=%d Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)",
                               i, ThisTask, (unsigned long long) P[i].ID, iter, P[i].Type, P[i].HsmlDM, P[i].DhsmlNgbFactorDM, LeftDM[i], RightDM[i],
                               (float) P[i].NumNgbDM, RightDM[i] - LeftDM[i], particle_set_to_maxhsmlDM_flag, particle_set_to_minhsmlDM_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
						PRINT_WARNING("SLOW CONVERGENCE IN PBHEF LOOP!");
                    }

                    /* need to redo this particle */
                    npleft++;

                    if(LeftDM[i] > 0 && RightDM[i] > 0)
                        if((RightDM[i] - LeftDM[i]) < 1.0e-3 * LeftDM[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            redo_particle = 0;
                            continue;
                        }

                    if((particle_set_to_maxhsmlDM_flag==0)&&(particle_set_to_minhsmlDM_flag==0))
                    {
                        if(P[i].NumNgbDM < (desnumngb - desnumngbdev)) {LeftDM[i] = DMAX(P[i].HsmlDM, LeftDM[i]);}
                        else
                        {
                            if(RightDM[i] != 0) {if(P[i].HsmlDM < RightDM[i]) {RightDM[i] = P[i].HsmlDM;}} else {RightDM[i] = P[i].HsmlDM;}
                        }

                        // right/left define upper/lower bounds from previous iterations
                        if(RightDM[i] > 0 && LeftDM[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(RightDM[i]/LeftDM[i]);}
                            if(P[i].NumNgbDM > 1)
                            {
                                double jumpvar = P[i].DhsmlNgbFactorDM * log( desnumngb / P[i].NumNgbDM ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                P[i].HsmlDM *= exp(jumpvar);
                            } else {
                                P[i].HsmlDM *= 2.0;
                            }
                            if((P[i].HsmlDM<RightDM[i])&&(P[i].HsmlDM>LeftDM[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(P[i].HsmlDM > RightDM[i] / hfac) {P[i].HsmlDM = RightDM[i] / hfac;}
                                    if(P[i].HsmlDM < LeftDM[i] * hfac) {P[i].HsmlDM = LeftDM[i] * hfac;}
                                }
                            } else {
                                if(P[i].HsmlDM>RightDM[i]) P[i].HsmlDM=RightDM[i];
                                if(P[i].HsmlDM<LeftDM[i]) P[i].HsmlDM=LeftDM[i];
                                P[i].HsmlDM = pow(P[i].HsmlDM * LeftDM[i] * RightDM[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(RightDM[i] == 0 && LeftDM[i] == 0)
                            {
                                char buf[1000]; sprintf(buf, "RightDM[i] == 0 && LeftDM[i] == 0 && P[i].HsmlDM=%g\n", P[i].HsmlDM); terminate(buf);
                            }

                            if(RightDM[i] == 0 && LeftDM[i] > 0)
                            {
                                if (P[i].NumNgbDM > 1)
                                    {fac_lim = log( desnumngb / P[i].NumNgbDM ) / NUMDIMS;} // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if((P[i].NumNgbDM < 2*desnumngb)&&(P[i].NumNgbDM > 0.1*desnumngb))
                                {
                                    double slope = P[i].DhsmlNgbFactorDM;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P[i].DhsmlNgbFactorDM==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac < fac_lim+0.231)
                                    {
                                        P[i].HsmlDM *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        P[i].HsmlDM *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    {P[i].HsmlDM *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }

                            if(RightDM[i] > 0 && LeftDM[i] == 0)
                            {
                                if(P[i].NumNgbDM > 1)
                                    {fac_lim = log( desnumngb / P[i].NumNgbDM ) / NUMDIMS;} // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if(fac_lim < -1.535) {fac_lim = -1.535;} // decreasing N_ngb by factor ~100

                                if((P[i].NumNgbDM < 2*desnumngb)&&(P[i].NumNgbDM > 0.1*desnumngb))
                                {
                                    double slope = P[i].DhsmlNgbFactorDM;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P[i].DhsmlNgbFactorDM==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac > fac_lim-0.231)
                                    {
                                        P[i].HsmlDM *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        {P[i].HsmlDM *= exp(fac_lim-0.231);} // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    {P[i].HsmlDM *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minhsml_flag], RightDM && LeftDM > 0
                    } // closes redo_particle, neither maxHsml or minHsml
                    /* resets for max/min values */
                    if(P[i].HsmlDM < minsoft) {P[i].HsmlDM = minsoft;}
                    if(particle_set_to_minhsmlDM_flag==1) {P[i].HsmlDM = minsoft;}
                    if(P[i].HsmlDM > maxsoft) {P[i].HsmlDM = maxsoft;}
                    if(particle_set_to_maxhsmlDM_flag==1) {P[i].HsmlDM = maxsoft;}
                } // redo particle
                else
                {
                    P[i].TimeBin = -P[i].TimeBin - 1; // Mark as inactive
                    redo_particle = 0;
                }
            } //  if(dm_density_isactive(i)), active particle
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i]), end DM loop

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 10) {PRINT_STATUS("PBHEF ngb iteration %d: need to repeat for %d%09d particles", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("PBHEF failed to converge in neighbour iteration in dm_density()\n"); fflush(stdout); endrun(1156);}
        }
    }
    while(ntot > 0);

    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(RightDM); myfree(LeftDM);

    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }


    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results
     ( any quantities that only need to be evaluated once, on the final iteration --
     won't save much b/c the real cost is in the neighbor loop for each particle, but it's something )
     -- also, some results (for example, viscosity suppression below) should not be calculated unless
     the quantities are 'stabilized' at their final values -- */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(dm_density_isactive(i))
        {
         /* finally, convert NGB to the more useful format, NumNgb^(1/NDIMS),
            which we can use to obtain the corrected particle sizes. Because of how this number is used above, we --must-- make
            sure that this operation is the last in the loop here */
            if(P[i].NumNgbDM > 0) {P[i].NumNgbDM=pow(P[i].NumNgbDM,1./NUMDIMS);} else {P[i].NumNgbDM=0;}

        } // dm_density_isactive(i)


    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])

    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_PBHEFDMDENSCOMPUTE] += timecomp; CPU_Step[CPU_PBHEFDMDENSWAIT] += timewait;
    CPU_Step[CPU_PBHEFDMDENSCOMM] += timecomm; CPU_Step[CPU_PBHEFDMDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

/*! This function computes the alpha coefficient for PBH evaporation, according to the Mosbech et al. (2022) analytical fit */
double calculate_alpha(double m_pbh_initial_grams)
{
    if (m_pbh_initial_grams < 1.0e18)
    {
        const double c1 = -0.3015;
        const double c2 = 0.3113;
        const double p_exponent = -0.0008;
        return c1 + c2 * pow(m_pbh_initial_grams, p_exponent);
    }
    else {return 2.011e-4;}
}


/*
 * Initialize the PBH mass evolution lookup table.
 * This function integrates the mass loss rate over the simulation time.
 */
void init_pbh_mass_evolution(void)
{
#ifndef PBH_EVAPORATION_FEEDBACK_NO_MASS_LOSS
    if(ThisTask == 0) printf("Initializing PBH mass evolution table...\n");

    int i;
    double current_mass = All.PBH_InitialMass;
    double current_mass3 = current_mass * current_mass * current_mass;
    double current_a = All.TimeBegin;

    // Determine integration limits and step
    double a_start = All.TimeBegin;
    double a_end = All.TimeMax;
    double da = (a_end - a_start) / (double)(PBH_TABLE_SIZE - 1);

    // Calculate Mass Loss Constant (M^3 rate)
    // The code calculates All.PBH_EvaporationConstant = (hbar * c^6 / G^2) in code units.
    // We need K_mass3 such that d(M^3)/dt = -3 * K_mass_loss * alpha
    // dM/dt = - (hbar * c^4 / G^2) * alpha / M^2
    // M^2 dM = - (hbar * c^4 / G^2) * alpha dt
    // integrated: M^3/3 = - (hbar * c^4 / G^2) * alpha * t + C
    // d(M^3)/dt = -3 * (hbar * c^4 / G^2) * alpha
    // So K_mass3 = 3 * (hbar * c^4 / G^2) * alpha

    // We need (hbar * c^4 / G^2) in code units.
    // All.PBH_EvaporationConstant has c^6. So divide by c^2.

    double K_BH_const_code = All.PBH_EvaporationConstant / (C_LIGHT_CODE * C_LIGHT_CODE);
    double decay_rate_M3 = 3.0 * K_BH_const_code * All.PBH_Alpha;

    for(i = 0; i < PBH_TABLE_SIZE; i++)
    {
        All.PBH_Table_ScaleFactor[i] = current_a;
        All.PBH_Table_Mass[i] = current_mass;

        if (i < PBH_TABLE_SIZE - 1)
        {
            double dt;
            double next_a = a_start + (i + 1) * da;

            if(All.ComovingIntegrationOn)
            {
                // Integrate dt = da / (a * H(a))
                // Use midpoint for better accuracy if needed, or just simple step
                double a_mid = 0.5 * (current_a + next_a);
                double obs_hubble = hubble_function(a_mid);
                dt = (next_a - current_a) / (a_mid * obs_hubble);
            }
            else
            {
                // Non-cosmological: All.Time is time. All.TimeBegin is start time.
                // So da is dt.
                dt = next_a - current_a;
            }

            // Evolve mass
            current_mass3 -= decay_rate_M3 * dt;
            if(current_mass3 < 0) current_mass3 = 0;
            current_mass = pow(current_mass3, 1.0/3.0);

            current_a = next_a;
        }
    }

    if(ThisTask == 0)
    {
        printf("PBH Mass Evolution Table Initialized.\n");
        printf("Initial Mass: %g, Final Mass (at a=%g): %g\n\n", All.PBH_Table_Mass[0], All.PBH_Table_ScaleFactor[PBH_TABLE_SIZE-1], All.PBH_Table_Mass[PBH_TABLE_SIZE-1]);
    }
#else
    if(ThisTask == 0) printf("PBH Mass Evolution not explicitly enabled. The PBH mass is fixed at the initial value.\n");
#endif
}

/*
 * Get the current PBH mass by interpolation from the lookup table.
 */
void get_current_pbh_mass(double a, double *mass_out)
{
#ifndef PBH_EVAPORATION_FEEDBACK_NO_MASS_LOSS
    // Handle out of bounds
    if (a <= All.PBH_Table_ScaleFactor[0])
    {
        *mass_out = All.PBH_Table_Mass[0];
        return;
    }
    if (a >= All.PBH_Table_ScaleFactor[PBH_TABLE_SIZE-1])
    {
        *mass_out = All.PBH_Table_Mass[PBH_TABLE_SIZE-1];
        return;
    }

    // Binary search or direct index (since uniform spacing in 'a')
    // Optimisation: uniform spacing allows O(1) lookup
    double da = (All.TimeMax - All.TimeBegin) / (double)(PBH_TABLE_SIZE - 1);
    int idx = (int)((a - All.PBH_Table_ScaleFactor[0]) / da);

    if (idx < 0) idx = 0;
    if (idx >= PBH_TABLE_SIZE - 1) idx = PBH_TABLE_SIZE - 2;

    // Check if idx is correct (safe check due to floating point)
    while (idx < PBH_TABLE_SIZE - 1 && a > All.PBH_Table_ScaleFactor[idx+1]) idx++;
    while (idx > 0 && a < All.PBH_Table_ScaleFactor[idx]) idx--;

    // Linear interpolation
    double a0 = All.PBH_Table_ScaleFactor[idx];
    double a1 = All.PBH_Table_ScaleFactor[idx+1];
    double m0 = All.PBH_Table_Mass[idx];
    double m1 = All.PBH_Table_Mass[idx+1];

    double f = (a - a0) / (a1 - a0);
    *mass_out = m0 + f * (m1 - m0);
#else
    *mass_out = All.PBH_InitialMass;
#endif
}

#endif /* #if defined(PBH_EVAPORATION_FEEDBACK) || defined(PBH_EVAPORATION_FEEDBACK_DM) */
