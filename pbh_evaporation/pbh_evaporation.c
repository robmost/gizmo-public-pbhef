#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "pbh_evaporation.h"

#if defined(PBH_EVAPORATION_FEEDBACK) || defined(PBH_EVAPORATION_FEEDBACK_DM)

/*
 * Calculates the alpha coefficient for PBH evaporation rate based on initial PBH mass. Following Mosbech et al. (2022).
 */
double pbh_evaporation_calculate_alpha_coefficient(double m_pbh_initial_grams)
{
    if (m_pbh_initial_grams < 1.0e18)
    {
        const double c1 = -0.3015;
        const double c2 = 0.3113;
        const double p_exponent = -0.0008;
        return c1 + c2 * pow(m_pbh_initial_grams, p_exponent);
    }
    else
    {
        return 2.011e-4;
    }
}
#endif // PBH_EVAPORATION_FEEDBACK || PBH_EVAPORATION_FEEDBACK_DM


#ifdef PBH_EVAPORATION_FEEDBACK

/*
 * Retrieves the local DM density at a gas particle. This function checks if the DM density has been computed for the current step.
 * If it has, it returns the computed value; otherwise, it returns 0.0.
 */
double pbh_evaporation_get_local_dm_density_at_gas_particle(int target_gas_particle_index)
{
    if(SphP[target_gas_particle_index].Flag_DM_Density_ForPBH_Computed_This_Step == All.NumCurrentTiStep)
    {
        return SphP[target_gas_particle_index].DM_Density_ForPBH;
    }
    else
    {
        return 0.0;
    }
}


/*
 * Applies the PBH evaporation feedback to gas particles. This function calculates the specific heating rate based on the local DM density
 * and updates the internal energy of the gas particles accordingly.
 */
void apply_receiver_pbh_evaporation_feedback(void)
{
    int i;
    double pbh_alpha;
    double m_pbh_initial_mass_grams_code_units_cubed;
    double rho_chi_local, rho_gas_local;
    double specific_heating_rate;

    pbh_alpha = pbh_evaporation_calculate_alpha_coefficient(All.PBH_InitialMass_grams);

    // If alpha is not strictly positive, then there's no heating from this mechanism.
    if (pbh_alpha <= 0.0)
    {
        if (ThisTask == 0 && All.Time == All.TimeBegin)
        {
            printf("PBH_EVAPORATION_FEEDBACK Info: Alpha coefficient calculated as %g (<=0) using PBH_InitialMass_grams = %g grams.\n", pbh_alpha, All.PBH_InitialMass_grams);
            printf("PBH evaporation heating will be zero for this configuration of PBH mass.\n");
        }
        return;
    }

    // Pre-calculate (M_PBH,0 [code units])^3 using All.PBH_InitialMass_grams
    double m_pbh_initial_mass_code_units = All.PBH_InitialMass_grams / UNIT_MASS_IN_CGS;
    m_pbh_initial_mass_grams_code_units_cubed = pow(m_pbh_initial_mass_code_units, 3.0);
    if (m_pbh_initial_mass_grams_code_units_cubed == 0.0) {
        // This case should ideally not be reached if All.PBH_InitialMass_grams is positive.
        // However, including it for robustness against potential floating point issues with very small masses,
        // or if pow() somehow returned zero.
        if (ThisTask == 0 && All.Time == All.TimeBegin) {
             printf("PBH_EVAPORATION_FEEDBACK Warning: m_pbh_initial_mass_grams_code_units_cubed is zero. No heating.\n");
        }
        return;
    }


    for(i = 0; i < NumGas; i++)
    {
        if(P[i].Type == 0 && P[i].Mass > 0) // Particle is gas and has mass
        {
            rho_gas_local = SphP[i].Density;
            if (rho_gas_local <= 0) continue;

            rho_chi_local = pbh_evaporation_get_local_dm_density_at_gas_particle(i);
            if (rho_chi_local <= 0) continue;

            // Calculate the specific heating rate (DtInternalEnergy)
            // Equation: (f * rho_chi / rho_gas) * (hbar*c^6/G^2) * (alpha / m_PBH,0^3)
            // Use All.PBH_MassFraction_f and All.PBH_EvaporationConstant
            specific_heating_rate = (All.PBH_MassFraction_f * rho_chi_local / rho_gas_local) * All.PBH_EvaporationConstant * (pbh_alpha / m_pbh_initial_mass_grams_code_units_cubed);


            if(specific_heating_rate < 0) specific_heating_rate = 0.0;

            SphP[i].DtInternalEnergy += specific_heating_rate;

            // Optional diagnostics
            // #ifdef PBH_OUTPUT_HEATING_RATE
            // SphP[i].PBH_HeatingRate_Gas = specific_heating_rate;
            // #endif
        }
    }
}

/*
 * Computes the DM density for gas particles. This function is called at the beginning of each time step.
 * It initializes the DM density for gas particles and sets the flag indicating that the DM density has been computed.
 */
void pbh_evaporation_calculate_dm_density_for_gas_particles(void)
{
    if(ThisTask == 0 && All.Time == All.TimeBegin)
    {
        printf("PBH_EVAPORATION_FEEDBACK: `pbh_evaporation_calculate_dm_density_for_gas_particles` needs full implementation!\n");
    }
    for(int i = 0; i < NumGas; i++)
    {
        if(P[i].Type == 0)
        {
            SphP[i].DM_Density_ForPBH = 0.0;
        }
    }
    // Actual computation logic for DM density to be added here.
    // This will fill SphP[i].DM_Density_ForPBH and set SphP[i].Flag_DM_Density_ForPBH_Computed_This_Step = All.NumCurrentTiStep;
}

#endif // PBH_EVAPORATION_FEEDBACK