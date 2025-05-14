#ifndef PBH_EVAPORATION_H
#define PBH_EVAPORATION_H

// This header file declares the functions for the Primordial Black Hole Evaporation module.

// Conditional compilation based on the GIZMO flags defined in Config.sh

// --- Common functions for any PBH evaporation feedback model ---
#if defined(PBH_EVAPORATION_FEEDBACK) || defined(PBH_EVAPORATION_FEEDBACK_DM)

/**
 * @brief Calculates the alpha coefficient for PBH evaporation rate.
 *
 * This coefficient depends on the initial mass of the Primordial Black Holes.
 * The formula is based on the user's provided specification.
 *
 * @param m_pbh_initial_grams The initial mass of a single PBH in grams.
 * @return The dimensionless alpha coefficient.
 */
double pbh_evaporation_calculate_alpha_coefficient(double m_pbh_initial_grams);

#endif // PBH_EVAPORATION_FEEDBACK || PBH_EVAPORATION_FEEDBACK_DM


// --- Functions specific to the Receiver-Based PBH Evaporation Feedback ---
// (Feedback energy calculated at the gas particle's position)
#ifdef PBH_EVAPORATION_FEEDBACK

/**
 * @brief Calculates the local Dark Matter density at each active gas particle's position.
 *
 * This function orchestrates the neighbor search and kernel summation required to
 * estimate rho_chi. The results are typically stored in a temporary field
 * (e.g., SphP[i].DM_Density_ForPBH) for use by the feedback application routine.
 * This function should be called once per relevant step before applying the feedback.
 */
void pbh_evaporation_calculate_dm_density_for_gas_particles(void);

/**
 * @brief Applies the PBH evaporation feedback energy to gas particles.
 *
 * This function iterates over active gas particles, calculates the specific heating
 * rate based on local gas density, pre-calculated local DM density, and PBH
 * parameters, and adds this to SphP[i].DtInternalEnergy.
 * It should be called after pbh_evaporation_calculate_dm_density_for_gas_particles().
 */
void apply_receiver_pbh_evaporation_feedback(void);

#endif // PBH_EVAPORATION_FEEDBACK


// --- Functions specific to the Donor-Based PBH Evaporation Feedback (Future Implementation) ---
#ifdef PBH_EVAPORATION_FEEDBACK_DM

// void apply_donor_pbh_evaporation_feedback(void); // Example for future use

#endif // PBH_EVAPORATION_FEEDBACK_DM


#endif // PBH_EVAPORATION_H