/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "MimickerBHScalarLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "AngularMomConservation.hpp"
#include "EnergyConservation.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FluxExtraction.hpp"
#include "ScalarExtraction.hpp"
#include "InitialScalarData.hpp"
#include "MimickerBH.hpp"
#include "MatterEvolution.hpp"
#include "ScalarField.hpp"
#include "SpatialDepScalarPotential.hpp"

// Initial data for field and metric variables
void MimickerBHScalarLevel::initialData()
{
    CH_TIME("MimickerBHScalarLevel::initialData");
    if (m_verbosity)
        pout() << "MimickerBHScalarLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for fields
    SetValue set_zero(0.0);
    MimickerBH mimicker_bh(m_p.bg_params, m_dx); // just calculates chi
    InitialScalarData initial_sf(m_p.initial_params, m_dx);
    auto compute_pack = make_compute_pack(set_zero, mimicker_bh);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(ExcisionEvolution<ScalarFieldWithPotential, MimickerBH>(
                       m_dx, m_p.center, mimicker_bh),
                   m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

// Things to do in RHS update, at each RK4 step
void MimickerBHScalarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                        const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    SpatialDepScalarPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    MimickerBH mimicker_bh(m_p.bg_params, m_dx);
    MatterEvolution<ScalarFieldWithPotential, MimickerBH> my_matter(
        scalar_field, mimicker_bh, m_p.sigma, m_dx, m_p.center);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(ExcisionEvolution<ScalarFieldWithPotential, MimickerBH>(
                       m_dx, m_p.center, mimicker_bh),
                   a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

void MimickerBHScalarLevel::specificPostTimeStep()
{
    // Check for nans on every level
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());

    // At any level, but after the min_level timestep
    int min_level = m_p.extraction_params.min_extraction_level();
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    if (calculate_diagnostics)
    {
        fillAllGhosts();
        SpatialDepScalarPotential potential(m_p.initial_params);
        ScalarFieldWithPotential scalar_field(potential);
        MimickerBH mimicker_bh(m_p.bg_params, m_dx);
        AngularMomConservation<ScalarFieldWithPotential, MimickerBH>
            angular_mom(scalar_field, mimicker_bh, m_dx, m_p.center);
        EnergyConservation<ScalarFieldWithPotential, MimickerBH> energies(
            scalar_field, mimicker_bh, m_dx, m_p.center);
        BoxLoops::loop(make_compute_pack(angular_mom, energies), m_state_new,
                       m_state_diagnostics, SKIP_GHOST_CELLS);

        // excise within/outside specified radii, no simd
        BoxLoops::loop(
            ExcisionDiagnostics<ScalarFieldWithPotential, MimickerBH>(
                m_dx, m_p.center, mimicker_bh, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
    }

    // write out the integral after each timestep on minimum level
    if (m_p.activate_extraction == 1 || m_p.activate_scalar_extraction == 1)
    {
        if (m_level == min_level)
        {
            bool first_step = (m_time == m_dt);
            // integrate the densities and write to a file
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double rhoEnergy_sum = amr_reductions.sum(c_rhoEnergy);
            double rhoAngMom_sum = amr_reductions.sum(c_rhoAngMom);

            SmallDataIO integral_file(m_p.data_path + "EnergyIntegrals", m_dt,
                                      m_time, m_restart_time,
                                      SmallDataIO::APPEND, first_step);
            // remove any duplicate data if this is post restart
            integral_file.remove_duplicate_time_data();

            std::vector<double> data_for_writing = {rhoEnergy_sum,
                                                    rhoAngMom_sum};

            // write data
            if (first_step)
            {
                integral_file.write_header_line(
                    {"Energy density", "Ang. Mom. density"});
            }
            integral_file.write_time_data_line(data_for_writing);

            // Now refresh the interpolator and do the interpolation
            bool fill_ghosts = false;
            m_gr_amr.m_interpolator->refresh(fill_ghosts);
            m_gr_amr.fill_multilevel_ghosts(
                VariableType::diagnostic, Interval(c_fluxEnergy, c_fluxAngMom));
            if (m_p.activate_extraction)
                {
                    FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                         m_restart_time);
                    my_extraction.execute_query(m_gr_amr.m_interpolator, m_p.data_path);
                }

                if (m_p.activate_scalar_extraction)
                {
                    ScalarExtraction phi_extraction(
                        m_p.scalar_extraction_params, m_dt, m_time, first_step,
                        m_restart_time);
                    phi_extraction.execute_query(m_gr_amr.m_interpolator);
                }
            
        }
    }
}

void MimickerBHScalarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
}
