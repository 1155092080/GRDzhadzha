/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXEXTRACTION_HPP_
#define FLUXEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the force components on
//!  spheroidal shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the force
   components over spheroidal shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class FluxExtraction : public SphericalExtraction
{
  private:
    static constexpr int m_scalar = 0;
  public:
    //! The constructor
    FluxExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_fluxEnergy, VariableType::diagnostic);
        add_var(c_fluxAngMom, VariableType::diagnostic);
        add_var(c_phi, VariableType::evolution);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    FluxExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                   double a_time, double a_restart_time = 0.0)
        : FluxExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    // the references of the vars as used in the integrator
    enum M_VARS
    {
        m_fluxEnergy,
        m_fluxAngMom,
        NUM_EXTRACTION_COMPS
    };

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator,
                       std::string a_datapath)
    {
        // extract the values of the Flux scalars on the spheres
        extract(a_interpolator);

        // this would write out the values at every point on the sphere
        if (m_params.write_extraction)
        {
            write_extraction(a_datapath + "FluxExtractionOut_");
        }
        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<double>, std::vector<double>>>
            mode_integrals(m_num_modes);
            
        // note that this is normalised by multiplying by radius
        auto normalised_rhoEnergy =
            [](std::vector<double> rhoEnergy_parts, double r, double, double)
        {
            // here the std::vector<double> passed will have the
            // values of scalar_Re = scalar and scalar_Im = 0.0 as its first two
            // components
            return std::make_pair(r * rhoEnergy_parts[m_scalar],
                                  0);
        };
				 // add the modes that will be integrated
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            constexpr int es = 0;
            // mode.first is l, mode.second is m
            add_mode_integrand(es, mode.first, mode.second,
                               normalised_rhoEnergy, mode_integrals[imode]);
        }
				
				 
        // Setup to integrate fluxes
        std::vector<std::vector<double>> force_integrals(NUM_EXTRACTION_COMPS);
        add_var_integrand(m_fluxEnergy, force_integrals[m_fluxEnergy],
                          IntegrationMethod::trapezium);
        add_var_integrand(m_fluxAngMom, force_integrals[m_fluxAngMom],
                          IntegrationMethod::trapezium);

        // do the integration over the surface
        // integrate() do integration on all integrands added by add_integrand
        // Since add_mode_integrand is implemented using add_integrand, modes can
        // also be integrated by calling integrate()
        integrate();
        // write the integrals
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            std::string integrals_filename = a_datapath + m_params.integral_file_prefix +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_for_writing = {
                std::move(mode_integrals[imode].first),
                std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            write_integrals(integrals_filename, integrals_for_writing, labels);
        }

        // write the integrals
        std::vector<std::string> labels(NUM_EXTRACTION_COMPS);
        labels[m_fluxEnergy] = "Energy Flux";
        labels[m_fluxAngMom] = "Ang. Mom. Flux";
        std::string filename = a_datapath + "FluxIntegrals";
        write_integrals(filename, force_integrals, labels);  
    }
};

#endif /* FLUXEXTRACTION_HPP_ */
