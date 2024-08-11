/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPATIALDEPSCALARPOTENTIAL_HPP_
#define SPATIALDEPSCALARPOTENTIAL_HPP_

#include "simd.hpp"

class SpatialDepScalarPotential
{
  protected:
    //    const double m_mu;
    const InitialScalarData::params_t m_initial_params;

  public:
    //! The constructor
    SpatialDepScalarPotential(const InitialScalarData::params_t a_initial_params)
        : m_initial_params(a_initial_params)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars,
                           const Coordinates<data_t> &coords) const
    {
        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        
        const data_t x_p = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        // the isotropic radius (boosted)
        const data_t r2 = x_p * x_p + y * y + z * z;
        const data_t r = sqrt(r2);
        // The potential value at phi
        // 1/2 m^2 phi^2
        const double mu = m_initial_params.mass;
        const double lambda = 0.02;
        // V_of_phi = 0.5 * mu * mu * pow(atan(1.0/200.0*r),2) * vars.phi * vars.phi + 1.0/2.0/3.0/4.0 * lambda * vars.phi * vars.phi * vars.phi * vars.phi;
        V_of_phi = 0.5 * mu * mu * vars.phi * vars.phi + 1.0/2.0/3.0/4.0 * lambda * vars.phi * vars.phi * vars.phi * vars.phi;

        // The potential gradient at phi wrt the field
        // m^2 phi
        // dVdphi = mu * mu * pow(atan(1.0/200.0*r),2) * vars.phi + 1.0/2.0/3.0 * lambda * vars.phi * vars.phi * vars.phi;
        dVdphi = mu * mu * vars.phi + 1.0/2.0/3.0 * lambda * vars.phi * vars.phi * vars.phi;
    }
};

#endif /* SCALARPOTENTIAL_HPP_ */
