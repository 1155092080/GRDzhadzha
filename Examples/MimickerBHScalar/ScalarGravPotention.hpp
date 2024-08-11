/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARGRAVPOTENTIAL_HPP_
#define SCALARGRAVPOTENTIAL_HPP_

#include "simd.hpp"
#include "Coordinates.hpp"

class ScalarGravPotential
{
  protected:
    //    const double m_mu;
    const InitialScalarData::params_t m_initial_params;

  public:
    //! The constructor
    ScalarPotential(const InitialScalarData::params_t a_initial_params)
        : m_initial_params(a_initial_params)
    {
    }

    template <class data_t>
    ALWAYS_INLINE data_t Grav(const vars_t<data_t> &vars, const Coordinates<data_t> &coords) const
    {
        return 0.5 * pow(X * X + Y * Y + Z * Z + 0.0000001, -0.5);
    } // V

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // 1/2 m^2 phi^2
        const double mu = m_initial_params.mass;
        const double lambda = 0.02;
        V_of_phi = 0.5 * mu * mu * vars.phi * vars.phi + 1.0/2.0/3.0/4.0 * lambda * vars.phi * vars.phi * vars.phi * vars.phi;

        // The potential gradient at phi wrt the field
        // m^2 phi
        dVdphi = mu * mu * vars.phi + 1.0/2.0/3.0 * lambda * vars.phi * vars.phi * vars.phi;
    }
};

#endif /* SCALARPOTENTIAL_HPP_ */
