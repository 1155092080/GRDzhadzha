/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "ScalarField.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <random>
#include <cmath>

//! Class which creates a constant scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    struct params_t
    {
        double mass;
        double amplitude;
        std::array<double, CH_SPACEDIM> center;
    };
    constexpr static int n_modes = 1000;
    double v_amp = 0.1;
    double fw_trunc_radius = 200.0;
    double fw_trunc_k = 0.2;
    double v_theta[n_modes];
    double v_phi[n_modes];
    double v_phase[n_modes];
    double omg[n_modes];

    //! The constructor for the class
    InitialScalarData(const params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
      // random number generator
      double pi = M_PI;
      std::random_device dev;
      std::mt19937 rng;
      std::seed_seq sseq{1, 2, 3};
      rng.seed(sseq);
      std::uniform_real_distribution<double> distributiont(0.0,pi);   //random number from 0 to pi
      std::uniform_real_distribution<double> distributionpi(0.0,2*pi);//random number from 0 to 2pi
      // Generate velocity and random phase
      double mu = m_params.mass; // scalar field mass
      for (int i = 0; i < n_modes; i++){
        v_theta[i] = distributiont(rng); // velocity diirection theta
        v_phi[i] = distributionpi(rng); // velocity diirection phi
        v_phase[i] = distributionpi(rng); // velocity phase
        omg[i] = std::sqrt(mu*mu + mu*mu*v_amp*v_amp);
      }
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        double pi = M_PI;
        

        ScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t radius = coords.get_radius();
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        double mu = m_params.mass; // scalar field mass
        data_t phi_of_r = 0.0;
        data_t Pi_of_r = 0.0;
        // set the field vars
        data_t r = sqrt(x*x + y*y + z*z);
        data_t trunc_coeff = 1.0/(1.0 + exp(fw_trunc_k*(r - fw_trunc_radius)));
        for (int i = 0; i < n_modes; i++){
          phi_of_r += m_params.amplitude/std::sqrt(n_modes) * cos(-mu*v_amp*(sin(v_theta[i])*(cos(v_phi[i])*x + sin(v_phi[i])*y) + cos(v_theta[i])*z) + v_phase[i]);
          Pi_of_r += - m_params.amplitude/std::sqrt(n_modes) * omg[i] * sin(-mu*v_amp*(sin(v_theta[i])*(cos(v_phi[i])*x + sin(v_phi[i])*y) + cos(v_theta[i])*z) + v_phase[i]);
        }
        vars.phi = trunc_coeff*phi_of_r;
        vars.Pi = trunc_coeff*Pi_of_r;
        current_cell.store_vars(vars);
    }

  protected:
    const double m_dx;
    const params_t m_params;
};

#endif /* INITIALSCALARDATA_HPP_ */
