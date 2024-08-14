/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MIMICKERBH_HPP_
#define MIMICKERBH_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Class which computes the initial conditions for a boosted BH Mimicker
// Undoosted elem. ds^2 = -Adt^2 + B dx^2
// Lorentz boost x = \gamma (x' - vt'), t = \gamma (t' - v*x')
// Addisional shift to fix coordinates to the BH
// x' = \tilde x + v \tilde t, t' = \tilde t
// BH appears fixed, but we're in the rest frame of the unboosted observer!
// See https://arxiv.org/pdf/2106.08280.pdf for details

class MimickerBH
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double rs = 2;  // Radius to smooth the Newtonian potential to avoid singularity
        double mu = 0.25;   // Mass of scalar field particles
        double pi = M_PI;
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double velocity = 0.0; //!< The boost velocity in the x direction
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    MimickerBH(params_t a_params, double a_dx) : m_params(a_params), m_dx(a_dx)
    {
        // check this boost param is sensible
        if ((m_params.velocity > 1.0) || (m_params.velocity < -1.0))
        {
            MayDay::Error("The boost velocity parameter must be in the range "
                          "-1.0 < velocity < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, coords);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);

        current_cell.store_vars(chi, c_chi);
    }

    /// Mimicker boosted solution as above
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {
        // Mimicker params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const double M = m_params.mass;
        const double v = m_params.velocity;
        const double pi = m_params.pi;
        const data_t rs = m_params.rs;
        const data_t mu = m_params.mu;
        // std::cout << "rs = " << rs << ", bh_mass = " << M << std::endl; 
        const double dens = M*3.0/(4.0*pi*pow(m_params.rs,3));
        const double v2 = v * v;
        const double boost2 = 1.0 / (1 - v2);
        const double boost = sqrt(boost2);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        const data_t x_p = coords.x * boost;
        const double y = coords.y;
        const double z = coords.z;

        // the radius (boosted)
        const data_t r2 = x_p * x_p + y * y + z * z;
        const data_t r = sqrt(r2);

        // find the H, A, B metric quantities, using weak field approximation, ds^2=-Adt^2+Bdx^2 = -(1+2H)dt^2+(1-2H)dx^2, where H is Newtonian potential
        const data_t H = simd_conditional(simd_compare_gt(r, rs), mu*2.0/3.0*pi*pow(rs,2)*dens + mu*M*(1/rs -1/r), mu*2.0/3.0*pi*pow(r,2)*dens); //(r>rs) ? mu*2.0/3.0*pi*pow(rs,2)*dens + mu*M*(1/rs -1/r) : mu*2.0/3.0*pi*pow(r,2)*dens
        const data_t A = (1+2*H);
        const data_t B = (1-2*H);

        // Calculate the gradients of A and B
        Tensor<1, data_t> dAdx, dBdx;
        get_metric_derivs(dAdx, dBdx, coords);

        // populate ADM vars
        // (Keep lapse negative within the horizon)
        vars.lapse = sqrt(A) / boost * sqrt(B / (B - A * v2));

        using namespace TensorAlgebra;
        FOR2(i, j) { vars.gamma[i][j] = delta(i, j) * B; }
        vars.gamma[0][0] = boost2 * (vars.gamma[0][0] - A * v2);
        const auto gamma_UU = compute_inverse_sym(vars.gamma);

        // this adjustment gives the shift which achieves x' = x - vt
        FOR1(i) { vars.shift[i] = delta(i, 0) * A * v / boost2 / (B - A * v2); }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k) { vars.d1_gamma[i][j][k] = delta(i, j) * dBdx[k]; }
        FOR1(i)
        {
            vars.d1_gamma[0][0][i] =
                (vars.d1_gamma[0][0][i] - dAdx[i] * v2) * boost2;
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = 0.5 * vars.lapse *
                               (dAdx[i] / A + dBdx[i] / B -
                                (dBdx[i] - v2 * dAdx[i]) / (B - v2 * A));
        }

        // use the fact that shift^i = lapse^2 * shift_i + v^i
        // and v^i is a constant vector
        FOR2(i, j)
        {
            vars.d1_shift[i][j] =
                delta(i, 0) * vars.shift[0] *
                (dAdx[j] / A - (dBdx[j] - dAdx[j] * v2) / (B - A * v2));
        }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_phys.ULL[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
    }

  protected:
    /// Work out the gradients of the quantities H, A and B appearing in the
    /// Isotropic Schwarzschild solution
    template <class data_t>
    void get_metric_derivs(Tensor<1, data_t> &dAdx, Tensor<1, data_t> &dBdx,
                           const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and boost v
        const double M = m_params.mass;
        const double v = m_params.velocity;
        const data_t mu = m_params.mu;
        const double pi = m_params.pi;
        const data_t rs = m_params.rs;
        const double dens = M*3.0/(4.0*pi*pow(m_params.rs,3));
        const double boost = pow(1 - v * v, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;
        const data_t x_p = coords.x * boost;

        // the coordinate radius (boosted), subject to a floor
        const data_t r2 = x_p * x_p + x[1] * x[1] + x[2] * x[2];
        const data_t r = sqrt(r2);

        // find H
        const data_t H = simd_conditional(simd_compare_gt(r, rs), mu*2.0/3.0*pi*pow(rs,2)*dens + mu*M*(1/rs -1/r), mu*2.0/3.0*pi*pow(r,2)*dens); //(r>rs) ? mu*2.0/3.0*pi*pow(rs,2)*dens + mu*M*(1/rs -1/r) : mu*2.0/3.0*pi*pow(r,2)*dens
        const data_t A = (1+2*H);
        const data_t B = (1-2*H);

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drdx;
        FOR1(i) { drdx[i] = x[i] / r; }
        drdx[0] *= boost * boost;

        // derivs of quantities
        Tensor<1, data_t> dHdx;
        FOR1(i)
        {
            dHdx[i] = 
            dHdx[i] = simd_conditional(simd_compare_gt(r, rs),  M * mu / r / r * drdx[i], mu*4.0/3.0*pi*r*dens*drdx[i]); //(r>rs) ? M * mu / r / r * drdx[i]: mu*4.0/3.0*pi*r*dens*drdx[i]
            dAdx[i] = 2 * dHdx[i];  
            dBdx[i] = -2 * dHdx[i];
        }
    }

  public:
    // used to decide when to excise - Since we are considering newtonian planet, there's no horizon, thus no need to excise
    // note that this is not templated over data_t
    bool check_if_excised(const Coordinates<double> &coords) const
    {
        bool is_excised = false;

        return is_excised;
    }
};

#endif /* MIMICKERBHFIXEDBG_HPP_ */
