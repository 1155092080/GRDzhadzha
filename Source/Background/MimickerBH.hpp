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

//! Class which computes the initial conditions for a BH Mimicker
// From https://arxiv.org/pdf/2106.08280
// and metric from https://arxiv.org/pdf/2110.06937

class MimickerBH
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double L0 = 5.0;  // Radius of compact star
        double pi = M_PI;
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    MimickerBH(params_t a_params, double a_dx) : m_params(a_params), m_dx(a_dx)
    {
        // check the L0 so that avoid horizon
        if ((m_params.L0 < 9.0*m_params.mass/4.0))
        {
            MayDay::Error("L0 should be larger than 9m/4 to avoid horizon");
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
        data_t chi = TensorAlgebra::compute_determinant(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
    }

    /// Mimicker solution as above
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {
        // Mimicker params - mass M and radius L0
        const double M = m_params.mass;
        const double L0 = m_params.L0;
        const double pi = m_params.pi;
        // const double dens = M*3.0/(4.0*pi*pow(m_params.L0,3));

        // work out coordinates
        const data_t x = coords.x;
        const data_t y = coords.y;
        const data_t z = coords.z;

        // the radius (boosted)
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const data_t L3 = pow(L0, 3);

        // Inner and outer metric components

        const data_t alpha_inner = 3.0/2.0 * sqrt(1 - 2*M/L0) - 1.0/2.0 * sqrt(1 - 2*M*r2/L3);
        const data_t alpha_outer = sqrt(1 - 2*M/r);
        const data_t H = simd_conditional(simd_compare_gt(r, L0), M/r2/(r - 2*M), M/(L3 - 2*M*r2)); // g_ab = eta_ab + 2H el_a el_b

        const Tensor<1, data_t> el = {x, y, z};

        // Calculate the gradients in el and H
        Tensor<1, data_t> dHdx;
        Tensor<1, data_t> drdx;
        Tensor<2, data_t> dldx;
        get_Mimiker_derivs(dHdx, dldx, drdx, H, coords);

        // populate ADM vars
        // (No horizon, so lapse always positive)
        vars.lapse = simd_conditional(simd_compare_gt(r, L0), alpha_outer, alpha_inner); // r>L0? alpha_outer:alpha_inner

        FOR2(i, j)
        {
            vars.gamma[i][j] =
                TensorAlgebra::delta(i, j) + 2.0 * H * el[i] * el[j];
        }

        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(vars.gamma);

        // For static metric, shift = 0
        FOR1(i) { vars.shift[i] = 0.0; }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                2.0 * (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                       H * el[j] * dldx[i][k]);
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = simd_conditional(simd_compare_gt(r, L0),
                                                M /r2 / alpha_outer * drdx[i], 
                                                M / L3 / sqrt(1 - 2*M*r2/L3)* r * drdx[i]
            );
        }

        // use the fact that shift^i = lapse^2 * shift_i = 0
        FOR2(i, j)
        {
            vars.d1_shift[i][j] = 0.0;
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
    /// Work out the gradients of the quantities H and el and r
    template <class data_t>
    void get_Mimiker_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                       Tensor<1, data_t> &drdx, const data_t &H,
                       const Coordinates<data_t> &coords) const
    {
        // black hole mimicker params - mass M and radium L0
        const double M = m_params.mass;
        const double L0 = m_params.L0;
        const double pi = m_params.pi;

        // work out where we are on the grid and useful quantities
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const data_t L3 = L0 * L0 * L0;

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        FOR1(i) { drdx[i] = x[i] / r; }

        // derivs of quantities H and el
        FOR1(i)
        {
            dHdx[i] = -pow(H, 2)/M;
            dHdx[i] *= simd_conditional(simd_compare_gt(r, L0),  3*pow(r, 2)-4*M*r, -4*M*r);
            dHdx[i] *= drdx[i];
        }

        FOR2(i, j)
        {
            dldx[i][j] = delta(i, j); // Because l = {x, y, z}
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
