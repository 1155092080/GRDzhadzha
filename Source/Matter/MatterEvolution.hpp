/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTEREVOLUTION_HPP_
#define MATTEREVOLUTION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS of matter variables only, metric vars assumed analytic
/*!
     The class calculates the RHS evolution for the matter variables.
     It does not assume a specific form of matter or background -
     it is templated over a matter class matter_t, and over a background metric,
     background_t
*/

template <class matter_t, class background_t> class MatterEvolution
{
  public:
    //! Inherit the variable definitions from the Matter vars
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    //  Need d2 of certain matter vars
    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    // This is used for the non evolved ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    //!  Constructor of class MatterEvolution
    MatterEvolution(matter_t a_matter, background_t a_background, double sigma,
                    double dx, std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_background(a_background), m_sigma(sigma),
          m_deriv(dx), m_dx(dx), m_center(a_center)
    {
    }

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy matter data from chombo gridpoint into local variable
        const auto matter_vars = current_cell.template load_vars<MatterVars>();

        // compute the background metric vars
        MetricVars<data_t> metric_vars;
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        m_background.compute_metric_background(metric_vars, coords);

        // compute derivs for matter grid vars
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterDiff2Vars>(current_cell);
        const auto advec = m_deriv.template advection<MatterVars>(
            current_cell, metric_vars.shift);

        // the RHS
        MatterVars<data_t> matter_rhs;
        VarsTools::assign(matter_rhs, 0.); // All components set to zero

        // add evolution of matter fields and dissipation
        m_matter.matter_rhs(matter_rhs, matter_vars, metric_vars, d1, d2,
                            advec, coords);
        m_deriv.add_dissipation(matter_rhs, current_cell, m_sigma);

        // Write the rhs into the output vars for this cell
        current_cell.store_vars(matter_rhs);
    }

  protected:
    const matter_t m_matter;              //!< The matter object
    const background_t m_background;      //!< The metric background
    const FourthOrderDerivatives m_deriv; //!< An object for calculating
                                          //!< derivatives of the vars
    const double m_sigma;                 //!< Sigma for dissipation
    const double m_dx;                    //!< Grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< Grid center
};

#endif /* MATTEREVOLUTION_HPP_ */
