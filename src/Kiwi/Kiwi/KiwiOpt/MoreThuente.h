/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_MORETHUENTE_H
#define KIWI_MORETHUENTE_H

#include <cstddef>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Optimization {

/*!
 *
 * Moré & Thuente Line Search Algortihm
 *
 * Original comment block:
 *
 * The purpose of cvsrch is to find a step which satisfies
 * a sufficient decrease condition and a curvature condition.
 * The user must provide a subroutine which calculates the
 * function and the gradient.
 *
 * At each stage the subroutine updates an interval of
 * uncertainty with endpoints stx and sty. The interval of
 * uncertainty is initially chosen so that it contains a
 * minimizer of the modified function
 *
 *      f(x + stp * s) - f(x) - f_tol * stp * (gradf(x)' * s).
 *
 * If a step is obtained for which the modified function
 * has a nonpositive function value and nonnegative derivative,
 * then the interval of uncertainty is chosen so that it
 * contains a minimizer of f(x + stp * s).
 *
 * The algorithm is designed to find a step which satisfies
 * the sufficient decrease condition
 *
 *       f(x + stp * s) <= f(x) + f_tol * stp * (gradf(x)' * s),
 *
 * and the curvature condition
 *
 *       abs(gradf(x + stp * s)' * s)) <= gtol * abs(gradf(x)' * s).
 *
 * If f_tol is less than gtol and if, for example, the function
 * is bounded below, then there is always a step which satisfies
 * both conditions. If no step can be found which satisfies both
 * conditions, then the algorithm usually stops when rounding
 * errors prevent further progress. In this case stp only
 * satisfies the sufficient decrease condition.
 *
 *
 * f_tol and gtol are nonnegative input variables. Termination
 *   occurs when the sufficient decrease condition and the
 *   directional derivative condition are satisfied.
 *
 * x_tol is a nonnegative input variable. Termination occurs
 *   when the relative width of the interval of uncertainty
 *   is at most x_tol.
 *
 * Sources:
 *
 * Moré, J. J., & Thuente, D. J. (1994).
 * Line search algorithms with guaranteed sufficient decrease.
 * ACM Transactions on Mathematical Software (TOMS), 20(3), 286-307.
 *
 * Argonne National Laboratory. MINPACK Project. June 1983
 *     Jorge J. More', David J. Thuente
 * https://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/dctep.f
 *
 *
 * https://github.com/JuliaNLSolvers/NLSolvers.jl
 * https://github.com/JuliaNLSolvers/LineSearches.jl/blob/master/src/morethuente.jl
 *
 *
 * @tparam ObjectiveFunction
 */
template<class ObjectiveFunction>
class MoreThuente {
 private:
  struct Function {
    double alpha;
    double value;
    double derivative;
  };

  ObjectiveFunction& objectiveFunction;

  int max_iterations = 10;
  int iteration;

  double f_tol = 1e-4;
  double gtol = 0.9;
  double x_tol = 1e-8;

 public:
  void setTol(double xTol) {
    x_tol = xTol;
  }

  void setAlphaMin(double alphaMin) {
    alphamin = alphaMin;
  }

  void setAlphaMax(double alphaMax) {
    alphamax = alphaMax;
  }

 private:
  double alpha_opt;

 public:
  double getAlphaOpt() const {
    return alpha_opt;
  }

 public:
  void setC1(double f_tol) {
    MoreThuente::f_tol = f_tol;
  }
  void setC2(double gtol) {
    MoreThuente::gtol = gtol;
  }

 private:
  /*
   *
   * phi(alpha) = f(R_k(alpha p_k)), with alpha>0
   *
   * Note that compared to the usual line search, here we have a `reaction` R_k;
   *
   * In the Euclidean case:       R_k (p) = x_k + p;
   * In the Grassmannian case:    R_k (p) = exp ( - ad_(x_k) p ) x_k exp (  ad_(x_k) p )
   *                              x_k in this case is the density matrix at the current iteration and
   *                              p is an element of the tangent space at the current iteration.
   */

  // phi(0) = f(x_k)

  // Maximum alpha
  double alphamax = 65536;
  // Minimum alpha
  double alphamin = 1e-16;

  struct Interval {
    double min;
    double max;
  };

 public:
  MoreThuente(ObjectiveFunction& objFunc) : objectiveFunction(objFunc) {
  }

  [[nodiscard]] auto getIterations() const -> int {
    return iteration;
  }

  void setMaxIterations(int maxIt) {
    max_iterations = maxIt;
  }

  auto evaluate_auxiliary_function_from_phi(const Function& phi_a, const double phi_0_val, const double phi_0_deriv) -> Function {
    Function psi_a;

    psi_a.alpha = phi_a.alpha;
    psi_a.value = phi_a.value - phi_0_val - f_tol * phi_a.alpha * phi_0_deriv;
    psi_a.derivative = phi_a.derivative - f_tol * phi_0_deriv;

    return psi_a;
  }

  auto evaluate_phi_from_auxiliary_function(const Function& psi_a, const double phi_0_val, const double phi_0_deriv) -> Function {
    Function phi_a;

    phi_a.alpha = psi_a.alpha;
    phi_a.value = psi_a.value - phi_0_val + f_tol * phi_a.alpha * phi_0_deriv;
    phi_a.derivative = psi_a.derivative + f_tol * phi_0_deriv;

    return phi_a;
  }

  /*!
   * This is the most important function. It calls the objective function and evaluates it and its gradient.
   * @param alpha
   * @return Function
   */
  auto evaluate_phi_from_objective_function(double alpha) -> Function {
    ++iteration;

    objectiveFunction.evaluate(alpha);

    Function phi;

    phi.value = objectiveFunction.getLineSearchValue();
    // std::cout << "a = " << alpha << std::endl;
    std::cout << "phi(" << std::fixed << std::setprecision(4) << alpha << ") = " << std::setprecision(10) << phi.value
              << std::endl;
    phi.derivative = objectiveFunction.getLineSearchDerivative();
    phi.alpha = alpha;

    return phi;
  }

  /**
   *
   * @param verbose
   */
  auto operator()(double alpha_0, double phi_0_val, double phi_0_deriv) -> void {
    // info is an integer output variable set as follows:
    //
    //   info = 0  Improper input parameters.
    //
    //   info = 1  The sufficient decrease condition and the
    //              directional derivative condition hold.
    //
    //   info = 2  Relative width of the interval of uncertainty
    //            is at most x_tol.
    //
    //   info = 3  Number of calls to df has reached maxfev.
    //
    //   info = 4  The step is at the lower bound alphamin.
    //
    //   info = 5  The step is at the upper bound alphamax.
    //
    //   info = 6  Rounding errors prevent further progress.
    //              There may not be a step which satisfies the
    //              sufficient decrease and curvature conditions.
    //              Tolerances may be too small.
    int info = 0;
    int info_cstep = 1;

    constexpr double zero{};

    if (alpha_0 <= zero || f_tol <= zero || gtol < zero || x_tol < zero || alphamin < zero || alphamax < alphamin ||
        max_iterations <= 0) {
      throw std::domain_error("Input to More-Thuente line search is out of bounds.");
    }

    // Attention: this must be changed in the future!!!
    // TODO: Take value and derivative as function argument!
    const Function phi_0{0, phi_0_val, phi_0_deriv};

    // if (phi_0.derivative >= zero) {
    //  throw std::runtime_error("Search direction is not a direction of descent.");
    //}

    //
    // Initialize local variables
    //

    iteration = 0;

    bool bracketed = false;
    bool stage1 = true;
    double finit = phi_0.value;
    double ginit = phi_0.derivative;
    double dgtest = f_tol * ginit;
    double width = alphamax - alphamin;
    double width1 = 2 * width;

    Function f_x = {zero, finit, ginit};
    Function f_y = {zero, finit, ginit};

    if (!std::isfinite(alpha_0)) {
      alpha_0 = 1;
    }

    // It's not clear where this comes from:
    Interval interval{f_x.alpha, alpha_0 + 4 * (alpha_0 - f_x.alpha)};
    alpha_0 = std::max(alpha_0, alphamin);
    alpha_0 = std::min(alpha_0, alphamax);

    // I think this is unnecessary
    Function phi;
    phi.alpha = alpha_0;

    // Function phi = evaluate_phi_from_objective_function(alpha_0);
    // if (!std::isfinite(phi.value) || !std::isfinite(phi.derivative)) {
    //  throw std::runtime_error("Initial alpha lead to infinite function values.");
    //}

    while (true) {
      //
      // Set the minimum and maximum steps to correspond
      // to the present interval of uncertainty.
      //

      if (bracketed) {
        interval.min = std::min(f_x.alpha, f_y.alpha);
        interval.max = std::max(f_x.alpha, f_y.alpha);
      }
      else {
        interval.min = f_x.alpha;
        interval.max = phi.alpha + 4 * (phi.alpha - f_x.alpha);
      }

      //
      // Ensure stmin and stmax (used in cstep) don't violate alphamin and alphamax
      // Not part of original FORTRAN translation
      //

      interval.min = std::max(alphamin, interval.min);
      interval.max = std::min(alphamax, interval.max);

      //
      // Force the step to be within the bounds alphamax and alphamin
      //

      phi.alpha = std::max(phi.alpha, alphamin);
      phi.alpha = std::min(phi.alpha, alphamax);

      //
      // If an unusual termination is to occur then let
      // alpha be the lowest point obtained so far.
      //

      if ((bracketed && (phi.alpha <= interval.min || phi.alpha >= alphamax)) || iteration >= max_iterations - 1 ||
          info_cstep == 0 || (bracketed && (interval.max - interval.min <= x_tol * interval.max))) {
        phi.alpha = f_x.alpha;
      }

      //
      // Evaluate the function and gradient at alpha
      // and compute the directional derivative.
      //

      phi = evaluate_phi_from_objective_function(phi.alpha);

      if (!std::isfinite(phi.value) || !std::isfinite(phi.derivative)) {
        throw std::runtime_error("Alpha lead to infinite function values.");
      }

      if (std::fabs(phi.derivative - 0) < std::numeric_limits<double>::epsilon()) {
        alpha_opt = phi.alpha;
        return;
      }

      double ftest1 = finit + phi.alpha * dgtest;

      //
      // Test for convergence.
      //

      // What does info_cstep stand for?

      if ((bracketed && (phi.alpha <= interval.min || phi.alpha >= interval.max)) || info_cstep == 0) {
        info = 6;
      }
      if (phi.alpha == alphamax && phi.value <= ftest1 && phi.derivative <= dgtest) {
        info = 5;
      }
      if (phi.alpha == alphamin && (phi.value > ftest1 || phi.derivative >= dgtest)) {
        info = 4;
      }
      if (iteration >= max_iterations) {
        info = 3;
      }
      if (bracketed && (interval.max - interval.min <= x_tol * interval.max)) {
        info = 2;
      }
      if (phi.value <= ftest1 && (std::fabs(phi.derivative) <= -gtol * ginit)) {
        info = 1;
      }

      //
      // Check for termination
      //
      if (info != 0) {
        alpha_opt = phi.alpha;
        return;
      }

      //
      // In the first stage we seek a step for which the modified
      // function has a nonpositive value and nonnegative derivative.
      //

      if (stage1 && (phi.value <= ftest1) && (phi.derivative >= std::min(f_tol, gtol) * ginit)) {
        stage1 = false;
      }

      //
      // A modified function is used to predict the step only if
      // we have not obtained a step for which the modified
      // function has a nonpositive function value and nonnegative
      // derivative, and if a lower function value has been
      // obtained but the decrease is not sufficient.
      //

      if (stage1 && (phi.value < f_x.value) && (phi.value > ftest1)) {
        //
        // Evaluate the modified/auxiliary function and derivatives
        //

        auto psi = evaluate_auxiliary_function_from_phi(phi, finit, ginit);
        auto psi_x = evaluate_auxiliary_function_from_phi(f_x, finit, ginit);
        auto psi_y = evaluate_auxiliary_function_from_phi(f_y, finit, ginit);

        //
        // Call cstep to update the interval of uncertainty
        // and to compute the new step
        //
        cstep(psi_x.alpha, psi_x.value, psi_x.derivative, psi_y.alpha, psi_y.value, psi_y.derivative, psi.alpha,
              psi.value, psi.derivative, bracketed, interval.min, interval.max, info_cstep);

        //
        // Reset the function and gradient values for f.
        //
        phi = evaluate_phi_from_auxiliary_function(psi, finit, ginit);
        f_x = evaluate_phi_from_auxiliary_function(psi_x, finit, ginit);
        f_y = evaluate_phi_from_auxiliary_function(psi_y, finit, ginit);
      }
      else {
        //
        // Call cstep to update the interval of uncertainty
        // and to compute the new step.
        //

        cstep(f_x.alpha, f_x.value, f_x.derivative, f_y.alpha, f_y.value, f_y.derivative, phi.alpha, phi.value,
              phi.derivative, bracketed, interval.min, interval.max, info_cstep);
      }

      // The code should (hopefully) be smart enough to handle this case.
      if (info_cstep == -1) {
        evaluate_phi_from_objective_function(0);
        alpha_opt = 0;
      }

      //
      // Force a sufficient decrease in the size of the
      // interval of uncertainty.
      //

      if (bracketed) {
        if (std::abs(f_y.alpha - f_x.alpha) >= 2. / 3. * width1) {
          phi.alpha = f_x.alpha + (f_y.alpha - f_x.alpha) / 2;
        }
        width1 = width;
        width = std::abs(f_y.alpha - f_x.alpha);
      }
    }
  }

 private:
  /*!
   *
   * Original comment block:
   *
   *  Subroutine cstep
   *
   * The purpose of cstep is to compute a safeguarded step for
   * a linesearch and to update an interval of uncertainty for
   * a minimizer of the function.
   *
   * The parameter stx contains the step with the least function
   * value. The parameter stp contains the current step. It is
   * assumed that the derivative at stx is negative in the
   * direction of the step. If bracketed is set true then a
   * minimizer has been bracketed in an interval of uncertainty
   * with endpoints stx and sty.
   *
   * The subroutine statement is
   *
   * subroutine cstep(stx, fx, dgx,
   *                  sty, fy, dgy,
   *                  stp, f, dg,
   *                  bracketed, alphamin, alphamax, info)
   *
   * where
   *
   * stx, fx, and dgx are variables which specify the step,
   *   the function, and the derivative at the best step obtained
   *   so far. The derivative must be negative in the direction
   *   of the step, that is, dgx and stp-stx must have opposite
   *   signs. On output these parameters are updated appropriately
   *
   * sty, fy, and dgy are variables which specify the step,
   *   the function, and the derivative at the other endpoint of
   *   the interval of uncertainty. On output these parameters are
   *   updated appropriately
   *
   * stp, f, and dg are variables which specify the step,
   *   the function, and the derivative at the current step.
   *   If bracketed is set true then on input stp must be
   *   between stx and sty. On output stp is set to the new step
   *
   * bracketed is a logical variable which specifies if a minimizer
   *   has been bracketed. If the minimizer has not been bracketed
   *   then on input bracketed must be set false. If the minimizer
   *   is bracketed then on output bracketed is set true
   *
   * alphamin and alphamax are input variables which specify lower
   *   and upper bounds for the step
   *
   * info is an integer output variable set as follows:
   *   If info = 1,2,3,4,5, then the step has been computed
   *   according to one of the five cases below. Otherwise
   *   info = 0, and this indicates improper input parameters
   *
   * Argonne National Laboratory. MINPACK Project. June 1983
   * Jorge J. More', David J. Thuente
   */
  static auto cstep(double& stx, double& fx, double& dgx, double& sty, double& fy, double& dgy, double& alpha,
                    double& f, double& dg, bool& bracketed, double& alphamin, double& alphamax, int& info) -> void {
    using std::abs;
    using std::max;

    info = 0;
    bool bound = false;
    double alphaf = -1;

    if (bracketed && (((alpha <= std::min(stx, sty)) || alpha >= std::max(stx, sty)) || dgx * (alpha - stx) >= 0 ||
                      alphamax < alphamin)) {
      std::cout << "Minimizer not bracketed!\n";
      info = -1;
      return;
    }

    // if (std::isnan(f_t.alpha) || std::isnan(f_x.alpha) || std::isnan(f_y.alpha)) {
    //  throw std::runtime_error("Got nan values in the step computation function.");
    //}

    double sgnd = dg * (dgx / std::fabs(dgx));

    // We cover here all the cases presented in the More-Thuente paper in section 4.
    //
    // First case. A higher function value.
    // The minimum is bracketed. If the cubic step is closer
    // to stx than the quadratic step, the cubic step is taken,
    // else the average of the cubic and quadratic steps is taken
    //
    if (f > fx) { // Case 1 from section 4.
      info = 1;
      bound = true;
      double theta = 3 * (fx - f) / (alpha - stx) + dgx + dg;
      // Use s to prevent overflow/underflow of theta^2 and dgx * dg
      double s = max({abs(theta), abs(dgx), abs(dg)});
      double gamma = s * std::sqrt((theta / s) * (theta / s) - (dgx / s) * (dg / s));
      if (alpha < stx) {
        gamma = -gamma;
      }
      double p = gamma - dgx + theta;
      double q = gamma - dgx + gamma + dg;
      double r = p / q;
      double alphac = stx + r * (alpha - stx);
      double alphaq = stx + (dgx / ((fx - f) / (alpha - stx) + dgx)) / 2 * (alpha - stx);

      if (std::fabs(alphac - stx) < std::fabs(alphaq - stx)) {
        alphaf = alphac;
      }
      else {
        alphaf = (alphac + alphaq) / 2;
      }
      bracketed = true;
    }

    //
    // Second case. A lower function value and derivatives of
    // opposite sign. The minimum is bracketed. If the cubic
    // step is closer to stx than the quadratic (secant) step,
    // the cubic step is taken, else the quadratic step is taken
    //
    else if (sgnd < 0.) { // Case 2 from section 4.
      info = 2;
      bound = false;
      double theta = 3 * (fx - f) / (alpha - stx) + dgx + dg;
      // Use s to prevent overflow/underflow of theta^2 and dgx * dg
      double s = max({abs(theta), abs(dgx), abs(dg)});
      double gamma = s * std::sqrt((theta / s) * (theta / s) - (dgx / s) * (dg / s));

      if (alpha > stx) {
        gamma = -gamma;
      }
      double p = gamma - dg + theta;
      double q = gamma - dg + gamma + dgx;
      double r = p / q;

      double alphac = alpha + r * (stx - alpha);
      double alphaq = alpha + (dg / (dg - dgx)) * (stx - alpha);

      if (std::fabs(alphac - alpha) < std::fabs(alphaq - alpha)) {
        alphaf = alphac;
      }
      else {
        alphaf = alphaq;
      }
      bracketed = true;
    }
    //
    // Third case. A lower function value, derivatives of the
    // same sign, and the magnitude of the derivative decreases.
    // The cubic step is only used if the cubic tends to infinity
    // in the direction of the step or if the minimum of the cubic
    // is beyond alpha. Otherwise the cubic step is defined to be
    // either alphamin or alphamax. The quadratic (secant) step is also
    // computed and if the minimum is bracketed then the the step
    // closest to stx is taken, else the step farthest away is taken
    //
    else if (abs(dg) < abs(dgx)) {
      // Case 3 from section 4.
      info = 3;
      bound = true;
      double theta = 3 * (fx - f) / (alpha - stx) + dgx + dg;
      // Use s to prevent overflow/underflow of theta^2 and dgx * dg
      double s = max({abs(theta), abs(dgx), abs(dg)});
      //
      // The case gamma = 0 only arises if the cubic does not tend
      // to infinity in the direction of the step
      //
      double gamma = s * std::sqrt(max(0., (theta / s) * (theta / s) - (dgx / s) * (dg / s)));

      if (alpha > stx) {
        gamma = -gamma;
      }
      double p = gamma - dg + theta;
      double q = gamma + dgx - dg + gamma;
      double r = p / q;
      double alphac;
      double alphaq;
      if (r < 0. && (gamma != 0)) {
        alphac = alpha + r * (stx - alpha);
      }
      else if (alpha > stx) {
        alphac = alphamax;
      }
      else {
        alphac = alphamin;
      }
      alphaq = alpha + (dg / (dg - dgx)) * (stx - alpha);

      if (bracketed) {
        if (abs(alpha - alphac) < abs(alpha - alphaq)) {
          alphaf = alphac;
        }
        else {
          alphaf = alphaq;
        }
      }
      else {
        if (abs(alpha - alphac) > abs(alpha - alphaq)) {
          // alphaf = std::min(alpha+0.66*(sty-alpha), alphaf);
          alphaf = alphac;
        }
        else {
          // alphaf = max(alpha+0.66*(sty-alpha), alphaf);
          alphaf = alphaq;
        }
      }
    }
    //
    // Fourth case. A lower function value, derivatives of the
    // same sign, and the magnitude of the derivative does
    // not decrease. If the minimum is not bracketed, the step
    // is either alphamin or alphamax, else the cubic step is taken
    //
    else {
      info = 4;
      bound = false;
      if (bracketed) {
        double theta = 3 * (f - fy) / (sty - alpha) + dgy + dg;
        // Use s to prevent overflow/underflow of theta^2 and dgy * dg;
        double s = max({abs(theta), abs(dgy), abs(dg)});
        double gamma = s * std::sqrt((theta / s) * (theta / s) - (dgy / s) * (dg / s));

        if (alpha > sty) {
          gamma = -gamma;
        }
        double p = gamma - dg + theta;
        double q = gamma - dg + gamma + dgy;
        double r = p / q;
        double alphac = alpha + r * (sty - alpha);
        alphaf = alphac;
      }
      else if (alpha > stx) {
        alphaf = alphamax;
      }
      else {
        alphaf = alphamin;
      }
    }

    //
    // Update the interval of uncertainty. This update does not
    // depend on the new step or the case analysis above
    //

    if (f > fx) {
      sty = alpha;
      fy = f;
      dgy = dg;
    }
    else {
      if (sgnd < 0) {
        sty = stx;
        fy = fx;
        dgx = dg;
      }
      stx = alpha;
      fx = f;
      dgx = dg;
    }

    alphaf = std::min(alphamax, alphaf);
    alphaf = std::max(alphamin, alphaf);
    alpha = alphaf;

    if (bracketed && bound) {
      if (sty > stx) {
        alpha = std::min(stx + 2. / 3. * (sty - stx), alpha);
      }
      else {
        alpha = std::max(stx + 2. / 3. * (sty - stx), alpha);
      }
    }
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_MORETHUENTE_H
