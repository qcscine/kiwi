/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_OPTIMIZER_H
#define KIWI_OPTIMIZER_H

namespace Scine {
namespace Kiwi {
namespace Optimization {

/*!
 * Abstract base class for an optimizer.
 */
class Optimizer {
 protected:
  double thresh_ = 1e-10;

  double phi_0_value;
  double phi_0_derivative;

  double current_alpha;
  double phi_alpha_value;
  double phi_alpha_derivative;

  // Optimizer() = 0;

 public:
  virtual auto evaluate(double alpha) -> void = 0;

  virtual auto evaluateDirection() -> void = 0;

  virtual auto getLineSearchValue() const -> double final {
    return phi_alpha_value;
  }

  virtual auto getLineSearchDerivative() const -> double final {
    return phi_alpha_derivative;
  }

  virtual auto applyUpdate() -> void = 0;

  virtual auto getError() -> double = 0;

  virtual auto getZeroValue() const -> double final {
    return phi_0_value;
  }

  virtual auto getZeroDerivative() const -> double final {
    return phi_0_derivative;
  }

  virtual auto printIteration() -> void = 0;
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_OPTIMIZER_H
