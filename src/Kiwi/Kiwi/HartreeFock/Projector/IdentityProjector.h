/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_IDENTITYPROJECTOR_H
#define KIWI_IDENTITYPROJECTOR_H

#include <Kiwi/HartreeFock/FockMatrix/BoHelper.h>
#include <Kiwi/HartreeFock/FockMatrix/PreBoHelper.h>
#include <Kiwi/HartreeFock/Projector/Projector.h>
#include <map>

namespace Scine {
namespace Kiwi {

// namespace FockMatrix {
// class BoHelper;
// class PreBoHelper;
//}

class IdentityProjector : public Projector {
  using Projector::Projector;

 public:
  auto init(bool verbose = true) -> void override;

  auto finalize() -> void override {
  }

  auto updateFockMatrices(bool incremental = false, std::size_t numFockRebuilds = 0) -> void override;

  auto evaluateEnergy() -> void override;

  auto getErrorMap() -> std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> override;

  auto getFockMatrixFormationTime() -> std::chrono::duration<double, std::milli> override;

 private:
  std::unique_ptr<std::map<Utils::ElementType, FockMatrix::BoHelper>> boHelper_;
  std::unique_ptr<FockMatrix::PreBoHelper> preBoHelper_;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_IDENTITYPROJECTOR_H
