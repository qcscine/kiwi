/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_PARTICLEDENSITYTASK_H
#define KIWI_PARTICLEDENSITYTASK_H

#include "Task.h"
#include <Kiwi/KiwiUtils/ParticleDensity/Evaluator.h>

namespace Scine {
namespace Kiwi {

class ParticleDensityTask : public Task {
 private:
  ParticleDensity::Settings settings_;
  bool useExternalRdm_ = false;
  bool useNaturalOrbitals_ = false;
  std::string rdmFileNameRestricted_;
  std::string rdmFileNameAlpha_;
  std::string rdmFileNameBeta_;

 public:
  ParticleDensityTask(YAML::Node& input);

  auto name() const -> const std::string final;

  auto run() -> void final;

  auto settings(YAML::Node& input) -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_PARTICLEDENSITYTASK_H
