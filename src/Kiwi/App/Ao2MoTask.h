/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_AO2MOTASK_H
#define KIWI_AO2MOTASK_H

#include "Task.h"

namespace Scine {
namespace Kiwi {

class Ao2MoTask : public Task {
 private:
  bool _write = false;
  double _integralThresh = 10e-16;
  bool useNaturalOrbitals_ = false;

 public:
  // ReferenceCalculationTask() = default;

  Ao2MoTask(YAML::Node& input);

  // ReferenceCalculationTask(const ReferenceCalculationTask&) = delete;
  // ReferenceCalculationTask& operator=(const ReferenceCalculationTask&) = delete;
  //~ReferenceCalculationTask() override = default;

  auto name() const -> const std::string final;

  auto run() -> void final;

  auto settings(YAML::Node& input) -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_AO2MOTASK_H
