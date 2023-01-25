/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_NATURALORBITALSTASK_H
#define KIWI_NATURALORBITALSTASK_H

#include "Task.h"
#include <Utils/Geometry/ElementTypes.h>
#include <map>
#include <string>

namespace Scine {
namespace Kiwi {

class NaturalOrbitalsTask : public Task {
  std::map<Utils::ElementType, std::string> elementType_LineMap_;

 public:
  NaturalOrbitalsTask(YAML::Node& input);

  auto name() const -> const std::string final;

  auto run() -> void final;

  auto settings(YAML::Node& input) -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_NATURALORBITALSTASK_H
