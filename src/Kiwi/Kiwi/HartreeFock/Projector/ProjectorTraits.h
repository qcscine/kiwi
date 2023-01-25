/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_PROJECTORTRAITS_H
#define KIWI_PROJECTORTRAITS_H

#include <Kiwi/HartreeFock/Projector/IdentityProjector.h>
#include <Kiwi/HartreeFock/SymmetryProjection/SymmetryTypes.h>

namespace Scine {
namespace Kiwi {

template<SymmetryType Symmetry>
class ProjectorTrait {};

template<>
class ProjectorTrait<SymmetryType::None> {
 public:
  using ProjectorType = IdentityProjector;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_PROJECTORTRAITS_H
