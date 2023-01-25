/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_NATURALORBITALS_H
#define KIWI_NATURALORBITALS_H

#include <Utils/Geometry/ElementTypes.h>
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace Scine {
namespace Kiwi {

class Data;
class Molecule;

/**
 * @class NaturalOrbitals @file NaturalOrbitals.h
 * @brief Given a map that contains the ElementType and a string with one or multiple rdm-files, this is used for
 * generating natural orbitals. If the type is unrestricted, first the alpha, then the beta rdm must be written in the
 * line.
 * Moreover, the NOs are automatically pruned if the natural occupation number is below 1e-16.
 */
class NaturalOrbitals {
 private:
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;
  std::map<Utils::ElementType, std::string> elementType_LineMap_;

 public:
  NaturalOrbitals(std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule,
                  const std::map<Utils::ElementType, std::string>& elementType_LineMap)
    : data_(std::move(data)), molecule_(std::move(molecule)), elementType_LineMap_(elementType_LineMap) {
  }

  auto generateNaturalOrbitals(bool prune = true) -> void;

 private:
  auto generateOrbitals(const Eigen::MatrixXd& rdm, bool prune = true) -> Eigen::MatrixXd;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_NATURALORBITALS_H
