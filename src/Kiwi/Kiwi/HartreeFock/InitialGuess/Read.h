/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_READ_H
#define KIWI_READ_H

#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Utils/Geometry/ElementTypes.h>
#include <memory>
#include <string>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;

/**
 * @class Read @file Read.h
 * @brief This class is a wrapper around the MolecularOrbitalsIO class which handles storing and reading the Molecular
 * Orbitals and writing them into a Data object.
 */
class Read {
 private:
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;

  const std::string extension_ = ".orbitals";

 public:
  Read() = default;

  Read(std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule)
    : data_(std::move(data)), molecule_(std::move(molecule)) {
  }

  auto readOrbitals() -> void;

  auto writeOrbitals() -> void;

  auto cleanUp() -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_READ_H
