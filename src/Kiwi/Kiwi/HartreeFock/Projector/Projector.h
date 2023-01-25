/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_PROJECTOR_H
#define KIWI_PROJECTOR_H

#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Eigen/Dense>
#include <memory>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;

struct ProjectionParameters {
  int N = 0;
  int parity = 1;
  std::vector<Utils::ElementType> typeVector;
};

class Projector {
 public:
  Projector(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data,
            ProjectionParameters projectionParameters = ProjectionParameters())
    : molecule_(std::move(molecule)), data_(std::move(data)), projectionParameters_(std::move(projectionParameters)) {
  }

  virtual ~Projector() = default;

 protected:
  std::shared_ptr<Molecule> molecule_;

  std::shared_ptr<Data> data_;

  ProjectionParameters projectionParameters_;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_;

  double energy_;

  std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> errorMap_;

  std::map<Utils::ElementType, Utils::DensityMatrix> incrementalDensityMatrix_;

  std::map<Utils::ElementType, Utils::DensityMatrix> lastDensityMatrix_;

  std::chrono::duration<double, std::milli> fockMatrixFormationTime_;

 public:
  virtual std::chrono::duration<double, std::milli> getFockMatrixFormationTime();

 protected:
  auto computeNaturalOrbitals(Utils::ElementType type) -> void;

 public:
  virtual auto init(bool verbose = true) -> void = 0;

  virtual auto finalize() -> void = 0;

  virtual auto updateDensity(Utils::ElementType type) -> void;

  virtual auto resetDensity(const std::map<Utils::ElementType, Utils::DensityMatrix>& D,
                            const std::map<Utils::ElementType, Utils::DensityMatrix>& D_OAO,
                            const std::map<Utils::ElementType, Utils::MolecularOrbitals>& C,
                            const std::map<Utils::ElementType, Utils::MolecularOrbitals>& C_OAO) -> void;

  virtual auto setLastDensity(const std::map<Utils::ElementType, Utils::DensityMatrix>& lastDensity) -> void final;

  virtual auto evaluateDensity(bool incremental = false) -> void final;

  virtual auto updateDensityFromRotation(Utils::ElementType type, const Utils::SpinAdaptedMatrix& orbitalRotation) -> void;

  virtual auto updateFockMatrices(bool incremental = false, std::size_t numFockRebuilds = 0) -> void = 0;

  virtual auto evaluateEnergy() -> void = 0;

  auto getEnergy() const -> double {
    return energy_;
  }

  virtual auto getErrorMap() -> std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> = 0;

  auto getProjectionParameters() const -> ProjectionParameters {
    return projectionParameters_;
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_PROJECTOR_H
