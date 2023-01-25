/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_EVALUATOR_H
#define KIWI_EVALUATOR_H

#include <Kiwi/KiwiUtils/Data.h>
#include <Eigen/Dense>
#include <iostream>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
enum class ElementType : unsigned;
}
namespace Kiwi {

class Molecule;
class Data;

namespace ParticleDensity {

class Engine;
enum class DensityMode { MO, PartDens };
enum class Format { Grid, Cube };

struct Settings {
  double xmin = 0;
  double xmax = 0;
  double ymin = 0;
  double ymax = 0;
  double zmin = 0;
  double zmax = 0;

  int N1D = 100;

  int N3DX = 100;
  int N3DY = 100;
  int N3DZ = 100;

  DensityMode densityMode;

  Format format;

  int moIndex;

  int msIndex = 0;

  std::string filename;

  Utils::ElementType type;
};

/**
 * @class Evaluator
 * @brief This class handles the calculation of the particle density on a grid, given a density matrix in the MO basis.
 * TODO change from MO to AO basis to speed up calculation.
 */
class Evaluator {
 public:
  Evaluator(Settings settings, std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule)
    : settings_(std::move(settings)), data_(std::move(data)), molecule_(std::move(molecule)) {
  }

  auto setDensityMatrix(const Utils::DensityMatrix& D) -> void {
    densityMatrix.alphaMatrix() = D.alphaMatrix();
    densityMatrix.betaMatrix() = D.betaMatrix();
    densityMatrix.restrictedMatrix() = D.restrictedMatrix();
  }

  auto setDensityMatrix(const Utils::SpinAdaptedMatrix& D) -> void {
    densityMatrix = D;
  }

  void constructGrid();

  void constructCube();

  void runCalculationGrid();

  void runCalculationCube();

  void writeGrid() const;

  void writeCube() const;

 private:
  Settings settings_;
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;

  // The density matrix in the MO basis. Is a SpinAdaptedMatrix here because not properties of the DensityMatrix are
  // required.
  Utils::SpinAdaptedMatrix densityMatrix;

  // `grid` stores the density as a one-dimensional array with the according coordinates.
  std::vector<std::pair<Eigen::Vector3d, double>> grid;
  // `cube` stores the density as required for a cube file. The coordinates are implicit here.
  std::vector<std::vector<std::vector<double>>> cube;

  Eigen::Vector3d step;

  void compute(double& destination, Engine& engine, double x, double y, double z);
};

} // namespace ParticleDensity
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_EVALUATOR_H
