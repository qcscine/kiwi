/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_DATA_H
#define KIWI_DATA_H

#include <Kiwi/KiwiUtils/Molecule.h>
#include <Kiwi/KiwiUtils/SymmetryEnums.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Geometry/ElementTypes.h>
#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>
#include <unordered_map>

namespace Scine {
namespace Kiwi {

using Matrix = Eigen::MatrixXd;
using ElementPair = std::pair<Utils::ElementType, Utils::ElementType>;
using Grid = std::pair<std::vector<Eigen::VectorXd>, std::vector<double>>;

// Note: < is defiend as >, such that electrons always come first.
// In the future, ElementType will be substituted by ParticleType, but for now, it is fine.
inline bool operator<(Utils::ElementType left, Utils::ElementType right) {
  return static_cast<std::underlying_type_t<Utils::ElementType>>(left) >
         static_cast<std::underlying_type_t<Utils::ElementType>>(right);
}

struct HartreeFockData {
  enum class FockMatrixContribution { F, G, J, K, L };

  /* Last density matrices --> for incremental Fock build*/
  std::map<Utils::ElementType, Utils::DensityMatrix> D_last;
  /* Fock matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> F;
  /* Orthogonalized Fock matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> F_OAO;
  /* G matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> G;
  /* I(nteraction) matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> I;
  /* L: specific interaction matrices */
  std::map<Utils::ElementType, std::map<Utils::ElementType, Matrix>> L;
  /* Coulomb matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> K;
  /* Exchange matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> J;
};

/**
 * @class Data
 * @brief Stores all important integrals and matrices.
 *
 * This class contains all essential data for all HF and post-HF methods.
 *
 */
class Data {
 public:
  Data() {
    hartreeFockData = std::make_shared<HartreeFockData>();
  }

  Data(const std::shared_ptr<Molecule>& molecule);

  ~Data() = default;

  /* Molecule: contains all physical information */
  std::shared_ptr<Molecule> molecule;
  /* Core matrices */
  std::map<Utils::ElementType, Matrix> H;
  /* Overlap matrices */
  std::map<Utils::ElementType, Matrix> S;
  /* Transformation matrices */
  std::map<Utils::ElementType, Matrix> X;
  /* MO-Coefficient matrices */
  std::map<Utils::ElementType, Utils::MolecularOrbitals> C;
  /* Orthogonal MO-Coefficient matrices */
  std::map<Utils::ElementType, Utils::MolecularOrbitals> C_OAO;
  /* Density matrices */
  std::map<Utils::ElementType, Utils::DensityMatrix> D;
  /* Density matrices */
  std::map<Utils::ElementType, Utils::DensityMatrix> D_OAO;
  /* Exchange matrices */
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> O_NO;
  /* Coulomb super-matrix */
  std::map<ElementPair, Matrix> Coulomb;
  /* Exchange super-matrix */
  std::map<Utils::ElementType, Matrix> Exchange;
  /* Unique pairs of elements. For easy looping */
  std::vector<ElementPair> uniquePairs;
  /* Hartree-Fock energy */
  double HartreeFockEnergy;
  /* Hartree-Fock specific matrices */
  std::shared_ptr<HartreeFockData> hartreeFockData;
  /* Projection grid */
  Grid grid;
  /* Two-integrals in memory? */
  bool integralDirect = true;
  /* Flag that shows whether the orbitals are natural orbitals */
  bool naturalOrbitals = false;

  /**
   * Generates vector of `ElementPair`s. Can be used for easy looping over unique pairs of particle types.
   * @param mol
   * @return
   */
  static auto generateUniquePairs(const std::shared_ptr<Molecule>& mol) -> std::vector<ElementPair>;

  void oneBodyIntegrals(bool prinTimings = true);

  void twoBodyIntegrals(bool prinTimings = true);

  void makeExchange(bool prinTimings = true);

  void makeLoewdinOrtho(double loewdinThresh = 10e-10, const bool verbose = true);

  static inline auto getUnique(Utils::ElementType left, Utils::ElementType right) -> ElementPair {
    ElementPair ret;
    if (left < right) {
      ret = {left, right};
    }
    else {
      ret = {right, left};
    }
    return ret;
  }

  /**
   * @brief Get the 2-body tensor index for a super matrix.
   */
  inline static auto index2(const int& i, const int& j, const int& dim) -> int {
    return i * dim + j;
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_DATA_H
