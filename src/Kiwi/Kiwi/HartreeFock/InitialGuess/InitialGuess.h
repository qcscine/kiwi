/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_INITIALGUESS_H
#define KIWI_INITIALGUESS_H

#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Geometry/ElementTypes.h>
#include <map>
#include <memory>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;

/**
 * @class InitialGuess
 * @brief This class is responsible for performing the initial guess for the density matrices of electrons and nuclei.
 */
class InitialGuess {
 private:
  Utils::DensityMatrix D_;
  Utils::MolecularOrbitals C_;
  bool _verbose = true;
  Utils::ElementType _type;

 public:
  InitialGuess(const InitialGuessSCF guess, const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
               const Utils::ElementType type, bool verbose = true);

  auto getDensity() const -> Utils::DensityMatrix;

  auto getMOs() const -> Utils::MolecularOrbitals;

  auto setVerbose(bool verbose) -> void;

  template<InitialGuessSCF>
  static auto perform(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data, const Utils::ElementType type,
                      Utils::DensityMatrix& D, Utils::MolecularOrbitals& C, bool verbose = true) -> void;

  static auto CoreGuessCoefficients(const Eigen::MatrixXd& H, const Eigen::MatrixXd& X) -> Eigen::MatrixXd;

  static auto HueckelGuessCoefficients(const Eigen::MatrixXd& coefficients, const Eigen::MatrixXd& S,
                                       const Eigen::MatrixXd& X, const Eigen::VectorXd& ionizationPotentials) -> Eigen::MatrixXd;

  /**
   * @brief This is a Hueckel guess, similar to the one implemented in Serenity. (Thanks Jan!)
   * @note It can only be used for electrons. The other Hueckel guess can be used for nuclei, too, although I am not
   *       sure if it makes sense.
   * @param mol
   * @param dataa
   * @return
   */
  static auto HueckelSTO3gGuessCoefficients(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& dataa,
                                            bool verbose = true) -> Eigen::MatrixXd;

  /**
   * @brief Function that returns the negative ionization potentials of 1st and 2nd row elements, i.e, the parameters
   *        for the extended Hueckel guess.
   */
  static auto EHTParameters() -> std::map<Utils::ElementType, std::vector<double>>;

  static auto finalize(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data, const Utils::ElementType type,
                       const Utils::DensityMatrix& D, const Utils::MolecularOrbitals& C) -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_INITIALGUESS_H
