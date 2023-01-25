/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_MOLECULE_H
#define KIWI_MOLECULE_H

#include <Utils/DataStructures/BasisSet.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Kiwi {

using ParticleType = Utils::Integrals::ParticleType;
using Spin = Utils::Integrals::Spin;

struct MoleculeSettings {
  bool isRestricted = true;
  bool useHighSpinApproximation = true;
  bool usePureSpherical = true;
  bool onlyElectrons = false;
  size_t multiplicity = 1;
  int charge = 0;
};

struct SpinadaptedInt {
  int restricted;
  int alpha;
  int beta;
};

struct ParticleTypeData {
  Utils::AtomCollection positions;
  Utils::Integrals::BasisSet basisSet;
  std::size_t N;
  std::size_t LAO;
  std::size_t LMO = 0;
  std::vector<Spin> msVector;
  ParticleType typeInfo;
  bool isRestricted = false;
  std::size_t index;
  SpinadaptedInt occ;
  SpinadaptedInt virt;
};

/**
 * @class Molecule
 * @file Molecule.h
 * @brief Class that stores the data required for performing calculations:
 *           - structure of the entire molecule_
 *           - point charges
 *           - the basis sets for all particle types
 *           - Spin, multiplicity, charge
 *           - ...
 */
class Molecule : public std::map<Utils::ElementType, ParticleTypeData, std::greater<>> { // greater guarantees that
                                                                                         // electrons come first -> will
                                                                                         // be fixed in the future
 public:
  Molecule() = default;

  ~Molecule() = default;

  Molecule(const Molecule&) = default;

  /**
   * @brief Most basic constructor, for an electrons-only calculation.
   *        Note: structureBoolVectorPair.second must be empty... Otherwise, it would not be electrons-only.
   * @param basisSetName
   * @param structureBoolVectorPair
   * @param multiplicity
   * @param charge
   * @param useHighSpinApprox
   */
  // Molecule(std::string basisSetName, const std::pair<Utils::AtomCollection, std::vector<bool>>&
  // structureBoolVectorPair,
  //         size_t multiplicity = 1, int charge = 0, bool useHighSpinApprox = true, bool usePureSpherical = true,
  //         bool isRest = true);

  // Molecule(const std::unordered_map<Utils::ElementType, std::string>& basisSetNameMap,
  //         const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, size_t multiplicity = 1,
  //         int charge = 0, bool useHighSpinApprox = true, bool usePureSpherical = true, bool isRest = true,
  //         bool onlyElectrons = false);

  // Molecule(const std::unordered_map<Utils::ElementType, Utils::Integrals::BasisSet>& basisSetMap,
  //         const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, size_t multiplicity =
  //         1, int charge = 0, bool useHighSpinApprox = true, bool isRest = true);

  // Molecule(const std::unordered_map<Utils::ElementType, std::string>& electronBasisSetNameMap,
  //         const std::unordered_map<Utils::ElementType, std::string>& nuclearBasisSetNameMap,
  //         const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, size_t multiplicity =
  //         1, int charge = 0, bool useHighSpinApprox = true, bool usePureSpherical = true, bool isRest = true);

  // Molecule(const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, size_t multiplicity =
  // 1,
  //         int charge = 0, bool useHighSpinApprox = true, bool isRest = true);

  Molecule(std::string basisSetName, const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair,
           MoleculeSettings molSettings);

  Molecule(const std::unordered_map<Utils::ElementType, std::string>& basisSetNameMap,
           const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, MoleculeSettings molSettings);

  Molecule(const std::unordered_map<Utils::ElementType, Utils::Integrals::BasisSet>& basisSetMap,
           const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, MoleculeSettings molSettings);

  Molecule(const std::unordered_map<Utils::ElementType, std::string>& electronBasisSetNameMap,
           const std::unordered_map<Utils::ElementType, std::string>& nuclearBasisSetNameMap,
           const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, MoleculeSettings molSettings);

  Molecule(const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, MoleculeSettings molSettings);

  const Utils::AtomCollection& getGeometry() const;
  const Utils::AtomCollection& getPointCharges() const;
  bool hasPointCharges() const;
  size_t getNumberOfTypes() const;
  size_t getNumberOfParticles() const;
  int getCharge() const;
  size_t getMultiplicity() const;
  bool useHighSpinApprox() const;
  bool isPurePreBO() const;
  bool isRestricted() const;
  double getTotalMass() const;

  auto addAdditionalNuclearCenters(const Utils::AtomCollection& atoms,
                                   const std::unordered_map<Utils::ElementType, std::string>& nuclearBasisSetNameMap) -> void;

  Molecule getBOMolecule() const;

  auto makeBO() -> void;

  auto uncontractElectronicBasis() -> void;

 private:
  // geometry collects the entire molecule_, i.e., point charges and quantum nuclei.
  Utils::AtomCollection _geometry;
  // pointCharges collects only the classical nuclei, i.e., the point charges.
  Utils::AtomCollection _pointCharges;
  bool _hasPointCharges;
  bool _isPurePreBO;
  size_t _numberOfTypes;
  size_t _numberOfParticles;
  // This Map stores all the relevant data about all particle types.
  //  ParticleDataMap _particleDataMap;
  // charge of the molecule_
  int _charge;
  // multiplicity of the electrons: 2S+1
  size_t _multiplicity;
  // Unpaired electrons
  bool _hasUnpairedElectrons;
  // By default, the nuclei of all types are assumed to be in the configuration with only a certain m_s (the highes one)
  // E.g., all protons are assumed to be in m_s = +1/2
  bool _useHighSpinApprox = true;
  double _totalMass = 0;
  // Restricted or unrestricted?
  bool _isRestricted;

  bool usePureSpherical_ = true;

  std::string inputFileName_;

 public:
  const std::string& getInputFileName() const;
  void setInputFileName(const std::string& inputFileName);

 private:
  /**
   * @brief Populate the electron-related data
   */
  void initElectrons();
  /**
   * @brief Populate the nuclei-related data
   */
  void initNuclei(const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair);
};
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_MOLECULE_H
