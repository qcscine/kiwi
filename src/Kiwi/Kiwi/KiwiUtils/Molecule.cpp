/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiUtils/Molecule.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/Settings.h>
#include <iostream>

using namespace Scine;
using namespace Kiwi;

Molecule::Molecule(const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair, MoleculeSettings molSettings)
  : _geometry(structureBoolVectorPair.first),
    _charge(molSettings.charge),
    _multiplicity(molSettings.multiplicity),
    _useHighSpinApprox(molSettings.useHighSpinApproximation),
    _isRestricted(molSettings.isRestricted) {
  if (molSettings.multiplicity < 1) {
    throw std::runtime_error("Multiplicity must be at least 1.");
  }
  if (molSettings.multiplicity != 1 && molSettings.isRestricted) {
    throw std::runtime_error("Restricted calculation is only possible for multiplicity=1.");
  }

  // Point charges:
  if (std::find(structureBoolVectorPair.second.begin(), structureBoolVectorPair.second.end(), false) !=
      structureBoolVectorPair.second.end()) {
    _isPurePreBO = false;
    _hasPointCharges = true;
    for (auto i = 0UL; i < structureBoolVectorPair.second.size(); ++i) {
      if (!structureBoolVectorPair.second[i]) {
        _pointCharges.push_back(structureBoolVectorPair.first.at(i));
      }
    }
  }
  else {
    _isPurePreBO = true;
    _hasPointCharges = false;
  }

  initElectrons();

  if (std::find(structureBoolVectorPair.second.begin(), structureBoolVectorPair.second.end(), true) !=
      structureBoolVectorPair.second.end()) {
    initNuclei(structureBoolVectorPair);
  }

  if (_isPurePreBO) {
    for (auto const& elem : *this) {
      _totalMass += elem.second.N * elem.second.typeInfo.mass;
    }
  }

  std::size_t count = 0;
  for (auto& elem : *this) {
    elem.second.index = count;
    ++count;
  }
}

Molecule::Molecule(std::string basisSetName, const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair,
                   MoleculeSettings molSettings)
  : Molecule(structureBoolVectorPair, molSettings) {
  usePureSpherical_ = molSettings.usePureSpherical;
  Integrals::LibintIntegrals eval;
  Utils::Settings& settings = eval.settings();
  settings.modifyBool("use_pure_spherical", usePureSpherical_);

  if (std::find(structureBoolVectorPair.second.begin(), structureBoolVectorPair.second.end(), true) !=
      structureBoolVectorPair.second.end()) {
    throw std::runtime_error("Electrons-only constructor was called, but I have found a quantum nucleus.");
  }

  (*this)[Utils::ElementType::E].basisSet = eval.initializeBasisSet(basisSetName, _geometry);
  (*this)[Utils::ElementType::E].LAO = (*this)[Utils::ElementType::E].basisSet.nbf();
}

Molecule::Molecule(const std::unordered_map<Utils::ElementType, std::string>& basisSetNameMap,
                   const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair,
                   MoleculeSettings molSettings)
  : Molecule(structureBoolVectorPair, molSettings) {
  usePureSpherical_ = molSettings.usePureSpherical;
  Integrals::LibintIntegrals eval;
  Utils::Settings& settings = eval.settings();
  settings.modifyBool("use_pure_spherical", usePureSpherical_);

  // Electrons + Nuclei
  if (molSettings.onlyElectrons) {
    (*this)[Utils::ElementType::E].basisSet = eval.initializeBasisSet(basisSetNameMap, _geometry);
    (*this)[Utils::ElementType::E].LAO = (*this)[Utils::ElementType::E].basisSet.nbf();
  }
  else {
    for (const auto& name : basisSetNameMap) {
      //    if(name.first!=Utils::ElementType::E && !(std::find(_geometry.getElements().begin(),
      //    _geometry.getElements().end(), name.first)!=_geometry.getElements().end())) {
      //      throw std::runtime_error("Basis set name for type '" + Utils::ElementInfo::symbol(name.first) +
      //                               "' does not match a nucleus, specified in the molcule.");
      //    }
      (*this)[name.first].basisSet = eval.initializeBasisSet(name.second, (*this)[name.first].positions);
      (*this)[name.first].LAO = (*this)[name.first].basisSet.nbf();
    }
  }
}

Molecule::Molecule(const std::unordered_map<Utils::ElementType, Utils::Integrals::BasisSet>& basisSetMap,
                   const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair,
                   MoleculeSettings molSettings)
  : Molecule(structureBoolVectorPair, molSettings) {
  for (const auto& basis : basisSetMap) {
    (*this)[basis.first].basisSet = basis.second;
    (*this)[basis.first].LAO = (*this)[basis.first].basisSet.nbf();
  }
}

Molecule::Molecule(const std::unordered_map<Utils::ElementType, std::string>& electronBasisSetNameMap,
                   const std::unordered_map<Utils::ElementType, std::string>& nuclearBasisSetNameMap,
                   const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair,
                   MoleculeSettings molSettings)
  : Molecule(structureBoolVectorPair, molSettings) {
  Integrals::LibintIntegrals eval;
  Utils::Settings& settings = eval.settings();
  settings.modifyBool("use_pure_spherical", molSettings.usePureSpherical);

  // Electrons:
  (*this)[Utils::ElementType::E].basisSet = eval.initializeBasisSet(electronBasisSetNameMap, _geometry);
  (*this)[Utils::ElementType::E].LAO = (*this)[Utils::ElementType::E].basisSet.nbf();

  // Nuclei:
  for (const auto& name : nuclearBasisSetNameMap) {
    if ((*this).find(name.first) == (*this).end()) {
      throw std::runtime_error("Basis set name for type '" + Utils::ElementInfo::symbol(name.first) +
                               "' does not match a nucleus, specified in the molcule");
    }
    (*this)[name.first].basisSet = eval.initializeBasisSet(name.second, (*this)[name.first].positions);
    (*this)[name.first].LAO = (*this)[name.first].basisSet.nbf();
  }
}

auto Molecule::addAdditionalNuclearCenters(const Utils::AtomCollection& atoms,
                                           const std::unordered_map<Utils::ElementType, std::string>& nuclearBasisSetNameMap)
    -> void {
  Integrals::LibintIntegrals eval;
  Utils::Settings& settings = eval.settings();
  settings.modifyBool("use_pure_spherical", usePureSpherical_);

  // const std::vector<ParticleType>& particleTypeInfo = Utils::Integrals::getParticleTypeInfo();

  for (const auto& atom : atoms) {
    Utils::AtomCollection atomc({atom.getElementType()}, atom.getPosition());
    auto type = atom.getElementType();
    (*this)[type].positions.push_back(atom);
    auto additionalBasis = eval.initializeBasisSet(nuclearBasisSetNameMap.at(type), atomc);
    (*this)[type].basisSet.append(additionalBasis);

    (*this)[type].basisSet.deleteShellPairs();
    eval.generateShellPairs((*this)[type].basisSet);

    (*this)[type].LAO = (*this)[type].basisSet.nbf();
  }
}

void Molecule::initNuclei(const std::pair<Utils::AtomCollection, std::vector<bool>>& structureBoolVectorPair) {
  // Stores the indicies in the AtomCollection, at which the quantum nuclei are.
  std::unordered_map<Utils::ElementType, std::vector<size_t>> quantumNucleiMap;

  _isPurePreBO = true;
  for (auto i = 0; i < _geometry.size(); ++i) {
    if (structureBoolVectorPair.second[i]) {
      quantumNucleiMap[_geometry.at(i).getElementType()].push_back(i);
    }
    else {
      _isPurePreBO = false;
    }
  }
  if (_isPurePreBO) {
    _totalMass = 0.;
  }

  const std::vector<ParticleType>& particleTypeInfo = Utils::Integrals::getParticleTypeInfo();

  // Quantum nuclei:
  for (auto const& type : quantumNucleiMap) {
    (*this)[type.first] = ParticleTypeData();
    Utils::ElementTypeCollection elements;
    for (auto i = 0UL; i < quantumNucleiMap[type.first].size(); ++i) {
      elements.push_back(type.first);
    }
    Utils::PositionCollection positions;
    positions.resize(type.second.size(), Eigen::NoChange_t());
    auto counter = 0;
    for (auto const& index : quantumNucleiMap[type.first]) {
      positions(counter, 0) = _geometry.getPosition(index)(0);
      positions(counter, 1) = _geometry.getPosition(index)(1);
      positions(counter, 2) = _geometry.getPosition(index)(2);
      ++counter;
    }
    Utils::AtomCollection atoms(elements, positions);
    (*this)[type.first].positions = atoms;
    (*this)[type.first].N = (*this)[type.first].positions.size();
    bool found = false;
    for (auto const& pt : particleTypeInfo) {
      if (pt.symbol == type.first) {
        (*this)[type.first].typeInfo = pt;
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error("Sorry, there is no data implemented for nuclei of type " +
                               Utils::ElementInfo::symbol(type.first));
    }
    // spin = 2 S, and 2S+1 is the dimension of the rerpesentation for spin S.
    (*this)[type.first].msVector.resize((*this)[type.first].typeInfo.spin + 1, 0);
    if (_useHighSpinApprox) {
      (*this)[type.first].msVector[0] = (*this)[type.first].N;
    }
    // Otherwise: choose lowest spin. No multiplicity implemented, yet.
    else {
      auto Ntmp = (*this)[type.first].N;
      while (Ntmp > 0) {
        for (auto& elem : (*this)[type.first].msVector) {
          elem += 1;
          Ntmp -= 1;
          if (Ntmp == 0) {
            break;
          }
        }
      }
    }
  }
}

void Molecule::initElectrons() {
  const std::vector<ParticleType>& particleTypeInfo = Utils::Integrals::getParticleTypeInfo();

  // Electrons:
  (*this)[Utils::ElementType::E] = ParticleTypeData();
  (*this)[Utils::ElementType::E].positions = _geometry;
  for (auto const& pt : particleTypeInfo) {
    if (pt.symbol == Utils::ElementType::E) {
      (*this)[Utils::ElementType::E].typeInfo = pt;
      break;
    }
  }
  if (_isRestricted) {
    (*this)[Utils::ElementType::E].isRestricted = true;
  }
  std::size_t Nel = 0;
  for (auto const& atom : _geometry) {
    Nel += Utils::ElementInfo::Z(atom.getElementType());
  }
  Nel -= _charge;
  if (_isRestricted && ((Nel % 2) != 0)) {
    throw std::runtime_error("Restricted calculation only possible if number of electrons is even.");
  }
  (*this)[Utils::ElementType::E].N = Nel;
  (*this)[Utils::ElementType::E].msVector.resize(2, 0);
  // Even number of electrons
  _hasUnpairedElectrons = true;
  auto unpairedElectrons = _multiplicity - 1;
  if (Nel % 2 == 0) {
    if ((unpairedElectrons % 2) != 0U) {
      throw std::runtime_error("Multiplicity cannot be even, if number of electrons is even.");
    }
    (*this)[Utils::ElementType::E].msVector[0] = Nel / 2;
    (*this)[Utils::ElementType::E].msVector[1] = Nel / 2;
    if (unpairedElectrons == 0) {
      _hasUnpairedElectrons = false;
    }
    else {
      (*this)[Utils::ElementType::E].msVector[0] += unpairedElectrons / 2;
      (*this)[Utils::ElementType::E].msVector[1] -= unpairedElectrons / 2;
    }
  }
  // Odd number of electrons
  else {
    if ((unpairedElectrons % 2) != 1) {
      throw std::runtime_error("Multiplicity does not match number of unpaired electrons!");
    }
    (*this)[Utils::ElementType::E].msVector[0] = Nel / 2 + 1;
    (*this)[Utils::ElementType::E].msVector[1] = Nel / 2;
    if (unpairedElectrons > 1) {
      (*this)[Utils::ElementType::E].msVector[0] += (unpairedElectrons - 1) / 2;
      (*this)[Utils::ElementType::E].msVector[1] -= (unpairedElectrons - 1) / 2;
    }
  }
}

const Utils::AtomCollection& Molecule::getGeometry() const {
  return _geometry;
}

const Utils::AtomCollection& Molecule::getPointCharges() const {
  return _pointCharges;
}

bool Molecule::hasPointCharges() const {
  return _hasPointCharges;
}

size_t Molecule::getNumberOfTypes() const {
  return _numberOfTypes;
}

size_t Molecule::getNumberOfParticles() const {
  return _numberOfParticles;
}

int Molecule::getCharge() const {
  return _charge;
}

size_t Molecule::getMultiplicity() const {
  return _multiplicity;
}

bool Molecule::useHighSpinApprox() const {
  return _useHighSpinApprox;
}

double Molecule::getTotalMass() const {
  return _totalMass;
}

bool Molecule::isPurePreBO() const {
  return _isPurePreBO;
}

bool Molecule::isRestricted() const {
  return _isRestricted;
}

auto Molecule::getBOMolecule() const -> Molecule {
  Molecule ret(*this);
  ret.makeBO();

  return ret;
}

auto Molecule::makeBO() -> void {
  _isPurePreBO = false;
  _hasPointCharges = true;
  _pointCharges = _geometry;
  _numberOfTypes = 1;
  _numberOfParticles = this->at(Utils::ElementType::E).N;

  std::vector<Utils::ElementType> eraseVector;

  for (const auto& elem : *this) {
    if (elem.first != Utils::ElementType::E) {
      eraseVector.push_back(elem.first);
    }
  }
  for (auto const& elem : eraseVector) {
    (*this).erase(elem);
  }
}

auto Molecule::uncontractElectronicBasis() -> void {
  this->at(Utils::ElementType::E).basisSet.uncontract();
  Integrals::LibintIntegrals::generateShellPairs(this->at(Utils::ElementType::E).basisSet);
  (*this)[Utils::ElementType::E].LAO = this->at(Utils::ElementType::E).basisSet.nbf();
}
const std::string& Molecule::getInputFileName() const {
  return inputFileName_;
}
void Molecule::setInputFileName(const std::string& inputFileName) {
  inputFileName_ = inputFileName;
}
