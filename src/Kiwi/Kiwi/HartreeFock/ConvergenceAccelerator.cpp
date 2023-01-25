/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/ConvergenceAccelerator.h>
#include <utility>

using namespace Scine;
using namespace Kiwi;

ConvergenceAccelerator::ConvergenceAccelerator(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data,
                                               Utils::scf_mixer_t mixer, bool verbose)
  : _molecule(std::move(molecule)), _data(std::move(data)), _mixer(mixer), verbose_(verbose) {
  switch (_mixer) {
    case Utils::scf_mixer_t::fock_diis:
      for (auto const& elem : *_molecule) {
        _fockDiisVector[elem.first] = std::make_unique<Utils::FockDiis>();
        _fockDiisCounter[elem.first] = 0;
        _ediisStarted[elem.first] = false;
        _fockDiisStarted[elem.first] = false;
      }
      break;
    case Utils::scf_mixer_t::ediis:
      for (auto const& elem : *_molecule) {
        _ediisVector[elem.first] = std::make_unique<Utils::Ediis>();
        _eDiisCounter[elem.first] = 0;
        _ediisStarted[elem.first] = false;
        _fockDiisStarted[elem.first] = false;
      }
      break;
    case Utils::scf_mixer_t::ediis_diis:
      for (auto const& elem : *_molecule) {
        _ediisVector[elem.first] = std::make_unique<Utils::Ediis>();
        _fockDiisVector[elem.first] = std::make_unique<Utils::FockDiis>();
        _eDiisCounter[elem.first] = 0;
        _ediisStarted[elem.first] = false;
        _fockDiisCounter[elem.first] = 0;
        _fockDiisStarted[elem.first] = false;
      }
    case Utils::scf_mixer_t::none:
      break;
    case Utils::scf_mixer_t::fock_simple:
      throw std::runtime_error("fock_simple and charge_simple not implemted.");
      break;
    case Utils::scf_mixer_t::charge_simple:
      throw std::runtime_error("fock_simple and charge_simple not implemted.");
      break;
  }
}

auto ConvergenceAccelerator::init() -> void {
  if (_mixer == Utils::scf_mixer_t::fock_diis) {
    for (auto const& elem : *_molecule) {
      _fockDiisVector[elem.first]->setNAOs(elem.second.LMO);
      _fockDiisVector[elem.first]->setUnrestricted(!elem.second.isRestricted);
      _fockDiisVector[elem.first]->setSubspaceSize(_subspaceSize);
      _fockDiisVector[elem.first]->setOrthogonal(true);
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
      _fockDiisCounter[elem.first] = 0;
      _fockDiisStarted[elem.first] = false;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis) {
    for (auto const& elem : *_molecule) {
      _ediisVector[elem.first]->setNAOs(elem.second.LMO);
      _ediisVector[elem.first]->setUnrestricted(!elem.second.isRestricted);
      _ediisVector[elem.first]->setSubspaceSize(_subspaceSize);
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis_diis) {
    for (auto const& elem : *_molecule) {
      _ediisVector[elem.first]->setNAOs(elem.second.LMO);
      _ediisVector[elem.first]->setUnrestricted(!elem.second.isRestricted);
      _ediisVector[elem.first]->setSubspaceSize(_subspaceSize);
      _fockDiisVector[elem.first]->setNAOs(elem.second.LMO);
      _fockDiisVector[elem.first]->setUnrestricted(!elem.second.isRestricted);
      _fockDiisVector[elem.first]->setSubspaceSize(_subspaceSize);
      _fockDiisVector[elem.first]->setOrthogonal(true);
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
      _fockDiisCounter[elem.first] = 0;
      _fockDiisStarted[elem.first] = false;
    }
  }
  for (auto const& elem : *_molecule) {
    _oldFockMatrix[elem.first].resize(elem.second.LMO);
    _oldFockMatrix[elem.first].alphaMatrix().setZero();
    _oldFockMatrix[elem.first].betaMatrix().setZero();
    _oldFockMatrix[elem.first].restrictedMatrix().setZero();
  }
}

auto ConvergenceAccelerator::reset() -> void {
  if (_mixer == Utils::scf_mixer_t::fock_diis) {
    for (auto const& elem : *_molecule) {
      _fockDiisVector[elem.first]->restart();
      _fockDiisCounter[elem.first] = 0;
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
      _fockDiisCounter[elem.first] = 0;
      _fockDiisStarted[elem.first] = false;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis) {
    for (auto const& elem : *_molecule) {
      _ediisVector[elem.first]->restart();
      _eDiisCounter[elem.first] = 0;
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis_diis) {
    for (auto const& elem : *_molecule) {
      _fockDiisVector[elem.first]->restart();
      _ediisVector[elem.first]->restart();
      _eDiisCounter[elem.first] = 0;
      _ediisStarted[elem.first] = false;
      _fockDiisCounter[elem.first] = 0;
      _fockDiisStarted[elem.first] = false;
    }
  }
  else {
    for (auto const& elem : *_molecule) {
      _oldFockMatrix[elem.first].alphaMatrix().setZero();
      _oldFockMatrix[elem.first].betaMatrix().setZero();
      _oldFockMatrix[elem.first].restrictedMatrix().setZero();
    }
  }
}

auto ConvergenceAccelerator::_check(Utils::ElementType type, double FockDiisThresh, double EDiisThresh, bool isRestricted)
    -> Utils::scf_mixer_t {
  if (_mixer == Utils::scf_mixer_t::fock_diis) {
    if (_checkType(type, isRestricted, FockDiisThresh)) {
      if (!_fockDiisStarted[type]) {
        print(type, "DIIS");
        _fockDiisStarted[type] = true;
      }
      return Utils::scf_mixer_t::fock_diis;
    }
    else {
      if (_fockDiisStarted[type]) {
        _fockDiisVector[type]->restart();
        _fockDiisStarted[type] = false;
      }
      return Utils::scf_mixer_t::none;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis) {
    if (_checkType(type, isRestricted, EDiisThresh)) {
      if (!_ediisStarted[type]) {
        print(type, "EDIIS");
        _ediisStarted[type] = true;
      }
      return Utils::scf_mixer_t::ediis;
    }
    else {
      if (_ediisStarted[type]) {
        _ediisVector[type]->restart();
        _ediisStarted[type] = false;
      }
      return Utils::scf_mixer_t::none;
    }
  }
  else if (_mixer == Utils::scf_mixer_t::ediis_diis) {
    if (_checkType(type, isRestricted, EDiisThresh)) {
      if (_checkType(type, isRestricted, FockDiisThresh)) {
        if (!_fockDiisStarted[type]) {
          print(type, "DIIS");
          _fockDiisStarted[type] = true;
        }
        return Utils::scf_mixer_t::fock_diis;
      }
      else if (_fockDiisStarted[type]) {
        _fockDiisVector[type]->restart();
        _fockDiisStarted[type] = false;
      }
      if (!_ediisStarted[type]) {
        print(type, "EDIIS");
        _ediisStarted[type] = true;
      }
      return Utils::scf_mixer_t::ediis;
    }
    else {
      if (_ediisStarted[type]) {
        _ediisVector[type]->restart();
        _ediisStarted[type] = false;
      }
      return Utils::scf_mixer_t::none;
    }
  }
  else {
    throw std::runtime_error("fock_simple and charge_simple not implemted.");
  }
}

auto ConvergenceAccelerator::update(const std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble>& errorMap,
                                    std::size_t restartAt) -> void {
  _errorMap = errorMap;
  _restartAt = restartAt;
  Utils::scf_mixer_t which;
  for (auto const& elem : *_molecule) {
    double Fock = 0;
    double E = 0;
    if (elem.first == Utils::ElementType::E) {
      Fock = _electronicFockDiisThresh;
      E = _electronicEDiisThresh;
    }
    else {
      Fock = _nuclearFockDiisThresh;
      E = _nuclearEDiisThresh;
    }
    which = _check(elem.first, Fock, E, elem.second.isRestricted);
    if (which == Utils::scf_mixer_t::fock_diis) {
      if (_fockDiisCounter[elem.first] >= _restartAt) {
        _fockDiisVector[elem.first]->restart();
        _fockDiisCounter[elem.first] = 0;
        _fockDiisStarted[elem.first] = false;
      }
      _fockDiisVector[elem.first]->addMatrices(_data->hartreeFockData->F_OAO[elem.first], _data->D_OAO[elem.first]);
      _data->hartreeFockData->F_OAO[elem.first] = _fockDiisVector[elem.first]->getMixedFockMatrix();
      _fockDiisCounter[elem.first] += 1;
    }
    else if (which == Utils::scf_mixer_t::ediis) {
      if (_eDiisCounter[elem.first] >= _restartAt) {
        _ediisVector[elem.first]->restart();
        _eDiisCounter[elem.first] = 0;
        _ediisStarted[elem.first] = false;
      }
      _ediisVector[elem.first]->addMatrices(_data->HartreeFockEnergy, _data->hartreeFockData->F_OAO[elem.first],
                                            _data->D_OAO[elem.first]);
      _data->hartreeFockData->F_OAO[elem.first] = _ediisVector[elem.first]->getMixedFockMatrix();
      _eDiisCounter[elem.first] += 1;
    }
  }
}
auto ConvergenceAccelerator::update(const std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble>& errorMap,
                                    Utils::ElementType type, std::size_t restartAt) -> void {
  _errorMap = errorMap;
  _restartAt = restartAt;
  Utils::scf_mixer_t which;
  double Fock = 0;
  double E = 0;
  if (type == Utils::ElementType::E) {
    Fock = _electronicFockDiisThresh;
    E = _electronicEDiisThresh;
  }
  else {
    Fock = _nuclearFockDiisThresh;
    E = _nuclearEDiisThresh;
  }
  which = _check(type, Fock, E, _molecule->at(type).isRestricted);
  if (which == Utils::scf_mixer_t::fock_diis) {
    if (_fockDiisCounter[type] >= _restartAt) {
      _fockDiisVector[type]->restart();
      _fockDiisCounter[type] = 0;
      _fockDiisStarted[type] = false;
    }
    _fockDiisVector[type]->addMatrices(_data->hartreeFockData->F_OAO[type], _data->D_OAO[type]);
    _data->hartreeFockData->F_OAO[type] = _fockDiisVector[type]->getMixedFockMatrix();
    _fockDiisCounter[type] += 1;
  }
  else if (which == Utils::scf_mixer_t::ediis) {
    if (_eDiisCounter[type] >= _restartAt) {
      _ediisVector[type]->restart();
      _eDiisCounter[type] = 0;
      _ediisStarted[type] = false;
    }
    _ediisVector[type]->addMatrices(_data->HartreeFockEnergy, _data->hartreeFockData->F_OAO[type], _data->D_OAO[type]);
    _data->hartreeFockData->F_OAO[type] = _ediisVector[type]->getMixedFockMatrix();
    _eDiisCounter[type] += 1;
  }
}
