/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/DataStructures/IntegralSpecifier.h>

using namespace Scine;
using namespace Kiwi;

Data::Data(const std::shared_ptr<Molecule>& mol) : molecule(mol), uniquePairs(generateUniquePairs(mol)) {
  hartreeFockData = std::make_shared<HartreeFockData>();
}

auto Data::generateUniquePairs(const std::shared_ptr<Molecule>& molecule) -> std::vector<ElementPair> {
  std::vector<ElementPair> ret;

  std::vector<Utils::ElementType> elements;

  for (auto const& elem : *molecule) {
    elements.push_back(elem.first);
  }

  for (auto i = 0UL; i < elements.size(); ++i) {
    for (auto j = 0UL; j <= i; ++j) {
      ret.push_back(getUnique(elements[i], elements[j]));
    }
  }

  return ret;
}

void Data::oneBodyIntegrals(const bool printTimings) {
  //
  // One-Body Integrals:
  //
  if (printTimings) {
    auto& clock = Clock::getInstance();
    std::cout << std::left << std::setw(30) << "Building one-body integrals" << std::right << std::setw(20);
    clock.time("1body");
  }
  for (auto const& typeData : *molecule) {
    Utils::Integrals::IntegralSpecifier specifier;
    specifier.typeVector = {typeData.second.typeInfo};
    // Overlap
    {
      specifier.op = Utils::Integrals::Operator::Overlap;
      auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, typeData.second.basisSet, typeData.second.basisSet);
      S[typeData.first] = resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
    }
    // H core
    // Pure pre-BO
    if (molecule->isPurePreBO()) {
      specifier.op = Utils::Integrals::Operator::KineticCOM;
      specifier.totalMass = molecule->getTotalMass();
      auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, typeData.second.basisSet, typeData.second.basisSet);
      H[typeData.first] = resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
    }
    else {
      assert(molecule->hasPointCharges());
      {
        specifier.op = Utils::Integrals::Operator::Kinetic;
        auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, typeData.second.basisSet, typeData.second.basisSet);
        H[typeData.first] = resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
      }
      {
        specifier.op = Utils::Integrals::Operator::PointCharges;
        specifier.atoms = molecule->getPointCharges();
        auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, typeData.second.basisSet, typeData.second.basisSet);
        H[typeData.first] += resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
      }
    }
  }
  if (printTimings) {
    auto& clock = Clock::getInstance();
    clock.time("1body");
  }
}

void Data::twoBodyIntegrals(const bool printTimings) {
  if (integralDirect) {
    return;
  }

  if (printTimings) {
    auto& clock = Clock::getInstance();
    std::cout << std::left << std::setw(30) << "Building two-body integrals" << std::right << std::setw(20);
    clock.time("2body");
  }
  Utils::Integrals::IntegralSpecifier specifier;
  for (auto const& pair : uniquePairs) {
    auto tp1 = molecule->at(pair.first);
    auto tp2 = molecule->at(pair.second);
    specifier.typeVector = {tp1.typeInfo, tp2.typeInfo};
    if (molecule->isPurePreBO()) {
      specifier.op = Utils::Integrals::Operator::CoulombCOM;
      specifier.totalMass = molecule->getTotalMass();
      auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, tp1.basisSet, tp2.basisSet);
      Coulomb[pair] = std::move(resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}]);
    }
    else {
      specifier.op = Utils::Integrals::Operator::Coulomb;
      auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, tp1.basisSet, tp2.basisSet);
      Coulomb[pair] = std::move(resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}]);
    }
  }
  if (printTimings) {
    auto& clock = Clock::getInstance();
    clock.time("2body");
  }
}

void Data::makeExchange(const bool printTimings) {
  if (integralDirect) {
    return;
  }

  if (printTimings) {
    auto& clock = Clock::getInstance();
    std::cout << std::left << std::setw(30) << "Transposing two-body integrals" << std::right << std::setw(20);
    clock.time("2bodyT");
  }

  for (auto const& typeData : *molecule) {
    auto const& type = typeData.first;
    Matrix const& coulomb = Coulomb[getUnique(type, type)];
    Exchange[type] = Matrix::Zero(coulomb.rows(), coulomb.cols());
    Matrix& exchange = Exchange[type];

    assert(coulomb.rows() == coulomb.cols());
    auto dim = static_cast<int>(std::sqrt(coulomb.rows()));

    for (auto p = 0; p < dim; ++p) {
      for (auto q = 0; q < dim; ++q) {
        exchange.block(p * dim, q * dim, dim, dim) = coulomb.block(p * dim, q * dim, dim, dim).transpose();
      }
    }
  }

  if (printTimings) {
    auto& clock = Clock::getInstance();
    clock.time("2bodyT");
  }
}

void Data::makeLoewdinOrtho(double loewdinThresh, const bool verbose) {
  for (auto const& elem : *molecule) {
    this->X[elem.first] = Loewdin(this->S[elem.first], loewdinThresh, verbose);
    molecule->at(elem.first).LMO = this->X[elem.first].cols();
  }
}
