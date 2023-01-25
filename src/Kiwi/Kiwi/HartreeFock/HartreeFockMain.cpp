/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kiwi/HartreeFock/InitialGuess/InitialGuess.h"
#include <Kiwi/HartreeFock/AlternatingScf.h>
#include <Kiwi/HartreeFock/BFGSScf.h>
#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/InitialGuess/Read.h>
#include <Kiwi/HartreeFock/SecondOrder.h>
#include <Kiwi/HartreeFock/SerialScf.h>
#include <Kiwi/HartreeFock/StabilityAnalysis.h>
#include <iostream>
#include <utility>

using namespace Scine;
using namespace Kiwi;

HartreeFockMain::HartreeFockMain(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data,
                                 HartreeFockSettings& settings, bool hasGuess, bool verbose)
  : molecule_(std::move(molecule)), data_(std::move(data)), settings_(settings), verbose_(verbose) {
  // Default settings. Can be overwritten later with `setTrahSettings`
  trahSettings_ = std::make_shared<TRAHSettings>();

  if (!hasGuess) {
    std::cout << "-----------------\n";
    std::cout << "  Initial guess\n";
    std::cout << "-----------------\n\n";
    makeGuess(molecule_, data_, settings_, verbose_);
  }

  std::cout << std::endl
            << std::left << std::setw(18) << "Particle type" << std::setw(18) << "Number of AOs" << std::setw(18)
            << "Number of MOs";
  std::cout << std::setw(5) << "N" << std::setw(9) << "N-alpha" << std::setw(9) << "N-beta" << std::endl;
  std::cout << std::string(78, '-') << std::endl;

  for (auto const& elem : *molecule_) {
    std::cout << std::right << std::setw(13) << Utils::ElementInfo::symbol(elem.first) << std::string(5, ' ')
              << std::setw(13) << elem.second.LAO << std::string(5, ' ') << std::setw(13) << elem.second.LMO
              << std::string(5, ' ');
    std::cout << std::setw(3) << elem.second.N << "  ";
    std::cout << std::setw(7) << elem.second.msVector[0] << "  ";
    std::cout << std::setw(7) << elem.second.msVector[1] << "  ";
    std::cout << std::endl;
  }
}

auto HartreeFockMain::makeGuess(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data,
                                HartreeFockSettings& settings, bool verbose) -> void {
  for (auto& elem : *molecule) {
    elem.second.occ.restricted = elem.second.N / 2;
    elem.second.occ.alpha = elem.second.msVector[0];
    elem.second.occ.beta = elem.second.msVector[1];
    elem.second.virt.restricted = elem.second.LMO - elem.second.N / 2;
    elem.second.virt.alpha = elem.second.LMO - elem.second.msVector[0];
    elem.second.virt.beta = elem.second.LMO - elem.second.msVector[1];
  }

  if (settings.guess == InitialGuessSCF::Read) {
    Read read(data, molecule);
    read.readOrbitals();
    std::cout << "Reading orbitals from files: " << molecule->getInputFileName() << std::endl;
  }
  else {
    for (auto const& elem : *molecule) {
      if (elem.first == Utils::ElementType::E) {
        InitialGuess guess(settings.guess, molecule, data, Utils::ElementType::E, verbose);
      }
      else {
        InitialGuess guess(settings.nuclearGuess, molecule, data, elem.first, verbose);
      }
    }
  }
}

auto HartreeFockMain::scf() -> void {
  if (settings_.guess == InitialGuessSCF::SAD || settings_.guess == InitialGuessSCF::SADNO ||
      (molecule_->size() > 1 && settings_.nuclearGuess == InitialGuessSCF::SND)) {
    Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.singleIteration();
  }

  if (settings_.guess == InitialGuessSCF::Read) {
    Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.singleIteration();
  }

  if (settings_.neScfType == NuclearElectronicSCF::Serial) {
    Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.run();
  }
  else if (settings_.neScfType == NuclearElectronicSCF::Alternating) {
    Kiwi::AlternatingScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.run();
  }
  else if (settings_.neScfType == NuclearElectronicSCF::TRAH) {
    if (settings_.guess == InitialGuessSCF::Hueckel) {
      Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
      scf.singleIteration();
    }
    Kiwi::SecondOrder<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.setTrahSettings(trahSettings_);
    scf.run();
  }
  else if (settings_.neScfType == NuclearElectronicSCF::BFGS) {
    if (settings_.guess == InitialGuessSCF::Hueckel) {
      Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
      scf.singleIteration();
    }
    Kiwi::BFGSScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
    scf.run();
  }

  if (settings_.stabilityAnalysis) {
    runStabilityAnalysis(molecule_, data_);
  }
}

void HartreeFockMain::setVerbose(bool verbose) {
  verbose_ = verbose;
}

auto HartreeFockMain::singleIteration() -> void {
  Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(molecule_, data_, settings_, ProjectionParameters(), verbose_);
  scf.singleIteration();
}

auto HartreeFockMain::writeOrbitals() -> void {
  Read read(this->data_, this->molecule_);
  read.writeOrbitals();
}
