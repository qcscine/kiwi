/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_ALTERNATINGSCF_H
#define KIWI_ALTERNATINGSCF_H

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/Scf.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

/**
 * @brief This is the alternating implementation of the Roothaan--Hall HF method which is only sensible for
 * nuclear-electronic calculations. With this method, the nuclear RH equations are converged, and afterwards the
 * electronic ones. Then the algortihm starts over until global convergence is reached.
 * @tparam Symmetry
 */
template<SymmetryType Symmetry>
class AlternatingScf : public Scf<Symmetry> {
  using Scf<Symmetry>::Scf;
  using Scf<Symmetry>::Scf::molecule_;
  using Scf<Symmetry>::Scf::data_;
  using Scf<Symmetry>::Scf::projector_;
  using Scf<Symmetry>::Scf::currentEnergy_;
  using Scf<Symmetry>::Scf::errorMap_;
  using Scf<Symmetry>::Scf::lastEnergy_;
  using Scf<Symmetry>::Scf::settings_;
  using Scf<Symmetry>::Scf::last_;
  using Scf<Symmetry>::Scf::convergenceAccelerator_;
  using Scf<Symmetry>::Scf::restartAt_;
  using Scf<Symmetry>::Scf::converged_;
  using Scf<Symmetry>::Scf::convergenceCheck;
  using Scf<Symmetry>::Scf::resetIncrementalCounter_;
  using Scf<Symmetry>::Scf::numFockRebuilds_;
  using Scf<Symmetry>::Scf::incremental_;
  using Scf<Symmetry>::Scf::direct_;
  using Scf<Symmetry>::Scf::verbose_;

 private:
  int _macroIter;
  int printCount_;

  auto runScf() -> void override;
  auto printHeader() -> void override;
  auto printSeparator() -> void;
  auto printIteration() -> void override;
  auto scf(Molecule::iterator elem) -> void;

  std::map<Utils::ElementType, std::size_t> _iterationMap;
  std::map<Utils::ElementType, bool> _firstMap;
};

template<SymmetryType Symmetry>
auto AlternatingScf<Symmetry>::scf(Molecule::iterator elem) -> void {
  std::size_t nestedCounter = 0;

  convergenceAccelerator_.reset();

  // Specialized optimization for nuclear-electronic ediis-diis:
  // In the first iteration, do not enable diis because the electronic guess is too bad.
  // if (_firstMap[elem->first] && elem->first != Utils::ElementType::E && settings_.accelerator ==
  // Utils::scf_mixer_t::ediis_diis) {
  //  // For nuclei: disable DIIS in the first iteration
  //  convergenceAccelerator_.setNuclearFockDiisThresh(0);
  //}
  // else {
  //  convergenceAccelerator_.setNuclearFockDiisThresh(settings_.nuclearDiisThresh);
  //}

  if (verbose_)
    printSeparator();
  while (true) {
    // Store the last density for incremental Fock matrix building
    projector_.setLastDensity(data_->D);

    projector_.updateDensity(elem->first);

    if (direct_ && !settings_.disableIncrementalFock) {
      if (!last_ && resetIncrementalCounter_ < settings_.resetIncremental) {
        incremental_ = true;
        ++resetIncrementalCounter_;
      }
      else {
        incremental_ = false;
        resetIncrementalCounter_ = 0;
        ++numFockRebuilds_;
        if (verbose_)
          std::cout << "                      *** Resetting incremental Fock-matrix formation. ***\n";
      }
    }
    if (_macroIter == 0) {
      incremental_ = false;
    }

    projector_.evaluateDensity(incremental_);

    projector_.updateFockMatrices(incremental_, numFockRebuilds_);

    errorMap_ = projector_.getErrorMap();
    projector_.evaluateEnergy();
    currentEnergy_ = projector_.getEnergy();

    data_->HartreeFockEnergy = currentEnergy_;

    ++_macroIter;
    ++_iterationMap[elem->first];
    ++nestedCounter;

    if (verbose_)
      printIteration();

    bool energyConverged = false;

    if (std::abs(lastEnergy_ - currentEnergy_) < settings_.energyThreshold) {
      energyConverged = true;
    }

    bool converged = true;

    if (elem->second.isRestricted) {
      if (errorMap_[elem->first].restricted > settings_.fockMatrixNormThresh) {
        converged = false;
      }
    }
    else {
      if (errorMap_[elem->first].alpha > settings_.fockMatrixNormThresh) {
        converged = false;
      }
      if (molecule_->at(elem->first).msVector[1] > 0) {
        if (errorMap_[elem->first].beta > settings_.fockMatrixNormThresh) {
          converged = false;
        }
      }
    }

    if (converged && energyConverged) {
      break;
    }

    // if (!last_ && !_firstMap[elem->first]) {
    if (!last_) {
      convergenceAccelerator_.update(errorMap_, elem->first, restartAt_);
      // convergenceAccelerator_.update(errorMap_, restartAt_);
    }

    lastEnergy_ = currentEnergy_;

    if (nestedCounter >= settings_.maxIterationsNested) {
      break;
    }
  }

  if (_firstMap[elem->first]) {
    _firstMap[elem->first] = false;
  }

  // nuclei first
  if (elem != molecule_->begin()) {
    scf(--elem);
  }

  // electrons first
  //++elem;
  // if (elem == molecule_->end()) {
  //  return;
  //}
  // scf(elem);
}

template<SymmetryType Symmetry>
auto AlternatingScf<Symmetry>::runScf() -> void {
  projector_.init(verbose_);
  lastEnergy_ = projector_.getEnergy();

  convergenceAccelerator_.init();
  convergenceAccelerator_.reset();

  const auto& maxIter = settings_.maxIterations;
  // const auto& eneThresh = settings_.energyThreshold;

  if (verbose_)
    printHeader();
  last_ = false;
  if (molecule_->size() > 1) {
    restartAt_ = 20;
  }
  else {
    restartAt_ = 16;
  }
  converged_ = false;
  _macroIter = 0;
  for (auto const& elem : *molecule_) {
    _iterationMap[elem.first] = 0;
    _firstMap[elem.first] = true;
  }

  // variables for incremental Fock matrix formation
  resetIncrementalCounter_ = 0;
  numFockRebuilds_ = 0;
  incremental_ = false;
  direct_ = data_->integralDirect;

  while (true) {
    // nuclei first
    scf(--molecule_->end());

    if (verbose_)
      printSeparator();

    /*
     *  SCF failed
     */
    if (_macroIter >= int(maxIter)) {
      if (verbose_)
        std::cout << "SCF did not converge after " << maxIter << " iterations." << std::endl;
      lastEnergy_ = currentEnergy_;
      break;
    }

    if (!last_) {
      /*
       *  SCF convergence criterium reached with DIIS
       */
      if (convergenceCheck()) {
        if (verbose_)
          std::cout << "                      *** Energy converged: Last iteration without DIIS ***" << std::endl;
        last_ = true;
      }
    }
    else {
      /*
       *  SCF convergence criterium reached without DIIS
       */
      if (convergenceCheck()) {
        if (verbose_)
          std::cout << "                      *** SCF converged ***" << std::endl;
        lastEnergy_ = currentEnergy_;
        break;
      }
      /*
       *  SCF convergence criterium **NOT** reached without DIIS
       */
      else {
        if (verbose_) {
          std::cout << "                      *** Convergence criterium was violated in last iteration ***" << std::endl;
          std::cout << "                      *** SCF will continue ***" << std::endl;
        }
        convergenceAccelerator_.reset();
        last_ = false;
      }
    }

    /*
     * DIIS optimization
     */
    // if (Symmetry!=SymmetryType::ParityAm) {
    // if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 10) {
    //  restartAt_ = 4;
    //}
    // if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 100) {
    //  restartAt_ = 8;
    //}
    // else {
    //  restartAt_ = 16;
    //}
    //}
    // if (!last_) {
    //  convergenceAccelerator_.update(errorMap_, restartAt_);
    //}

    lastEnergy_ = currentEnergy_;
  }

  projector_.finalize();
}
template<SymmetryType Symmetry>
auto AlternatingScf<Symmetry>::printSeparator() -> void {
  using std::cout;
  using std::endl;
  // cout << "Print count = " << printCount_ << endl;
  cout << std::string(15 + molecule_->size() * 10 + 21 + 19 + printCount_ * 18 + 6, '-') << endl;
}

template<SymmetryType Symmetry>
auto AlternatingScf<Symmetry>::printHeader() -> void {
  using std::cout;
  using std::endl;
  using std::setw;

  printCount_ = 0;
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      ++printCount_; // [F,D]
    }
    else {
      if (elem.second.msVector[0] > 0) {
        ++printCount_; // [F,D]
      }
      if (elem.second.msVector[1] > 0) {
        ++printCount_; // [F,D]
      }
    }
  }
  // cout << "Print count = " << printCount_ << endl;
  cout << std::right << endl;
  cout << std::string(15 + molecule_->size() * 10 + 21 + 19 + printCount_ * 18 + 6, '-') << endl;
  cout << setw(15) << "Iteration";
  for (auto const& elem : *molecule_) {
    cout << setw(8) << "Iter-" << setw(2) << Utils::ElementInfo::symbol(elem.first);
  }
  cout << setw(21) << "Energy/Ha" << setw(19) << "Delta E/Ha";
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      cout << setw(16) << "[F,D]:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
    }
    else {
      if (elem.second.msVector[0] > 0) {
        cout << setw(16) << "[F,D]-alpha:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
      }
      if (elem.second.msVector[1] > 0) {
        cout << setw(16) << "[F,D]-beta:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
      }
    }
  }
  cout << endl;
  cout << std::string(15 + molecule_->size() * 10 + 21 + 19 + printCount_ * 18 + 6, '-') << endl;
}

template<SymmetryType Symmetry>
auto AlternatingScf<Symmetry>::printIteration() -> void {
  using std::cout;
  using std::endl;
  using std::setprecision;
  using std::setw;

  cout << std::fixed << std::right;
  cout << setw(15) << _macroIter;
  for (auto const& elem : *molecule_) {
    cout << setw(10) << _iterationMap[elem.first];
  }
  cout << setw(21) << setprecision(12) << currentEnergy_;
  cout << setw(19) << setprecision(12) << currentEnergy_ - lastEnergy_;
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      cout << setw(18) << setprecision(10) << errorMap_[elem.first].restricted;
    }
    else {
      if (elem.second.msVector[0] > 0) {
        cout << setw(18) << setprecision(10) << errorMap_[elem.first].alpha;
      }
      if (elem.second.msVector[1] > 0) {
        cout << setw(18) << setprecision(10) << errorMap_[elem.first].beta;
      }
    }
  }
  cout << endl;
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_ALTERNATINGSCF_H
