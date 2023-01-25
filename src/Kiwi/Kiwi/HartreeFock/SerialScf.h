/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SERIALSCF_H
#define KIWI_SERIALSCF_H

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/Scf.h>
#include <Kiwi/HartreeFock/SimultaneousDiis/SimultaneousDiis.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

/**
 * @class SerialScf
 * @brief This is the serial implementation of the Roothaan-Hall method. This is the only way for BO methods.
 * For pre-BO, it means that nuclear and electronic Fock matrices are diagonalized in a serial way and the density
 * matrices are simultaneously updated.
 * @tparam Symmetry
 */
template<SymmetryType Symmetry>
class SerialScf : public Scf<Symmetry> {
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

  auto runScf() -> void override;
  auto printHeader() -> void override;
  auto printIteration() -> void override;
};

template<SymmetryType Symmetry>
auto SerialScf<Symmetry>::runScf() -> void {
  projector_.init(verbose_);
  lastEnergy_ = projector_.getEnergy();

  const auto& maxIter = settings_.maxIterations;
  const auto& eneThresh = settings_.energyThreshold;

  std::shared_ptr<SimultaneousDiis> ptr_diis;

  if (verbose_)
    printHeader();
  last_ = false;
  restartAt_ = 16;
  converged_ = false;
  _macroIter = 0;

  // variables for incremental Fock matrix formation
  resetIncrementalCounter_ = 0;
  numFockRebuilds_ = 0;
  incremental_ = false;
  direct_ = data_->integralDirect;

  bool useSimultaneousDiis = true;

  if (molecule_->size() == 1) {
    useSimultaneousDiis = false;
  }

  if (!useSimultaneousDiis) {
    convergenceAccelerator_.init();
    convergenceAccelerator_.reset();
  }
  else {
    std::cout << "                      *** Enabling Simultaneous DIIS. ***\n";
    ptr_diis = std::make_shared<SimultaneousDiis>(molecule_, data_);
    ptr_diis->setSubspaceSize(16);
    ptr_diis->setOrthonormal(true);
  }

  while (true) {
    // Store the last density for incremental Fock matrix building
    projector_.setLastDensity(data_->D);

    for (auto const& elem : *molecule_) {
      projector_.updateDensity(elem.first);
    }

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
    projector_.evaluateEnergy();
    currentEnergy_ = projector_.getEnergy();

    ++_macroIter;
    errorMap_ = projector_.getErrorMap();
    if (verbose_)
      printIteration();

    /*
     *  SCF failed
     */
    if (_macroIter == int(maxIter)) {
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
        converged_ = true;
        lastEnergy_ = currentEnergy_;
        break;
      }
      /*
       *  SCF convergence criterium **NOT** reached without DIIS
       */
      if (verbose_) {
        std::cout << "                      *** Convergence criterium was violated in last iteration ***" << std::endl;
        std::cout << "                      *** SCF will continue ***" << std::endl;
      }
      if (!useSimultaneousDiis) {
        convergenceAccelerator_.reset();
      }
      else {
        ptr_diis->restart();
      }
      last_ = false;
    }

    /*
     * DIIS optimization
     */
    // if (Symmetry!=SymmetryType::ParityAm) {
    if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 100) {
      restartAt_ = 4;
    }
    else if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 1000) {
      restartAt_ = 8;
    }
    else {
      restartAt_ = 16;
    }
    lastEnergy_ = currentEnergy_;

    if (!last_) {
      if (!useSimultaneousDiis) {
        convergenceAccelerator_.update(errorMap_, restartAt_);
      }
      else {
        ptr_diis->addMatrices();
        ptr_diis->setMixedFockMatrix();
      }
    }
  }

  projector_.finalize();
}

template<SymmetryType Symmetry>
auto SerialScf<Symmetry>::printHeader() -> void {
  using std::cout;
  using std::endl;
  using std::setw;

  int count = 0;
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      ++count; // [F,D]
    }
    else {
      if (elem.second.msVector[0] > 0) {
        ++count; // [F,D]
      }
      if (elem.second.msVector[1] > 0) {
        ++count; // [F,D]
      }
    }
  }
  cout << std::right << endl;
  cout << std::string(15 + 21 + 19 + count * 18 + 6, '-') << endl;
  cout << setw(15) << "Iteration";
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
  cout << std::string(15 + 21 + 19 + count * 18 + 6, '-') << endl;
}

template<SymmetryType Symmetry>
auto SerialScf<Symmetry>::printIteration() -> void {
  using std::cout;
  using std::endl;
  using std::setprecision;
  using std::setw;

  cout << std::fixed << std::right;
  cout << setw(15) << _macroIter;
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

#endif // KIWI_SERIALSCF_H
