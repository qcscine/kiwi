/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_BFGSSCF_H
#define KIWI_BFGSSCF_H

#include <Kiwi/HartreeFock/BFGS/Optimizer.h>
#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/Scf.h>
#include <Kiwi/KiwiOpt/MoreThuente.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

/**
 * @class BFGS
 * @brief This is the serial implementation of the Roothaan-Hall method. This is the only way for BO methods.
 * For pre-BO, it means that nuclear and electronic Fock matrices are diagonalized in a serial way and the density
 * matrices are simultaneously updated.
 * @tparam Symmetry
 */
template<SymmetryType Symmetry>
class BFGSScf : public Scf<Symmetry> {
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
auto BFGSScf<Symmetry>::runScf() -> void {
  // projector_.init(verbose_);
  // lastEnergy_ = projector_.getEnergy();

  convergenceAccelerator_.init();
  convergenceAccelerator_.reset();

  auto interface = std::make_shared<Kiwi::BFGS::Interface<Symmetry>>(molecule_, data_, projector_.getProjectionParameters());
  auto optimizer = std::make_shared<Kiwi::BFGS::Optimizer<Symmetry>>(interface);
  Kiwi::Optimization::MoreThuente<Kiwi::BFGS::Optimizer<Symmetry>> moreThuente(*optimizer);
  moreThuente.setMaxIterations(5);

  // double currentThresh = 0.01;
  std::size_t total_micro_iterations = 0;
  double alpha_0 = 1.;

  const auto& maxIter = settings_.maxIterations;
  // const auto& eneThresh = settings_.energyThreshold;

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

  while (true) {
    // p_i = - B_i g_i
    optimizer->evaluateDirection();

    // optimize:
    // f(a) := E( Exp ( -a*p_i) )
    // moreThuente.setTol(currentThresh);
    moreThuente(alpha_0, optimizer->getZeroValue(), optimizer->getZeroDerivative());
    // a = a_opt
    // g_i+1 = grad E ( Exp (-a_opt * p_i) )
    // D_i+1 =  Exp (-a_opt * p_i)

    // Build B_i+1 from B_i, g_i+1, g_i, a*p_i
    optimizer->applyUpdate();

    total_micro_iterations += moreThuente.getIterations();

    currentEnergy_ = optimizer->getLineSearchValue();

    ++_macroIter;
    errorMap_ = interface->getPtrData()->projector.getErrorMap();
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
      convergenceAccelerator_.reset();
      last_ = false;
    }

    ///*
    // * DIIS optimization
    // */
    //// if (Symmetry!=SymmetryType::ParityAm) {
    // if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 100) {
    //  restartAt_ = 4;
    //}
    // else if (std::abs(lastEnergy_ - currentEnergy_) < eneThresh * 1000) {
    //  restartAt_ = 8;
    //}
    // else {
    //  restartAt_ = 16;
    //}
    lastEnergy_ = currentEnergy_;

    // if (!last_) {
    //  convergenceAccelerator_.update(errorMap_, restartAt_);
    //}
  }

  interface->getPtrData()->projector.finalize();
  // projector_ = ptrData->projector;
  // projector_.finalize();
}

template<SymmetryType Symmetry>
auto BFGSScf<Symmetry>::printHeader() -> void {
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
auto BFGSScf<Symmetry>::printIteration() -> void {
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

#endif // KIWI_BFGSSCF_H
