/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/NaturalOrbitals.h>
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Kiwi {

auto NaturalOrbitals::generateNaturalOrbitals(bool prune) -> void {
  std::vector<std::string> fileNames; // Store file names here
  // Trim leading and final spaces in the string.
  for (const auto& input : elementType_LineMap_) {
    auto type = input.first;
    std::string line = input.second;

    line.erase(line.begin(), std::find_if(line.begin(), line.end(), [&](int ch) { return !std::isspace(ch); }));
    line.erase(std::find_if(line.rbegin(), line.rend(), [&](int ch) { return !std::isspace(ch); }).base(), line.end());
    // Split the string
    boost::algorithm::split(fileNames, line, boost::is_any_of(" "), boost::token_compress_on);
    std::vector<Eigen::MatrixXd> rdms;
    for (auto const& fname : fileNames) {
      rdms.push_back(csvFileReader(fname));
    }
    if (!molecule_->at(type).isRestricted) {
      if (rdms.size() == 1 && (molecule_->at(type).msVector[1] > 0)) {
        throw std::runtime_error("Provide RDM file for spin alpha and beta!");
      }
    }

    if (molecule_->at(type).isRestricted) {
      data_->C.at(type).restrictedMatrix() = data_->C.at(type).restrictedMatrix() * generateOrbitals(rdms[0], prune);
    }
    else {
      data_->C.at(type).alphaMatrix() = data_->C.at(type).alphaMatrix() * generateOrbitals(rdms[0], prune);
      if (molecule_->at(type).msVector[1] > 0) {
        data_->C.at(type).betaMatrix() = data_->C.at(type).betaMatrix() * generateOrbitals(rdms[1], prune);
      }
      else {
        data_->C.at(type).betaMatrix().resizeLike(data_->C.at(type).alphaMatrix());
        data_->C.at(type).betaMatrix().setZero();
      }
    }
  }

  for (const auto& elem : *molecule_) {
    if (elem.second.isRestricted) {
      molecule_->at(elem.first).LMO = data_->C.at(elem.first).restrictedMatrix().cols();
    }
    else {
      molecule_->at(elem.first).LMO = data_->C.at(elem.first).alphaMatrix().cols();
      if (molecule_->at(elem.first).msVector[1] > 0) {
        // The alpha and beta rdms can obviously have different eigenvalues.
        // However, we would not want to have alpha and beta orbitals of a different size.
        // Therefore, we adjust the orbitals such that they are of equal size where we insert zeros in the
        // NOs that do not contribute.
        if (data_->C.at(elem.first).alphaMatrix().cols() > data_->C.at(elem.first).betaMatrix().cols()) {
          int sizeBeta = data_->C.at(elem.first).betaMatrix().cols();
          int sizeAlpha = data_->C.at(elem.first).alphaMatrix().cols();
          data_->C.at(elem.first).betaMatrix().conservativeResizeLike(data_->C.at(elem.first).alphaMatrix());
          for (int i = sizeBeta; i < sizeAlpha; ++i) {
            data_->C.at(elem.first).betaMatrix().col(i).setZero();
          }
        }
        if (data_->C.at(elem.first).betaMatrix().cols() > data_->C.at(elem.first).alphaMatrix().cols()) {
          int sizeBeta = data_->C.at(elem.first).alphaMatrix().cols();
          int sizeAlpha = data_->C.at(elem.first).betaMatrix().cols();
          data_->C.at(elem.first).alphaMatrix().conservativeResizeLike(data_->C.at(elem.first).betaMatrix());
          for (int i = sizeBeta; i < sizeAlpha; ++i) {
            data_->C.at(elem.first).alphaMatrix().col(i).setZero();
          }
        }
      }
    }
  }

  std::cout << "\n-------------------------------\n";
  std::cout << "  Natural orbital composition\n";
  std::cout << "-------------------------------\n";

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

auto NaturalOrbitals::generateOrbitals(const Eigen::MatrixXd& rdm, bool prune) -> Eigen::MatrixXd {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(rdm);
  Eigen::VectorXd eigvalues = es.eigenvalues().reverse();
  Eigen::MatrixXd eigvecs = es.eigenvectors().rowwise().reverse();
  Eigen::MatrixXd ret;

  int rows = eigvecs.rows();

  if (prune) {
    int cols = 0;
    for (auto i = 0; i < eigvalues.size(); ++i) {
      if (std::abs(eigvalues(i)) > 1e-16) {
        ++cols;
        ret.conservativeResize(rows, cols);
        ret.col(cols - 1) = eigvecs.col(i);
        std::cout << cols << std::endl;
      }
    }
  }
  else {
    ret = eigvecs;
  }

  std::cout << std::endl;

  return ret;
}

} // namespace Kiwi
} // namespace Scine
