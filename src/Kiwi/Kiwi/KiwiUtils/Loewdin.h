/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_LOEWDIN_H
#define KIWI_LOEWDIN_H

#include <Eigen/Eigenvalues>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {

/**
 * @brief Performs the canonical Loewdin orthonormalization.
 *
 * General procedure for the canonical Loewdin orthonormalization:
 *
 *              X = U * z^(-1/2)
 *
 * Important -> keep track of the dimensions:
 *
 *         ini x fin  = ini x fin * fin x fin
 *
 *  U is the matrix that diagonalizes S, and s is the diagonal matrix, where all orbitals, below a given
 *  threshold were discarded.
 *
 *
 *
 * The orthonormalization is done iteratively, such that the least amount of orbitals will be discarded in case of
 * a nearly singular overlap matrix.
 *
 * @param overlap
 * @param loewdinThresh
 * @param verbose
 * @return transformation matrix X
 */
inline auto Loewdin(const Eigen::MatrixXd& overlap, double loewdinThresh = 10e-10, const bool verbose = true) -> Eigen::MatrixXd {
  const std::size_t initialDim = overlap.rows(); // Dimension of the uncontracted tensors
  auto rank = 0UL;
  double thresh = 1.0e-16;   // Starting Loewdin threshold.
  double loewdinError = NAN; // Eigenvalues of the overlap matrix below this threshold will be discarded.

  Eigen::MatrixXd X;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap.selfadjointView<Eigen::Lower>());
  std::size_t nZeroEigenvalues = 0;

  if (verbose) {
    std::cout << "Initial dimension: " << initialDim << std::endl;
    std::cout << "Desired maxixmum Loewdin error: " << std::scientific << std::setprecision(2) << loewdinThresh << std::endl;
    std::cout << std::string(30, '-') << std::endl;
    std::cout << std::setw(10) << std::left << "Iteration" << std::string(4, ' ') << std::setw(12) << "Error" << std::endl;
    std::cout << std::string(30, '-') << std::endl;
  }
  auto iteration = 0;
  while (true) {
    ++iteration;
    nZeroEigenvalues = 0;
    // discarding eigenvalues if less than Loewdin threshold
    for (auto i = 0L; i < es.eigenvalues().size(); i++) {
      if (std::abs(es.eigenvalues()[i]) < thresh) {
        ++nZeroEigenvalues;
      }
    }
    rank = es.eigenvalues().size() - nZeroEigenvalues;

    // Take only eigenpairs for which the eigenvalue is != 0
    Eigen::VectorXd nonSingularEigenvalues = es.eigenvalues().tail(rank);
    Eigen::MatrixXd nonSingularEigenvectors = es.eigenvectors().rightCols(rank);

    X = nonSingularEigenvectors * nonSingularEigenvalues.array().sqrt().inverse().matrix().asDiagonal();

    Eigen::MatrixXd TEST = X.transpose() * overlap.selfadjointView<Eigen::Lower>() * X;
    // Generate TEST matrix
    for (auto i = 0UL; i < rank; i++) {
      // Diminish by one -> Matrix should be zero.
      TEST(i, i) -= 1.0;
    }

    loewdinError = TEST.norm(); // Search error in transformation of the overlap

    if (verbose) {
      std::cout << std::setw(10) << std::right << iteration++ << std::string(4, ' ') << std::setw(12) << std::scientific
                << std::setprecision(4) << loewdinError << std::endl;
    }
    // Success: contracting hamiltonian terms
    if (loewdinError < loewdinThresh) {
      if (verbose) {
        std::cout << std::string(30, '-') << std::endl;
        std::cout << "Loewdin orthonormalization done." << std::endl;
      }
      break;
    }
    // The 'thresh' must be increased slowly so as to discard as few orbitals, as possible.
    thresh *= 5;
  }

  if (verbose) {
    if (rank == initialDim) {
      std::cout << "No orbitals were discarded." << std::endl;
    }
    else {
      if ((initialDim - rank) == 1)
        std::cout << "One orbital was discarded." << std::endl;
      else
        std::cout << (initialDim - rank) << " orbitals were discarded." << std::endl;
      std::cout << "Final dimension: " << rank << std::endl;
    }
  }
  return X;
}

/**
 * @brief Performs the canonical Loewdin orthonormalization.
 *
 * General procedure for the canonical Loewdin orthonormalization:
 *
 *              X = U * z^(-1/2)
 *
 * Important -> keep track of the dimensions:
 *
 *         ini x fin  = ini x fin * fin x fin
 *
 *  U is the matrix that diagonalizes S, and s is the diagonal matrix, where all orbitals, below a given
 *  threshold were discarded.
 *
 *
 *
 * The orthonormalization is done in a single step.
 *
 * @param overlap
 * @param loewdinThresh
 * @param verbose
 * @return transformation matrix X
 */
inline auto SingleStepLoewdin(const Eigen::MatrixXd& S, const double epsilon = 1.0e-10, const bool verbose = false)
    -> Eigen::MatrixXd {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S.selfadjointView<Eigen::Lower>());

  auto rank = es.eigenvalues().size();
  auto nZero = 0;
  for (auto i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) > epsilon) {
      nZero = i;
      break;
    }
  }
  rank -= nZero;
  Eigen::MatrixXd nonSingularEigenvalues = es.eigenvalues().tail(rank).array().sqrt().inverse().matrix().asDiagonal();
  Eigen::MatrixXd nonSingularEigenvectors = es.eigenvectors().rightCols(nonSingularEigenvalues.rows());

  if (verbose) {
    auto initialDim = S.cols();
    assert(rank == nonSingularEigenvalues.cols());
    if (rank == initialDim) {
      std::cout << "No orbitals were discarded." << std::endl;
    }
    else {
      if ((initialDim - rank) == 1) {
        std::cout << "One orbital was discarded." << std::endl;
      }
      else {
        std::cout << (initialDim - rank) << " orbitals were discarded." << std::endl;
      }
      std::cout << "Final dimension: " << rank << std::endl;
    }
  }
  return nonSingularEigenvectors * nonSingularEigenvalues;
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_LOEWDIN_H
