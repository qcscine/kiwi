/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiOpt/Bisection.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <Kiwi/KiwiOpt/DirectConjugateGradient.h>
#include <Kiwi/KiwiOpt/LineSearchInterface.h>
#include <Kiwi/KiwiOpt/MoreThuente.h>
#include <Kiwi/KiwiOpt/ObjectiveFunction.h>
#include <Kiwi/KiwiOpt/Optimization.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <Kiwi/KiwiOpt/PreconditionedDirectConjugateGradient.h>
#include <gmock/gmock.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <utility>

namespace Scine {

using namespace testing;

class OptimizationTest : public Test {};

class Rosenbrock : public Kiwi::Optimization::ObjectiveFunction {
 public:
  Rosenbrock() {
    gradient.resize(2);
    parameters.resize(2);
    Hessian.resize(2, 2);
  }

  // Rosenbrock(double param_a, double param_b) : a(param_a), b(param_b) {
  //  gradient.resize(2);
  //  parameters.resize(2);
  //}

  auto evaluate(const Eigen::VectorXd& x) -> void final {
    parameters = x;

    value = std::pow((1 - parameters[0]), 2) + 100 * std::pow(parameters[1] - parameters[0] * parameters[0], 2);

    gradient[0] = 2 * (parameters[0] - 1) + 400 * parameters[0] * parameters[0] * parameters[0] -
                  400 * parameters[0] * parameters[1];
    gradient[1] = 200 * (parameters[1] - parameters[0] * parameters[0]);

    Hessian(0, 0) = 1200 * parameters[0] * parameters[0] - 400 * parameters[1] + 2;
    Hessian(0, 1) = -400 * parameters[0];
    Hessian(1, 0) = -400 * parameters[0];
    Hessian(1, 1) = 200;
  }
};

class Wood : public Kiwi::Optimization::ObjectiveFunction {
 public:
  Wood() {
    gradient.resize(4);
    parameters.resize(4);
  }

  auto evaluate(const Eigen::VectorXd& x) -> void final {
    double f0 = 10 * (x[1] - pow(x[0], 2));
    double f1 = 1 - x[0];
    double f2 = std::sqrt(90) * (x[3] - pow(x[2], 2));
    double f3 = 1 - x[2];
    double f4 = std::sqrt(10) * (x[1] + x[3] - 2);
    double f5 = (x[1] - x[3]) / std::sqrt(10.);
    value = f0 * f0 + f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4 + f5 * f5;

    double dx0f = -400 * (x[1] - pow(x[0], 2)) * x[0] - 2 * (1 - x[0]);
    double dx1f = 200 * (x[1] - pow(x[0], 2)) + 20 * (x[1] + x[3] - 2) + 2 / 10.0 * (x[1] - x[3]);
    double dx2f = -4 * 90 * (x[3] - pow(x[2], 2)) * x[2] - 2 * (1 - x[2]);
    double dx3f = 2 * 90 * (x[3] - pow(x[2], 2)) + 20 * (x[1] + x[3] - 2) - 2 / 10.0 * (x[1] - x[3]);

    parameters = x;
    gradient.setZero();
    gradient << dx0f, dx1f, dx2f, dx3f;
  }
};

class BiggsEXP6 : public Kiwi::Optimization::ObjectiveFunction {
 public:
  BiggsEXP6() {
    gradient.resize(6);
    parameters.resize(6);
  }

  auto evaluate(const Eigen::VectorXd& x) -> void final {
    value = 0;
    gradient.setZero();

    for (int i = 1; i <= 13; i++) {
      double z = i / 10.;
      double y = std::exp(-z) - 5 * std::exp(-10 * z) + 3 * std::exp(-4 * z);
      double t = x[2] * std::exp(-x[0] * z) - x[3] * std::exp(-x[1] * z) + x[5] * std::exp(-x[4] * z) - y;

      double dfdx0 = -z * x[2] * std::exp(-x[0] * z);
      double dfdx1 = z * x[3] * std::exp(-x[1] * z);
      double dfdx2 = std::exp(-x[0] * z);
      double dfdx3 = -std::exp(-x[1] * z);
      double dfdx4 = -z * x[5] * std::exp(-x[4] * z);
      double dfdx5 = std::exp(-x[4] * z);

      value += t * t;
      gradient[0] += 2 * t * dfdx0;
      gradient[1] += 2 * t * dfdx1;
      gradient[2] += 2 * t * dfdx2;
      gradient[3] += 2 * t * dfdx3;
      gradient[4] += 2 * t * dfdx4;
      gradient[5] += 2 * t * dfdx5;
    }

    parameters = x;
  }
};

TEST_F(OptimizationTest, ConjugateGradient) {
  int dim = 100;
  int maxIt = 1000;
  double sparsity = 1.;
  Eigen::MatrixXd A;
  A = Eigen::MatrixXd::Random(dim, dim);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i, j) *= sparsity;
      A(j, i) = A(i, j);
    }
  }

  class Evaluator {
    Eigen::MatrixXd A_;

   public:
    Evaluator(Eigen::MatrixXd A) : A_(std::move(A)) {
    }
    auto evaluate(const Eigen::VectorXd& x) -> Eigen::VectorXd {
      return A_ * x;
    }
  };

  using namespace std::chrono;

  Eigen::VectorXd b = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd eigenSolution(dim);
  Eigen::VectorXd x0 = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd bcp = b;
  Eigen::VectorXd x0cp = x0;

  auto start1 = high_resolution_clock::now();
  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> eigenCg;
  eigenCg.setMaxIterations(maxIt);
  eigenCg.compute(A);
  Eigen::VectorXd solCg = eigenCg.solve(b);
  auto stop1 = high_resolution_clock::now();
  auto duration1 = duration_cast<microseconds>(stop1 - start1);

  std::cout << "Eigen solution\n";
  std::cout << "#iterations:     " << eigenCg.iterations() << std::endl;
  std::cout << "estimated error: " << eigenCg.error() << std::endl;
  std::cout << "duration: " << duration1.count() << std::endl;

  auto start2 = high_resolution_clock::now();
  Kiwi::Optimization::DirectConjugateGradient<double> cg(std::move(x0));
  cg.setMaxIterations(maxIt);
  cg.updateImage(std::move(b));
  Evaluator evaluator(A);
  cg.optimize(evaluator);
  auto stop2 = high_resolution_clock::now();
  auto duration2 = duration_cast<microseconds>(stop2 - start2);
  std::cout << "Kiwi CG\n";
  std::cout << "#iterations:     " << cg.getIterations() << std::endl;
  std::cout << "estimated error: " << cg.getError() << std::endl;
  std::cout << "duration: " << duration2.count() << std::endl;
  Eigen::VectorXd cgSol = cg.getSolution();

  auto start3 = high_resolution_clock::now();
  Kiwi::Optimization::PreconditionedDirectConjugateGradient<double> pcg(std::move(x0cp));
  Eigen::MatrixXd precond;
  precond.resizeLike(A);
  precond.setZero();
  for (int i = 0; i < dim - 1; ++i) {
    precond(i, i) = A(i, i);
    precond(i, i + 1) = A(i, i + 1);
    precond(i + 1, i) = A(i + 1, i);
  }
  precond(dim - 1, dim - 1) = A(dim - 1, dim - 1);
  precond = precond.inverse();
  pcg.setMaxIterations(maxIt);
  pcg.updateImage(std::move(bcp));
  pcg.optimize(evaluator, precond);
  auto stop3 = high_resolution_clock::now();
  auto duration3 = duration_cast<microseconds>(stop3 - start3);
  std::cout << "Kiwi PCG\n";
  std::cout << "#iterations:     " << pcg.getIterations() << std::endl;
  std::cout << "estimated error: " << pcg.getError() << std::endl;
  std::cout << "duration: " << duration3.count() << std::endl;
  Eigen::VectorXd pcgSol = pcg.getSolution();

  EXPECT_THAT((cgSol - solCg).norm(), DoubleNear(0.0, 1e-10));
  EXPECT_THAT((pcgSol - solCg).norm(), DoubleNear(0.0, 1e-10));
}

TEST_F(OptimizationTest, Davidson) {
  class Evaluator {
    Eigen::MatrixXd A_;

   public:
    Evaluator(Eigen::MatrixXd A) : A_(std::move(A)) {
    }
    auto evaluate(const Eigen::MatrixXd& P) -> Eigen::MatrixXd {
      return A_ * P;
    }
  };

  int dim = 50;
  int maxIt = 1000;
  double sparsity = 0.001;
  Eigen::MatrixXd A;
  A = Eigen::MatrixXd::Random(dim, dim);
  // A = (A + A.transpose()).eval();

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      A(i, j) *= sparsity;
      A(j, i) = A(i, j);
    }
  }

  int numberOfEigenpairs = 2;
  int maxSubspaceDim = 10;

  using namespace std::chrono;

  Kiwi::Optimization::Davidson davidson(numberOfEigenpairs, maxSubspaceDim, dim);

  davidson.setMaxIterations(maxIt);

  davidson.setDiagonalPreconditioner(A.diagonal().asDiagonal());

  auto start2 = high_resolution_clock::now();

  Evaluator evaluator(A);
  davidson.compute(evaluator, true);
  auto stop2 = high_resolution_clock::now();
  auto duration2 = duration_cast<seconds>(stop2 - start2);
  std::cout << "Kiwi\n";
  std::cout << "#iterations:     " << davidson.getIterations() << std::endl;
  std::cout << "estimated error: " << davidson.getError() << std::endl;
  std::cout << "duration: " << duration2.count() << std::endl;

  Eigen::MatrixXd eigvecs = davidson.getEigvecs();
  Eigen::VectorXd eigvals = davidson.getEigvals();

  std::cout << "Eigenvalues\n" << eigvals.block(0, 0, numberOfEigenpairs, 1) << std::endl;

  auto start1 = high_resolution_clock::now();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
  eig.compute(A);
  auto stop1 = high_resolution_clock::now();
  auto duration1 = duration_cast<seconds>(stop1 - start1);
  std::cout << "Eigen\n";
  std::cout << "duration: " << duration1.count() << std::endl;

  Eigen::MatrixXd refVec = eig.eigenvectors();
  Eigen::VectorXd refVal = eig.eigenvalues();

  std::cout << "Reference eigenvalues\n" << refVal.block(0, 0, numberOfEigenpairs, 1) << std::endl;

  EXPECT_THAT((eigvals.block(0, 0, numberOfEigenpairs, 1) - refVal.block(0, 0, numberOfEigenpairs, 1)).norm(),
              DoubleNear(0.0, 1e-10));
}

TEST_F(OptimizationTest, Bisection) {
  double reference = 1.52137971;

  class Evaluator {
   public:
    Evaluator() = default;

    auto evaluate(double x) -> double {
      return std::pow(x, 3) - x - 2;
    }
  };

  Evaluator eval;

  Kiwi::Optimization::Bisection bisection;

  bisection.setLowerBound(1);
  bisection.setUpperBound(2);
  bisection.setTolerance(1e-7);

  auto result = bisection.compute(eval);

  if (bisection.wasSuccessful()) {
    std::cout << "Bisection was successful after " << bisection.getIterations() << " iterations.\n";
  }
  else {
    std::cout << "Bisection was not successful after " << bisection.getIterations() << " iterations.\n";
  }
  std::cout << "Error = " << std::scientific << std::abs(eval.evaluate(result) - eval.evaluate(reference)) << std::endl;

  EXPECT_NEAR(result, reference, 1e-6);
}

TEST_F(OptimizationTest, LineSearch) {
  double c1 = 0.001;
  double c2 = 0.1;

  std::vector<double> alpha_0{1e-3, 1e-1, 1e1, 1e3};
  std::vector<double> reference1{1.365, 1.4413720790892741, 10., 36.88760696396662};
  std::vector<double> reference2{1.596000000186075, 1.5960000000049346, 1.5959999997572036, 1.595999998872531};

  class ObjectiveFunction {
   protected:
    double beta_;
    double value_;
    double argument_;
    double derivative_;

   public:
    double getLineSearchValue() const {
      return value_;
    }
    double getLineSearchDerivative() const {
      return derivative_;
    }

    virtual auto evaluate(double alpha) -> void = 0;
  };

  class Function1 : public ObjectiveFunction {
   public:
    virtual ~Function1() = default;

    Function1() {
      beta_ = 2;
    }

    auto evaluate(double alpha) -> void override {
      argument_ = alpha;
      value_ = -alpha / (std::pow(alpha, 2) + 2);
      derivative_ = (std::pow(alpha, 2) - 2) / std::pow((std::pow(alpha, 2) + 2), 2);
    }
  };

  class Function2 : public ObjectiveFunction {
   public:
    virtual ~Function2() = default;

    Function2() {
      beta_ = 0.004;
    }

    auto evaluate(double alpha) -> void override {
      argument_ = alpha;
      value_ = std::pow(alpha + beta_, 5) - 2 * std::pow(alpha + beta_, 4);
      derivative_ = 5 * std::pow(alpha + beta_, 4) - 8 * std::pow(alpha + beta_, 3);
    }
  };

  Function1 function1;
  Function2 function2;

  Kiwi::Optimization::MoreThuente<Function1> lineSearch1(function1);
  Kiwi::Optimization::MoreThuente<Function2> lineSearch2(function2);

  lineSearch1.setC1(c1);
  lineSearch1.setC2(c2);

  lineSearch2.setC1(0.1);
  lineSearch2.setC2(0.1);

  for (int i = 0; i < int(alpha_0.size()); ++i) {
    {
      function1.evaluate(0);
      lineSearch1(alpha_0[i], function1.getLineSearchValue(), function1.getLineSearchDerivative());
      double alpha_star = lineSearch1.getAlphaOpt();
      std::cout << "Iterations - 1 = " << lineSearch1.getIterations() << std::endl;
      EXPECT_NEAR(alpha_star, reference1[i], 1e-8);
    }
    //{
    //  function2.evaluate(0);
    //  lineSearch2(alpha_0[i], function2.getLineSearchValue(), function2.getLineSearchDerivative());
    //  double alpha_star = lineSearch2.getAlphaOpt();
    //  std::cout << "Iterations - 2 = " << lineSearch2.getIterations() << std::endl;
    //  EXPECT_NEAR(alpha_star, reference2[i], 1e-8);
    //}
  }
}

TEST_F(OptimizationTest, Optimization) {
  Eigen::VectorXd guess;
  guess.resize(2);
  guess.setZero();

  auto rosenbrock = std::make_shared<Rosenbrock>();

  auto interface = std::make_shared<Kiwi::Optimization::LineSearchInterface>(rosenbrock, guess);

  Kiwi::Optimization::Optimization gradientDescent(interface);
  gradientDescent.setThresh(1e-10);
  gradientDescent.setMaxIterations(10000);
  gradientDescent.optimize();

  Eigen::Vector2d result = interface->getResult();

  EXPECT_NEAR(result[0], 1.0, 1e-7);
  EXPECT_NEAR(result[1], 1.0, 1e-7);
}

TEST_F(OptimizationTest, OptimizationWood) {
  Eigen::VectorXd guess;
  guess.resize(4);
  guess << -3, -1, -3, -1;

  auto wood = std::make_shared<Wood>();

  auto interface = std::make_shared<Kiwi::Optimization::LineSearchInterface>(wood, guess);

  Kiwi::Optimization::Optimization gradientDescent(interface);
  gradientDescent.setThresh(1e-5);
  gradientDescent.setMaxIterations(10000);
  gradientDescent.optimize();

  Eigen::VectorXd result = interface->getResult();

  EXPECT_NEAR(result[0], 1.0, 1e-5);
  EXPECT_NEAR(result[1], 1.0, 1e-5);
  EXPECT_NEAR(result[2], 1.0, 1e-5);
  EXPECT_NEAR(result[3], 1.0, 1e-5);
}

} // namespace Scine
