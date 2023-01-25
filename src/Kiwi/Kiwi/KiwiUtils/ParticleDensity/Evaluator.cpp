/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/ParticleDensity/Engine.h>
#include <Kiwi/KiwiUtils/ParticleDensity/Evaluator.h>
#include <Utils/Constants.h>
#include <fstream>

namespace Scine {
namespace Kiwi {
namespace ParticleDensity {

void Evaluator::constructGrid() {
  Eigen::Vector3d start;
  start << settings_.xmin, settings_.ymin, settings_.zmin;
  Eigen::Vector3d end;
  end << settings_.xmax, settings_.ymax, settings_.zmax;

  step = (end - start) * (1. / (settings_.N1D - 1));
  std::cout << "Step: " << step.transpose() << std::endl;
  std::cout << "N = " << settings_.N1D << std::endl;
  std::cout << "1/(N-1) " << (1. / (settings_.N1D - 1)) << std::endl;

  grid.reserve(settings_.N1D);

  for (auto n = 0; n < settings_.N1D; ++n) {
    grid.emplace_back(start + n * step, 0.0);
  }
}

void Evaluator::runCalculationGrid() {
  std::cout << std::endl << std::endl;
  std::cout << "- Start evaluation" << std::flush << std::endl;

  auto& clock = Clock::getInstance();
  clock.time("density");
#pragma omp parallel
  {
    Engine engine(molecule_->at(settings_.type).basisSet);
#pragma omp for schedule(dynamic)
    for (auto i = 0; i < int(grid.size()); ++i) {
      compute(grid[i].second, engine, grid[i].first(0), grid[i].first(1), grid[i].first(2));
      // if (omp_get_thread_num() == 0) {
      //  std::cout << " " << i + 1 << " / " << grid.size() / omp_get_max_threads() << std::endl;
      //}
    }
  }
  std::cout << "- Finished evaluation...      \t\t" << std::flush;
  clock.time("density");
}

void Evaluator::constructCube() {
  const auto& NX = settings_.N3DX;
  const auto& NY = settings_.N3DY;
  const auto& NZ = settings_.N3DZ;

  step(0) = (settings_.xmax - settings_.xmin) / (NX - 1);
  step(1) = (settings_.ymax - settings_.ymin) / (NY - 1);
  step(2) = (settings_.zmax - settings_.zmin) / (NZ - 1);

  using std::cout;
  using std::endl;
  using std::flush;
  using std::setprecision;
  using std::setw;

  cout << endl << "- I will be looking at this portion of space:" << endl << "\tx\t\ty\t\tz" << setprecision(3);
  cout << endl << setw(13) << settings_.xmin << setw(16) << settings_.ymin << setw(16) << settings_.zmin;
  cout << endl
       << setw(13) << settings_.xmin + (NX - 1) * step(0) << setw(16) << settings_.ymin + (NY - 1) * step(1) << setw(16)
       << settings_.zmin + (NZ - 1) * step(2) << endl
       << setprecision(14);

  auto& clock = Clock::getInstance();
  clock.time("allocateCube");
  cout << "- Constructing the cube... \t\t" << flush;
  clock.time("allocateCube");
  cube = std::vector<std::vector<std::vector<double>>>(NX, std::vector<std::vector<double>>(NY, std::vector<double>(NZ, 0.0)));
  clock.time("allocateCube");
  cout << "- num points to be evaluated: " << NX * NY * NZ << endl << endl;
}

void Evaluator::runCalculationCube() {
  std::cout << std::endl << std::endl;
  std::cout << "- Start evaluation" << std::flush << std::endl;

  auto& clock = Clock::getInstance();
  clock.time("density");

  const auto& NX = settings_.N3DX;
  const auto& NY = settings_.N3DY;
  const auto& NZ = settings_.N3DZ;

#pragma omp parallel
  {
    Engine engine(molecule_->at(settings_.type).basisSet);

#pragma omp parallel for schedule(dynamic) default(shared) collapse(3)
    for (auto ix = 0; ix < NX; ix++) {
      for (auto iy = 0; iy < NY; iy++) {
        for (auto iz = 0; iz < NZ; iz++) {
          auto x = settings_.xmin + ix * step(0);
          auto y = settings_.ymin + iy * step(1);
          auto z = settings_.zmin + iz * step(2);
          compute(cube[ix][iy][iz], engine, x, y, z);
        }
      }
    }
  }
  std::cout << "- Finished evaluation...      \t\t" << std::flush;
  clock.time("density");
}

void Evaluator::compute(double& destination, Engine& engine, double x, double y, double z) {
  switch (settings_.densityMode) {
    case DensityMode::MO: {
      auto result = engine.evaluateOrbitals(x, y, z);
      if (molecule_->at(settings_.type).isRestricted) {
        destination = data_->C[settings_.type].restrictedMatrix().col(settings_.moIndex).transpose() * result;
      }
      else {
        if (settings_.msIndex == 0) {
          destination = data_->C[settings_.type].alphaMatrix().col(settings_.moIndex).transpose() * result;
        }
        else {
          destination = data_->C[settings_.type].betaMatrix().col(settings_.moIndex).transpose() * result;
        }
      }
    } break;
    case DensityMode::PartDens: {
      auto result = engine.evaluateDensity(x, y, z);
      if (molecule_->at(settings_.type).isRestricted) {
        auto resultMO =
            data_->C[settings_.type].restrictedMatrix().transpose() * result * data_->C[settings_.type].restrictedMatrix();
        destination = (resultMO * densityMatrix.restrictedMatrix()).trace();
      }
      else {
        destination = (data_->C[settings_.type].alphaMatrix().transpose() * result *
                       data_->C[settings_.type].alphaMatrix() * densityMatrix.alphaMatrix())
                          .trace();
        if (molecule_->at(settings_.type).msVector[1] > 0) {
          destination += (data_->C[settings_.type].betaMatrix().transpose() * result *
                          data_->C[settings_.type].betaMatrix() * densityMatrix.betaMatrix())
                             .trace();
        }
      }
    } break;
  }
}

void Evaluator::writeGrid() const {
  std::cout << "- Writing on disk    \t\t";
  std::ofstream file(settings_.filename.c_str(), std::ios::out);

  for (const auto& point : grid) {
    file << std::setprecision(4) << std::fixed << Utils::toAngstrom(Utils::Bohr(point.first[0])) << ","
         << Utils::toAngstrom(Utils::Bohr(point.first[1])) << "," << Utils::toAngstrom(Utils::Bohr(point.first[2]))
         << "," << std::scientific << std::setprecision(14) << point.second << "\n";
  }
  file.close();
}

void Evaluator::writeCube() const {
  std::cout << "- Writing on disk    \t\t";

  const auto& NX = settings_.N3DX;
  const auto& NY = settings_.N3DY;
  const auto& NZ = settings_.N3DZ;
  std::ofstream file(settings_.filename.c_str(), std::ios::out);

  file << "Cube file generated with Kiwi." << std::endl;
  file << "For information on the format, visit http://paulbourke.net/dataformats/cube/." << std::endl;
  // Number of atoms -- origin of volumetric data: x y z
  file << std::setw(5) << molecule_->getGeometry().size() << std::setw(12) << std::setprecision(6) << std::fixed
       << settings_.xmin << std::setw(12) << settings_.ymin << std::setw(12) << settings_.zmin << std::endl;
  // Next three lines give number of vosxels followed by axis vector.
  file << std::setw(5) << NX << std::setw(12) << std::setprecision(6) << std::fixed << step(0) << std::setw(12) << 0.0
       << std::setw(12) << 0.0 << "\n";
  file << std::setw(5) << NY << std::setw(12) << std::setprecision(6) << std::fixed << 0.0 << std::setw(12) << step(1)
       << std::setw(12) << 0.0 << "\n";
  file << std::setw(5) << NZ << std::setw(12) << std::setprecision(6) << std::fixed << 0.0 << std::setw(12) << 0.0
       << std::setw(12) << step(2) << "\n";
  // For each atom:
  // atomic number, charge, x y z
  for (const auto& atom : molecule_->getGeometry()) {
    file << std::setw(5) << Utils::ElementInfo::Z(atom.getElementType()) << std::setw(12) << std::setprecision(6)
         << std::fixed << 0.0 << std::setw(12) << atom.getPosition()(0) << std::setw(12) << atom.getPosition()(1)
         << std::setw(12) << atom.getPosition()(2) << std::endl;
  }

  for (auto ix = 0; ix < NX; ix++) {
    for (auto iy = 0; iy < NY; iy++) {
      for (auto iz = 0; iz < NZ; iz++) {
        file << std::setw(12) << std::setprecision(6) << std::fixed << std::scientific << cube[ix][iy][iz] << "    ";
        if (iz % 6 == 5) {
          file << "\n";
        }
      }
      file << "\n";
    }
  }
  file.close();
}

} // namespace ParticleDensity
} // namespace Kiwi
} // namespace Scine