/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/KiwiUtils/BoFciDumper.h>
#include <Utils/IO/Regex.h>
#include <Utils/Technical/ScopedLocale.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

namespace Scine {
namespace Kiwi {

double BoFciDumper::integralThreshold = 1e-16;

void BoFciDumper::write(const std::string& filename, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                        const BoFciDumpData& fciDumpData, format f) {
  std::ofstream fout;
  if (f == format::fci) {
    fout.open(filename);
  }
  if (!fout.is_open()) {
    throw std::runtime_error("Problem when opening/creating file " + filename);
  }
  return write(fout, hCore, eris, fciDumpData, f);
}

void BoFciDumper::write(std::ostream& out, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                        const BoFciDumpData& fciDumpData, format f) {
  if (f == format::fci) {
    writeFci(out, hCore, eris, fciDumpData);
  }
  else {
    throw std::runtime_error("Unsupported format to write FciDump files to.");
  }
}

void BoFciDumper::writeFci(std::ostream& out, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                           const BoFciDumpData& fciDumpData) {
  auto scopedLocale = Utils::ScopedLocale::cLocale();
  writeHeaderFci(out, fciDumpData);
  out << std::scientific << std::left << std::setprecision(15);
  writeTwoElectronIntegralsFci(out, eris);
  writeOneElectronIntegralsFci(out, hCore);
  // Write core energy
  out << std::setw(25) << fciDumpData.coreEnergy << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 0
      << std::setw(10) << 0;
}
void BoFciDumper::writeHeaderFci(std::ostream& out, const Scine::Kiwi::BoFciDumpData& fciDumpData) {
  out << "&FCI NORB=" << fciDumpData.nOrbitals << ",NELEC=" << fciDumpData.nElectrons
      << ",MS2=" << fciDumpData.spinPolarization << "," << std::endl;
  out << " ORBSYM=";
  for (auto sym : fciDumpData.orbitalSymmetry)
    out << sym << ",";
  out << std::endl;
  out << " ISYM=" << fciDumpData.wfSymIfLowestOccupied << std::endl;
  if (fciDumpData.unrestrictedReference)
    out << "UHF" << std::endl;
  out << "&END" << std::endl;
}

void BoFciDumper::writeOneElectronIntegralsFci(std::ostream& out, const Eigen::MatrixXd& hCore) {
  for (int row = 0; row < hCore.rows(); ++row) {
    for (int col = 0; col < row + 1; ++col) {
      if (std::abs(hCore(row, col)) > integralThreshold) {
        out << std::setw(25) << hCore(row, col) << std::setw(10) << row + 1 << std::setw(10) << col + 1 << std::setw(10)
            << 0 << std::setw(10) << 0 << std::endl;
      }
    }
  }
}

void BoFciDumper::writeTwoElectronIntegralsFci(std::ostream& out, const Eigen::MatrixXd& eris) {
  auto N = int(std::sqrt(eris.rows()));

  for (auto i = 0; i != N; ++i) {
    for (auto j = 0; j <= i; ++j) {
      for (auto k = 0; k <= i; ++k) {
        const auto l_max = (i == k) ? j : k;
        for (auto l = 0; l <= l_max; ++l) {
          double value = eris(i * N + j, k * N + l);
          if (std::abs(value) > integralThreshold) {
            auto index = Integrals::TwoBody::getMappedIndex<Integrals::TwoBody::IntegralSymmetry::eightfold>({i, j, k, l});
            out << std::setw(25) << value << std::setw(10) << index[0] + 1 << std::setw(10) << index[1] + 1
                << std::setw(10) << index[2] + 1 << std::setw(10) << index[3] + 1 << std::endl;
          }
        }
      }
    }
  }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> BoFciDumper::read(const std::string& filename,
                                                                              BoFciDumper::format f) {
  std::ifstream fin;
  fin.open(filename);
  if (!fin.is_open()) {
    throw std::runtime_error("Problem when opening file " + filename);
  }
  return read(fin, f);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> BoFciDumper::read(std::istream& in, BoFciDumper::format f) {
  if (f == format::fci)
    return readFci(in);
  throw std::runtime_error("Unsupported format to read FciDump files from.");
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> BoFciDumper::readFci(std::istream& in) {
  Eigen::MatrixXd hCoreMo;
  Eigen::MatrixXd eris;
  auto fciDumpData = parseFciHeader(in);
  std::tie(hCoreMo, eris) = parseFciIntegrals(in, fciDumpData);
  return std::make_tuple(hCoreMo.selfadjointView<Eigen::Lower>(), eris, fciDumpData);
}

BoFciDumpData BoFciDumper::parseFciHeader(std::istream& in) {
  auto scopedLocale = Utils::ScopedLocale::cLocale();
  std::string number = Utils::Regex::capturingIntegerNumber();
  std::regex r(number);
  BoFciDumpData fciDumpData{};
  std::string line;
  // skip lines until &FCI is found or eof.
  do {
    std::getline(in, line);
  } while ((line.find("&FCI") == std::string::npos) && !in.eof());
  do {
    std::string field;
    auto spaceStrippedLine = std::stringstream(line);
    if (spaceStrippedLine.str().find("ORBSYM") != std::string::npos) {
      std::smatch match;
      while (std::regex_search(line, match, r)) {
        fciDumpData.orbitalSymmetry.push_back(std::stoi(match[0]));
        line = match.suffix();
      }
    }
    else {
      while (std::getline(spaceStrippedLine, field, ',')) {
        if (field.find("NORB") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.nOrbitals = std::stoi(match[0]);
        }
        else if (field.find("NELEC") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.nElectrons = std::stoi(match[0]);
        }
        else if (field.find("MS2") != std::string::npos) {
          auto pos = field.find("MS2") + 3;
          auto substring = field.substr(pos);
          std::smatch match;
          std::regex_search(substring, match, r);
          fciDumpData.spinPolarization = std::stoi(match[0]);
        }
        else if (field.find("ISYM") != std::string::npos) {
          std::smatch match;
          std::regex_search(field, match, r);
          fciDumpData.wfSymIfLowestOccupied = std::stoi(match[0]);
        }
        else if (field.find("UHF") != std::string::npos) {
          fciDumpData.unrestrictedReference = true;
        }
      }
    }
    std::getline(in, line);
  } while (line.find("&END") == std::string::npos && !in.eof());
  return fciDumpData;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> BoFciDumper::parseFciIntegrals(std::istream& in, BoFciDumpData& fciDumpData) {
  int N = fciDumpData.nOrbitals;
  Eigen::MatrixXd hCoreMo = Eigen::MatrixXd::Zero(N, N);
  Eigen::MatrixXd eris = Eigen::MatrixXd::Zero(N * N, N * N);
  std::string line;
  while (std::getline(in, line)) {
    std::stringstream buffer(line);
    double value = 0.0;
    int i = 0, j = 0, k = 0, l = 0;
    buffer >> value >> i >> j >> k >> l;
    if (i != 0 && j != 0 && k != 0 && l != 0) {
      auto indices = Integrals::TwoBody::getSymmetricIndices<Integrals::TwoBody::IntegralSymmetry::eightfold>(
          {i - 1, j - 1, k - 1, l - 1});
      for (auto const& idx : indices) {
        eris(idx[0] * N + idx[1], idx[2] * N + idx[3]) = value;
      }
    }
    else if (i == 0 && j == 0 && k == 0 && l == 0)
      fciDumpData.coreEnergy = value;
    else if (k == 0 && l == 0) {
      hCoreMo(i - 1, j - 1) = value;
      if (i != j) {
        hCoreMo(j - 1, i - 1) = value;
      }
    }
  }
  return std::make_tuple(hCoreMo, eris);
}

} // namespace Kiwi
} // namespace Scine
