// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include <netcdf>
#include "MAT_Matrix.h"
#include "SignatureSymmetric.h"
// clang-format on


template<typename T>
MyMatrix<T> ReadMatrix(std::vector<int16_t> const& V) {
  std::cerr << " V=[";
  for (size_t i=0; i<V.size(); i++) {
    int16_t val = V[i];
    T val_T = UniversalScalarConversion<T,int16_t>(val);
    std::cerr << " " << val_T;
  }
  std::cerr << " ]\n";
  //
  int dim = 9;
  MyMatrix<T> M(dim, dim);
  size_t pos = 0;
  for (int u=0; u<dim; u++) {
    for (int v=0; v<=u; v++) {
      int16_t val = V[pos];
      T val_T = UniversalScalarConversion<T,int16_t>(val);
      M(u, v) = val_T;
      M(v, u) = val_T;
      pos += 1;
    }
  }
  bool test = IsPositiveDefinite(M, std::cerr);
  if (!test) {
    std::cerr << "The matrix is not positive definite\n";
    throw TerminalException{1};
  }
  return M;
}





void ReadMatrices(std::string const& FileName, size_t& index, int n_short) {
  netCDF::NcFile dataFileI(FileName, netCDF::NcFile::read);
  std::cerr << "ReadMatrices, FileName=" << FileName << "\n";
  using T = mpq_class;
  //
  netCDF::NcVar var1 = dataFileI.getVar("minQ");
  size_t number_perfect = var1.getDim(0).getSize();
  std::vector<int> V_minQ(number_perfect);
  var1.getVar(V_minQ.data());
  //
  size_t dim = 45;
  netCDF::NcVar var2 = dataFileI.getVar("form");
  std::vector<size_t> count{1, dim};
  std::vector<size_t> start{0, 0};
  std::vector<int16_t> V_form(dim);
  //
  for (size_t i=0; i<number_perfect; i++) {
    if (V_minQ[i] == n_short) {
      start[0] = i;
      var2.getVar(start, count, V_form.data());
      MyMatrix<T> M = ReadMatrix<T>(V_form);
      std::string MatFile = "Matrix_" + std::to_string(n_short) + "_" + std::to_string(index);
      WriteMatrixFile(MatFile, M);
      index += 1;
      std::cerr << "Now index=" << index << "\n";
    }
  }

}




int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "NC_AccessForms [n_short]\n";
      std::cerr << "---\n";
      std::cerr << "With n_shoft the number of shortest vector divided by 2 (so 136 is the maximal value)\n";
      throw TerminalException{1};
    }
    std::string n_short_str = argv[1];
    int n_short = ParseScalar<int>(n_short_str);
    //
    // Now reading the
    //
    size_t index = 0;
    size_t n_files = 64;
    for (size_t i_file=0; i_file<n_files; i_file++) {
      std::string FileName = "perfect_" + StringNumber(i_file, 4) + ".nc";
      ReadMatrices(FileName, index, n_short);
    }
    std::cerr << "We have written " << index << " perfect forms to file\n";
    std::cerr << "Normal termination of NC_AccessForms\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in NC_ComputeAverage\n";
    exit(e.eVal);
  }
}
