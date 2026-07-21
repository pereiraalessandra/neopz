/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrixAccelerate methods.
 */

#ifdef MACOSX
#include "TPZSYSMPAccelerate.h"
#include "pzfmatrix.h"

#include <Accelerate/Accelerate.h>

#include <complex>
#include <limits>
#include <type_traits>
#include <vector>

template<class TVar>
void TPZSYsmpMatrixAccelerate<TVar>::CopyFrom(const TPZMatrix<TVar> *  mat)
{
  auto *from = dynamic_cast<const TPZSYsmpMatrixAccelerate<TVar> *>(mat);
  if (from) {
    *this = *from;
  }
  else
  {
    auto *from2 = dynamic_cast<const TPZSYsmpMatrix<TVar> *>(mat);
    if (from2) {
      *this = *from2;
      return;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Called with incompatible type\n.";
    PZError<<"Aborting...\n";
    DebugStop();
  }
}

template<class TVar>
int TPZSYsmpMatrixAccelerate<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrixAccelerate") ^ TPZSYsmpMatrix<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZSYsmpMatrixAccelerate<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                              TPZFMatrix<TVar> &z,
                                              const TVar alpha,const TVar beta,const int opt) const {
  // computes z = beta * y + alpha * this*x  (symmetric: opt is irrelevant)
  this->MultAddChecks(x,y,z,alpha,beta,opt);

  if constexpr ((std::is_same_v<TVar,float>) || (std::is_same_v<TVar,double>)){
    const int64_t r = this->Rows();

    if(r < 0 || static_cast<std::size_t>(r+1) > this->fIA.size() ||
       this->fIA[0] != 0){
      PZError << __PRETTY_FUNCTION__ << "\nERROR: invalid fIA\n";
      DebugStop();
    }

    const int64_t nnz = this->fIA[r];
    if(nnz < 0 || static_cast<std::size_t>(nnz) > this->fJA.size() ||
       static_cast<std::size_t>(nnz) > this->fA.size()){
      PZError << __PRETTY_FUNCTION__ << "\nERROR: invalid CSR structure\n";
      DebugStop();
    }

    const int64_t x_cols = x.Cols();

    // When beta==0, MultAddChecks skips its y/x column check, but PrepareZ
    // still sizes z from y.Cols() -- without this check, z could end up
    // with the wrong column count for the SparseMultiplyAdd call.
    if(y.Cols() != x_cols){
      PZError << __PRETTY_FUNCTION__
              << "\nERROR: y and x must have the same number of columns\n";
      DebugStop();
    }

    // Accelerate's structures use 32-bit row/column counts.
    // Fall back to NeoPZ for larger matrices.
    constexpr int64_t maxIdx = std::numeric_limits<int>::max();
    if(r > maxIdx || x_cols > maxIdx){
      TPZSYsmpMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
      return;
    }

    // Nothing to multiply: z = beta*y (or 0).
    if(r == 0 || nnz == 0 || x_cols == 0 || alpha == (TVar)0){
      this->PrepareZ(y,z,beta,opt);
      return;
    }

    // NeoPZ's upper triangle (CSR), read as CSC, is Accelerate's lower
    // triangle (Apple's convention for symmetric matrices).
    // A == A^T needs no transpose.
    //
    // fIA/fJA (int64_t) can't be reinterpret_cast to `long`/`int` (same
    // size, distinct types), so they are copied.
    std::vector<long> columnStarts(r+1);
    for(int64_t i = 0; i <= r; i++){
      const int64_t value = this->fIA[i];
      if(value < 0 || value > nnz || (i > 0 && value < this->fIA[i-1])){
        PZError << __PRETTY_FUNCTION__ << "\nERROR: fIA is not monotonic\n";
        DebugStop();
      }
      columnStarts[i] = static_cast<long>(value);
    }
    std::vector<int> rowIndices(nnz);
    for(int64_t k = 0; k < nnz; k++){
      const int64_t col = this->fJA[k];
      if(col < 0 || col >= r){
        PZError << __PRETTY_FUNCTION__
                << "\nERROR: invalid fJA[" << k << "] = " << col << '\n';
        DebugStop();
      }
      rowIndices[k] = static_cast<int>(col);
    }

    this->PrepareZ(y,z,beta,opt);

    SparseAttributes_t attr = {};
    attr.kind = SparseSymmetric;
    attr.triangle = SparseLowerTriangle;

    SparseMatrixStructure structure = {};
    structure.rowCount = static_cast<int>(r);
    structure.columnCount = static_cast<int>(r);
    structure.columnStarts = columnStarts.data();
    structure.rowIndices = rowIndices.data();
    structure.attributes = attr;
    structure.blockSize = 1;

    TVar *dataPtr = const_cast<TVar*>(this->fA.begin());

    const int64_t x_rows = x.Rows();
    const int64_t z_rows = z.Rows();
    const int64_t z_cols = z.Cols();

    if constexpr (std::is_same_v<TVar,double>){
      SparseMatrix_Double A{structure, dataPtr};
      DenseMatrix_Double X{static_cast<int>(x_rows), static_cast<int>(x_cols),
                            static_cast<int>(x_rows), SparseAttributes_t{},
                            const_cast<double*>(x.Elem())};
      DenseMatrix_Double Z{static_cast<int>(z_rows), static_cast<int>(z_cols),
                            static_cast<int>(z_rows), SparseAttributes_t{},
                            z.Elem()};
      SparseMultiplyAdd(static_cast<double>(alpha), A, X, Z);
    } else {
      SparseMatrix_Float A{structure, dataPtr};
      DenseMatrix_Float X{static_cast<int>(x_rows), static_cast<int>(x_cols),
                           static_cast<int>(x_rows), SparseAttributes_t{},
                           const_cast<float*>(x.Elem())};
      DenseMatrix_Float Z{static_cast<int>(z_rows), static_cast<int>(z_cols),
                           static_cast<int>(z_rows), SparseAttributes_t{},
                           z.Elem()};
      SparseMultiplyAdd(static_cast<float>(alpha), A, X, Z);
    }
  } else {
    // long double and complex are not supported by this Accelerate API.
    TPZSYsmpMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
  }
}

template class TPZSYsmpMatrixAccelerate<double>;
template class TPZSYsmpMatrixAccelerate<float>;
template class TPZSYsmpMatrixAccelerate<long double>;
template class TPZSYsmpMatrixAccelerate<std::complex<float>>;
template class TPZSYsmpMatrixAccelerate<std::complex<double>>;
template class TPZSYsmpMatrixAccelerate<std::complex<long double>>;
#endif
