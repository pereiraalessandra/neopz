/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrixAccelerate methods.
 */

#ifdef MACOSX
#include "TPZYSMPAccelerate.h"
#include "pzfmatrix.h"

#include <Accelerate/Accelerate.h>

#include <complex>
#include <limits>
#include <type_traits>
#include <vector>

template<class TVar>
void TPZFYsmpMatrixAccelerate<TVar>::CopyFrom(const TPZMatrix<TVar> *  mat)
{
  auto *from = dynamic_cast<const TPZFYsmpMatrixAccelerate<TVar> *>(mat);
  if (from) {
    *this = *from;
  }
  else
  {
    auto *from2 = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);
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
int TPZFYsmpMatrixAccelerate<TVar>::ClassId() const{
    return Hash("TPZFYsmpMatrixAccelerate") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZFYsmpMatrixAccelerate<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                              TPZFMatrix<TVar> &z,
                                              const TVar alpha,const TVar beta,const int opt) const {
  // computes z = beta * y + alpha * opt(this)*x
  this->MultAddChecks(x,y,z,alpha,beta,opt);

  if constexpr ((std::is_same_v<TVar,float>) || (std::is_same_v<TVar,double>)){
    const int64_t m_rows = this->Rows();
    const int64_t m_cols = this->Cols();

    if(m_rows < 0 || static_cast<std::size_t>(m_rows+1) > this->fIA.size() ||
       this->fIA[0] != 0){
      PZError << __PRETTY_FUNCTION__ << "\nERROR: invalid fIA\n";
      DebugStop();
    }

    const int64_t nnz = this->fIA[m_rows];
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
    if(m_rows > maxIdx || m_cols > maxIdx || x_cols > maxIdx){
      TPZFYsmpMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
      return;
    }

    const int64_t r = (opt) ? m_cols : m_rows;
    // Nothing to multiply: z = beta*y (or 0).
    if(r == 0 || nnz == 0 || x_cols == 0 || alpha == (TVar)0){
      this->PrepareZ(y,z,beta,opt);
      return;
    }

    // Reading NeoPZ's CSR as CSC (Accelerate's format) describes A^T, not A. 
    // Row/column are relabeled, not reordered. SparseGetTranspose recovers A.
    //
    // fIA/fJA (int64_t) can't be reinterpret_cast to `long`/`int` (same
    // size, distinct types), so they are copied.
    //
    // fJA[k]==-1 is a reserved sentinel. The generic MultAdd only handles
    // it for opt!=0. For opt==0 there's no safe fallback.
    std::vector<long> columnStarts(m_rows+1);
    for(int64_t i = 0; i <= m_rows; i++){
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
      if(col == -1 && opt != 0){
        TPZFYsmpMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
        return;
      }
      if(col < 0 || col >= m_cols){
        PZError << __PRETTY_FUNCTION__
                << "\nERROR: invalid fJA[" << k << "] = " << col << '\n';
        DebugStop();
      }
      rowIndices[k] = static_cast<int>(col);
    }

    this->PrepareZ(y,z,beta,opt);

    SparseAttributes_t attr = {};
    attr.kind = SparseOrdinary;

    SparseMatrixStructure structure = {};
    // this structure describes A^T (m_cols x m_rows), see comment above
    structure.rowCount = static_cast<int>(m_cols);
    structure.columnCount = static_cast<int>(m_rows);
    structure.columnStarts = columnStarts.data();
    structure.rowIndices = rowIndices.data();
    structure.attributes = attr;
    structure.blockSize = 1;

    TVar *dataPtr = const_cast<TVar*>(this->fA.begin());

    const int64_t x_rows = x.Rows();
    const int64_t z_rows = z.Rows();
    const int64_t z_cols = z.Cols();

    // SparseMultiplyAdd computes Z += alpha*A*X; z already holds beta*y.
    if constexpr (std::is_same_v<TVar,double>){
      SparseMatrix_Double At{structure, dataPtr}; // structurally A^T
      SparseMatrix_Double A = opt ? At : SparseGetTranspose(At);
      DenseMatrix_Double X{static_cast<int>(x_rows), static_cast<int>(x_cols),
                            static_cast<int>(x_rows), SparseAttributes_t{},
                            const_cast<double*>(x.Elem())};
      DenseMatrix_Double Z{static_cast<int>(z_rows), static_cast<int>(z_cols),
                            static_cast<int>(z_rows), SparseAttributes_t{},
                            z.Elem()};
      SparseMultiplyAdd(static_cast<double>(alpha), A, X, Z);
    } else {
      SparseMatrix_Float At{structure, dataPtr};
      SparseMatrix_Float A = opt ? At : SparseGetTranspose(At);
      DenseMatrix_Float X{static_cast<int>(x_rows), static_cast<int>(x_cols),
                           static_cast<int>(x_rows), SparseAttributes_t{},
                           const_cast<float*>(x.Elem())};
      DenseMatrix_Float Z{static_cast<int>(z_rows), static_cast<int>(z_cols),
                           static_cast<int>(z_rows), SparseAttributes_t{},
                           z.Elem()};
      SparseMultiplyAdd(static_cast<float>(alpha), A, X, Z);
    }
  } else {
    // long double and complex: this Accelerate API only has Double/Float variants.
    TPZFYsmpMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
  }
}

template class TPZFYsmpMatrixAccelerate<double>;
template class TPZFYsmpMatrixAccelerate<float>;
template class TPZFYsmpMatrixAccelerate<long double>;
template class TPZFYsmpMatrixAccelerate<std::complex<float>>;
template class TPZFYsmpMatrixAccelerate<std::complex<double>>;
template class TPZFYsmpMatrixAccelerate<std::complex<long double>>;
#endif
