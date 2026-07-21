/**
 * @file
 * @brief Contains TPZSYsmpMatrixAccelerate class which implements a symmetric sparse matrix
 * whose MultAdd is accelerated using Apple's Accelerate framework (Sparse Solvers API).
 * Only available when compiling for macOS without USING_MKL (see CMakeLists.txt).
 */

#ifndef SYSMPMATACCELERATE_H
#define SYSMPMATACCELERATE_H

#include "TPZSYSMPMatrix.h"

/**
 * @brief Implements a symmetric sparse matrix whose MultAdd uses
 * Apple's Accelerate framework Sparse Solvers API. \ref matrix "Matrix"
 * @ingroup matrix
 * @note Accelerate's Sparse Solvers API does not support complex types,
 * so this class only accelerates float and double instantiations. Other
 * types fall back to DebugStop.
 */
template<class TVar>
class TPZSYsmpMatrixAccelerate : public TPZSYsmpMatrix<TVar>{
public :

  TPZSYsmpMatrixAccelerate() : TPZRegisterClassId(&TPZSYsmpMatrixAccelerate::ClassId),
  TPZSYsmpMatrix<TVar>() {}
  /** @brief Constructors from parent class*/
  using TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix;
  /** @brief Copy constructor */
  TPZSYsmpMatrixAccelerate(const TPZSYsmpMatrixAccelerate<TVar> &cp) = default;
  /** @brief Move constructor*/
  TPZSYsmpMatrixAccelerate(TPZSYsmpMatrixAccelerate<TVar> &&cp) = default;
  /** @brief Copy-assignment operator*/
  TPZSYsmpMatrixAccelerate &operator=(const TPZSYsmpMatrixAccelerate<TVar> &copy) = default;
  /** @brief Move-assignment operator*/
  TPZSYsmpMatrixAccelerate &operator=(TPZSYsmpMatrixAccelerate<TVar> &&copy) = default;

  /** @brief Copy constructor from generic sparse matrix*/
  TPZSYsmpMatrixAccelerate(const TPZSYsmpMatrix<TVar> &cp)
    : TPZSYsmpMatrix<TVar>(cp) {}
  /** @brief Move constructor from generic sparse matrix*/
  TPZSYsmpMatrixAccelerate(TPZSYsmpMatrix<TVar> &&rval)
    : TPZSYsmpMatrix<TVar>(rval) {}
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixAccelerate &operator=(const TPZSYsmpMatrix<TVar> &cp)
  { TPZSYsmpMatrix<TVar>::operator=(cp); return *this;}
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixAccelerate &operator=(TPZSYsmpMatrix<TVar> &&rval)
  { TPZSYsmpMatrix<TVar>::operator=(rval); return *this;}

  inline TPZSYsmpMatrixAccelerate<TVar>*NewMatrix() const override {return new TPZSYsmpMatrixAccelerate<TVar>{};}
  CLONEDEF(TPZSYsmpMatrixAccelerate)
  /** @brief Destructor */
  ~TPZSYsmpMatrixAccelerate() = default;

  /** @brief Creates a copy from another sparse matrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override;

  void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
               const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;

  int ClassId() const override;
};

#endif
