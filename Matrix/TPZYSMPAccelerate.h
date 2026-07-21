/**
 * @file
 * @brief Contains TPZFYsmpMatrixAccelerate class which implements a non symmetric sparse matrix
 * whose MultAdd is accelerated using Apple's Accelerate framework (Sparse Solvers API).
 * Only available when compiling for macOS without USING_MKL.
 */

#ifndef YSMPMATACCELERATE_H
#define YSMPMATACCELERATE_H

#include "TPZYSMPMatrix.h"

/**
 * @brief Implements a non symmetric sparse matrix whose MultAdd uses
 * Apple's Accelerate framework Sparse Solvers API. \ref matrix "Matrix"
 * @ingroup matrix
 * @note Accelerate's Sparse Solvers API does not support complex types,
 * so this class only accelerates float and double instantiations. Other
 * types fall back to DebugStop.
 */
template<class TVar>
class TPZFYsmpMatrixAccelerate : public TPZFYsmpMatrix<TVar>{
public :

  TPZFYsmpMatrixAccelerate() : TPZRegisterClassId(&TPZFYsmpMatrixAccelerate::ClassId),
  TPZFYsmpMatrix<TVar>() {}
  /** @brief Constructors from parent class*/
  using TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix;
  /** @brief Copy constructor */
  TPZFYsmpMatrixAccelerate(const TPZFYsmpMatrixAccelerate<TVar> &cp) = default;
  /** @brief Move constructor*/
  TPZFYsmpMatrixAccelerate(TPZFYsmpMatrixAccelerate<TVar> &&cp) = default;
  /** @brief Copy-assignment operator*/
  TPZFYsmpMatrixAccelerate &operator=(const TPZFYsmpMatrixAccelerate<TVar> &copy) = default;
  /** @brief Move-assignment operator*/
  TPZFYsmpMatrixAccelerate &operator=(TPZFYsmpMatrixAccelerate<TVar> &&copy) = default;

  /** @brief Copy constructor from generic sparse matrix*/
  TPZFYsmpMatrixAccelerate(const TPZFYsmpMatrix<TVar> &cp)
    : TPZFYsmpMatrix<TVar>(cp) {}
  /** @brief Move constructor from generic sparse matrix*/
  TPZFYsmpMatrixAccelerate(TPZFYsmpMatrix<TVar> &&rval)
    : TPZFYsmpMatrix<TVar>(rval) {}
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixAccelerate &operator=(const TPZFYsmpMatrix<TVar> &cp)
  { TPZFYsmpMatrix<TVar>::operator=(cp); return *this;}
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixAccelerate &operator=(TPZFYsmpMatrix<TVar> &&rval)
  { TPZFYsmpMatrix<TVar>::operator=(rval); return *this;}

  inline TPZFYsmpMatrixAccelerate<TVar>*NewMatrix() const override {return new TPZFYsmpMatrixAccelerate<TVar>{};}
  CLONEDEF(TPZFYsmpMatrixAccelerate)
  /** @brief Destructor */
  ~TPZFYsmpMatrixAccelerate() = default;

  /** @brief Creates a copy from another sparse matrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override;

  void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
               const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;

  int ClassId() const override;
};

#endif
