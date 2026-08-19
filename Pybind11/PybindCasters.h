/// @file PybindCasters.h
/// @brief Type casters converting Python sequences to the NeoPZ containers
/// @note conversion is by copy, so output parameters never reach Python.
/// Must be included by every translation unit, which Bindings.h ensures

#ifndef PZ_PYBIND_CASTERS_H
#define PZ_PYBIND_CASTERS_H

#include <pybind11/pybind11.h>

#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzvec.h"

namespace pybind11 {
namespace detail {

template <typename T>
struct type_caster<TPZVec<T>> {
  using value_conv = make_caster<T>;

  PYBIND11_TYPE_CASTER(TPZVec<T>, const_name("list[") + value_conv::name + const_name("]"));

  bool load(handle src, bool convert) {
    if (!isinstance<sequence>(src) || isinstance<str>(src) || isinstance<bytes>(src))
      return false;
    auto seq = reinterpret_borrow<sequence>(src);
    value.Resize(static_cast<int64_t>(seq.size()));
    int64_t i = 0;
    for (auto item : seq) {
      value_conv conv;
      if (!conv.load(item, convert))
        return false;
      value[i++] = cast_op<T &&>(std::move(conv));
    }
    return true;
  }

  static handle cast(const TPZVec<T> &src, return_value_policy policy, handle parent) {
    list out(src.size());
    for (int64_t i = 0; i < src.size(); i++) {
      auto item = reinterpret_steal<object>(value_conv::cast(src[i], policy, parent));
      if (!item)
        return handle();
      PyList_SET_ITEM(out.ptr(), static_cast<ssize_t>(i), item.release().ptr());
    }
    return out.release();
  }
};

template <typename T, int N>
struct type_caster<TPZManVector<T, N>> {
  using value_conv = make_caster<T>;
  using ManVec = TPZManVector<T, N>;

  PYBIND11_TYPE_CASTER(ManVec, const_name("list[") + value_conv::name + const_name("]"));

  bool load(handle src, bool convert) {
    if (!isinstance<sequence>(src) || isinstance<str>(src) || isinstance<bytes>(src))
      return false;
    auto seq = reinterpret_borrow<sequence>(src);
    value.Resize(static_cast<int64_t>(seq.size()));
    int64_t i = 0;
    for (auto item : seq) {
      value_conv conv;
      if (!conv.load(item, convert))
        return false;
      value[i++] = cast_op<T &&>(std::move(conv));
    }
    return true;
  }

  static handle cast(const ManVec &src, return_value_policy policy, handle parent) {
    list out(src.size());
    for (int64_t i = 0; i < src.size(); i++) {
      auto item = reinterpret_steal<object>(value_conv::cast(src[i], policy, parent));
      if (!item)
        return handle();
      PyList_SET_ITEM(out.ptr(), static_cast<ssize_t>(i), item.release().ptr());
    }
    return out.release();
  }
};

template <typename T>
struct type_caster<TPZFMatrix<T>> {
  using value_conv = make_caster<T>;

  PYBIND11_TYPE_CASTER(TPZFMatrix<T>, const_name("list[list[") + value_conv::name + const_name("]]"));

  bool load(handle src, bool convert) {
    if (!isinstance<sequence>(src) || isinstance<str>(src) || isinstance<bytes>(src))
      return false;
    auto rows = reinterpret_borrow<sequence>(src);
    const int64_t nrows = static_cast<int64_t>(rows.size());
    int64_t ncols = 0;
    if (nrows > 0) {
      if (!isinstance<sequence>(rows[0]))
        return false;
      ncols = static_cast<int64_t>(reinterpret_borrow<sequence>(rows[0]).size());
    }
    value.Resize(nrows, ncols);
    for (int64_t i = 0; i < nrows; i++) {
      if (!isinstance<sequence>(rows[i]))
        return false;
      auto row = reinterpret_borrow<sequence>(rows[i]);
      if (static_cast<int64_t>(row.size()) != ncols)
        return false;
      for (int64_t j = 0; j < ncols; j++) {
        value_conv conv;
        if (!conv.load(row[j], convert))
          return false;
        value(i, j) = cast_op<T &&>(std::move(conv));
      }
    }
    return true;
  }

  static handle cast(const TPZFMatrix<T> &src, return_value_policy policy, handle parent) {
    list out(src.Rows());
    for (int64_t i = 0; i < src.Rows(); i++) {
      list row(src.Cols());
      for (int64_t j = 0; j < src.Cols(); j++) {
        auto item = reinterpret_steal<object>(
            value_conv::cast(src.GetVal(i, j), policy, parent));
        if (!item)
          return handle();
        PyList_SET_ITEM(row.ptr(), static_cast<ssize_t>(j), item.release().ptr());
      }
      PyList_SET_ITEM(out.ptr(), static_cast<ssize_t>(i), row.release().ptr());
    }
    return out.release();
  }
};

} // namespace detail
} // namespace pybind11

#endif // PZ_PYBIND_CASTERS_H
