/**
 *  @file    Utils.hpp
 *  @brief   Collections of funtions.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef UTILS_H
#define UTILS_H
#include <chrono>
#include "spdlog/spdlog.h"

namespace utils
{

/**
* Returns True if two Eigen::DenseBase are element-wise equal within a 
* tolerance.
* The tolerance values are positive, typically very small numbers. The 
* relative difference (rtol * abs(b)) and the absolute difference atol are 
* added together to compare against the absolute difference between a and b.
*
* @param a
* @param b
* @param rtol
* @param atol
* @return True if the two matrices are almost equal.
*/
template<typename DerivedA, typename DerivedB>
bool all_close(
	const Eigen::DenseBase<DerivedA>& a,
    const Eigen::DenseBase<DerivedB>& b,
    const typename DerivedA::RealScalar& rtol
        = Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
    const typename DerivedA::RealScalar& atol
        = Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
{
  return ((a.derived() - b.derived()).array().abs()
          <= (atol + rtol * b.derived().array().abs())).all();
}

} // namespace utils

#endif  // UTILS_H
