/**
 *  @file    Performances.hpp
 *  @brief   Check performances for different methods.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef PERFORMANCES_H
#define PERFORMANCES_H

#include <vector>
#include <Eigen/SparseCore>
#include "spdlog/spdlog.h"

namespace bench {

namespace performances {

const unsigned int n = 10000;
const Eigen::VectorXd u = Eigen::VectorXd::Random(n + 1);

/**
* Check time taken to compute a simple looped difference.
* @return The computed vector.
*/
Eigen::VectorXd loop(){
	Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; i++) {
        d[i] = u[i+1] - u[i];       
    }
    return d;
}
/**
* check time taken to compute a vectorized looped difference.
* @return The computed vector.
*/
Eigen::VectorXd vectorized_loop(){
	// Eigen::VectorXd u = Eigen::VectorXd::Constant(n + 1, 2.4);
	Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
	d = u.segment(1, u.size()-1) - u.segment(0, u.size()-1);
    return d;
}

}
}

#endif  // PERFORMANCES_H