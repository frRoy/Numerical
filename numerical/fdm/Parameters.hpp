/**
 *  @file    Parameters.hpp
 *  @brief   Definition of the diffusion problem.
 *  @author  Francois Roy
 *  @date    12/04/2019
 */
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include <Eigen/Core>
#include "spdlog/spdlog.h"


namespace numerical 
{
namespace fdm 
{

/**
* Defines the parameters for the scalar finite difference problem.
*/
template<typename T>
struct Parameters
{
    T alpha, theta;
    std::vector<std::vector<T>> lengths;
    T t0, tend;
    int nt;
    std::vector<int> n;
    Parameters()
    {
    	// default parameters
    	alpha = 1.0;
        theta = 0.5;
        // unit cube
    	lengths = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
        // 10 divisions per dimensions
    	n = {10, 10, 10};
        // initial time
    	t0 = 0.0;
        // final time
    	tend = 1.0;
        // number of uniform time step
    	nt = 10;
    }
};

} // end namespace fdm
}  // end namesapce numerical

#endif  // PARAMETERS_H
