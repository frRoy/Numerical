/**
 *  @file    FVSolver.hpp
 *  @brief   The finite volume solver.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef FVSOLVER_H
#define FVSOLVER_H

#include "numerical/common/Solver.hpp"
#include "FVProblem.hpp"

namespace numerical {

namespace fvm {

/**
*
*/
template <typename T>
class FVSolver : public numerical::Solver<T> {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
public:
    FVSolver(numerical::fvm::FVProblem<T>* problem)
        : numerical::Solver<T>(problem){
        this->m_problem;
    }
};

}  // namespace fvm

} // namespace numerical

#endif  // FVSOLVER_H
