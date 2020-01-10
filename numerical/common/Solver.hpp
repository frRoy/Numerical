/**
 *  @file    Solver.hpp
 *  @brief   The solver base class.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef SOLVER_H
#define SOLVER_H

#include "Problem.hpp"


namespace numerical {

/**
*
* @param problem The problem definition.
*/
template <typename T>
class Solver {
protected:
	Problem<T>* m_problem;
public:
	Solver(Problem<T>* problem)
    : m_problem(problem) {
	}
};

} // namespace numerical

#endif  // SOLVER_H
