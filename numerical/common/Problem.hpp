/**
 *  @file    Problem.hpp
 *  @brief   The problem base class.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef PROBLEM_H
#define PROBLEM_H

#include "Mesh.hpp"


namespace numerical {

/**
*
*/
template <typename T>
class Problem {
protected:
	Mesh<T>* m_mesh;  // mesh
public:
	Problem() {
		std::vector<double> lengths = {0.0, 1.0};
        double t0 = 0.0, tend = 1.0;
        int nx = 10, nt = 10;
		m_mesh = new Mesh<T>(lengths, t0, tend, nx, nt);
	}
	virtual ~Problem(){
        delete m_mesh;
    }
};
    
} // namespace numerical

#endif  // PROBLEM_H
