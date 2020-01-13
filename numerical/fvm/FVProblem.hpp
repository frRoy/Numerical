/**
 *  @file    FVProblem.hpp
 *  @brief   The finite volume problem definition.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef FVPROBLEM_H
#define FVPROBLEM_H

#include "numerical/common/Problem.hpp"
#include "FVMesh.hpp"

namespace numerical {

namespace fvm {

/**
*
*/
template <typename T>
class FVProblem : public numerical::Problem<T> {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
    FVMesh<T>* m_mesh;  // mesh
public:
    FVProblem()
        : numerical::Problem<T>(){
        std::vector<double> lengths = {0.0, 1.0};
        double t0 = 0.0, tend = 1.0;
        int nx = 10, nt = 10;
        m_mesh = new FVMesh<T>(lengths, t0, tend, nx, nt);
        // Vec v = this->m_mesh->vertices();
        Vec v = m_mesh->vertices();
        for(int i=0; i<v.size(); i++){
            // spdlog::info("{}", v[i]);
        }
    }
    ~FVProblem(){
        delete m_mesh;
    }
};

}  // namespace fvm

} // namespace numerical

#endif  // FVPROBLEM_H
