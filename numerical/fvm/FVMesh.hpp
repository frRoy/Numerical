/**
 *  @file    FVMesh.hpp
 *  @brief   The finite volume element mesh for a line.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef FVMESH_H
#define FVMESH_H

#include "numerical/common/Mesh.hpp"

namespace numerical {

namespace fvm {

/**
* Holds data structures for a uniform mesh on a line in space, plus a 
* uniform mesh in time.
* 
* The finite volume mesh has n elements. Two discretization schemes are
* permitted, cell centered and node centered.
*
* Cell-centered discretizations assume solutions are defined at the centers of 
* the primal grid cells with the primal cells serving as the control volumes. 
* The cell center coordinates are typically defined as the averages of the 
* coordinates of the cell's vertexes.
*
* Node-centered discretizations assume solutions are defined at the primal 
* mesh nodes. For node-centered schemes, control volumes are constructed
* around the mesh nodes by the median-dual partition: the centers
* of primal cells are connected with the midpoints of the surrounding faces. 
* These non-overlapping control volumes cover the entire computational domain 
* and compose a mesh that is dual to the primal mesh.
*
* @param lengths A 2-lists of min and max coordinates in the spatial direction.
* @param t0 The initial time in time mesh.
* @param tend The final time in time mesh.
* @param nt Number of cells in time mesh.
* @param n Number of cells in the spatial direction.
* @param disc The discretization scheme, default = 0 (node centered).
*/
template <typename T>
class FVMesh : public numerical::Mesh<T> {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
    const int m_discretization; // 0 for node centered, 1 for cell centered
public:
    FVMesh(const std::vector<T>& lengths, T t0, T tend, int nx, int nt, 
        int disc=0)
        : numerical::Mesh<T>(lengths, t0, tend, nx, nt), 
        m_discretization(disc){
        if(m_discretization != 0 && m_discretization != 1){
              throw std::invalid_argument(
                "Not a valid discretization scheme.");
        }
    }
    // TODO: Use DOF MAP instead
    /**
    * @return the list of global coordinates.
    */
    Vec vertices() {
        Vec v;
        if(m_discretization == 0){
            v = Vec::LinSpaced(this->m_nx, 
                this->m_lengths[0] + this->m_dx/2.0, 
                this->m_lengths[1] - this->m_dx/2.0);   
        } else{
            v = Vec::LinSpaced(this->m_nx + 1, 
                this->m_lengths[0], 
                this->m_lengths[1]);
        }
        return v;
    }
};

}  // namespace fvm

} // namespace numerical

#endif  // FVMESH_H
