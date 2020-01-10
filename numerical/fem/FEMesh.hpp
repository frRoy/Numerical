/**
 *  @file    FEMesh.hpp
 *  @brief   The finite element mesh for a line.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#ifndef FEMESH_H
#define FEMESH_H

#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "utils/Utils.hpp"

namespace numerical {

namespace fem {

/**
* Holds data structures for a uniform mesh on a line in space, plus a 
* uniform mesh in time.
* 
* The finite element mesh has n elements and an associated polynomial of 
* degree d.
*
* @param lengths A 2-lists of min and max coordinates in the spatial direction.
* @param t0 The initial time in time mesh.
* @param tend The final time in time mesh.
* @param nt Number of cells in time mesh.
* @param n Number of cells in the spatial direction.
* @param d The polynomial degree, default = 1.
*/
template <typename T>
class FEMesh {  // TODO inherit form common/Mesh
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
	const std::vector<T> m_lengths;
	const T m_t0, m_tend;
    const int m_nt, m_nx, m_d;
	T m_dt, m_dx;
	Vec m_x;
	Vec m_t;
public:
	FEMesh(const std::vector<T>& lengths, T t0, T tend, int nx, int nt, 
        int d=1) 
	  : m_lengths(lengths),
	  m_t0(t0),
      m_tend(tend),
      m_nx(nx),
      m_nt(nt),
      m_d(d)
	{
      	if(m_lengths.size() != 2){
      		throw std::invalid_argument("length is not a 2-list.");
       	}
       	m_dx = (m_lengths[1] - m_lengths[0])/ T(m_nx);
        if(m_dx <= 0){
            throw std::invalid_argument("not a valid interval.");
        }
        m_dt = m_tend / T(m_nt);
        m_x = Vec::LinSpaced(m_nx + 1, m_lengths[0], m_lengths[1]);
        m_t = Vec::LinSpaced(m_nt + 1, m_t0, m_tend);
	}
    /**
    *  Local vertex to global vertex mapping (cells).
    */
    std::vector<std::vector<int>> cells(){
        std::vector<std::vector<int>> out;
        for(int i=0; i<m_nx; i++){
            out.push_back({i, i+1});
        }
        return out;
    }
    /**
    *  Local to global degree of freedom mapping (dof_map).
    */
    std::vector<std::vector<int>> dof_map(){
        std::vector<std::vector<int>> out;
        std::vector<int> m;
        if(m_d == 0){
            for(int i; i<m_nx; i++){
                out.push_back({i});
            }
        } else {
            for(int i; i<m_nx; i++){
                for(int j; j<m_d+1; j++){
                    m.push_back(i*m_d + j);
                }
                out.push_back(m);
            }
        }
        return out;
    }
    /**
    * @return the constant time increment.
    */
    T dt() {
        return m_dt;
    }
    /**
    * @return the space increment.
    */
    T dx() {
        return m_dx;
    }
    /**
    *  The indices of the nodes on the left boundary.
    *
    * @return The indices of the nodes on the left boundary. 
    */
    std::vector<int> left(){
        return {0};
    }
    /**
    *  The indices of the nodes on the right boundary.
    *
    * @return The indices of the nodes on the right boundary, 
    */
    std::vector<int> right(){
        std::vector<int> out = left();
        for(int& i : out){
            i += m_nx;
        }
    	return out;
    }
    /**
    * @return the time list.
    */
    Vec t() {
        return m_t;
    }
	/**
	* @return the list of coordinates.
	*/
	Vec vertices() {
		return m_x;
	}
};

}  // namespace fem

} // namespace numerical

#endif  // MESH_H
