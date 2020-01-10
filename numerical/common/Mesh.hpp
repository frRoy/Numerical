/**
 *  @file    Mesh.hpp
 *  @brief   The mesh base class for a line.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#ifndef MESH_H
#define MESH_H

#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "utils/Utils.hpp"

namespace numerical {

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
*/
template <typename T>
class Mesh {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
typedef  T (Mesh<T>::*fctptr)(T x);
typedef std::vector<fctptr> Basis;
protected:
	const std::vector<T> m_lengths;
	const T m_t0, m_tend;
    const int m_nt, m_nx;
	T m_dt, m_dx;
	Vec m_x;
	Vec m_t;
public:
	Mesh(const std::vector<T>& lengths, T t0, T tend, int nx, int nt) 
	  : m_lengths(lengths),
	  m_t0(t0),
      m_tend(tend),
      m_nx(nx),
      m_nt(nt)
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
    * Local vertex to global vertex mapping (cells).
    */
    virtual std::vector<std::vector<int>> cells(){
        std::vector<std::vector<int>> out;
        for(int i=0; i<m_nx; i++){
            out.push_back({i, i+1});
        }
        return out;
    }
    /**
    * Defines a mapping \f$m(i, j, k)\f$ from a mesh point with indices 
    * \f$(i, j, k)\f$ to the corresponding unknown index \f$p\f$ in the 
    * equation system:
    *
    * \f[
    *  p = m(i, j, k) = i+ j(n_x + 1) + k((n_x + 1)(n_y+1)))
    * \f]
    *
    * where \f$n_x\f$  and \f$n_y\f$ are the number of division along the 
    * \f$x\f$- and \f$y\f$-directions.
    *
    * We number the points along the x axis, starting with y = y_0, and 
    * z = z_0. We then progress one mesh line at a time 
    * (from y = y_0 to y = y_end) until the slice 
    * located at z = z_0 has been processed. We then continue on the next 
    * slice until the entire meshed cube has been mapped.
    *
    * The coordinates are returned in a vector of 3D Eigen::Vector.
    * @return The mesh node coordinates.
    */
    virtual std::vector<Eigen::Matrix<T, 3, 1>> coordinates() {
        // get the number of divisions per axis
        int n_x = m_nx;
        int n_y = 0;
        int n_z = 0;
        Vec y = Vec::Zero(n_x + 1);
        Vec z = Vec::Zero(n_x + 1); 
        // generate the vector of coordinates
        Eigen::Matrix<T, 3, 1> node;
        std::vector<Eigen::Matrix<T, 3, 1>> coords;
        for (int k=0; k<n_z + 1; ++k){
            for (int j=0; j<n_y + 1; ++j) {
                for (int i=0; i<n_x + 1; ++i){
                    node = {m_x[i], y[j], z[k]};
                    // spdlog::info("x: {} y: {} z: {}", m_x[i], y[j], z[k]);
                    coords.push_back(node);
                }
            }
        }
        return coords;
    }
    /**
    * @return the constant time increment.
    */
    virtual T dt() {
        return m_dt;
    }
    /**
    * @return the space increment.
    */
    virtual T dx() {
        return m_dx;
    }
    /**
    * The indices of the nodes on the left boundary.
    *
    * @return The indices of the nodes on the left boundary. 
    */
    virtual std::vector<int> left(T u=0){
        return {0};
    }
    /**
    * The indices of the nodes on the right boundary.
    *
    * @return The indices of the nodes on the right boundary, 
    */
    virtual std::vector<int> right(T u=0){
        std::vector<int> out = left();
        for(int& i : out){
            i += m_nx;
        }
    	return out;
    }
    /**
    * @return the time list.
    */
    virtual Vec t() {
        return m_t;
    }
	/**
	* @return the list of global coordinates.
	*/
	virtual Vec vertices() {
		return m_x;
	}
};

} // namespace numerical

#endif  // MESH_H
