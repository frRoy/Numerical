/**
 *  @file    Mesh.hpp
 *  @brief   Definition of the finite difference mesh.
 *  @author  Francois Roy
 *  @date    12/04/2019
 */
#ifndef MESH_H
#define MESH_H

#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include "spdlog/spdlog.h"

namespace numerical {

namespace fdm {

typedef Eigen::VectorXd Vec;

/**
* Holds data structures for a uniform mesh on a hypercube in
* space, plus a uniform mesh in time.
*
* @param lengths List of 2-lists of min and max coordinates in each spatial 
* direction.
* @param t0 The initial time in time mesh.
* @param tend The final time in time mesh.
* @param nt Number of cells in time mesh.
* @param n List of number of cells in the spatial directions.
*/
template <typename T>
class Mesh {
private:
	const std::vector<std::vector<T>> m_lengths;
	const T m_t0, m_tend;
	T m_dt;
	const std::vector<int> m_n;
	std::vector<T> m_dx;
	const int m_nt;
	std::vector<Vec> m_x;
	Vec m_t;
	std::size_t m_dim; 
public:
	Mesh(const std::vector<std::vector<T>>& lengths, T t0, T tend, 
		const std::vector<int>& n, int nt) 
	  : m_lengths(lengths),
	  m_t0(t0),
      m_tend(tend),
      m_n(n),
      m_nt(nt)
	{
        spdlog::info("param tend: {}", m_tend);
        m_dim = m_lengths.size(); 
        if (m_dim < 1){
        	throw std::invalid_argument(
        		"At least one 2-lists of min and max coords must be given."
        		);
        }
        if (m_n.size() != m_lengths.size()){
        	throw std::invalid_argument("n and lengths have different size.");
        }
        for(int i=0; i<m_dim; ++i){
        	if(m_lengths[i].size() != 2){
        		throw std::invalid_argument("length item is not a 2-list.");
        	}
        	m_dx.push_back(T (m_lengths[i][1] - m_lengths[i][0])/m_n[i]);
        }
        m_dt = m_tend / T(nt);
        for(int i=0; i<m_dim; ++i){
        	Vec temp(m_n[i] + 1);
        	for(int j=0; j<m_n[i]; ++j){
        		temp[j] = j * m_dx[i];
        	}
            m_x.push_back(temp);
        }
        m_t = Vec(m_nt + 1);
        for(int i=0; i<m_nt; ++i){
        	m_t[i] = i * m_nt;
        }
	}
    
    /**
    *  The indices of the nodes on the left boundary (1D, 2D, 3D).
    */
    std::vector<int> left_indices(){
    	return {0};
    }

    /**
    *  The indices of the nodes on the right boundary (1D, 2D, 3D.
    */
    std::vector<int> right_indices(){
    	return {0};
    }

    /**
    *  The indices of the nodes on the bottom boundary (2D, 3D.
    */
    std::vector<int> bottom_indices(){
    	return {0};
    }

    /**
    *  The indices of the nodes on the top boundary (2D, 3D).
    */
    std::vector<int> top_indices(){
    	return {0};
    }

    /**
    *  The indices of the nodes on the front boundary (3D).
    */
    std::vector<int> front_indices(){
    	return {0};
    }

    /**
    *  The indices of the nodes on the back boundary (3D).
    */
    std::vector<int> back_indices(){
    	return {0};
    }

	/**
	* @return the space dimension.
	*/
	std::size_t dim(){
		return m_dim;
	}
	/**
	* @return the space increment for each dimension.
	*/
	std::vector<T> dx() {
		return m_dx;
	}
	/**
	* @return the time increment.
	*/
	T dt() {
		return m_dt;
	}
	/**
	* @return the list of coordinates in each dimensions.
	*/
	std::vector<Vec> x() {
		return m_x;
	}
	/**
	* @return the time list.
	*/
	Vec t() {
		return m_t;
	}
};


}  // namespace fdm

} // namespace numerical

#endif  // MESH_H
