/**
 *  @file    FDMesh.hpp
 *  @brief   The finite difference mesh for the hypercube.
 *  @author  Francois Roy
 *  @date    12/04/2019
 */
#ifndef FDMESH_H
#define FDMESH_H

#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "utils/Utils.hpp"

namespace numerical {

namespace fdm {

/**
* Holds data structures for a uniform mesh on a hypercube in
* space, plus a uniform mesh in time.
*
*  \anchor fig_mesh
*  \image html bloc.png <"Figure 1: A typical mesh for a bloc. The scene shows 
*  the bottom (orange), left (magenta) and front (blue) boundaries.">
*  \image latex bloc.eps "bloc" width=10cm
*
* @param lengths List of 2-lists of min and max coordinates in each spatial 
* direction.
* @param t0 The initial time in time mesh.
* @param tend The final time in time mesh.
* @param nt Number of cells in time mesh.
* @param n List of number of cells in the spatial directions.
*/
template <typename T>
class FDMesh {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
	const std::vector<std::vector<T>> m_lengths;
	const T m_t0, m_tend;
	T m_dt;
	const std::vector<int> m_nx;
	std::vector<T> m_dx;
	const int m_nt;
	std::vector<Vec> m_x;
	Vec m_t;
	std::size_t m_dim; 
public:
	FDMesh(const std::vector<std::vector<T>>& lengths, T t0, T tend, 
		const std::vector<int>& nx, int nt) 
	  : m_lengths(lengths),
	  m_t0(t0),
      m_tend(tend),
      m_nx(nx),
      m_nt(nt)
	{
        m_dim = m_lengths.size(); 
        if (m_dim < 1){
        	throw std::invalid_argument(
        		"At least one 2-lists of min and max coords must be given."
        		);
        }
        if (m_nx.size() != m_lengths.size()){
        	throw std::invalid_argument("n and lengths have different size.");
        }
        for(int i=0; i<m_dim; ++i){
        	if(m_lengths[i].size() != 2){
        		throw std::invalid_argument("length item is not a 2-list.");
        	}
        	m_dx.push_back(T (m_lengths[i][1] - m_lengths[i][0])/m_nx[i]);
        }
        m_dt = m_tend / T(nt);
        for(int i=0; i<m_dim; ++i){
        	Vec temp(m_nx[i] + 1);
        	for(int j=0; j<m_nx[i] + 1; ++j){
        		temp[j] = j * m_dx[i];
        	}
            m_x.push_back(temp);
        }
        m_t = Vec(m_nt + 1);
        for(int i=0; i<m_nt; ++i){
        	m_t[i] = i * m_dt;
        }
	}
    /**
    *  The indices of the nodes on the back boundary, i.e. at 
    *  z = length[2][1]. Only for 3D models.
    *
    * @return The indices of the nodes on the back boundary, 
    */
    std::vector<int> back(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> temp;
        for (int i=0; i<nx[0] + 1; ++i){
                for (int j=0; j<nx[1] + 1; ++j){
                    temp.push_back(i + j * (nx[0] + 1));
                }
            }
        return temp;
    }
    /**
    *  The indices of the nodes on the bottom boundary, i.e. at 
    *  y = lengths[1][0]. Only for 2D and 3D models.
    *
    * @return The indices of the nodes on the bottom boundary, 
    */
    std::vector<int> bottom(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> temp;
        if (m_dim == 2){
            for (int i=0; i<nx[0] + 1; ++i){
                temp.push_back(i);
            }
        }
        else {  // m_dim = 3
            for (int i=0; i<nx[0] + 1; ++i){
                for (int k=0; k<nx[2] + 1; ++k){
                    temp.push_back(i + k * (nx[0] + 1) * (nx[1] + 1));
                }
            }
        }
        return temp;
    }
    /**
    * @return The indices of the external boundaries.
    */
    std::vector<int> boundaries(){
        // points in 1D
        // lines in 2D
        // surfaces in 3D
        std::vector<int> nx = division_per_axis();
        std::vector<int> temp;
        if (m_dim == 1){
            temp = {0, nx[0]};
        }
        else if (m_dim == 2){
            for (int j=0; j<nx[1] + 1; ++j){
                temp.push_back(j * (nx[0] + 1)); // left
                temp.push_back(j * (nx[0] + 1) + nx[0]); // right
            }
            for (int i=0; i<nx[0] + 1; ++i){
                temp.push_back(i);  // bottom 
                temp.push_back(i + nx[1] * (nx[0] + 1));  // top 
            }
            // make unique
            std::sort(temp.begin(), temp.end());
            auto last = std::unique(temp.begin(), temp.end());
            temp.erase(last, temp.end()); 
        }
        else {  // m_dim = 3
            for (int j=0; j<nx[1] + 1; ++j){
                for (int k=0; k<nx[2] + 1; ++k){
                    temp.push_back(j * (nx[0] + 1) + 
                        k * ((nx[0] + 1)*(nx[1] + 1)));  // left
                    temp.push_back(j * (nx[0] + 1) + 
                        k * ((nx[0] + 1)*(nx[1] + 1)) + nx[0]); // right
                }
            }
            for (int i=0; i<nx[0] + 1; ++i){
                for (int k=0; k<nx[2] + 1; ++k){
                    temp.push_back(i + k * (nx[0] + 1) * 
                        (nx[1] + 1));  // bottom
                    temp.push_back(i + k * (nx[0] + 1) * 
                        (nx[1] + 1) + (nx[0] + 1) * (nx[0] + 1));  // top
                }
            }
            for (int i=0; i<nx[0] + 1; ++i){
                for (int j=0; j<nx[1] + 1; ++j){
                    temp.push_back(i + j * (nx[0] + 1));  // back
                    temp.push_back(i + j * (nx[0] + 1) + 
                        nx[2] * (nx[0] + 1) * (nx[1] + 1));   // front
                }
            }
            // make unique
            std::sort(temp.begin(), temp.end());
            auto last = std::unique(temp.begin(), temp.end());
            temp.erase(last, temp.end()); 
        }
        return temp;
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
    std::vector<Eigen::Matrix<T, 3, 1>> coordinates() {
        // get the number of divisions per axis
        int n_x = m_nx[0];
        int n_y = 0;
        int n_z = 0;
        Vec y = Vec::Zero(n_x + 1);
        Vec z = Vec::Zero(n_x + 1); 
        if (m_dim >= 2){
            n_y = m_nx[1];
            y = m_x[1];
            z = Vec::Zero((n_x + 1)*(n_y + 1)); 
        } 
        if (m_dim == 3) {
            n_z = m_nx[2];
            z = m_x[2];
        }
        // generate the vector of coordinates
        Eigen::Matrix<T, 3, 1> node;
        std::vector<Eigen::Matrix<T, 3, 1>> coords;
        for (int k=0; k<n_z + 1; ++k){
            for (int j=0; j<n_y + 1; ++j) {
                for (int i=0; i<n_x + 1; ++i){
                    node = {m_x[0][i], y[j], z[k]};
                    // spdlog::info("x: {} y: {} z: {}", m_x[0][i], y[j], z[k]);
                    coords.push_back(node);
                }
            }
        }
        return coords;
    }
    /**
    *  1D, 2D, 3D
    * @return The indices of the corners. 
    */
    std::vector<int> corners(){
        // TODO define in 3D
        // 2D: 0, nx, ny * (nx+1), (ny+1) * (nx+1) -1
        return {0};
    }
    /**
    * @return the space dimension.
    */
    std::size_t dim(){
        return m_dim;
    }
    /**
    * If 1D, \c ny\c and \c nz\c = 0, if 2D \c nz\c = 0, else all are > 1. See 
    * \ref fig_mesh
    *
    * @return The number of divisions per axis.
    */
    std::vector<int> division_per_axis(){
        int n_x = m_nx[0];
        int n_y = 0;
        int n_z = 0;
        if (m_dim >= 2){
            n_y = m_nx[1];
        }
        if (m_dim == 3) {
            n_z = m_nx[2];
        }
        return {n_x, n_y, n_z};
    }
    /**
    * @return The indices of the domain (all). 
    */
    std::vector<int> domain(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> all;
        if (m_dim == 1){
            all == utils::linear_spaced<int>(0, nx[0], nx[0] + 1);
        }
        else if (m_dim ==2){
            all == utils::linear_spaced<int>(0, (nx[0] + 1) * (nx[1] + 1) - 1, 
                (nx[0] + 1) * (nx[1] + 1));
        } else {
            all == utils::linear_spaced<int>(0, 
                (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1) - 1, 
                (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1));
        }
        return all;
    }
    /**
    * @return the constant time increment.
    */
    T dt() {
        return m_dt;
    }
    /**
    * @return the space increment for each dimension.
    */
    std::vector<T> dx() {
        return m_dx;
    }
    /**
    *  2D, 3D 
    * @return The nodes of the edges.
    */
    std::vector<int> edges(){
        // TODO define in 3D
        return {0};
    }
    /**
    *  The indices of the nodes on the front, i.e. at 
    *  z = length[2][0]. Only for 3D models.
    *
    * @return The indices of the nodes on the front boundary,  
    */
    std::vector<int> front(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> out = back();
        for(int& i : out){
            i += nx[2] * (nx[0] + 1) * (nx[1] + 1);
        }
        return out;
    }
    /**
    *  @return The interior nodes.
    */
    std::vector<int> interior(){
        // TODO
        return {0};
    }
    /**
    *  The indices of the nodes on the left boundary, 
    *  i.e. at x = lengths[0][0].
    *
    * @return The indices of the nodes on the left boundary. 
    */
    std::vector<int> left(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> temp;
        if (m_dim == 1){
            temp = {0};
        }
        else if (m_dim == 2){
            for (int j=0; j<nx[1] + 1; ++j){
                temp.push_back(j * (nx[0] + 1));
            }
        }
        else {  // m_dim = 3
            for (int j=0; j<nx[1] + 1; ++j){
                for (int k=0; k<nx[2] + 1; ++k){
                    temp.push_back(j * (nx[0] + 1) + 
                        k * ((nx[0] + 1)*(nx[1] + 1)));
                }
            }
        }
        return temp;
    }
    /**
    *  The indices of the nodes on the right boundary, i.e. at 
    *  x = lengths[0][1]. It is just left + n_x.
    *
    * @return The indices of the nodes on the right boundary, 
    */
    std::vector<int> right(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> out = left();
        for(int& i : out){
            i += nx[0];
        }
    	return out;
    }
    /**
    *  Returns the boundary nodes in 3D.
    *  @return The surfaces. 
    */
    std::vector<int> surfaces(){
        std::vector<int> out;
        if (m_dim == 3){
            out = boundaries();
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
    *  The indices of the nodes on the top boundary i.e. at 
    *  y = lengths[1][1]. Only for 2D and 3D models.
    *
    * @return The indices of the nodes on the top boundary, 
    */
    std::vector<int> top(){
        std::vector<int> nx = division_per_axis();
        std::vector<int> out = bottom();
        if (m_dim == 2){
            for(int& i : out){
                i += nx[1] * (nx[0] + 1);
            }
        }
        else { // 3D
            for(int& i : out){
                i += (nx[0] + 1) * (nx[0] + 1);
            }
        }
        return out;
    }
	/**
	* @return the list of coordinates in each dimensions.
	*/
	std::vector<Vec> x() {
		return m_x;
	}
};

}  // namespace fdm

} // namespace numerical

#endif  // MESH_H
