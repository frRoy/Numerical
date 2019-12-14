/**
 *  @file    Problem.hpp
 *  @brief   Define a finite difference problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef PROBLEM_H
#define PROBLEM_H

#include <map>
#include <string>
#include <vector>
#include <Eigen/SparseCore>
#include "spdlog/spdlog.h"
#include "Parameters.hpp"
#include "Mesh.hpp"

namespace numerical {

namespace fdm {

/*
 * This class defines the finite difference problem of the diffusion type:
 *
 * /f[
 *    \frac{\partial u}{\partial t} = \nabla\left(\alpha\nabla u) + f
 * /f]
 *
 * over the hypercube.
 */
template <typename T>
class Problem {
typedef Eigen::SparseMatrix<T> SpMat;
typedef Eigen::Triplet<T> Trip;
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
private:
    int m_dim;
    
protected:
    Parameters<T>* m_params;
    Vec m_u;  // The solution vector
public:
    Problem(Parameters<T>* params): 
      m_params(params) {
        // define other variables variables
        m_dim = m_params->lengths.size();
  	    // define mesh
  	    Mesh<T> mesh = Mesh<T>(m_params->lengths, m_params->t0, m_params->tend, 
          m_params->n, m_params->nt);
    }
    virtual ~Problem(){
    }
    /**
    * Diffusion coefficient value. The default is a constant obtained from 
    * m_params.
    *
    * @param x The x-coordinate of the mesh node.
    * @param y The y-coordinate of the mesh node.
    * @param z The z-coordinate of the mesh node.
    * @param t The discrete time.
    * @return The diffusion coefficient value at a specified mesh location.
    */
    virtual T alpha(T x, T y, T z, T t){
  	    return m_params->alpha;
    }

    /**
    * @return The spatial dimension of the problem.
    */
    int dim(){
    	return m_dim;
    }

    /**
    * Left boundary value, i.e. for x = lengths[0][0].
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T left(int type, T x1, T x2, T t){
  	    return 0.0;
    }

    /**
    * Right boundary value, i.e. for x = lengths[0][1].
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T right(int type, T x1, T x2, T t){
        return 0.0;
    }

    /**
    * Bottom boundary value, i.e. at y = lengths[1][0]. Only for 2D and 
    * 3D models.
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T bottom(int type, T x1, T x2, T t){
        return 0.0;
    }

    /**
    * Top boundary value, i.e. at y = lengths[1][1]. Only for 2D and 3D models.
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T top(int type, T x1, T x2, T t){
        return 0.0;
    }

    /**
    * Front boundary value, i.e. at cz = length[2][0]. Only for 3D models.
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T front(int type, T x1, T x2, T t){
        return 0.0;
    }

    /**
    * Back boundary value, i.e. at z = length[2][0]. Only for 3D models.
    *
    * @param type Dirichlet = 0, Neumann = 1.
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T back(int type, T x1, T x2, T t){
        return 0.0;
    }

    /**
    * Initial  value.
    *
    * @param x The x-coordinate of the mesh node.
    * @param y The y-coordinate of the mesh node.
    * @param z The z-coordinate of the mesh node.
    * @return The initial value at a specified mesh location.
    */
    virtual T initial_value(T x, T y, T z){
  	    return 0.0;
    }

    /**
    *  Get the reference solution at a specified mesh location.
    *
    * @param x The x-coordinate of the mesh node.
    * @param y The y-coordinate of the mesh node.
    * @param z The z-coordinate of the mesh node.
    * @param t The discrete time.
    * @return The reference solution at a specified mesh location.
    */
    virtual T reference(T x, T y, T z, T t){
    	return 0.0;
    }

    /**
    *  Get the computed solution at a specified mesh location.
    *
    * @param x The x-coordinate of the mesh node.
    * @param y The y-coordinate of the mesh node.
    * @param z The z-coordinate of the mesh node.
    * @param t The discrete time.
    * @return The computed solution at a specified mesh location.
    */
    virtual T solution(T x, T y, T z, T t){
    	return 0.0;
    }

    /**
    * Source term.
    *
    * @param x The x-coordinate of the mesh node.
    * @param y The y-coordinate of the mesh node.
    * @param z The z-coordinate of the mesh node.
    * @param t The discrete time.
    * @return The source term at a specified mesh location.
    */
    virtual T source(T x, T y, T z, T t){
  	    return 0.0;
    }

};

}  // namespace fdm

}  // namespace numerical

#endif  // PROBLEM_H
