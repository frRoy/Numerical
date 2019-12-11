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
typedef Eigen::ArrayXd Arr;
private:
    int m_dim;
protected:
    Parameters<T>* m_params;
    Vec m_u;  // The solution vector
public:
    Problem(Parameters<T>* params): m_params(params) {
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
    * Left boundary value, on the left node (1D) or on a line parallel to the
    * y-axis (2D) or on a plane parallel to y0z (3D).
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
    * Right boundary value, on the right node (1D) or on a line parallel to the
    * y-axis (2D) or on a plane parallel to y0z (3D).
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
    * Bottom boundary value, on a line parallel to the x-axis (2D) or on a 
    * plane parallel to x0z (3D).
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
    * Top boundary value, on a line parallel to the x-axis (2D) or on a 
    * plane parallel to x0z (3D).
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
    * Front boundary value, on a plane parallel to x0y (3D).
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
    * Back boundary value, on a plane parallel to x0y (3D).
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
    * Neumann boundary value.
    *
    * @param x1 The first coordinate of the mesh node.
    * @param x2 The second coordinate of the mesh node.
    * @param t The discrete time.
    * @return The Neumann boundary value at a specified mesh location.
    */
    virtual T neumann(T x1, T x2, T t){
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
    virtual T reference(T x1, T x2, T t){
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
    virtual T solution(T x1, T x2, T t){
    	return 0.0;
    }

    /**
    *  Solve the finite difference problem.
    */
    virtual void solve(){
    	// sets m_u here
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
