/**
 *  @file    FEProblem.hpp
 *  @brief   Defines basis functions.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#ifndef FEPROBLEM_H
#define FEPROBLEM_H

#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "FEMesh.hpp"

namespace numerical {
namespace fem {
/**
* Defines the finite element problem over a line.
*
* This class only define the 1D diffusion problem with Dirichlet/Neumann 
* boundary conditions and heterogenous diffusion coefficient (for now).
*
* Support for hyperbolic equations (wave and convection diffusion) will be 
* added later.
*/
template <typename T>
class FEProblem {  // TODO inherit from common/Problem
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
typedef std::vector<Eigen::Matrix<T, 3, 1>> Coord;
// function pointer to member functions
typedef  T (FEProblem<T>::*fctptr)(T x1, T x2, T time);
private:
  Vec m_t;
  std::vector<T> m_dx;
  T m_dt;
protected:
  FEMesh<T>* m_mesh;  // mesh
  // boundary types (0=Dirichlet, 1=Neumann) for left and right boundaries.
  int m_bc_types[2];   
public:
  FEProblem():
    m_bc_types{0, 0} {
	// define mesh
	// m_mesh = new Mesh<T>();
  }
  virtual ~FEProblem(){
    // delete m_mesh;
  }
  /**
  *  Domain LHS (bilinear form)
  */
  virtual void dlhs(T e, T phi, T r, T s, T X, T x, T h){

  }
  /**
  *  Domain RHS (linear form)
  */
  virtual void drhs(T e, T phi, T r, T X, T x, T h){
  	
  }
  /**
  *  Boundary RHS (bilinear form)
  */
  virtual void blhs(T e, T phi, T r, T s, T X, T x, T h){
  	
  }
  /**
  *  Boundary LHS (linear form)
  */
  virtual void brhs(T e, T phi, T r, T X, T x, T h){
  	
  }
  /**
  *  Get the reference solution at a specified mesh location.
  *
  * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
  *    node.
  * @param t The discrete time.
  * @return The reference solution at a specified mesh location.
  */
  virtual T reference(const Eigen::Matrix<T, 3, 1>& x, T t){
    return 0.0;
  }
  /**
  * User defined function: right boundary.
  *
  * Right boundary value, i.e. for x = lengths[0][1].
  *
  * @param u The state variable.
  * @param t The discrete time.
  * @return The boundary value at a specified mesh location.
  */
  virtual T right(T u, T t){
    return 0.0;
  }
  /**
  * User defined function: left boundary.
  *
  * Left boundary value, i.e. for x = lengths[0][1].
  *
  * @param u The state variable.
  * @param t The discrete time.
  * @return The boundary value at a specified mesh location.
  */
  virtual T left(T u, T t){
    return 0.0;
  }
  /**
  * User defined function: source term.
  *
  * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
  *    node.
  * @param t The discrete time.
  * @param u The state variable for nonlinear problems.
  * @return The source term at a specified mesh location.
  */
  virtual T source(const Eigen::Matrix<T, 3, 1>& x, T t, T u=0){
    return 0.0;
  }
};


}  // namespace fem
} // namespace numerical

#endif  // FEMPROBLEM_H