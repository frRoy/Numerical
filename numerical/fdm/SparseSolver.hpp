/**
 *  @file    SparseSolver.hpp
 *  @brief   Solves a finite difference problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef SPARSESOLVER_H
#define SPARSESOLVER_H

#include <vector>
#include <math.h> 
#include <Eigen/SparseCore>
#include<Eigen/SparseCholesky>
#include "spdlog/spdlog.h"
#include "Problem.hpp"
#include <iostream>

namespace numerical {

namespace fdm {

/**
 * This class only computes the diffusion problem with Dirichlet boundary 
 * conditions and constant diffusion coefficient (for now).
 */
template <typename T>
class SparseSolver {
typedef Eigen::SparseMatrix<T> SpMat;
typedef Eigen::Triplet<T> Trip;
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
typedef std::vector<Eigen::Matrix<T, 3, 1>> Coord;
private:
  SpMat m_A;
  Problem<T>* m_problem;
  int m_n, m_n_x, m_n_y, m_n_z, m_n_t, m_dim;
  T m_dt, m_dx, m_dy, m_dz, m_theta;
  Vec m_alpha, m_x, m_y, m_z, m_b, m_u, m_f, m_f_n;
public:
  SparseSolver(Problem<T>* problem)
    : m_problem(problem)
    {  
        m_n = m_problem->n();
        m_dim = m_problem->dim();
        m_A = SpMat(m_n, m_n);
        m_b = Vec::Zero(m_n);
        m_u = m_problem->u_0();
        T dx, dy, dz, dt, theta;
        theta = m_problem->theta();
        m_n_x = m_problem->n_x();
        m_n_y = m_problem->n_y();
        m_n_z = m_problem->n_z();
        m_n_t = m_problem->n_t();
        m_dt = m_problem->dt();
        m_dx = m_problem->dx()[0];
        m_dy = m_problem->dx()[1];
        m_dz = m_problem->dx()[2];
        m_theta = m_problem->theta();
        m_alpha = Vec::Zero(m_n);
        m_f = Vec::Zero(m_n);
        m_f_n = Vec::Zero(m_n);
        m_x = Vec::Zero(m_n);
        m_y = Vec::Zero(m_n);
        m_z = Vec::Zero(m_n);
        const Coord& coords = m_problem->coordinates();
        // define x, y, z, and alpha
        for(int i=0; i<m_n; i++) {
            m_x[i] = coords[i][0];
            m_y[i] = coords[i][1];
            m_z[i] = coords[i][2];
            m_alpha[i] = m_problem->alpha(coords[i], 0.0);
        }
  }
  ~SparseSolver(){
  }
  /**
  * Assembles the sparse coefficient matrix \f$\mathbf{A}\f$. The matrix is 
  * defined by its diagonals. The number of diagonals is related to the
  * number of direct neighbors for interior mesh nodes. In 1D, the matrix has 
  * 3 non-zero diagonals, in 2D it has 5, and in 3D, 7. 
  *
  * The diagonals of the matrix \f$\mathbf{A}\f$ is are filled by vectorization
  * of the loops for efficiency.
  */
  void assemble_a(){
    std::vector<Trip> trp;
    Vec diagonal = Vec::Constant(m_n, 1.0);
    Vec lower = Vec::Zero(m_n - 1);
    Vec upper = Vec::Zero(m_n - 1);
    // The loops are vectorized for efficiency -- see bench/performances
    if (m_dim == 1){
        // spdlog::info("1D");
        T d = m_dt/m_dx/m_dx*m_theta/2.0;
        spdlog::info("dx: {}, dt: {}, theta: {}, alpha: {}", 
                      m_dx, m_dt, m_theta, m_alpha[0]);
        spdlog::info("Fx: {}", d*m_alpha[0]*2.0);
        diagonal[0] = 0.0;
        diagonal[m_n - 1] = 0.0;
        diagonal.segment(1, m_n-2) += d * m_alpha.segment(2, m_n-2);
        diagonal.segment(1, m_n-2) += d * 2.0 * m_alpha.segment(1, m_n-2);
        diagonal.segment(1, m_n-2) += d * m_alpha.segment(0, m_n-2);
        lower.segment(0, (m_n-1)-1) += -d * m_alpha.segment(1, m_n-2);
        lower.segment(0, (m_n-1)-1) += -d * m_alpha.segment(0, m_n-2);
        upper.segment(1, (m_n-1)-1) += -d * m_alpha.segment(2, m_n-2);
        upper.segment(1, (m_n-1)-1) += -d * m_alpha.segment(1, m_n-2);
        // boundary conditions
        if(m_problem->bc_type(0) == 0){ // left Dirichlet
            diagonal[0] = 1.0;
            upper[0] = 0.0;
          } else{  // left Neumann
            diagonal[0] = d * 2.0 * m_alpha[0];
            upper[0] = -d * 2.0 * m_alpha[0];
          }
        if(m_problem->bc_type(1) == 0){ // right Dirichlet
            diagonal[m_n-1] = 1.0;
            lower[(m_n - 1)-1] = 0.0;
          } else {  // right Neumann
            diagonal[m_n-1] = d * 2.0 * m_alpha[m_n-1];
            lower[(m_n - 1)-1] = -d * 2.0 * m_alpha[m_n-1];
          }
        // insert diagonals in A
        for(int i=0; i<m_n; i++){
            trp.push_back(Trip(i,i,diagonal[i]));    
        }
        for(int i=1; i<m_n - 1; i++){
            trp.push_back(Trip(i,i-1,lower[i-1]));    
        }
        for(int i=1; i<m_n - 1; i++){
            trp.push_back(Trip(i,i+1,upper[i]));    
        }
        // create sparse matrix
        m_A.setFromTriplets(trp.begin(), trp.end());
    } else if (m_dim == 2){
        // spdlog::info("2D");
        // TODO
    } else {  // 3D
        // spdlog::info("3D");
        // TODO
    }
  }
  /**
  * Assembles the RHS vector \f$\mathbf{b}\f$.
  *
  * \f[
  *    b_i = u_i^n + F\left(1-\Theta\right)u_{i+1}^n-2u_i^n+u_{i-1}^n +
  *        \Delta t \Theta f_i^{n+1} + \Delta t \left(1-\Theta\right)f_i^n
  * \f]
  *
  * using vectorization we get:
  *
  * \f[
  *    b[1:n_x-1] = u_n[1:n_x-1] + \left(1-\Theta\right)F
  *        \left(u_n[2:n_x]-2u_n[1:n_x-1]+u_n[0:n_x-2]\right) + 
  *        \Theta\Delta t f[1:n_x-1](n+1) + 
  *        \left(1-\Theta\right)\Delta t f[1:n_x-1](n)
  * \f]
  *
  */
  void assemble_b(T t, Vec& u_n){
      const Coord& coords = m_problem->coordinates();
      // TODO define the diffusion coefficient and source term only if they
      // depend on time
      for(int i=0; i<m_n; i++) {
          m_alpha[i] = m_problem->alpha(coords[i], t);
          m_f_n[i] = m_problem->source(coords[i], t);
          m_f[i] = m_problem->source(coords[i], t + m_n_t * m_dt);
      }
      if (m_dim == 1){
          // spdlog::info("1D");
          T d = m_dt/m_dx/m_dx*(1.0 - m_theta) / 2.0;
          // spdlog::info("d: {}", d);
          m_b.segment(1, m_n-2) = u_n.segment(1, m_n-2);
          m_b.segment(1, m_n-2) += d * ((m_alpha.segment(2, m_n-2) + 
            m_alpha.segment(1, m_n-2)).array() * (u_n.segment(2, m_n-2) - 
            u_n.segment(1, m_n-2)).array()).matrix();
          m_b.segment(1, m_n-2) -= d * ((m_alpha.segment(1, m_n-2) + 
            m_alpha.segment(0, m_n-2)).array() * (u_n.segment(1, m_n-2) - 
            u_n.segment(0, m_n-2)).array()).matrix();
          m_b.segment(1, m_n-2) += m_dt * m_theta * m_f.segment(1, m_n-2);
          m_b.segment(1, m_n-2) += m_dt * (1.0 -m_theta) * 
            m_f_n.segment(1, m_n-2);
          // Boundary conditions
          m_problem->rhs_bc(m_b, t);
      } else if (m_dim == 2){
        // spdlog::info("2D");
        // TODO
      } else {  // 3D
        // spdlog::info("3D");
        // TODO
      }
  }
  /*
  * Solve the time dependent problem.
  * TODO: Create a VTKFile class to store the solution in a vtu file at each 
  * time steps
  */
  virtual void solve(){
      // Set initial condition
      Vec u_n = m_problem->u_0();
      Vec u = Vec::Zero(m_problem->n());
      int n_t = m_problem->n_t();
      T t, l2norm=0.0;
      Vec t_list = m_problem->t();
      assemble_a();
      // std::cout << m_A;
      Eigen::SimplicialLDLT<SpMat> solver;
      // Time loop
      T e=0.0;
      for(int n=0; n<n_t; n++){
        t = t_list[n]; 
        assemble_b(t, u_n);
        // Solve
        solver.compute(m_A);
        if(solver.info()!=Eigen::Success) {
            spdlog::error("decomposition failed");
            spdlog::error("{}", solver.info());
            return;
        }
        u = solver.solve(m_b);
        if(solver.info()!=Eigen::Success) {
            spdlog::error("solving failed");
            spdlog::error("{}", solver.info());
            return;
        }
        // spdlog::info("b[1]: {}, b[n-2]: {}", m_b[1], m_b[m_n-2]);
        // TODO Save result in file here.
        const Coord& coords = m_problem->coordinates();
        //spdlog::info("exact solution: {}, computed solution: {}", 
        //    m_problem->reference(coords[4], t+m_dt), u[4]);
        
        for(int i=0;i<m_n;i++){
            e += pow(m_problem->reference(coords[i], t+m_dt)-u[i], 2.0); 
        }
        u_n = u;
      }
      l2norm += pow(m_dx*m_dt*e, 0.5);
      spdlog::info("L2-norm: {}", l2norm);
  }
};

}  // namespace fdm

}  // namespace numerical

#endif  // SPARSESOLVER_H
