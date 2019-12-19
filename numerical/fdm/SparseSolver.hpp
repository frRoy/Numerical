/**
 *  @file    SparseSolver.hpp
 *  @brief   Solves a finite difference problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef SPARSESOLVER_H
#define SPARSESOLVER_H

#include <vector>
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
  Vec m_b;
  Vec m_u;
  Problem<T>* m_problem;
public:
  SparseSolver(Problem<T>* problem)
    : m_problem(problem)
    {
        int n = m_problem->n();
        m_A = SpMat(n, n);
        m_b = Vec(n);
        m_b = Eigen::VectorXd::Zero(n);
        m_u = m_problem->u_0();
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
  void assemble_a(T t=0.0){
    // define variables
    int n, n_x, n_y, n_z, n_t;
    T dx, dy, dz, dt, theta;
    theta = m_problem->theta();
    n_x = m_problem->n_x();
    n_y = m_problem->n_y();
    n_z = m_problem->n_z();
    n = m_problem->n();
    n_t = m_problem->n_t();
    dt = m_problem->dt();
    dx = m_problem->dx()[0];
    dy = m_problem->dx()[1];
    dz = m_problem->dx()[2];
    Vec alpha(n), x(n), y(n), z(n), ix(n_x  +1), iy(n_y + 1), iz(n_z + 1),
        it(n_t + 1);
    // index sets
    ix.setLinSpaced(n_x + 1, 0, n_x);
    Coord coords = m_problem->coordinates();
    // define the diffusion coefficient, x, y and z
    for(int i=0; i<alpha.size(); i++) {
        x[i] = coords[i][0];
        y[i] = coords[i][1];
        z[i] = coords[i][2];
        alpha[i] = m_problem->alpha(coords[i], t);
    }
    std::vector<Trip> trp;
    Vec diagonal = Eigen::VectorXd::Constant(n, 1.0);
    Vec lower = Eigen::VectorXd::Zero(n - 1);
    Vec upper = Eigen::VectorXd::Zero(n - 1);
    // The loops are vectorized for efficiency -- see bench/performances
    if (m_problem->dim() == 1){
        spdlog::info("1D");
        T d = dt/dx/dx*theta/2.0;
        diagonal[0] = 0.0;
        diagonal[n - 1] = 0.0;
        diagonal.segment(1, n-2) += d * alpha.segment(2, n-2);
        diagonal.segment(1, n-2) += d * 2.0 * alpha.segment(1, n-2);
        diagonal.segment(1, n-2) += d * alpha.segment(0, n-2);
        lower.segment(0, (n-1)-1) += -d * alpha.segment(1, n-2);
        lower.segment(0, (n-1)-1) += -d * alpha.segment(0, n-2);
        upper.segment(1, (n-1)-1) += -d * alpha.segment(2, n-2);
        upper.segment(1, (n-1)-1) += -d * alpha.segment(1, n-2);
        // boundary conditions
        if(m_problem->bc_type(0) == 0){ // left Dirichlet
            diagonal[0] = 1.0;  // TODO get indices from mesh
            upper[0] = 0.0;
          } else{  // Neumann
            diagonal[0] = d * 2.0 * alpha[0];  // TODO get indices from mesh
            upper[0] = -d * 2.0 * alpha[0];
          }
        if(m_problem->bc_type(1) == 0){ // right Dirichlet
            diagonal[n-1] = 1.0;  // TODO get indices from mesh
            lower[(n - 1)-1] = 0.0;
          } else {  // Neumann
            diagonal[n-1] = d * 2.0 * alpha[n-1];  // TODO get indices from mesh
            lower[(n - 1)-1] = -d * 2.0 * alpha[n-1];
          }
        for(int i=0; i<n; i++){
            trp.push_back(Trip(i,i,diagonal[i]));    
        }
        // std::cout << lower << "\n";
        for(int i=1; i<n - 1; i++){
            trp.push_back(Trip(i,i-1,lower[i-1]));    
        }
        // std::cout << upper << "\n";
        for(int i=1; i<n - 1; i++){
            trp.push_back(Trip(i,i+1,upper[i]));    
        }
        // create sparse matrix
        m_A.setFromTriplets(trp.begin(), trp.end());
    } else if (m_problem->dim() == 2){
        spdlog::info("2D");
    } else {  // 3D
        spdlog::info("3D");
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
      // define variables
      int n, n_x, n_y, n_z, n_t;
      T dx, dy, dz, dt, theta;
      theta = m_problem->theta();
      n_x = m_problem->n_x();
      n_y = m_problem->n_y();
      n_z = m_problem->n_z();
      n = m_problem->n();
      n_t = m_problem->n_t();
      dt = m_problem->dt();
      dx = m_problem->dx()[0];
      dy = m_problem->dx()[1];
      dz = m_problem->dx()[2];
      Vec alpha(n), f_n(n), f(n), x(n), y(n), z(n), ix(n_x  +1), iy(n_y + 1), 
          iz(n_z + 1), it(n_t + 1);
      // index sets
      ix.setLinSpaced(n_x + 1, 0, n_x);
      Coord coords = m_problem->coordinates();
      // define the diffusion coefficient and source term only they depend
      // on time
      for(int i=0; i<alpha.size(); i++) {
          x[i] = coords[i][0];
          y[i] = coords[i][1];
          z[i] = coords[i][2];
          alpha[i] = m_problem->alpha(coords[i], t);
          f_n[i] = m_problem->source(coords[i], t);
          f[i] = m_problem->source(coords[i], t + n_t * dt);
      }
      if (m_problem->dim() == 1){
          spdlog::info("1D");
          T d = dt/dx/dx*(1.0 - theta) / 2.0;
          m_b.segment(1, n-2) = u_n.segment(1, n-2);
          m_b.segment(1, n-2) += d * ((alpha.segment(2, n-2) + 
            alpha.segment(1, n-2)).array() * (u_n.segment(2, n-2) - 
            u_n.segment(1, n-2)).array()).matrix();
          m_b.segment(1, n-2) -= d * ((alpha.segment(1, n-2) + 
            alpha.segment(0, n-2)).array() * (u_n.segment(1, n-2) - 
            u_n.segment(0, n-2)).array()).matrix();
          m_b.segment(1, n-2) += dt * theta * f.segment(1, n-2);
          m_b.segment(1, n-2) += dt * (1.0 -theta) * f_n.segment(1, n-2);
          // Boundary conditions
          // TODO use m_problem->rhs_bc(bnd, b)
          if(m_problem->bc_type(0) == 0){ // left Dirichlet
              m_b[0] = m_problem->left(0., 0., t + n_t * dt);
            } else {  // Neumann
              // TODO
            }
          if(m_problem->bc_type(1) == 0){ // right Dirichlet
              m_b[n-1] = m_problem->right(0., 0., t + n_t * dt);
            } else {  // Neumann
              // TODO
            }
      } else if (m_problem->dim() == 2){
        spdlog::info("2D");
      } else {  // 3D
        spdlog::info("3D");
      }
  }
  /*
  * Solve the time dependent problem.
  */
  virtual void solve(){
      // Set initial condition
      Vec u_n = m_problem->u_0();
      Vec u = Vec::Zero(m_problem->n());
      int n_t = m_problem->n_t();
      T t;
      Vec t_list = m_problem->t();
      assemble_a();
      std::cout << m_A;
      Eigen::SimplicialLDLT<SpMat> solver;
      // Time loop
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
        spdlog::info("{}", u[3]);
        u_n = u;
      }
  }
};

}  // namespace fdm

}  // namespace numerical

#endif  // SPARSESOLVER_H
