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
#include "spdlog/spdlog.h"
#include "Parameters.hpp"
#include "Mesh.hpp"

namespace numerical {

namespace fdm {

/*
 * This class is used to solve the diffusion problem in 1D with uniform and
 * constant diffusion coefficient \f$\alpha\f$ and constant Dirichlet boundary 
 * conditions.
 */
template <typename T>
class SparseSolver {
typedef Eigen::SparseMatrix<T> SpMat;
typedef Eigen::Triplet<T> Trip;
typedef Eigen::VectorXd Vec;
typedef Eigen::ArrayXd Arr;
typedef numerical::fdm::Mesh<T> Mesh;
private:
  SpMat m_A;
  Vec m_b;
  Parameters<T> m_params;
  Vec m_x, m_t;  // spatial and temporal mesh nodes
  T m_alpha, m_theta;  // diffusion constant, theta value
  T m_d0, m_dl;  // Dirichlet boundary conditions
public:
  SparseSolver(Parameters<T> params): m_params(params) {
      // define variables
      T dx, dt;
      // define 1D mesh
      Mesh mesh = Mesh(m_params.lengths, m_params.t0, m_params.tend, 
        m_params.n, m_params.nt);
      dx = mesh.dx()[0];
      dt = mesh.dt();
      m_x = mesh.x()[0];
      m_t = mesh.t();
      // set boundary conditions here
  }
  ~SparseSolver(){
      spdlog::info("parameters destroyed");
      // delete m_A;
      // delete m_b;
  }

  /**
  * Assembles sparse coefficient matrix A.
  *
  * \f[
  *    A_{i,i-1}=-F\Theta,~A{i, i}=1+2F\Theta,~A_{i, i+1}=-F\Theta
  * \f].
  *
  */
  void assemble_a(){
    int nx = m_params.n[0];
    /*
    // The loops are vectorized for efficiency -- see bench/performances
    SpMat A(nx + 1, nx + 1);
    std::vector<Trip> trp;
    // diagonal terms
    // Eigen::VectorXd val = Eigen::VectorXd::Zero(nx + 1);  // initialized to zero
    Vec diagonal = Eigen::VectorXd::Constant(nx + 1, 1.0);
    diagonal[0] = 0.0;
    diagonal[nx] = 0.0;
    // segment(pos, n) the n coeffs in the range [pos : pos + n - 1]
    diagonal.segment(1, diagonal.size()-2) += Dl * m_a.segment(2, m_a.size()-2);
    diagonal.segment(1, diagonal.size()-2) += Dl * 2.0 * m_a.segment(1, m_a.size()-2);
    diagonal.segment(1, diagonal.size()-2) += Dl * m_a.segment(0, m_a.size()-2);
    // lower terms
    Vec lower = Eigen::VectorXd::Zero(nx);
    lower.segment(0, lower.size()-1) += -Dl * m_a.segment(1, m_a.size()-2);
    lower.segment(0, lower.size()-1) += -Dl * m_a.segment(0, m_a.size()-2);
    // upper terms
    Vec upper = Eigen::VectorXd::Zero(nx);
    upper.segment(1, upper.size()-1) += -Dl * m_a.segment(2, m_a.size()-2);
    upper.segment(1, upper.size()-1) += -Dl * m_a.segment(1, m_a.size()-2);

    // boundary conditions
    diagonal[0] = 1.0;
    upper[0] = 0.0;
    diagonal[nx] = 1.0;
    lower[nx-1] = 0.0;
    
    // std::cout << diagonal << "\n";
    for(int i=1; i<m_x.size() - 1; i++){
        trp.push_back(Trip(i,i,diagonal[i]));    
    }
    // std::cout << lower << "\n";
    for(int i=1; i<m_x.size() - 1; i++){
         trp.push_back(Trip(i,i-1,lower[i-1]));    
    }
    // std::cout << upper << "\n";
    for(int i=1; i<m_x.size() - 1; i++){
        trp.push_back(Trip(i,i+1,upper[i]));    
    }
    // create sparse matrix
    A.setFromTriplets(trp.begin(), trp.end());
    m_A = A;
    */
  }

  /**
  * Assembles RHS vector b.
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
  void assemble_b(T t){

  }

  /*
  * Solve the time dependent problem.
  */
  virtual Eigen::VectorXf solve(){
      int nx = m_params.n[0];
      // assemble();
      Vec u = Eigen::VectorXd::Zero(nx + 1);  // solution array at t[n+1]
      Vec u_n = Eigen::VectorXd::Zero(nx + 1);  // solution at t[n]
      Vec b = Eigen::VectorXd::Zero(nx + 1);  // right-hand side
      Arr arr_u = Eigen::ArrayXd::Zero(nx + 1);

      //  Set initial condition
      for(int i=0; i<u_n.size(); i++){
          u_n[i] = m_params.init(m_x[i], 0., 0.);    
      }
      // std::cout << u_n << "\n";
      // Time loop
      for(int n=0; n<m_params.nt; n++){
          b.segment(1, b.size()-2) = u_n.segment(1, u_n.size()-2);
      }
      Eigen::VectorXf solution = Eigen::VectorXf::Unit(4,1);
      return solution;
  }

};

/*
 * For a sparse matrix, return a vector of triplets, such that we can
 * reconstruct the matrix using setFromTriplet function
 * @param matrix A sparse matrix
 * @return A triplet with the row, column and value of the non-zero entries.
 * See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html for more
 * information on the triplet
 */
template <typename Derived>
std::vector<Eigen::Triplet<typename Derived::Scalar>> SparseMatrixToTriplets(
    const Derived& matrix) {
  using Scalar = typename Derived::Scalar;
  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(matrix.nonZeros());
  for (int i = 0; i < matrix.outerSize(); i++) {
    for (typename Derived::InnerIterator it(matrix, i); it; ++it) {
      triplets.push_back(
          Eigen::Triplet<Scalar>(it.row(), it.col(), it.value()));
    }
  }
  return triplets;
}

}  // namespace fdm

}  // namespace numerical

#endif  // SPARSESOLVER_H
