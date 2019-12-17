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
#include "Problem.hpp"

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
  SpMat* m_A;
  Vec* m_b;
  Vec* m_u;
  Problem<T>* m_problem;
public:
  SparseSolver(Problem<T>* problem)
    : m_problem(problem)
    {
        int n = m_problem->n();
        m_A = new SpMat(n, n);
        m_b = new Vec(n);
        m_u = new Vec(n);
  }
  ~SparseSolver(){
      delete m_A;
      delete m_b;
      delete m_u;
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
    // define variables
    int n, n_x, n_y, n_z, n_t;
    T dx, dy, dz, dt, t=0.0;
    n_x = m_problem->n_x();
    n_y = m_problem->n_y();
    n_z = m_problem->n_z();
    n = m_problem->n();
    n_t = m_problem->nt();
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
        T d = dt/dx/dx;
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
        // spdlog::info("before: {}", diagonal[0]);
        m_problem->diagonal_bc(0, diagonal, t); // left
        // spdlog::info("after: {}", diagonal[0]);
        m_problem->diagonal_bc(1, diagonal, t); // right
        upper[0] = 0.0;
        lower[(n - 1)-1] = 0.0;
        for(int i=1; i<n - 1; i++){
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
        m_A->setFromTriplets(trp.begin(), trp.end());
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
  void assemble_b(int n){

  }
  /*
  * Solve the time dependent problem.
  */
  virtual Eigen::VectorXf solve(){
      //  Set initial condition
      //for(int i=0; i<u_n.size(); i++){
          // u_n[i] = m_params.init(m_x[i], 0., 0.);    
      //}
      // std::cout << u_n << "\n";

      // spdlog::info("{}", m_problem->left(1.0, 1.0, 1.0));
      // Time loop

      Eigen::VectorXf solution = Eigen::VectorXf::Unit(4,1);
      return solution;
  }
};

}  // namespace fdm

}  // namespace numerical

#endif  // SPARSESOLVER_H
