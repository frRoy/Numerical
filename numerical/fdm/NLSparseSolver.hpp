/**
 *  @file    NLSparseSolver.hpp
 *  @brief   Nonlinear sparse solver.
 *  @author  Francois Roy
 *  @date    01/03/2020
 */
#ifndef NL_SPARSESOLVER_H
#define NL_SPARSESOLVER_H

#include "SparseSolver.hpp"

namespace numerical {

namespace fdm {

/**
 * Nonlinear sparse solver.
 */
template <typename T>
class NLSparseSolver : public numerical::fdm::SparseSolver<T> {
typedef Eigen::SparseMatrix<T> SpMat;
typedef Eigen::Triplet<T> Trip;
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
typedef std::vector<Eigen::Matrix<T, 3, 1>> Coord;
typedef numerical::fdm::SparseSolver<T> Base;
private:
public:
  NLSparseSolver(Problem<T>* problem)
    : numerical::fdm::SparseSolver<T>(problem)
    {
    }
  //void assemble_a(){
  //}
  void assemble_b(T t, Vec& u_n){
      const Coord& coords = this->m_problem->coordinates();
      // TODO Calculate defivative of f and alpha with respect to u
      for(int i=0; i<this->m_n; i++) {
          Base::m_alpha[i] = this->m_problem->alpha(coords[i], t, u_n[i]);
          this->m_f_n[i] = this->m_problem->source(coords[i], t, u_n[i]);
          this->m_f[i] = this->m_problem->source(coords[i], t + this->m_dt, 
            u_n[i]);
      }
  }
  T solve(){
    // TODO add nonlinear loop inside time loop
    spdlog::info("NL solve");
    this->assemble_a();
    assemble_b(0.0, this->m_u);
    std::cout << Base::m_A << std::endl;
    return 999.0;
  }
};

}  // namespace fdm

}  // namespace numerical

#endif  // NL_SPARSESOLVER_H