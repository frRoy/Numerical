/**
 *  @file    SparseSolverTest.cpp
 *  @brief   Test theSparseSolve class @see numerical::fdm::SparseSolver.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/SparseSolver.hpp"


namespace
{

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef Eigen::Triplet<double> T;

/**
* Run the solver for a 1D problem.
* 
* @return True.
*/
bool test_solver(){
    numerical::fdm::Parameters<double>* params = 
        new numerical::fdm::Parameters<double>();
    params->lengths = {{0.0, 1.0}};
    params->n = {10};
    numerical::fdm::FDProblem<double>* p = 
        new numerical::fdm::FDProblem<double>(params);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(p);
    s.solve();
    delete params;
    delete p;
    return true;
}

}  // end namespace

TEST_CASE( "SPARSESOLVER: 1D tests.", "[SparseSolver]" )
{
  CHECK(test_solver());
}
