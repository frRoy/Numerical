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
typedef Eigen::Triplet<double> T;

T trip = T(1,1,3.0);

bool test_solver(){
    numerical::fdm::Parameters<double>* params = 
        new numerical::fdm::Parameters<double>();
    params->lengths = {{0.0, 1.0}};
    params->n = {10};
    numerical::fdm::Problem<double>* p = 
        new numerical::fdm::Problem<double>(params);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(p);
    s.solve();
    s.assemble_a();
    delete params;
    delete p;
    return true;
}

}  // end namespace

TEST_CASE( "SparseSolver are computed", "[sparse solver]" )
{
  REQUIRE( trip.value() == Approx(3.0) );
  CHECK(test_solver());
}
