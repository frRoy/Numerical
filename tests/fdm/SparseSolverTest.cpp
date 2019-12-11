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

}  // end namespace

TEST_CASE( "SparseSolver are computed", "[sparse solver]" )
{
  REQUIRE( trip.value() == Approx(3.0) );
}
