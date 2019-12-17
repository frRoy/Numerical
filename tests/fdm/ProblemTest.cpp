/**
 *  @file    ProblemTest.cpp
 *  @brief   Test the Problem class @see numerical::fdm::Problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/Problem.hpp"


namespace
{

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
typedef numerical::fdm::Problem<double> Problem;
typedef numerical::fdm::Parameters<double> Parameters;

template<typename T>
T sum(T a, T b, T c, T d){
    return a + b + c +d;
}

/**
*/
double test_default() {
	Parameters* params = new Parameters();
    Problem p = Problem(params);
    delete params;
    return p.dim();
}

}  // end namespace

TEST_CASE( "Finite difference problem tests", "[fdm problem]" )
{
  REQUIRE( test_default() == 2 );
}