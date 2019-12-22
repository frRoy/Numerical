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
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

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

Vec test_rhs(){
	Parameters* params = new Parameters();
    Problem p = Problem(params);
    Vec rhs = Vec::Constant(p.n(), 1.0);
    Vec u_n = Vec::Constant(p.n(), 1.0);
    Vec alpha = Vec::Constant(p.n(), 1.0);
    Vec f_n = Vec::Constant(p.n(), 1.0);
    Vec f = Vec::Constant(p.n(), 1.0);
    delete params;
    p.rhs_bc(rhs, u_n, alpha, f_n, f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    return rhs;
}

}  // end namespace

TEST_CASE( "Finite difference problem tests", "[fdm problem]" )
{
  CHECK( test_default() == 2 );
  CHECK(test_rhs()[2] == Approx(0.0));
}