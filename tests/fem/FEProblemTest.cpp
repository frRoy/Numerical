/**
 *  @file    FEProblemTest.cpp
 *  @brief   Test the FEMProblem class.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fem/FEProblem.hpp"

namespace
{
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

}  // namespace

TEST_CASE( "FEM: Problem tests", "[FEProblem]" )
{
    CHECK(1.0 == Approx(1.0));
}
