/**
 *  @file    UtilsTest.cpp
 *  @brief   Test the utility functions.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include <Eigen/Core>
#include "utils/Utils.hpp"


namespace
{

/**
* true test
*/
bool all_close_true(){
    Eigen::VectorXd u = Eigen::VectorXd::Constant(11, 1.2);
    Eigen::VectorXd d = Eigen::VectorXd::Constant(11, 1.2);
    return utils::all_close(u, d);
}

/**
* false test
*/
bool all_close_false(){
    Eigen::VectorXd u = Eigen::VectorXd::Constant(11, 1.2);
    Eigen::VectorXd d = Eigen::VectorXd::Constant(11, 1.3);
    return utils::all_close(u, d);
}

}  // namespace

TEST_CASE( "UTILS: all_close tests.", "[Utils, all_close]" )
{
    CHECK( all_close_true() );
    CHECK_FALSE( all_close_false() );
}
