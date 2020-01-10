/**
 *  @file    ExecutionTimerTest.cpp
 *  @brief   Test the @c ExecutionTimer class.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include <Eigen/Core>
#include "utils/ExecutionTimer.hpp"


namespace
{

/**
* A loop that takes some time to execute.
* @param n The number of iterations.
*/
void loop(unsigned int n) {
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < u.size(); i++) {
        u[i] = 1.0;       
    }
}

/**
* Dummy test
*/
bool exec_time_ns(){ 
    utils::ExecutionTimer<std::chrono::nanoseconds> timer;
    unsigned int n = 100;
    loop(n);
    timer.stop();
    return true;
}

/**
* Dummy test
*/
bool exec_time_ms(){ 
    utils::ExecutionTimer<std::chrono::milliseconds> timer;
    unsigned int n = 100;
    loop(n);
    timer.stop();
    return true;
}

}  // namespace

TEST_CASE( "EXECUTION_TIMER: Execution time tests.", "[ExecutionTimer]" )
{
  REQUIRE( exec_time_ms() == true );
  REQUIRE( exec_time_ns() == true );
}
    