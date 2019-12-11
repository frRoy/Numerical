/**
 *  @file run.cpp
 *  @brief Run the performances test.
 *  @author Francois Roy
 *  @date 12/01/2019
 */
#include <cmath>
#include <string>
#include "spdlog/spdlog.h"
#include "Performances.hpp"
#include "utils/ExecutionTimer.hpp"
#include "utils/Utils.hpp"
#include <cassert>

int main(int argc, char* argv[])
{
	// unsigned int n = 10000;
	double outl, outv;
	Eigen::VectorXd ul, uv;
	{
	    utils::ExecutionTimer<std::chrono::nanoseconds> timer;
	    ul = bench::performances::loop();
	    outl = timer.stop() * 1.0e-6;
	}
	{
		utils::ExecutionTimer<std::chrono::nanoseconds> timer;
	    uv = bench::performances::vectorized_loop();
	    outv = timer.stop() * 1.0e-6;
	}
    assert(utils::all_close(ul, uv));
    spdlog::info("outl {:8f} ms, outv {:8f} ms", outl, outv);
}
