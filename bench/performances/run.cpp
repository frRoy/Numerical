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
	{  // vectorization
	double outl, outv;
	Eigen::VectorXd ul, uv;
	{  // scope for standard loop
	    utils::ExecutionTimer<std::chrono::nanoseconds> timer;
	    ul = bench::performances::loop();
	    outl = timer.stop() * 1.0e-6;
	}
	{  // scope for vectorized loop
		utils::ExecutionTimer<std::chrono::nanoseconds> timer;
	    uv = bench::performances::vectorized_loop();
	    outv = timer.stop() * 1.0e-6;
	}
    // comparison
    assert(utils::all_close(ul, uv));
    spdlog::info("vectorization test:");
    spdlog::info("outl {:8f} ms, outv {:8f} ms\n", outl, outv);
    }

    {  // sets vs vectors for indexing
    double outs, outv;
    Eigen::VectorXd usl, uvl;
    {  // scope for sets
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        Eigen::VectorXd  usl = bench::performances::uleft_s(400, 300);
        outs = timer.stop() * 1.0e-6;
    }
    {  // scope for vectors
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        Eigen::VectorXd  uvl = bench::performances::uleft_v(400, 300);
        outv = timer.stop() * 1.0e-6;
    }
    assert(utils::all_close(usl, uvl));
    spdlog::info("sets vs vectors indexing test:");
    spdlog::info("outs {:8f} ms, outv {:8f} ms\n", outs, outv);
    }

    {  // sets vs vectors for intersection
    double outs, outv, outvo;
    int ps, pv, pvo;
    {
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        ps = bench::performances::bottom_left_corner_s(400, 300);
        outs = timer.stop() * 1.0e-6;
    }
    {
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        pv = bench::performances::bottom_left_corner_v(400, 300);
        outv = timer.stop() * 1.0e-6;
    }
    {
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        pvo = bench::performances::bottom_left_corner_v_1(400, 300);
        outvo = timer.stop() * 1.0e-6;
    }
    assert(ps == pv);
    assert(pv == pvo);
    spdlog::info("sets vs vectors intersection test:");
    spdlog::info("outs {:8f} ms, outv {:8f} ms, outvo {:8f} ms\n", 
        outs, outv, outvo);
    }
    {  // sets vs vectors for union
    double outs, outv;
    std::vector<int> us, uv;
    {
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        us = bench::performances::union_s(400, 300);
        outs = timer.stop() * 1.0e-6;
    }
    {
        utils::ExecutionTimer<std::chrono::nanoseconds> timer;
        uv = bench::performances::union_v(400, 300);
        outv = timer.stop() * 1.0e-6;
    }
    for(int i=0;i<uv.size(); i++){
        // spdlog::info("{}", uv[i]);
        assert(us[i] == us[i]);
    }
    spdlog::info("sets vs vectors union test:");
    spdlog::info("outs {:8f} ms, outv {:8f} ms\n", outs, outv);
    }
}
