/**
 *  @file run.cpp
 *  @brief Run a diffusion problem.
 *  @author Francois Roy
 *  @date 12/01/2019
 */
#include <string>
#include "spdlog/spdlog.h"
#include "diffusion.hpp"


int main(int argc, char* argv[])
{
    std::string problem = "b";
    double out = 1.0;
    try{
    	if(problem == "a"){
    		spdlog::info("{}: problem {}", "Running diffusion", problem);
            // TODO return the L2-norm and assert smaller than 1e-6
    	    out = bench::diffusion::fdm_diffusion_a();
        }
        if(problem == "b"){
            // TODO return the L2-norm and assert smaller than 1e-6
            spdlog::info("{}: problem {}", "Running diffusion", problem);
            out = bench::diffusion::fdm_diffusion_b();
        }
        spdlog::info("The problem {}, computed with a L2-norm error of {}",
            problem, out);
    }
    catch (const std::invalid_argument& error){
        spdlog::error(error.what());
    }
    // auto out = numerical::fdm::init_value<double> (0.2, 1.0, 0.0, 1.0);
    // implicit instantiation of the function template
    // spdlog::info("initial value at x= 0.2: {:03.2f}", 
    //               numerical::fdm::init_value(0.2, 1.0, 0.0, 1.0));
    //std::string out1 = numerical::fdm::solver(
    // spdlog::info("{}", out1);
	// spdlog::error("Some error message with arg: {}", 1);
    // spdlog::warn("Easy padding in numbers like {:08d}", 12);
    // spdlog::critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
    // spdlog::info("Support for floats {:03.2f}", 1.23456);
    // spdlog::info("Positional args are {1} {0}..", "too", "supported");
    // spdlog::set_level(spdlog::level::debug); // Set global log level to debug
    // spdlog::debug("This message should be displayed..");    
    // change log pattern
    // spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
    // spdlog::debug("This message should be displayed..");
	return 0;

}