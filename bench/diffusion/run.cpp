/**
 *  @file run.cpp
 *  @brief Run a diffusion problem.
 *  @author Francois Roy
 *  @date 12/01/2019
 */
#include <string>
#include <cassert>
#include <map>
#include "spdlog/spdlog.h"
#include "diffusion.hpp"

// Set global log level to debug 
// spdlog::set_level(spdlog::level::debug);  
// change log pattern
// spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");

int main(int argc, char* argv[])
{   
    // choose problem here
    std::string problem = "d";
    // define problems here
    std::string problems[4] = {"a", "b", "c", "d"};
    int len = sizeof(problems)/sizeof(problems[0]);
    std::map<const std::string, double (*)()> problem_map = {
            {"a", bench::diffusion::fdm_diffusion_a},
            {"b", bench::diffusion::fdm_diffusion_b},
            {"c", bench::diffusion::fdm_diffusion_c},
            {"d", bench::diffusion::fdm_diffusion_d} 
        };
    // run problem
    double out = 1.0;
    for(std::string p : problems){
        if(p == problem){
            try{
                spdlog::info("{}: problem {}", "Running diffusion", p);
                out = problem_map[p]();
            }
            catch (const std::invalid_argument& error){
                spdlog::error(error.what());
            }     
        }
    }
	return 0;
}