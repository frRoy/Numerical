/**
 *  @file fdm_diffusion_a.cpp
 *  @brief Run a 1D diffusion problem with Dirichlet boundary conditions.
 *  @author Francois Roy
 *  @date 12/01/2019
 */
#include "spdlog/spdlog.h"
#include "numerical/fdm/Parameters.hpp"
#include "diffusion.hpp"

namespace bench {
namespace diffusion {
typedef numerical::fdm::Parameters<double> Params;
typedef numerical::fdm::Problem<double> Problem;

// double FDMDiffusionA::alpha(double x, double y, double z, double t){
//     return m_params.alpha;
// }


bool fdm_diffusion_a(){
    Params* p = new Params();
    spdlog::info("{}", p->alpha);
    auto out = reference_diffusion_a<double>(0.1, 0.01, 1.0);
    spdlog::info("{}", out);
    FDMDiffusionA *prob = new FDMDiffusionA(p);
    spdlog::info("{}\n", prob->dim());
    double alpha = prob->alpha(0.0, 0.0, 0.0, 0.0);
    spdlog::info("{}\n", alpha);
    delete p;
    delete prob;
	return true;
}

}  // namespace diffusion

}  // namespace bench
