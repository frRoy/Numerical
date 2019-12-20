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

double FDMDiffusionA::left(double y, double z, double t){
     return 0.0;
}
double FDMDiffusionA::right(double y, double z, double t){
     return 0.0;
}
double FDMDiffusionA::source(const Eigen::Matrix<double, 3, 1>& x, double t){
    std::vector<double> lx = FDMDiffusionA::m_params->lengths[0];
    double alpha = FDMDiffusionA::alpha(x, t);
    double length = lx[1] - lx[0];
    return 5*x[0]*(length - x[0]) + 10.0*alpha*t;
}
double FDMDiffusionA::reference(const Eigen::Matrix<double, 3, 1>& x, 
        double t){
    std::vector<double> lx = FDMDiffusionA::m_params->lengths[0];
    double length = lx[1] - lx[0];
    // spdlog::info("length: {}, x: {}, t: {}", length, x[0], t);
    return 5.0 * t * x[0] * (length - x[0]);
}

bool fdm_diffusion_a(){
    Params* p = new Params();
    p->lengths = {{0.0, 1.0}};
    p->n = {10};
    p->nt = {5};
    p->tend = 0.01;
    p->alpha = 0.01;
    FDMDiffusionA *prob = new FDMDiffusionA(p);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(prob);
    s.solve();
    delete p;
    delete prob;
	return true;
}

}  // namespace diffusion

}  // namespace bench
