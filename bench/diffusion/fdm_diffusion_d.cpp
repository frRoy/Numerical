/**
 *  @file fdm_diffusion_d.cpp
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

double FDMDiffusionD::left(double y, double z, double t){
    std::vector<double> lx = FDMDiffusionD::m_params->lengths[0];
    double length = lx[1] - lx[0];
    return 5.0*t*length*length/4.0;
}
double FDMDiffusionD::right(double y, double z, double t){
    std::vector<double> lx = FDMDiffusionD::m_params->lengths[0];
    double length = lx[1] - lx[0];
    return 5.0*t*length*length/4.0;
}
double FDMDiffusionD::source(const Eigen::Matrix<double, 3, 1>& x, double t){
    std::vector<double> lx = FDMDiffusionD::m_params->lengths[0];
    double alpha = FDMDiffusionD::alpha(x, t);
    double length = lx[1] - lx[0];
    return 5*(x[0]-length/2.0)*(x[0]-length/2.0) - 10.0*alpha*t;
}
double FDMDiffusionD::reference(const Eigen::Matrix<double, 3, 1>& x, 
        double t){
    std::vector<double> lx = FDMDiffusionD::m_params->lengths[0];
    double length = lx[1] - lx[0];
    // spdlog::info("length: {}, x: {}, t: {}", length, x[0], t);
    return 5.0 * t * (x[0] - length/2.0)* (x[0]-length/2.0);
}

double fdm_diffusion_d(){
    Params* p = new Params();
    p->lengths = {{0.0, 1.0}};
    p->n = {10};
    p->nt = {5};
    p->tend = 0.01;
    p->alpha = 0.01;
    FDMDiffusionD *prob = new FDMDiffusionD(p);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(prob);
    double l2norm = s.solve();
    delete p;
    delete prob;
	return l2norm;
}

}  // namespace diffusion

}  // namespace bench
