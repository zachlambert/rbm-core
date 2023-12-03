
#include "rbm/probability/gaussian.h"
#include <iostream>


int main() {

    rbm::Vector3d mean = rbm::Vector3d(1, 2, 3);
    rbm::Matrix3d covariance = rbm::Vector3d(1, 1, 1).asDiagonal();

    std::default_random_engine rng;
    rng.seed(0);

    rbm::GaussianDistribution<rbm::Vector3d> distribution;
    distribution.set_mean(mean);
    distribution.set_covariance(covariance);

    rbm::Vector3d x = rbm::Vector3d(3, 2, 1);
    std::cout << "p([" << x.transpose() << "]) = " << distribution.evaluate(x) << std::endl;

    for (std::size_t i = 0; i < 5; i++) {
        rbm::Vector3d yi = distribution.sample(rng);
        std::cout << "Sampled y[" << i << "] = [" << yi.transpose() << "]" << std::endl;
    }

    rbm::Distribution<rbm::Vector3d> dist_generic = distribution;
    std::cout << "Evaluating with type erasure: p([" << x.transpose() << "]) = " << dist_generic.evaluate(x) << std::endl;

    return 0;
}
