
#include "owl/probability/gaussian.h"
#include <iostream>


int main() {

    owl::Vector3d mean = owl::Vector3d(1, 2, 3);
    owl::Matrix3d covariance = owl::Vector3d(1, 1, 1).asDiagonal();

    std::default_random_engine rng;
    rng.seed(0);

    owl::GaussianDistribution<owl::Vector3d> distribution;
    distribution.set_mean(mean);
    distribution.set_covariance(covariance);

    owl::Vector3d x = owl::Vector3d(3, 2, 1);
    std::cout << "p([" << x.transpose() << "]) = " << distribution.evaluate(x) << std::endl;

    for (std::size_t i = 0; i < 5; i++) {
        owl::Vector3d yi = distribution.sample(rng);
        std::cout << "Sampled y[" << i << "] = [" << yi.transpose() << "]" << std::endl;
    }

    owl::Distribution<owl::Vector3d> dist_generic = distribution;
    std::cout << "Evaluating with type erasure: p([" << x.transpose() << "]) = " << dist_generic.evaluate(x) << std::endl;

    return 0;
}
