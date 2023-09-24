
#include "math/probability/gaussian.h"
#include <iostream>


int main() {

    math::Vector3d mean = math::Vector3d(1, 2, 3);
    math::Matrix3d covariance = math::Vector3d(1, 1, 1).asDiagonal();

    std::default_random_engine rng;
    rng.seed(0);

    math::GaussianDistribution<math::Vector3d> distribution;
    distribution.set_mean(mean);
    distribution.set_covariance(covariance);

    math::Vector3d x = math::Vector3d(3, 2, 1);
    std::cout << "p([" << x.transpose() << "]) = " << distribution.evaluate(x) << std::endl;

    for (std::size_t i = 0; i < 5; i++) {
        math::Vector3d yi = distribution.sample(rng);
        std::cout << "Sampled y[" << i << "] = [" << yi.transpose() << "]" << std::endl;
    }

    math::Distribution<math::Vector3d> dist_generic = distribution;
    std::cout << "Evaluating with type erasure: p([" << x.transpose() << "]) = " << dist_generic.evaluate(x) << std::endl;

    return 0;
}
