
#include "mbox/probability/gaussian.h"
#include <iostream>


int main() {

    mbox::Vector3d mean = mbox::Vector3d(1, 2, 3);
    mbox::Matrix3d covariance = mbox::Vector3d(1, 1, 1).asDiagonal();

    std::default_random_engine rng;
    rng.seed(0);

    mbox::GaussianDistribution<mbox::Vector3d> distribution;
    distribution.set_mean(mean);
    distribution.set_covariance(covariance);

    mbox::Vector3d x = mbox::Vector3d(3, 2, 1);
    std::cout << "p([" << x.transpose() << "]) = " << distribution.evaluate(x) << std::endl;

    for (std::size_t i = 0; i < 5; i++) {
        mbox::Vector3d yi = distribution.sample(rng);
        std::cout << "Sampled y[" << i << "] = [" << yi.transpose() << "]" << std::endl;
    }

    mbox::Distribution<mbox::Vector3d> dist_generic = distribution;
    std::cout << "Evaluating with type erasure: p([" << x.transpose() << "]) = " << dist_generic.evaluate(x) << std::endl;

    return 0;
}
