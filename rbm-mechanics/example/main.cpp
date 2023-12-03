
#include "rbm/mechanics/inertia.h"
#include <iostream>


int main()
{
    rbm::Cone3d cone;
    cone.pose.setIdentity();
    cone.radius = 1;
    cone.length = 1;
    double density = 1;

    rbm::SpatialInertia3d inertia = rbm::primitive_to_inertia(cone, density);
    std::cout << "Inertia:\n" << inertia.matrix() << std::endl;
}
