
#include "math/mechanics/inertia.h"
#include <iostream>


int main()
{
    math::Cone3d cone;
    cone.pose.setIdentity();
    cone.radius = 1;
    cone.length = 1;
    double density = 1;

    math::SpatialInertia3d inertia = math::primitive_to_inertia(cone, density);
    std::cout << "Inertia:\n" << inertia.matrix() << std::endl;
}
