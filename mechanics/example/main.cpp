
#include "owl/mechanics/inertia.h"
#include <iostream>


int main()
{
    owl::Cone3d cone;
    cone.pose.setIdentity();
    cone.radius = 1;
    cone.length = 1;
    double density = 1;

    owl::SpatialInertia3d inertia = owl::primitive_to_inertia(cone, density);
    std::cout << "Inertia:\n" << inertia.matrix() << std::endl;
}
