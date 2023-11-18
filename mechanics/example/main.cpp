
#include "mathbox/mechanics/inertia.h"
#include <iostream>


int main()
{
    mbox::Cone3d cone;
    cone.pose.setIdentity();
    cone.radius = 1;
    cone.length = 1;
    double density = 1;

    mbox::SpatialInertia3d inertia = mbox::primitive_to_inertia(cone, density);
    std::cout << "Inertia:\n" << inertia.matrix() << std::endl;
}
