
#include "rbm/geometry/primitive_collisions.h"
#include <iostream>

int test_edge_triangle() {
    rbm::Edge3d a;
    a.vertices[0] = rbm::Vector3d(0, 0, -1);
    a.vertices[1] = rbm::Vector3d(0, 0, 1);
    rbm::Triangle3d b;
    b.vertices[0] = rbm::Vector3d(-1, -1, 0);
    b.vertices[1] = rbm::Vector3d(1, -1, 0);
    b.vertices[2] = rbm::Vector3d(0, 1, 0);

    std::optional<rbm::Vector3d> point = rbm::intersection(a, b);
    if (!point.has_value()) {
        std::cout << "No intersection" << std::endl;
    }
    else {
        std::cout << "Intersection: " << point.value().transpose() << std::endl;
    }

    return 0;
}

int test_triangle_triangle() {
    rbm::Triangle3d a;
    a.vertices[0] = rbm::Vector3d(0, 0, -1);
    a.vertices[1] = rbm::Vector3d(0, 0, 1);
    a.vertices[2] = rbm::Vector3d(2, 0, 0);

    rbm::Triangle3d b;
    b.vertices[0] = rbm::Vector3d(-1, -1, 0);
    b.vertices[1] = rbm::Vector3d(1, -1, 0);
    b.vertices[2] = rbm::Vector3d(0, 1, 0);

    std::optional<rbm::Edge3d> edge = rbm::intersection(a, b);
    if (!edge.has_value()) {
        std::cout << "No intersection" << std::endl;
    }
    else {
        std::cout << "Intersection: (" << edge.value().vertices[0].transpose() << ", " << edge.value().vertices[1].transpose() << ")" << std::endl;
    }

    return 0;
}

int main() {
    test_edge_triangle();
    test_triangle_triangle();
}
