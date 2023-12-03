
#include "mbox/geometry/primitive_collisions.h"
#include <iostream>

int test_edge_triangle() {
    mbox::Edge3d a;
    a.vertices[0] = mbox::Vector3d(0, 0, -1);
    a.vertices[1] = mbox::Vector3d(0, 0, 1);
    mbox::Triangle3d b;
    b.vertices[0] = mbox::Vector3d(-1, -1, 0);
    b.vertices[1] = mbox::Vector3d(1, -1, 0);
    b.vertices[2] = mbox::Vector3d(0, 1, 0);

    std::optional<mbox::Vector3d> point = mbox::intersection(a, b);
    if (!point.has_value()) {
        std::cout << "No intersection" << std::endl;
    }
    else {
        std::cout << "Intersection: " << point.value().transpose() << std::endl;
    }

    return 0;
}

int test_triangle_triangle() {
    mbox::Triangle3d a;
    a.vertices[0] = mbox::Vector3d(0, 0, -1);
    a.vertices[1] = mbox::Vector3d(0, 0, 1);
    a.vertices[2] = mbox::Vector3d(2, 0, 0);

    mbox::Triangle3d b;
    b.vertices[0] = mbox::Vector3d(-1, -1, 0);
    b.vertices[1] = mbox::Vector3d(1, -1, 0);
    b.vertices[2] = mbox::Vector3d(0, 1, 0);

    std::optional<mbox::Edge3d> edge = mbox::intersection(a, b);
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
