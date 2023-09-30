
#include "owl/geometry/primitive_collisions.h"
#include <iostream>

int test_edge_triangle() {
    owl::Edge3d a;
    a.vertices[0] = owl::Vector3d(0, 0, -1);
    a.vertices[1] = owl::Vector3d(0, 0, 1);
    owl::Triangle3d b;
    b.vertices[0] = owl::Vector3d(-1, -1, 0);
    b.vertices[1] = owl::Vector3d(1, -1, 0);
    b.vertices[2] = owl::Vector3d(0, 1, 0);

    std::optional<owl::Vector3d> point = owl::intersection(a, b);
    if (!point.has_value()) {
        std::cout << "No intersection" << std::endl;
    }
    else {
        std::cout << "Intersection: " << point.value().transpose() << std::endl;
    }

    return 0;
}

int test_triangle_triangle() {
    owl::Triangle3d a;
    a.vertices[0] = owl::Vector3d(0, 0, -1);
    a.vertices[1] = owl::Vector3d(0, 0, 1);
    a.vertices[2] = owl::Vector3d(2, 0, 0);

    owl::Triangle3d b;
    b.vertices[0] = owl::Vector3d(-1, -1, 0);
    b.vertices[1] = owl::Vector3d(1, -1, 0);
    b.vertices[2] = owl::Vector3d(0, 1, 0);

    std::optional<owl::Edge3d> edge = owl::intersection(a, b);
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
