#pragma once

#include "parrot/json.h"
#include "owl/types/serialize.h"
#include "owl/transform/serialize.h"

#include "owl/geometry/vertex.h"
#include "owl/geometry/mesh.h"
#include "owl/geometry/primitive.h"
#include "owl/transform/euler.h"


template <typename Scalar, int Dim>
struct parrot::json_functions<owl::Box<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::Box<Scalar, Dim>& box) {
        bool valid = true;
        owl::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        box.pose = pose.toTransform();
        valid &= node["size"].read(box.size);
        return valid;
    }
    static void to_json(const owl::Box<Scalar, Dim>& box, parrot::JsonNode node) {
        node.set_object();
        node["pose"] = owl::EulerTransform<Scalar, Dim>(box.pose);
        node["size"] = box.size;
    }
};

template <typename Scalar, int Dim>
struct parrot::json_functions<owl::Sphere<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::Sphere<Scalar, Dim>& sphere) {
        bool valid = true;
        valid &= node["position"].read(sphere.position);
        valid &= node["radius"].read(sphere.radius);
        return valid;
    }
    static void to_json(const owl::Sphere<Scalar, Dim>& sphere, parrot::JsonNode node) {
        node.set_object();
        node["position"] = sphere.position;
        node["radius"] = sphere.radius;
    }
};

template <typename Scalar, int Dim>
struct parrot::json_functions<owl::Cylinder<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::Cylinder<Scalar, Dim>& cylinder) {
        bool valid = true;
        owl::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cylinder.pose = pose.toTransform();
        valid &= node["radius"].read(cylinder.radius);
        valid &= node["length"].read(cylinder.length);
        return valid;
    }
    static void to_json(const owl::Cylinder<Scalar, Dim>& cylinder, parrot::JsonNode node) {
        node.set_object();
        node["pose"] = owl::EulerTransform<Scalar, Dim>(cylinder.pose);
        node["radius"] = cylinder.radius;
        node["length"] = cylinder.length;
    }
};

template <typename Scalar, int Dim>
struct parrot::json_functions<owl::Cone<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::Cone<Scalar, Dim>& cone) {
        bool valid = true;
        owl::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cone.pose = pose.toTransform();
        valid &= node["radius"].read(cone.radius);
        valid &= node["length"].read(cone.length);
        return valid;
    }
    static void to_json(const owl::Cone<Scalar, Dim>& cone, parrot::JsonNode node) {
        node.set_object();
        node["pose"] = owl::EulerTransform<Scalar, Dim>(cone.pose);
        node["radius"] = cone.radius;
        node["length"] = cone.length;
    }
};
