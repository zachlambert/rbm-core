#pragma once

#include "core/std/json.h"
#include "math/types/json.h"
#include "math/transform/json.h"

#include "math/geometry/vertex.h"
#include "math/geometry/mesh.h"
#include "math/geometry/primitive.h"
#include "math/transform/euler.h"


template <typename Scalar, int Dim>
struct core::json_functions<math::Box<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::Box<Scalar, Dim>& box) {
        bool valid = true;
        math::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        box.pose = pose.toTransform();
        valid &= node["size"].read(box.size);
        return valid;
    }
    static void to_json(const math::Box<Scalar, Dim>& box, core::JsonNode node) {
        node.set_object();
        node["pose"] = math::EulerTransform<Scalar, Dim>(box.pose);
        node["size"] = box.size;
    }
};

template <typename Scalar, int Dim>
struct core::json_functions<math::Sphere<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::Sphere<Scalar, Dim>& sphere) {
        bool valid = true;
        valid &= node["position"].read(sphere.position);
        valid &= node["radius"].read(sphere.radius);
        return valid;
    }
    static void to_json(const math::Sphere<Scalar, Dim>& sphere, core::JsonNode node) {
        node.set_object();
        node["position"] = sphere.position;
        node["radius"] = sphere.radius;
    }
};

template <typename Scalar, int Dim>
struct core::json_functions<math::Cylinder<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::Cylinder<Scalar, Dim>& cylinder) {
        bool valid = true;
        math::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cylinder.pose = pose.toTransform();
        valid &= node["radius"].read(cylinder.radius);
        valid &= node["length"].read(cylinder.length);
        return valid;
    }
    static void to_json(const math::Cylinder<Scalar, Dim>& cylinder, core::JsonNode node) {
        node.set_object();
        node["pose"] = math::EulerTransform<Scalar, Dim>(cylinder.pose);
        node["radius"] = cylinder.radius;
        node["length"] = cylinder.length;
    }
};

template <typename Scalar, int Dim>
struct core::json_functions<math::Cone<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::Cone<Scalar, Dim>& cone) {
        bool valid = true;
        math::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cone.pose = pose.toTransform();
        valid &= node["radius"].read(cone.radius);
        valid &= node["length"].read(cone.length);
        return valid;
    }
    static void to_json(const math::Cone<Scalar, Dim>& cone, core::JsonNode node) {
        node.set_object();
        node["pose"] = math::EulerTransform<Scalar, Dim>(cone.pose);
        node["radius"] = cone.radius;
        node["length"] = cone.length;
    }
};
