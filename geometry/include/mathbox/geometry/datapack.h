#pragma once

#include "datapack/json.h"
#include "mbox/types/datapack.h"
#include "mbox/transform/datapack.h"

#include "mbox/geometry/vertex.h"
#include "mbox/geometry/mesh.h"
#include "mbox/geometry/primitive.h"
#include "mbox/transform/euler.h"


template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::Box<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::Box<Scalar, Dim>& box) {
        bool valid = true;
        mbox::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        box.pose = pose.toTransform();
        valid &= node["size"].read(box.size);
        return valid;
    }
    static void to_json(const mbox::Box<Scalar, Dim>& box, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = mbox::EulerTransform<Scalar, Dim>(box.pose);
        node["size"] = box.size;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::Sphere<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::Sphere<Scalar, Dim>& sphere) {
        bool valid = true;
        valid &= node["position"].read(sphere.position);
        valid &= node["radius"].read(sphere.radius);
        return valid;
    }
    static void to_json(const mbox::Sphere<Scalar, Dim>& sphere, datapack::JsonNode node) {
        node.set_object();
        node["position"] = sphere.position;
        node["radius"] = sphere.radius;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::Cylinder<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::Cylinder<Scalar, Dim>& cylinder) {
        bool valid = true;
        mbox::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cylinder.pose = pose.toTransform();
        valid &= node["radius"].read(cylinder.radius);
        valid &= node["length"].read(cylinder.length);
        return valid;
    }
    static void to_json(const mbox::Cylinder<Scalar, Dim>& cylinder, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = mbox::EulerTransform<Scalar, Dim>(cylinder.pose);
        node["radius"] = cylinder.radius;
        node["length"] = cylinder.length;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::Cone<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::Cone<Scalar, Dim>& cone) {
        bool valid = true;
        mbox::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cone.pose = pose.toTransform();
        valid &= node["radius"].read(cone.radius);
        valid &= node["length"].read(cone.length);
        return valid;
    }
    static void to_json(const mbox::Cone<Scalar, Dim>& cone, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = mbox::EulerTransform<Scalar, Dim>(cone.pose);
        node["radius"] = cone.radius;
        node["length"] = cone.length;
    }
};
