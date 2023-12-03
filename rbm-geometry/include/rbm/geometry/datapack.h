#pragma once

#include "datapack/json.h"
#include "rbm/types/datapack.h"
#include "rbm/transform/datapack.h"

#include "rbm/geometry/vertex.h"
#include "rbm/geometry/mesh.h"
#include "rbm/geometry/primitive.h"
#include "rbm/transform/euler.h"


template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::Box<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Box<Scalar, Dim>& box) {
        bool valid = true;
        rbm::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        box.pose = pose.toTransform();
        valid &= node["size"].read(box.size);
        return valid;
    }
    static void to_json(const rbm::Box<Scalar, Dim>& box, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = rbm::EulerTransform<Scalar, Dim>(box.pose);
        node["size"] = box.size;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::Sphere<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Sphere<Scalar, Dim>& sphere) {
        bool valid = true;
        valid &= node["position"].read(sphere.position);
        valid &= node["radius"].read(sphere.radius);
        return valid;
    }
    static void to_json(const rbm::Sphere<Scalar, Dim>& sphere, datapack::JsonNode node) {
        node.set_object();
        node["position"] = sphere.position;
        node["radius"] = sphere.radius;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::Cylinder<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Cylinder<Scalar, Dim>& cylinder) {
        bool valid = true;
        rbm::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cylinder.pose = pose.toTransform();
        valid &= node["radius"].read(cylinder.radius);
        valid &= node["length"].read(cylinder.length);
        return valid;
    }
    static void to_json(const rbm::Cylinder<Scalar, Dim>& cylinder, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = rbm::EulerTransform<Scalar, Dim>(cylinder.pose);
        node["radius"] = cylinder.radius;
        node["length"] = cylinder.length;
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::Cone<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Cone<Scalar, Dim>& cone) {
        bool valid = true;
        rbm::EulerTransform<Scalar, Dim> pose;
        valid &= node["pose"].read(pose);
        cone.pose = pose.toTransform();
        valid &= node["radius"].read(cone.radius);
        valid &= node["length"].read(cone.length);
        return valid;
    }
    static void to_json(const rbm::Cone<Scalar, Dim>& cone, datapack::JsonNode node) {
        node.set_object();
        node["pose"] = rbm::EulerTransform<Scalar, Dim>(cone.pose);
        node["radius"] = cone.radius;
        node["length"] = cone.length;
    }
};
