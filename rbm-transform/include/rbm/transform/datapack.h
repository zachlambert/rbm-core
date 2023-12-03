#pragma once

#include <datapack/json.h>
#include <rbm/types/datapack.h>
#include "rbm/transform/euler.h"


template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::EulerRotation<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::EulerRotation<Scalar, Dim>& rotation) {
        return node.read(rotation.coords());
    }
    static void to_json(const rbm::EulerRotation<Scalar, Dim>& rotation, datapack::JsonNode node) {
        node = rotation.coords();
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<rbm::EulerTransform<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, rbm::EulerTransform<Scalar, Dim>& transform) {
        if (!node.is_object()) return false;
        if (!node["rotation"].read(transform.rotation())) {
            return false;
        }
        if (!node["translation"].read(transform.translation())) {
            return false;
        }
        return true;
    }
    static void to_json(const rbm::EulerTransform<Scalar, Dim>& transform, datapack::JsonNode node) {
        node.set_object();
        node["rotation"] = transform.rotation();
        node["translation"] = transform.translation();
    }
};
