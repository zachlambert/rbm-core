#pragma once

#include "mathbox/transform/euler.h"
#include "datapack/json.h"
#include "mathbox/types/datapack.h"


template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::EulerRotation<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::EulerRotation<Scalar, Dim>& rotation) {
        return node.read(rotation.coords());
    }
    static void to_json(const mbox::EulerRotation<Scalar, Dim>& rotation, datapack::JsonNode node) {
        node = rotation.coords();
    }
};

template <typename Scalar, int Dim>
struct datapack::json_functions<mbox::EulerTransform<Scalar, Dim>> {
    static bool from_json(datapack::JsonConstNode node, mbox::EulerTransform<Scalar, Dim>& transform) {
        if (!node.is_object()) return false;
        if (!node["rotation"].read(transform.rotation())) {
            return false;
        }
        if (!node["translation"].read(transform.translation())) {
            return false;
        }
        return true;
    }
    static void to_json(const mbox::EulerTransform<Scalar, Dim>& transform, datapack::JsonNode node) {
        node.set_object();
        node["rotation"] = transform.rotation();
        node["translation"] = transform.translation();
    }
};
