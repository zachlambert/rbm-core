#pragma once

#include "owl/transform/euler.h"
#include "parrot/json.h"
#include "owl/types/serialize.h"


template <typename Scalar, int Dim>
struct parrot::json_functions<owl::EulerRotation<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::EulerRotation<Scalar, Dim>& rotation) {
        return node.read(rotation.coords());
    }
    static void to_json(const owl::EulerRotation<Scalar, Dim>& rotation, parrot::JsonNode node) {
        node = rotation.coords();
    }
};

template <typename Scalar, int Dim>
struct parrot::json_functions<owl::EulerTransform<Scalar, Dim>> {
    static bool from_json(parrot::JsonConstNode node, owl::EulerTransform<Scalar, Dim>& transform) {
        if (!node.is_object()) return false;
        if (!node["rotation"].read(transform.rotation())) {
            return false;
        }
        if (!node["translation"].read(transform.translation())) {
            return false;
        }
        return true;
    }
    static void to_json(const owl::EulerTransform<Scalar, Dim>& transform, parrot::JsonNode node) {
        node.set_object();
        node["rotation"] = transform.rotation();
        node["translation"] = transform.translation();
    }
};
