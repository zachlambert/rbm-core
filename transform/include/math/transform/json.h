#pragma once

#include "math/transform/euler.h"
#include "core/serialize/json.h"
#include "math/types/json.h"


template <typename Scalar, int Dim>
struct core::json_functions<math::EulerRotation<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::EulerRotation<Scalar, Dim>& rotation) {
        return node.read(rotation.coords());
    }
    static void to_json(const math::EulerRotation<Scalar, Dim>& rotation, core::JsonNode node) {
        node = rotation.coords();
    }
};

template <typename Scalar, int Dim>
struct core::json_functions<math::EulerTransform<Scalar, Dim>> {
    static bool from_json(core::JsonConstNode node, math::EulerTransform<Scalar, Dim>& transform) {
        if (!node.is_object()) return false;
        if (!node["rotation"].read(transform.rotation())) {
            return false;
        }
        if (!node["translation"].read(transform.translation())) {
            return false;
        }
        return true;
    }
    static void to_json(const math::EulerTransform<Scalar, Dim>& transform, core::JsonNode node) {
        node.set_object();
        node["rotation"] = transform.rotation();
        node["translation"] = transform.translation();
    }
};
