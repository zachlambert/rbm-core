#pragma once

#include <datapack/json.h>
#include "rbm/types/matrix.h"
#include "rbm/types/color.h"
#include "rbm/types/interval.h"


template <typename Scalar, int Rows, int Cols>
struct datapack::json_functions<rbm::Matrix<Scalar, Rows, Cols>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Matrix<Scalar, Rows, Cols>& matrix) {
        if (!node.is_array()) return false;
        for (std::size_t i = 0; i < matrix.rows(); i++) {
            if (matrix.cols() == 1) {
                if (!node[i].read(matrix(i, 0))) return false;
            }
            else {
                auto row = node[i];
                if (!row.is_array()) return false;
                for (std::size_t j = 0; j < matrix.cols(); j++) {
                    if (!row[j].read(matrix(i, j))) return false;
                }
            }
        }
        return true;
    }
    static void to_json(const rbm::Matrix<Scalar, Rows, Cols>& matrix, datapack::JsonNode node) {
        node.set_array();
        for (std::size_t i = 0; i < matrix.rows(); i++) {
            if (matrix.cols() == 1) {
                node[i] = matrix(i, 0);
            }
            else {
                auto row = node[i];
                row.set_array();
                for (std::size_t j = 0; j < matrix.cols(); j++) {
                    row[j] = matrix(i, j);
                }
            }
        }
    }
};

template <typename Scalar>
struct datapack::json_functions<rbm::ColorRGB<Scalar>> {
    static bool from_json(datapack::JsonConstNode node, rbm::ColorRGB<Scalar>& color) {
        return node.read(color.data);
    }
    static void to_json(const rbm::ColorRGB<Scalar>& color, datapack::JsonNode node) {
        node = color.data;
    }
};

template <typename Scalar>
struct datapack::json_functions<rbm::Interval<Scalar>> {
    static bool from_json(datapack::JsonConstNode node, rbm::Interval<Scalar>& interval) {
        if (!node.is_array()) return false;
        if (node.size() != 2) return false;
        node[0].read(interval.lower);
        node[1].read(interval.upper);
        return true;
    }
    static void to_json(const rbm::Interval<Scalar>& interval, datapack::JsonNode node) {
        node.set_array();
        node[0] = interval.lower;
        node[1] = interval.upper;
    }
};
