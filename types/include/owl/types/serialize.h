#pragma once

#include "owl/types/matrix.h"
#include "owl/types/color.h"
#include "owl/types/interval.h"
#include "parrot/json.h"


template <typename Scalar, int Rows, int Cols>
struct parrot::json_functions<owl::Matrix<Scalar, Rows, Cols>> {
    static bool from_json(parrot::JsonConstNode node, owl::Matrix<Scalar, Rows, Cols>& matrix) {
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
    static void to_json(const owl::Matrix<Scalar, Rows, Cols>& matrix, parrot::JsonNode node) {
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
struct parrot::json_functions<owl::ColorRGB<Scalar>> {
    static bool from_json(parrot::JsonConstNode node, owl::ColorRGB<Scalar>& color) {
        return node.read(color.data);
    }
    static void to_json(const owl::ColorRGB<Scalar>& color, parrot::JsonNode node) {
        node = color.data;
    }
};

template <typename Scalar>
struct parrot::json_functions<owl::Interval<Scalar>> {
    static bool from_json(parrot::JsonConstNode node, owl::Interval<Scalar>& interval) {
        if (!node.is_array()) return false;
        if (node.size() != 2) return false;
        node[0].read(interval.lower);
        node[1].read(interval.upper);
        return true;
    }
    static void to_json(const owl::Interval<Scalar>& interval, parrot::JsonNode node) {
        node.set_array();
        node[0] = interval.lower;
        node[1] = interval.upper;
    }
};
