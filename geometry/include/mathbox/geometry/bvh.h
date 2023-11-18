#pragma once

#include "owl/geometry/bounding.h"
#include "owl/geometry/primitive.h"
#include "owl/geometry/collision.h"
#include "owl/geometry/primitive_collisions.h"
#include <algorithm>
#include <iostream>
#include <stack>


namespace owl {

// Bounded volume hierarchy
template <typename Primitive>
requires is_primitive<Primitive>
class Bvh {
    using Scalar = typename Primitive::Scalar;
    static constexpr int Dim = Primitive::Dim;

    struct Node {
        BoundingBox<Scalar, Dim> bounding_box;
        int partition_dim;
        bool leaf;
        std::size_t below_id;
        std::size_t above_id;
        std::size_t primitive_index;
    };

    struct PrimitiveData {
        Primitive primitive;
        BoundingBox<Scalar, Dim> bounding_box;
        Vector<Scalar, Dim> centre;
        PrimitiveData(const Primitive& primitive):
            primitive(primitive),
            bounding_box(primitive.bounding_box()),
            centre(bounding_box.centre())
        {}
    };
public:
    void construct(const std::vector<Primitive>& primitives_in) {
        std::vector<PrimitiveData> primitive_data;
        typedef typename std::vector<PrimitiveData>::iterator iter_t;

        for (const auto& primitive: primitives_in) {
            primitive_data.emplace_back(primitive);
        }

        nodes.clear();
        root_id = nodes.size();
        nodes.emplace_back();

        struct State {
            std::size_t node_id;
            iter_t begin;
            iter_t end;
        };
        State root;
        root.node_id = root_id;
        root.begin = primitive_data.begin();
        root.end = primitive_data.end();

        std::stack<State> stack;
        stack.push(root);
        while (!stack.empty()) {
            State state = stack.top();
            stack.pop();

            Node& node = nodes[state.node_id];

            for (auto iter = state.begin; iter != state.end; iter++) {
                node.bounding_box.merge(iter->bounding_box);
            }

            if (state.end - state.begin == 1) {
                node.leaf = true;
                node.primitive_index = primitives_.size();
                primitives_.push_back(state.begin->primitive);
                continue;
            }
            node.leaf = false;

            // Simple method for now, partition by centre

            Scalar partition_extent = 0;
            for (int dim = 0; dim < Dim; dim++) {
                Scalar extent = node.bounding_box.size()(dim);
                if (extent > partition_extent) {
                    partition_extent = extent;
                    node.partition_dim = dim;
                }
            }

            std::sort(state.begin, state.end, [&node](auto& lhs, auto& rhs){
                return lhs.centre(node.partition_dim) < rhs.centre(node.partition_dim);
            });
            auto middle = state.begin + (state.end - state.begin) / 2;

            // Don't create new nodes until the current node is fully modified
            node.below_id = nodes.size();
            node.above_id = nodes.size() + 1;

            State below;
            below.node_id = node.below_id;
            below.begin = state.begin;
            below.end = middle;
            State above;
            above.node_id = node.above_id;
            above.begin = middle;
            above.end = state.end;

            stack.push(below);
            stack.push(above);
            nodes.resize(nodes.size() + 2);
        }
    }
    void query(const BoundingBox<Scalar, Dim>& bounding_box, std::vector<const Primitive*>& candidates) const {
        candidates.clear();
        std::stack<std::size_t> stack;
        stack.push(root_id);

        while (stack.size()) {
            std::size_t node_id = stack.top();
            stack.pop();
            const Node& node = nodes[node_id];

            if (!node.bounding_box.intersects(bounding_box)) {
                continue;
            }

            if (node.leaf) {
                candidates.push_back(&primitives_[node.primitive_index]);
                continue;
            }

            stack.push(node.below_id);
            stack.push(node.above_id);
        }
    }
    const std::vector<Primitive> primitives() const { return primitives_; }

    Transform<Scalar, Dim>& pose() { return pose_; }
    const Transform<Scalar, Dim>& pose() const { return pose_; }

    // For debuggging
    std::size_t node_count() const { return nodes.size(); }
    const BoundingBox<Scalar, Dim>& node_bounding_box(std::size_t i) const { return nodes[i].bounding_box; }

private:
    Transform<Scalar, Dim> pose_;
    std::size_t root_id;
    std::vector<Primitive> primitives_;
    std::vector<Node> nodes;
};

template <typename A, typename B>
requires primitives_compatible<A, B>
struct collision_details<Bvh<A>, Bvh<B>> {
    typedef typename A::Scalar Scalar;
    static constexpr int Dim = A::Dim;
    typedef std::vector<::owl::intersection_type<A, B>> intersection_type;

    static bool intersects(const Bvh<A>& a, const Bvh<B>& b) {
        std::vector<B*> b_candidates;
        for (auto a_candidate: a.primitives()) {
            a_candidate.transform(a.pose() * b.pose().inverse());
            BoundingBox<Scalar, Dim> a_box = a_candidate.bounding_box();
            b.query(a_box, b_candidates);
            for (auto b_candidate: b_candidates) {
                if (query(a_candidate, *b_candidate)) {
                    return true;
                }
            }
        }
        return false;
    }

    static std::optional<intersection_type> intersection(const Bvh<A>& a, const Bvh<B>& b) {
        intersection_type result;
        std::vector<const B*> b_candidates;

        for (auto a_candidate: a.primitives()) {
            a_candidate.transform(b.pose().inverse() * a.pose());

            BoundingBox<Scalar, Dim> a_box = a_candidate.bounding_box();
            b.query(a_box, b_candidates);

            for (auto b_candidate: b_candidates) {
                auto intersection = ::owl::intersection(a_candidate, *b_candidate);
                if (intersection.has_value()) {
                    intersection.value().transform(b.pose());
                    result.push_back(intersection.value());
                }
            }
        }
        if (result.empty()) {
            return std::nullopt;
        }
        return result;
    }
};

} // namespace owl
