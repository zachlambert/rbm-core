#pragma once

#include <unordered_map>
#include <list>
#include <algorithm>
#include <stack>
#include "math/geometry/bounding.h"

// KdTree:
// - Recursively subdivides a vector space
// - Each node has two children, subdividing a particular dimension
//   at a particular position (called partition) along that dimension.
// - Calling down(point / box) will continue to traverse up the tree
//   until either:
//   - It reaches the leaf node that contains the point
//   - It reaches the smallest node that contains the box

#if 0 // Add back when required

namespace math
{

template <typename Scalar, int Dim, typename Element>
class KdTree {
    struct NodeChildren {
        Scalar partition
        int lower;
        int upper;
        NodeChildren(Scalar partition, int lower, int upper):
            partition(partition), lower(lower), upper(upper)
        {}
    };
    struct Node {
        BoundingBox<Scalar, Dim> box; // When axis-aligned
        std::list<Element> elements;
        int dim;
        std::optional<NodeChildren> children;
        Node(const BoundingBox<Scalar, Dim>& box, int dim):
            box(box),
            dim(dim),
        {}
    };
public:
    KdTree():
        next_id(0)
    {}
    template <bool IsConst>
    class Handle_ {
    public:
        bool down(const BoundingBox<Scalar, Dim>& box) {
            if (!node().children.has_value()) {
                return false;
            }
            const NodeChildren& children = node().children.value();
            if (box.upper(dim) < children.partition) {
                id = children.lower;
                return true;
            }
            if (box.lower(dim) > children.partition) {
                id = children.upper;
                return true;
            }
            return false;
        }
        bool down(const Vector<Scalar, Dim>& point) {
            if (!node().children.has_value()) {
                return false;
            }
            const NodeChildren& children = node().children.value();
            if (point(dim) < children.partition) {
                id = children.lower;
                return true;
            }
            else {
                id = children.upper;
                return true;
            }
        }
        int dim() const {
            return node().dim;
        }
        Handle_ lower() const {
            return Handle_(parent, node().children.value().lower);
        }
        Handle_ upper() const {
            return Handle_(parent, node().children.value().upper);
        }
        void create_children(Scalar partition) const {
            int next_dim = (dim + 1) % Dim;
            BoundingBox<Scalar, Dim> lower_box, upper_box;
            node().box.split(dim, partition, lower_box, upper_box);
            node().children = NodeChildren(
                partition,
                parent->create_node(lower_box, next_dim),
                parent->create_node(upper_box, next_dim)
            );
        }
        std::conditional_t<IsConst, const Data&, Data&> elements() const {
            return node().elements;
        }
    private:
        typedef std::conditional_t<IsConst, const KdTree*, KdTree*> parent_t;
        Handle_(parent_t parent, int id):
            parent(parent),
            id(id)
        {}
        std::conditional_t<IsConst, const Node&, Node&> node() const { return parent->nodes[id]; }
        parent_t parent;
        int id;
    };
    typedef Handle_<false> Handle;
    typedef Handle_<true> ConstHandle;

    void create_root(const BoundingBox<Scalar, Dim>& box) {
        assert(!root.has_value());
        root = create_node(box, 0);
    }
    Handle root() {
        return Handle(this, root);
    }
private:
    int create_node(const BoundingBox<Scalar, Dim>& box, int dim) {
        int id;
        do {
            id = next_id++;
        } while (!nodes.contains(id));
        nodes.try_emplace(id, box, dim);
        return id;
    }

    std::unordered_map<int, Node> nodes;
    int next_id;
    std::optional<int> root;
};

template <typename Scalar, int Dim, typename Element>
struct ElementPoint {
    Element element;
    BoundingBox<Scalar, Dim> box;
};

template <typename Scalar, int Dim, typename Element>
struct ElementPoint {
    Element element;
    Vector<Scalar, Dim> point;
};

template <typename Scalar, int Dim, typename Element>
void construct(KdTree<Scalar, Dim, Element>& tree, const std::vector<ElementBox<Scalar, Dim, Element>>& elements) {
    struct State {
        typename decltype(elements)::iterator begin;
        typename decltype(elements)::iterator end;
        BoundingBox<Scalar, Dim> box;
        int dim;
        KdTree<Scalar, Dim, Element>::Handle handle;
    };

    auto get_partition = [](const State& state) -> Scalar {
        std::sort(state.begin, state.end, [](auto& lhs, auto& rhs) {
            return lhs.box.lower(state.dim) < rhs.box.lower(state.dim);
        });
        std::size_t bottom = 0;
        std::size_t top = state.end - state.begin;
        std::size_t i = (top + bottom) / 2;
        while (true) {
            auto iter = state.begin + i;
            Scalar partition = iter->box.lower(state.dim);

            long below = 0;
            for (auto iter2 = state.begin; iter2 != iter; iter2++) {
                if (iter2->box.upper(state.dim) < partition) {
                    below++;
                }
            }
            long above = state.end - iter;
            long diff = above - below;
            if (std::abs(diff) <= 1) {
                return partition;
            }
            if (diff > 0) {
                bottom = i;
            }
            else {
                top = i;
            }
            i = (bottom + top) / 2;
        }
        return partition;
    };

    auto get_box = [&](decltype(elements)::const_iterator begin, decltype(elements)::const_iterator end) {
        BoundingBox<Scalar, Dim> box;
        for (auto iter = begin; iter != end; iter++) {
            box.merge(iter->box);
        }
        return box;
    };

    std::stack<State> stack;
    {
        State state;
        state.begin = elements.begin();
        state.end = elements.end();
        state.box = get_box(state.begin, state.end);

        tree.create_root();
        state.handle = tree.root();
        state.dim = state.handle.dim();

        stack.push(state);
    }

    while (!stack.empty()) {
        State state = stack.top();
        stack.pop();
        Scalar partition = get_partition(state);

        auto upper_end = std::partition(state.begin, state.end, [&](const auto& element){
            return element.box.upper(state.dim) < partition || element.box.lower(state.dim) > partition;
        });
        state.handle.elements().clear();
        for (auto iter = upper_end; iter != state.end; iter++) {
            state.handle.elements().insert_back(*iter);
        }

        if (upper_end == state.begin) {
            continue;
        }
        auto upper_begin = std::partition(state.begin, state.end, [&](const auto& element){
            return element.box.upper(state.dim) < partition;
        });

        state.handle.create_children(partition);

        State new_lower;
        new_lower.begin = state.begin;
        new_lower.end = upper_begin-1;
        new_lower.handle = state.handle.lower();
        new_lower.dim = new_lower.handle.dim();

        State new_upper;
        new_upper.begin = upper_begin;
        new_upper.end = upper_end;
        new_upper.handle = state.handle.upper();
        new_upper.dim = new_upper.handle.dim();

        state.box.split(state.dim, partition, new_lower.box, new_upper.box);

        stack.push(new_lower);
        stack.push(new_upper);
    }
}

template <typename Scalar, int Dim, typename Element>
void construct(KdTree<Scalar, Dim, Element>& tree, const std::vector<ElementPoint<Scalar, Dim, Element>>& elements) {

}

template <typename, Scalar, int Dim, typename Element>
std::list<Element>& find(KdTree<Scalar, Dim, Element>& tree, const BoundingBox<Scalar, Dim>& box) {
    auto handle = tree.root();
    while (handle.down(box)) {}
    return handle.elements();
}
template <typname Scalar, int Dim, typename Data>
const std::list<Element>& find(const KdTree<Scalar, Dim, Element>& tree, const BoundingBox<Scalar, Dim>& box) {
    auto handle = tree.root();
    while (handle.down(box)) {}
    return handle.elements();
}

template <typename Scalar, int Dim, typename Element>
const std::list<Element>& find(KdTree<Scalar, Dim, Element>& tree, const Vector<Scalar, Dim>& point) {
    auto handle = tree.root();
    while (handle.down(point)) {}
    return handle.elements();
}
template <typename Scalar, int Dim, typename Element>
const std::list<Element>& find(const KdTree<Scalar, Dim, Element>& tree, const Vector<Scalar, Dim>& point) {
    auto handle = tree.root();
    while (handle.down(point)) {}
    return handle.elements();
}

} // namespace math

#endif
