#pragma once

#include <unordered_map>
#include <unordered_set>
#include <concepts>
#include <span>
#include <assert.h>

// Probably won't end up using this, but leave here for now since it was tricky to do

#if 0

namespace geometry::topology {

template <typename T>
class ListNode {
public:
    ListNode():
        prev(nullptr),
        next(nullptr)
    {}
private:
    T* prev;
    T* next;
    template <typename T_, typename GetNode_, bool IsConst>
    friend class ListIterator_;
    template <typename T_, typename GetNode_>
    friend class WeakList;
    template <typename T_, typename GetNode_>
    friend class List;
};

template <typename T, typename GetNode, bool IsConst>
class ListIterator_ {
public:
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = std::conditional_t<IsConst, const T*, T*>;
    using reference = std::conditional_t<IsConst, const T&, T&>;

    ListIterator_(pointer ptr): ptr(ptr) {}

    reference operator*() const {
        return *ptr;
    }
    pointer operator->() const {
        return ptr;
    }
    operator pointer() const {
        return ptr;
    }

    ListIterator_& operator++() {
        ptr = GetNode()(ptr)->next;
        return *this;
    }
    ListIterator_ operator++(int) {
        ListIterator_ temp = (*this);
        ++(*this);
        return temp;
    }

    ListIterator_& operator--() {
        ptr = GetNode()(ptr)->prev;
        return *this;
    }
    ListIterator_ operator--(int) {
        ListIterator_ temp = (*this);
        --(*this);
        return temp;
    }

    friend bool operator==(const ListIterator_& lhs, const ListIterator_& rhs) { return lhs.ptr == rhs.ptr; }
    friend bool operator!=(const ListIterator_& lhs, const ListIterator_& rhs) { return !(lhs == rhs); }

private:
    pointer ptr;
};

template <typename T, typename GetNode>
using ConstListIterator = ListIterator_<T, GetNode, true>;
template <typename T, typename GetNode>
using ListIterator = ListIterator_<T, GetNode, false>;

template <typename T, typename GetNode>
class WeakList {
public:
    typedef ConstListIterator<T, GetNode> const_iterator;
    typedef ListIterator<T, GetNode> iterator;

    WeakList():
        head(nullptr),
        tail(nullptr),
        size_(0)
    {}

    void add(T* element)
    {
        ListNode<T>* node = GetNode()(element);
        if (!head) {
            assert(!tail);
            assert(size_ == 0);
            head = element;
            tail = element;
        }
        else {
            assert(size_ > 0);
            node->prev = tail;
            GetNode()(tail)->next = element;
        }
        size_++;
    }

    void remove(T* element) {
        ListNode<T>* node = GetNode()(element);
        assert(size > 0);
        if (node->prev) {
            GetNode()(node->prev)->next = node->next;
        }
        else {
            head = node->next;
        }
        if (node->next) {
            GetNode()(node->next)->prev = node->prev;
        }
        else {
            tail = node->prev;
        }
        size_--;
    }

    ListIterator<T, GetNode> begin() { return iterator(head); }
    ListIterator<T, GetNode> end() { return iterator(nullptr); }
    ListIterator<T, GetNode> begin() const { return const_iterator(head); }
    ListIterator<T, GetNode> end() const { return const_iterator(nullptr); }
    ListIterator<T, GetNode> cbegin() { return const_iterator(head); }
    ListIterator<T, GetNode> cend() { return const_iterator(nullptr); }

    ListIterator<T, GetNode> from(T* start) { return iterator(start); }
    ListIterator<T, GetNode> from(T* start) const { return const_iterator(start); }
    ListIterator<T, GetNode> cfrom(T* start) { return const_iterator(start); }

    std::size_t size() const { return size; }
    bool empty() const { return size == 0; }

private:
    T* head;
    T* tail;
    std::size_t size_;

    template <typename T_, typename GetNode_>
    friend class List;
};

template <typename T, typename GetNode>
class List {
public:
    ~List() {
        T* node = list.head;
        T* next;
        while (node) {
            next = GetNode()(node)->next;
            delete node;
            node = next;
        }
    }

    T* allocate() {
        T* node = new T();
        list.add(node);
        return node;
    }
    void free(T* node) {
        list.remove(node);
        delete node;
    }

    ListIterator<T, GetNode> begin() { return list.begin(); }
    ListIterator<T, GetNode> end() { return list.end(); }
    ListIterator<T, GetNode> begin() const { return list.begin(); }
    ListIterator<T, GetNode> end() const { return list.end(); }
    ListIterator<T, GetNode> cbegin() { return list.cbegin(); }
    ListIterator<T, GetNode> cend() { return list.cend(); }

    ListIterator<T, GetNode> from(T* start) { return iterator(start); }
    ListIterator<T, GetNode> from(T* start) const { return const_iterator(start); }
    ListIterator<T, GetNode> cfrom(T* start) { return const_iterator(start); }

    std::size_t size() const { return size; }
    bool empty() const { return size == 0; }

private:
    WeakList<T, GetNode> list;
};

struct Facet;

struct Node {
    std::size_t dim;
    Facet* facet;
    Node* parent;
    Node* other;

    ListNode<Node> facet_nodes_node;
    ListNode<Node> ridges_node;
    ListNode<Node> nodes_node;

    struct GetFacetNodesNode { ListNode<Node>* operator()(Node* node) { return &node->facet_nodes_node; } };
    struct GetRidgesNode { ListNode<Node>* operator()(Node* node) { return &node->ridges_node; } };
    struct GetNodesNode { ListNode<Node>* operator()(Node* node) { return &node->nodes_node; } };

    WeakList<Node, GetRidgesNode> ridges;

    Node():
        dim(0),
        facet(nullptr),
        parent(nullptr),
        other(nullptr)
    {}
};

struct Facet {
    std::size_t dim;
    WeakList<Node, Node::GetFacetNodesNode> nodes;
    ListNode<Facet> facets_node;
    struct GetFacetsNode { ListNode<Facet>* operator()(Facet* facet) { return &facet->facets_node; } };
};

typedef List<Facet, Facet::GetFacetsNode> FacetsList;
typedef List<Node, Node::GetNodesNode> NodesList;

template <std::size_t Dim>
class PolytopeRef {
    static_assert(Dim >= 1);
    static constexpr std::size_t VERTEX = 0;
    static constexpr std::size_t EDGE = 1;
    static constexpr std::size_t FACE = 2;
    static constexpr std::size_t CELL = 3;
public:
    PolytopeRef(
        std::span<FacetsList, Dim + 1> facets_array,
        std::span<NodesList, Dim + 1> nodes_array,
        Node& polytope
    ):
        facets_array(facets_array),
        nodes_array(nodes_array),
        polytope_(polytope)
    {}

    Node* add_node(Facet* facet, Node* parent_node = nullptr) {
        Node* node = nodes_array[facet->dim].allocate();
        node->dim = facet->dim;
        node->facet = facet;
        node->parent = parent_node;

        node->facet->nodes.add(node);
        if (node->parent) {
            node->parent->ridges.add(node);
        }
        if (parent_node) {
            Facet& facet = *node->facet;
            for (auto other = facet.nodes.begin(); other != facet.nodes.end(); other++) {
                if ((Node*)other == node) continue;
                assert(other->parent);
                if (node->facet == other->facet && node->parent->parent == other->parent->parent) {
                    assert(!other->other);
                    node->other = other;
                    other->other = node;
                    break;
                }
            }
        }

        return node;
    }

    void remove_node(Node* node) {
        auto ridge = node->ridges.begin();
        while (ridge != node->ridges.end()) {
            auto to_remove = ridge;
            ridge++;
            remove_node(to_remove);
        }
        assert(node->ridges.empty());

        node->facet->nodes.remove(node);
        if (node->parent) {
            node->parent->ridges.remove(node);
        }
        if (node->other) {
            node->other->other = nullptr;
        }
        nodes_array[node->dim].free(node);
    }

    Facet* add_facet(std::size_t dim) {
        assert(dim <= Dim);
        return facets_array[dim].allocate();
    }
    void remove_facet(Facet* facet) {
        assert(facet->nodes.empty());
        facets_array[facet->dim].free(facet);
    }

    Facet* add_vertex() { return add_facet(VERTEX); }
    Facet* add_edge() { return add_facet(EDGE); }
    Facet* add_face() { return add_facet(FACE); }
    Facet* add_cell() { return add_facet(CELL); }

    template <std::size_t Dim_>
    std::conditional_t<Dim_==0, Node&, PolytopeRef<Dim_>> scope_node(Node* node) {
        static_assert(Dim_ <= Dim);
        if constexpr (Dim_ == 0) {
            return *node;
        }
        if constexpr (Dim_ > 0) {
            return PolytopeRef<Dim_ - 1>(
                facets_array.template subspan<0, Dim_>(),
                nodes_array.template subspan<0, Dim_>(),
                *node
            );
        }
    }

    Node* polytope() { return &polytope_; }

private:
    std::span<FacetsList, Dim + 1> facets_array;
    std::span<NodesList, Dim + 1> nodes_array;
    Node& polytope_;
};

template <std::size_t Dim>
class Polytope: public PolytopeRef<Dim> {
public:
    Polytope():
        PolytopeRef<Dim>(std::span{facets_array}, std::span{nodes_array}, polytope)
    {}
private:
    std::array<FacetsList, Dim + 1> facets_array;
    std::array<NodesList, Dim + 1> nodes_array;
    Node polytope;
};

typedef Polytope<1> PointPair;
typedef PolytopeRef<1> PointPairRef;
typedef Polytope<2> Polygon;
typedef PolytopeRef<2> PolygonRef;
typedef Polytope<3> Polyhedron;
typedef PolytopeRef<3> PolyhedronRef;
typedef Polytope<4> Hypervolume;
typedef PolytopeRef<4> HypervolumeRef;


} // namespace geometry::topology

#endif

// NOTES

// An n-dimenional mesh, is a subdivision of a n-dimensional surface
// If this surface is "manifold" - it is "closed" with no holes, it represents the surface
// or a 3d volume, and thus can represent an (n+1) dimensional volume.
// Alternatively, an n-dimensional volume can be represented by an n-dimensional mesh

// For an n-dimensional mesh, whether this is a surface in (n+1)-dimensional space, or a
// "volume" in n-dimensional space, depends on the position data attaches to the vertices.
// eg: 2d mesh:
// - 2D surface in 3D space is 3D positions are attached to the vertices
// - 2D volume in 2D space if 2D positions are attached to the verties
// Independent of the vertex attributes, these have the same data structure, however
// the 2D surface may "wrap-around", so represents a topological space that isn't valid
// for a 2D volume in 2D space.

// An n-dimensional mesh subdivides n-dimensional space.
// It is made of up of n-dimensional facets
// eg: Facets of a 2D surface in 3D space are 2D polygons
//     Facets of a 3D volume in 3D space are 3D polyhedrons
//     Facets of a 1D surface in 2D space are 1D line segments
//     Facets of a 2D volume in 2D space are 2D polygons
//     Facets of a 1D volume in 1D space are line segments
//     Facets of a 0D surface in 1D space are points
// Note: The dimension of the facet is the dimension of the surface/volume being subdivided

// For an n-dimensional facet, this repesents a volume in n-dimensional space
// and is defined by a closed surface.
// Therefore an n-dimensional facet is defined by an (n-1)-dimensional mesh
// eg: A 2D polygon is defined by a 1D mesh made up of line segments (edges)
//     A 3D polyhedron is defined by a 2D mesh made of of faces
//     A 1D line segment is a defined by a 0D "mesh" made up of (two) points
// The "facets" of the mesh of a facet, are called ridges.
// - In order for a collection of facets to be a valid mesh, the union of these must represent the whole surface/volume
//   This requires that each ridge is shared by two facets
// - In general, a facet may be defined by any number of ridges. A linked list is suitable for defining this.
//   Each ridge can act as a node for the linked-list of ridges for BOTH facets.
//   ie: It stores the following information:
//   - Facet 1, and the adjacent ridges relative to that facet
//   - Facet 2, and the adjacent ridges relative to that facet
// - A ridge is "adjacent" to another, relative to a given facet if:
//   - Both ridges share a lower-dimensional facet (eg: edges share a vertex, or faces share an edge)
//   - Both ridges are adjacent to the same higher-dimensional facet
// - Therefore, the number of adjacent ridges is defined by the number of lower-dimensional facets
//   - An Edge will only have two adjacent edges, since an edge only has two vertices
//   - A face can have any number of adjacent edges, defined by the number of edges
// eg:
// - 2D mesh: An edge will store the adjacent edges for both connected faces, in this case always having 2 adjacent edges
// - 3D mesh: A face will store the adjacent faces for both connected cells, in this case the number of adjacent faces being dependent on the number of edges
// To make this feasible the following restriction is applied:
// - For a 3D mesh, all faces will have a fixed number of edges, eg: 3 or 4
// The mesh can be simplified further if the highest-dimension facet is restricted to a fixed number of ridges:
// eg:
// - a 2D mesh is only made up of triangles
// - a 3D mesh is only made up of cubes or tetrahedrons
// However, to keep it more generic, the highest-dimension facets will be allowed arbitrary number of ridges
