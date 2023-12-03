#pragma once

#include "rbm/optional_tuple.h"
#include "rbm/diff/manifold.h"


template <typename Tuple, std::size_t Remaining>
struct optional_manifold_iterator {
    static int dynamic_dim(const Tuple& tuple, int dim) {
        static_assert(Remaining > 0);
        static constexpr int Index = Tuple::capacity() - Remaining;

        if (!tuple.template has_value<Index>()) return dim;
        dim += rbm::manifold_dynamic_dim(tuple.template value<Index>());

        return optional_manifold_iterator<Tuple, Remaining-1>::dynamic_dim(tuple, dim);
    }
    static void add(Tuple& tuple, const rbm::VectorXd& delta, int dim) {
        static_assert(Remaining > 0);
        static constexpr int Index = Tuple::capacity() - Remaining;

        if (!tuple.template has_value<Index>()) return;
        int element_dim = rbm::manifold_dynamic_dim(tuple.template value<Index>());
        rbm::manifold_add(tuple.template value<Index>(), delta.block(dim, 0, element_dim, 1));
        dim += element_dim;

        optional_manifold_iterator<Tuple, Remaining-1>::add(tuple, delta, dim);
    }
    static rbm::VectorXd difference(const Tuple& a, const Tuple& b, const rbm::VectorXd& delta, int dim) {
        static_assert(Remaining > 0);
        static constexpr int Index = Tuple::capacity() - Remaining;

        if (!a.template has_value<Index>()) return delta;
        assert(b.template has_value<Index>());

        int element_dim = rbm::manifold_dynamic_dim(b.template value<Index>());
        assert(rbm::manifold_dynamic_dim(b.template value<Index>()) == element_dim);

        delta.block(dim, 0, element_dim, 1) = rbm::manifold_difference(
            a.template value<Index>(),
            b.template value<Index>());
        dim += element_dim;

        return optional_manifold_iterator<Tuple, Remaining-1>::difference(a, b, delta, dim);
    }
};

template <typename Tuple>
struct optional_manifold_iterator<Tuple, 0> {
    static int dynamic_dim(const Tuple& tuple, int dim) {
        return dim;
    }
    static void add(Tuple& tuple, const rbm::VectorXd& delta, int dim) {
        // Do nothing
    }
    rbm::VectorXd difference(const Tuple& a, const Tuple& b, const rbm::VectorXd& delta, int dim) {
        return delta;
    }
};

template <typename ...Args>
struct rbm::manifold_details<rbm::OptionalTuple<Args...>> {
    typedef rbm::OptionalTuple<Args...> X;
    static constexpr int dim = Eigen::Dynamic;
    static int dynamic_dim(const X& x) {
        return optional_manifold_iterator<X, X::capacity()>::dynamic_dim(x, 0);
    }
    static void add(X& x, const rbm::VectorXd& delta) {
        optional_manifold_iterator<X, X::capacity()>::add(x, delta, 0);
    }
    static rbm::VectorXd difference(const X& a, const X& b) {
        rbm::VectorXd delta = rbm::VectorXd::Zero(dynamic_dim(a));
        return optional_manifold_iterator<X, X::capacity()>::difference(a, b, delta, 0);
    }
};

// Standard tuple

template <typename Tuple, std::size_t Remaining, int PrevDim>
struct tuple_manifold_iterator {

    static constexpr std::size_t Index = std::tuple_size_v<Tuple> - Remaining;

    using element_t = std::tuple_element_t<Index, Tuple>;
    static constexpr int element_dim = rbm::manifold_dim<element_t>;

    static constexpr int remaining_dim =
        element_dim == -1 ? -1 :
            element_dim + tuple_manifold_iterator<Tuple, Remaining-1, PrevDim + element_dim>::remaining_dim;

    static constexpr int dim =
        PrevDim == -1 || remaining_dim == -1 ? -1 :
            PrevDim + remaining_dim;

    static constexpr int new_prev_dim =
        PrevDim == -1 || element_dim == -1 ? -1 :
            PrevDim + element_dim;

    static int dynamic_dim(const Tuple& tuple, int current_dim)
    {
        static_assert(Remaining > 0);
        if constexpr(element_dim != -1) {
            current_dim += element_dim;
        }
        if constexpr(element_dim == -1) {
            current_dim += rbm::manifold_dynamic_dim(std::get<Index>(tuple));
        }
        return tuple_manifold_iterator<Tuple, Remaining-1, new_prev_dim>::dynamic_dim(tuple, current_dim);
    }

    static void add(Tuple& tuple, const rbm::Vectord<dim>& delta, int delta_index)
    {
        static_assert(Remaining > 0);
        if constexpr(dim != -1) {
            rbm::manifold_add(std::get<Index>(tuple), delta.template block<element_dim, 1>(delta_index, 0));
            delta_index += element_dim;
        }
        if constexpr(dim == -1) {
            int dynamic_dim = rbm::manifold_dynamic_dim(std::get<Index>(tuple));
            rbm::manifold_add(std::get<Index>(tuple), delta.block(delta_index, 0, dynamic_dim, 1));
            delta_index += dynamic_dim;
        }
        tuple_manifold_iterator<Tuple, Remaining-1, new_prev_dim>::add(tuple, delta, delta_index);
    }

    static rbm::Vectord<dim> difference(const Tuple& a, const Tuple& b, const rbm::Vectord<dim>& delta, int delta_index)
    {
        static_assert(Remaining > 0);

        if constexpr(dim != -1) {
            delta.template block<element_dim, 1>(delta_index, 0) = rbm::manifold_difference(
                a.template value<Index>,
                b.template value<Index>
            );
            delta_index += element_dim;
        }
        if constexpr(dim == -1) {
            int dynamic_dim_a = rbm::manifold_dynamic_dim(a.template value<Index>());
            int dynamic_dim_b = rbm::manifold_dynamic_dim(b.template value<Index>());
            assert(dynamic_dim_a == dynamic_dim_b);
            int dynamic_dim = dynamic_dim_a;

            delta.block(delta_index, 0, dynamic_dim, 1) = rbm::manifold_difference(
                a.template value<Index>,
                b.template value<Index>
            );
            delta_index += dynamic_dim;
        }

        return tuple_manifold_iterator<Tuple, Remaining-1, new_prev_dim>::difference(a, b, delta, delta_index);
    }
};

template <typename Tuple, int PrevDim>
struct tuple_manifold_iterator<Tuple, 0, PrevDim> {
    static constexpr int remaining_dim = 0;
    static constexpr int dim = PrevDim;
    static int dynamic_dim(const Tuple& tuple, int total_dim) {
        return total_dim;
    }
    static void add(Tuple& tuple, const rbm::VectorXd& delta, int delta_index) {
        // Do nothing
    }
    rbm::VectorXd difference(const Tuple& a, const Tuple& b, const rbm::VectorXd& delta, int delta_index) {
        return delta;
    }
};

template <typename ...Args>
struct rbm::manifold_details<std::tuple<Args...>> {
    typedef std::tuple<Args...> X;
    static constexpr int dim = tuple_manifold_iterator<X, std::tuple_size_v<X>, 0>::dim;
    static int dynamic_dim(const X& x) {
        return tuple_manifold_iterator<X, std::tuple_size_v<X>, 0>::dynamic_dim(x, 0);
    }
    static void add(X& x, const rbm::Vectord<dim>& delta) {
        tuple_manifold_iterator<X, std::tuple_size_v<X>, 0>::add(x, delta, 0);
    }
    static rbm::Vectord<dim> difference(const X& a, const X& b) {
        rbm::Vectord<dim> delta = rbm::Vectord<dim>::Zero(dynamic_dim(a));
        return tuple_manifold_iterator<X, std::tuple_size_v<X>, 0>::difference(a, b, delta, 0);
    }
};
