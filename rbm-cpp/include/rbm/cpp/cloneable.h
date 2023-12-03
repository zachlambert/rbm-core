#pragma once

#include <memory>
#include <concepts>


namespace rbm {

template <typename Interface>
concept is_cloneable_interface = requires(const Interface& interface) {
    { interface.clone() } -> std::convertible_to<std::unique_ptr<Interface>>;
};

template <typename Interface>
class Cloneable {
public:
    template <typename Impl>
    Cloneable(const Impl& impl):
        interface(std::make_unique<Impl>(impl))
    {}

    template <typename Impl>
    Cloneable(Impl&& impl):
        interface(std::make_unique<Impl>(std::move(impl)))
    {}

    Cloneable(const Cloneable& other):
        interface(other.interface->clone())
    {}
    Cloneable(Cloneable&& other):
        interface(std::move(other.interface))
    {}

    Cloneable& operator=(const Cloneable& other) {
        interface = other.interface->clone();
        return *this;
    }
    Cloneable& operator=(Cloneable&& other) {
        interface = std::move(other.interface);
        return *this;
    }

protected:
    std::unique_ptr<Interface> interface;
};

template <typename Interface>
class OptionalCloneable {
public:
    OptionalCloneable():
        interface(nullptr)
    {}

    template <typename Impl>
    OptionalCloneable(const Impl& impl):
        interface(std::make_unique<Impl>(impl))
    {}

    template <typename Impl>
    OptionalCloneable(Impl&& impl):
        interface(std::make_unique<Impl>(std::move(impl)))
    {}

    OptionalCloneable(const OptionalCloneable& other):
        interface(other ? other.interface->clone() : nullptr)
    {}
    OptionalCloneable(OptionalCloneable&& other):
        interface(std::move(other.interface))
    {}

    OptionalCloneable& operator=(const OptionalCloneable& other) {
        if (other) {
            interface = other.interface->clone();
        }
        return *this;
    }
    OptionalCloneable& operator=(OptionalCloneable&& other) {
        interface = std::move(other.interface);
        return *this;
    }

    operator bool() const {
        return (bool)interface;
    }

protected:
    std::unique_ptr<Interface> interface;
};

} // namespace rbm
