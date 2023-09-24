#pragma once

#include "math/types/matrix.h"
#include "core/std/cloneable.h"
#include <memory>

namespace math {

class HeightMapInterface {
public:
    virtual bool contains(const math::Vector2d& position) const = 0;
    virtual double height(const math::Vector2d& position) const = 0;
    virtual math::Vector3d normal(const math::Vector2d& position) const {
        assert(false);
        // Not implemented
        return math::Vector3d::UnitZ();
    }
    virtual bool has_normal() const {
        return false;
    }
    virtual std::unique_ptr<HeightMapInterface> clone() const = 0;
};

class HeightMap: public core::OptionalCloneable<HeightMapInterface> {
public:
    HeightMap() {}
    template <typename Impl>
    HeightMap(const Impl& impl):
        core::OptionalCloneable<HeightMapInterface>(impl)
    {}
    bool contains(const math::Vector2d& position) const {
        return interface->contains(position);
    }
    double height(const math::Vector2d& position) const {
        return interface->height(position);
    }
    math::Vector3d normal(const math::Vector2d& position) const {
        return interface->normal(position);
    }
    bool has_normal() const {
        return interface->has_normal();
    }
    operator const HeightMapInterface& () const {
        return *interface;
    }
};

} // namespace math
