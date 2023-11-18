#pragma once

#include "mbox/types/matrix.h"
#include "cpp_utils/cloneable.h"
#include <memory>

namespace mbox {

class HeightMapInterface {
public:
    virtual bool contains(const Vector2d& position) const = 0;
    virtual double height(const Vector2d& position) const = 0;
    virtual Vector3d normal(const Vector2d& position) const {
        assert(false);
        // Not implemented
        return Vector3d::UnitZ();
    }
    virtual bool has_normal() const {
        return false;
    }
    virtual std::unique_ptr<HeightMapInterface> clone() const = 0;
};

class HeightMap: public cpp_utils::OptionalCloneable<HeightMapInterface> {
public:
    HeightMap() {}
    template <typename Impl>
    HeightMap(const Impl& impl):
        cpp_utils::OptionalCloneable<HeightMapInterface>(impl)
    {}
    bool contains(const Vector2d& position) const {
        return interface->contains(position);
    }
    double height(const Vector2d& position) const {
        return interface->height(position);
    }
    Vector3d normal(const Vector2d& position) const {
        return interface->normal(position);
    }
    bool has_normal() const {
        return interface->has_normal();
    }
    operator const HeightMapInterface& () const {
        return *interface;
    }
};

} // namespace mbox
