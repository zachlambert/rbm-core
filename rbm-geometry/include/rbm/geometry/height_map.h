#pragma once

#include <memory>
#include <rbm/cpp/cloneable.h>
#include <rbm/types/matrix.h>

namespace rbm {

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

class HeightMap: public rbm::OptionalCloneable<HeightMapInterface> {
public:
    HeightMap() {}
    template <typename Impl>
    HeightMap(const Impl& impl):
        rbm::OptionalCloneable<HeightMapInterface>(impl)
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

} // namespace rbm
