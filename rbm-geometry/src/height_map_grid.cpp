
#include "rbm/geometry/height_map_grid.h"


namespace rbm {

double HeightMapGrid::height(const Vector2d& position) const {
    auto result = height_data_.query(dimensions_.to_normalized(position));
    if (!height_data_.contains_index(result.a) || !height_data_.contains_index(result.b)) {
        return 0;
    }
    double z00 = height_data_.get(result.a);
    double z01 = height_data_.get(result.a + Vector2i(0, 1));
    double z10 = height_data_.get(result.a + Vector2i(1, 1));
    double z11 = height_data_.get(result.b);
    double z0 = z00 * (1 - result.fraction(1)) + z01 * result.fraction(1);
    double z1 = z10 * (1 - result.fraction(1)) + z11 * result.fraction(1);
    double z = z0 * (1 - result.fraction(0)) + z1 * result.fraction(0);
    return z;
}

Vector3d HeightMapGrid::normal(const Vector2d& position) const {
    Vector2i index = normal_data_.query_closest(dimensions_.to_normalized(position));
    if (!normal_data_.contains_index(index)) {
        return Vector3d::UnitZ();
    }
    return normal_data_.get(index);
}

} // namespace rbm
