#pragma once

#include <Eigen/Core>
#include <array>
#include <concepts>
#include <string>
#include <cstring>
#include <optional>
#include <variant>
#include "math/types/matrix.h"
#include "math/types/color.h"

// TODO: Get this working, and add serialisation/deserialisation to .ply files

#if 0

namespace math  {

using scalar_t = std::variant<
    char,
    unsigned char,
    short,
    unsigned short,
    int,
    unsigned int,
    float,
    double
>;
std::size_t scalar_size(scalar_t scalar) {
    return std::visit([](auto scalar){ return sizeof(decltype(scalar)); }, scalar);
}
bool scalar_equal(scalar_t first, scalar_t second) {
    return std::visit([](auto first, auto second){
        return std::is_same_v<decltype(first), decltype(second)>;
    }, first, second);
}
bool scalar_convertible(scalar_t from, scalar_t to) {
    return std::visit([](auto from, auto to){
        typedef decltype(from) From;
        typedef decltype(to) To;
        if (std::is_floating_point_v<From>) {
            if (std::is_floating_point_v<To>) return false;
            // Allow narrowing conversion for floats
            return true;
        }
        else if (std::is_integral_v<From>) {
            if (!std::is_integral_v<To>) return false;
            // Don't allow narrowing conversion for integer types
            if (sizeof(To) < sizeof(From)) return false;
            if (std::is_signed_v<From>) {
                if (!std::is_signed_v<To>) return false;
            }
            else if (std::is_unsigned_v<From>) {
                if (!std::is_unsigned_v<To>) return false;
            }
            else {
                assert(false);
                return false;
            }
            return true;
        }
        else {
            assert(false);
            return false;
        }
    }, from, to);
}

struct PropertyDescription {
    std::string name;
    scalar_t data_type;
    std::optional<scalar_t> size_type;

    bool convertible_to(const PropertyDescription& property) const {
        if (!property.name == name) return false;
        if (!scalar_convertible(data_type, property.data_type));
        if (size_type.has_value()) {
            if (!property.size_type.has_value()) return false;
            if (!scalar_convertible(size_type.value(), property.size_type.value())) return false;
        }
        else { // if (!size_type.has_value())
            if (property.size_type.has_value()) return false;
        }
        return true;
    }
};
struct ElementDescription {
    std::string name;
    std::vector<PropertyDescription> properties;

    bool convertible_to(const ElementDescription& element) const {
        if (!element.name == name) return false;
        for (const auto& to: element.properties) {
            if (!supports_property(to)) return false;
        }
        return true;
    }
    bool contains(const PropertyDescription& property) const {
        for (const auto& from: properties) {
            if (from.convertible_to(property)) {
                return true;
            }
        }
        return false;
    }
};

namespace standard_elements {
    static constexpr std::string vertex = "vertex";
    static constexpr std::string face = "face";
    static constexpr std::string cell = "cell";
    static constexpr std::string edge = "edge";
    static constexpr std::string material = "material";
}
namespace standard_properties {
    static constexpr std::string vertex_index = "vertex_index";
}

class PolygonList {
public:
    class PropertyData;

    template <bool IsConst>
    class DataIterator {
    public:
        bool is_scalar() const {
            return !data->size_stride.has_value();
        }
        std::size_t size() const {
            assert(!is_scalar());
            return std::visit([this](auto size_scalar) -> std::size_t {
                typedef decltype(size_scalar) SizeType;
                return *(SizeType*)&data->data[i];
            }, description->size_type.value());
        }
        template <typename T>
        T get() const {
            assert(is_scalar());
            assert(scalar_convertible(description->data_type, T()));
            return std::visit([this](auto data_scalar) -> T {
                typedef decltype(data_scalar) DataType;
                return *(DataType*)&data->data[i];
            }, description->data_type);
        }
        template <typename T>
        T get(std::size_t index) const {
            assert(!is_scalar());
            assert(scalar_convertible(description->data_type, T()));
            return std::visit([this, index](auto data_scalar) -> T {
                typedef decltype(data_scalar) DataType;
                return *(DataType*)(&data->data[i] + data->size_stride.value() + index * data->data_stride);
            }, description->data_type);
        }
        DataIterator& operator++() {
            if (data->size_stride.has_value()) {
                // Must get size() before changing i
                std::size_t size_ = size();
                i += data->size_stride.value();
                i += size_ * data->data_stride;
            }
            else {
                i += data->data_stride;
            }
        }
        // TODO: Set current value, for scalar types
        // TODO: Create new value if at end, for scalar or list, so long as property_count < element_count
        friend bool operator==(const DataIterator& lhs, const DataIterator& rhs) {
            return lhs.data == rhs.data && lhs.i == rhs.i;
        }
    private:
        DataIterator(PropertyDescription* description, PropertyData* data, std::size_t i):
            description(description),
            data(data),
            i(i),
            size_type(size_type)
        {}
        PropertyDescription* description;
        PropertyData* data;
        std::size_t i;
        friend class PropertyHandle;
    };
    template <bool IsConst>
    class PropertyHandle {
    public:
        const PropertyDescription& description() const {
            return parent->elements[element_i].properties[i];
        }
        DataIterator begin() {
            const PropertyData& data = parent->property_datas[parent->element_datas[element_i].property_indices[i]];
            return DataIterator(&description(), &data, 0);
        }
        DataIterator end() {
            const PropertyData& data = parent->property_datas[parent->element_datas[element_i].property_indices[i]];
            return DataIterator(&description(), &data, data->bytes.size());
        }
        std::size_t property_count() const {
            return parent->property_datas[parent->element_datas[element_i].property_indices[i]].property_count
        }
        std::size_t element_count() const {
            return parent->element_datas[element_i].element_count;
        }
        bool is_complete() const {
            return property_count() == element_count();
        }
    private:
        ProperyHandle(const PolygonList* parent, std::size_t element_i, std::size_t i):
            parent(parent),
            element_i(element_i),
            i(i)
        {}
        const PolygonList* parent;
        std::size_t element_i;
        std::size_t i;
        friend class ElementHandle;
    }

    template <bool IsConst>
    class ElementHandle {
    public:
        ElementHandle& operator++() {
            i++;
            return *this;
        }
        ElementHandle operator++(int) {
            ElementHandle temp = *this;
            i++;
            return temp;
        }
        const ElementDescription& description() const {
            return parent->elements[i];
        }
        PropertyHandle begin() const {
            return PropertyHandle(parent, i, 0);
        }
        PropertyHandle end() const {
            return PropertyHandle(parent, i, description().properties.size());
        }
        PropertyHandle find(const std::string& name) const {
            auto iter = begin();
            while (iter != end() && iter.description().name != name) {}
            return iter;
        }
        PropertyHandle match(const PropertyDescription& property) const {
            auto iter = begin();
            while (iter != end() && !iter.description().convertible_to(property)) {}
            return iter;
        }

        std::size_t element_count() const {
            return parent->element_data[i].element_count;
        }
        void resize(std::size_t element_count) const {
            static_assert(!IsConst);
            parent->element_datas[i].element_count = element_count;
            // Don't resize properties here, since we can't resize property lists
            // Leave it to external code to fully populate data following a resize
            // or when creating a new property
        }

        PropertyHandle create_property(const Propertydescription& description) const {
            static_assert(!IsConst);
            parent->elements[i].properties.push_back(description);

            std::size_t property_index = parent->property_datas.size();
            parent->element_datas[i].property_indices.push_back(property_index);

            PropertyData property_data;
            property_data.data_stride = scalar_size(description.data_type);
            if (description.size_type.has_value()) {
                property_data.size_stride = scalar_size(description.size_type.value());
            }
            parent->property_datas.push_back(property_data);

            return PropertyHandle(parent, i, parent->elemets[i].properties.size() - 1);
        }
    private:
        typedef std::conditional_t<IsConst, const PolygonList*, PolygonList*> parent_t;
        ElementHandle(parent_t parent, std::size_t i):
            parent(parent),
            i(i)
        {}
        parent_t parent;
        std::size_t i;
        friend class PolygonList;
    };

    ElementHandle begin() const {
        return ElementHandle(this, 0);
    }
    ElementHandle end() const {
        return ElementHandle(this, elements.size());
    }
    ElementHandle find(const std::string& name) const {
        auto handle = begin();
        while (handle != end() && handle.description().name != name) {}
        return handle;
    }
    ElementHandle match(const ElementDescription& element) const {
        auto handle = begin();
        while (handle != end() && !handle.description().convertible_to(element)) {}
        return handle;
    }
private:
    std::vector<ElementDescription> elements;

    struct ElementData {
        std::vector<std::size_t> property_indices;
        std::size_t element_count;
        ElementData():
            element_count(0)
        {}
    };
    std::vector<ElementData> element_datas;
    struct PropertyData {
        std::size_t data_stride;
        std::optional<std::size_t> size_stride;
        std::vector<char> data;
        std::size_t property_count;
        PropertyData():
            data_stride(0),
            property_count(0)
        {}
    };
    std::vector<PropertyData> property_datas;
};

template <typaname Struct>
bool export_struct(const PolygonList& polygon_list, const ElementDescription& description, std::vector<Struct>& output) {
    auto element = polygon_list.match(description);
    if (element == polygon_list.end()) return false;

    // Scalar only
    {
        std::size_t offset = 0;
        for (const auto& property: description.properties) {
            if (property.size_type.has_value()) return false;
            offset += property.data_type;
        }
        if (offset != sizeof(Struct)) return false;
    }

    output.resize(element.element_count());

    std::size_t offset = 0;
    const std::size_t stride = sizeof(Struct);
    for (const auto& output_property: description.properties) {
        auto property = element.match(output_property);
        if (property == element.end()) {
            assert(false);
        }
        std::visit([&](auto output_value) {
            typedef decltype(output_value) OutputType;

            char* data = (char*)output.data() + offset;
            for (auto iter = property.begin(); iter != property.end(); iter++) {
                *(OutputType*)data = iter.get<OutputType>();
                data += stride;
            }

            offset += sizeof(OutputType);
        }, output_property.data_type);
    }
}

template <int Dim, typename IndexType>
struct Facet {
    std::array<IndexType, Dim> indices;
};

template <int Dim, typename IndexType>
ElementDescription facet_description() {
    static_assert(Dim >= 1 && Dim <= 3);

    ElementDescription element;
    switch (Dim) {
    case 1:
        element.name = standard_elements::edge;
    case 2:
        element.name = standard_elements::face;
    case 3:
        element.name = standard_elements::cell;
    default:
        break;
    }

    PropertyDescription property;
    property.name = standard_properties::vertex_index;
    property.data_type = IndexType();
    property.size_type = (unsigned char)0;
    element.properties.push_back(property);

    return element;
}

template <int Dim, typename IndexType>
bool export_facet(const PolygonList& polygon_list, std::vector<Facet<Dim, IndexType>>& output) {
    ElementDescription description = facet_description<Dim, IndexType>();
    auto element = polygon_list.match(description);
    if (element == polygon_list.end()) return false;

    output.resize(element.element_count());

    const std::size_t index_count = Dim + 1;

    auto property = element.begin().begin();
    auto iter = property.begin();
    for (auto& facet: output) {
        assert(iter != property.end());
        assert(iter.size() == index_count);
        for (std::size_t i = 0; i < index_count; i++) {
            facet.indices[i] = iter.get<IndexType>(i);
        }
        iter++;
    }
}

} // namespace math

#endif
