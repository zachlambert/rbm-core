#pragma once

namespace opt {

// TODO: Add back in if I need it for something

#if 0

template <int XDim, int HDim>
class Lagrangian: public ScalarFunction<XDim + HDim> {
public:
    static constexpr int ZDim = XDim + HDim;

    Lagrangian(const ScalarFunction<XDim>& f, const VectorFunction<XDim, HDim>& h):
        f(f),
        h(h)
    {}

    double operator()(const Vectord<ZDim>& z) const override {
        auto x = z.template head<XDim>();
        auto lambda = z.template tail<HDim>();
        return f(x) + lambda.dot(h(x));
    }

    Matrixd<1, ZDim> dx(const Vectord<ZDim>& z) const override {
        auto x = z.template head<XDim>();
        auto lambda = z.template tail<HDim>();
        Matrixd<1, ZDim> dz;
        dz.template block<1, XDim>(0, 0) = f.dx(x) + lambda.transpose() * h.dx(x);
        dz.template block<1, HDim>(0, XDim) = h(x).transpose();
        return dz;
    }

    Matrixd<ZDim, ZDim> dx2(const Vectord<ZDim>& z) const override{
        Vectord<XDim> x = z.template head<XDim>();
        Vectord<HDim> lambda = z.template tail<HDim>();
        Matrixd<HDim, XDim> dhdx = h.dx(x);

        Matrixd<ZDim, ZDim> dz2;

        dz2.template block<XDim, XDim>(0, 0) = f.dx2(x);
        for (std::size_t i = 0; i < HDim; i++) {
            dz2.template block<XDim, XDim>(0, 0) += lambda(i, 0) * h.dx2_i(i, x);
        }

        dz2.template block<HDim, XDim>(XDim, 0) = dhdx;
        dz2.template block<XDim, HDim>(0, XDim) = dhdx.transpose();

        dz2.template block<HDim, HDim>(XDim, XDim).setZero();

        return dz2;
    }
    
private:
    const ScalarFunction<XDim>& f;
    const VectorFunction<XDim, HDim>& h;
};

#endif

} // namespace opt