#pragma once

#include "math/probability/distribution.h"
#include "math/types/matrix.h"
#include "math/algebra/manifold.h"


namespace math {

template <typename T>
class GaussianDistribution: public Distribution<T> {
    static constexpr int Dim = math::manifold_dim<T>;
public:
    GaussianDistribution():
        dist_(0.0, 1.0)
    {
        set_mean(Vectord<Dim>::Zero());
        set_covariance(Matrixd<Dim, Dim>::Identity());
    }
    GaussianDistribution(const Vectord<Dim>& mean, const Matrixd<Dim, Dim>& covariance):
        dist_(0.0, 1.0)
    {
        set_mean(mean);
        set_covariance(covariance);
    }

    T sample(std::default_random_engine& rng) const override {
        VectorXd delta(manifold_dynamic_dim(mean_));
        for (std::size_t i = 0; i < delta.size(); i++) {
            delta(i) = dist_(rng);
        }
        return math::manifold_add(mean_, covariance_L_ * delta);
    }
    virtual double evaluate(const T& x) const override {
        Vectord<Dim> delta = math::manifold_difference(mean_, x);
        return normalization_ * std::exp(-0.5 * (delta.transpose() * covariance_inv_ * delta).value());
    }
    virtual std::unique_ptr<DistributionImpl<T>> clone() const override {
        return std::make_unique<GaussianDistribution<Vectord<Dim>>>(*this);
    }

    const T& mean() const { return mean_; }
    const Matrixd<Dim, Dim>& covariance() const { return covariance_; }
    const Matrixd<Dim, Dim>& covariance_inv() const { return covariance_inv_; }
    const Matrixd<Dim, Dim>& covariance_L() const { return covariance_L_; }
    double normalization() const { return normalization_; }

    void set_mean(const T& mean) {
        mean_ = mean;
    }
    void set_covariance(const Matrixd<Dim, Dim>& covariance) {
        covariance_ = covariance;
        covariance_L_ = Eigen::LLT<Matrixd<Dim, Dim>>(covariance_).matrixL();
        covariance_inv_ = covariance.inverse();
        normalization_ = 1.0 / std::sqrt(std::pow(2 * M_PI, Dim) * covariance_.determinant());
    }

private:
    mutable std::normal_distribution<double> dist_;
    T mean_;
    Matrixd<Dim, Dim> covariance_;
    Matrixd<Dim, Dim> covariance_inv_;
    Matrixd<Dim, Dim> covariance_L_; // From llt decomposition: cov = LL^T
    double normalization_;
};


// Override if using double

template <int Dim>
class GaussianDistribution<Vectord<Dim>>: public DistributionImpl<Vectord<Dim>> {
    typedef Vectord<Dim> T;
public:
    GaussianDistribution():
        dist_(0.0, 1.0)
    {
        if constexpr(Dim != -1) {
            set_mean(Vectord<Dim>::Zero());
            set_covariance(Matrixd<Dim, Dim>::Identity());
        }
    }
    GaussianDistribution(const Vectord<Dim>& mean, const Matrixd<Dim, Dim>& covariance):
        dist_(0.0, 1.0)
    {
        set_mean(mean);
        set_covariance(covariance);
    }

    T sample(std::default_random_engine& rng) const override {
        Vectord<Dim> delta(mean_.size());
        for (std::size_t i = 0; i < delta.size(); i++) {
            delta(i) = dist_(rng);
        }
        return mean_ + covariance_L_ * delta;
    }
    virtual double evaluate(const T& x) const override {
        Vectord<Dim> delta = x - mean_;
        return normalization_ * std::exp(-0.5 * (delta.transpose() * covariance_inv_ * delta).value());
    }
    virtual std::unique_ptr<DistributionImpl<T>> clone() const override {
        return std::make_unique<GaussianDistribution<Vectord<Dim>>>(*this);
    }

    const T& mean() const { return mean_; }
    const Matrixd<Dim, Dim>& covariance() const { return covariance_; }
    const Matrixd<Dim, Dim>& covariance_inv() const { return covariance_inv_; }
    const Matrixd<Dim, Dim>& covariance_L() const { return covariance_L_; }
    double normalization() const { return normalization_; }

    void set_mean(const T& mean) {
        mean_ = mean;
    }
    void set_covariance(const Matrixd<Dim, Dim>& covariance) {
        covariance_ = covariance;
        covariance_L_ = Eigen::LLT<Matrixd<Dim, Dim>>(covariance_).matrixL();
        covariance_inv_ = covariance.inverse();
        normalization_ = 1.0 / std::sqrt(std::pow(2 * M_PI, Dim) * covariance_.determinant());
    }

private:
    mutable std::normal_distribution<double> dist_;
    T mean_;
    Matrixd<Dim, Dim> covariance_;
    Matrixd<Dim, Dim> covariance_inv_;
    Matrixd<Dim, Dim> covariance_L_; // From llt decomposition: cov = LL^T
    double normalization_;
};


template <>
class GaussianDistribution<double>: public DistributionImpl<double> {
public:
    GaussianDistribution():
        dist_(0.0, 1.0)
    {
        set_mean(0);
        set_covariance(1);
    }
    GaussianDistribution(double mean, double covariance):
        dist_(0.0, 1.0)
    {
        set_mean(mean);
        set_covariance(covariance);
    }

    double sample(std::default_random_engine& rng) const override {
        return mean_ + covariance_L_ * dist_(rng);
    }
    virtual double evaluate(const double& x) const override {
        double delta = x - mean_;
        return normalization_ * std::exp(-0.5 * delta * covariance_inv_ * delta);
    }
    virtual std::unique_ptr<DistributionImpl<double>> clone() const override {
        return std::make_unique<GaussianDistribution<double>>(*this);
    }

    const double& mean() const { return mean_; }
    const double& covariance() const { return covariance_; }
    const double& covariance_inv() const { return covariance_inv_; }
    const double& covariance_L() const { return covariance_L_; }
    double normalization() const { return normalization_; }

    void set_mean(double mean) {
        mean_ = mean;
    }
    void set_covariance(double covariance) {
        covariance_ = covariance;
        covariance_L_ = std::sqrt(covariance);
        covariance_inv_ = 1.0 / covariance;
        normalization_ = 1.0 / std::sqrt(2 * M_PI * covariance_);
    }

private:
    mutable std::normal_distribution<double> dist_;
    double mean_;
    double covariance_;
    double covariance_inv_;
    double covariance_L_;
    double normalization_;
};


} // namespace math
