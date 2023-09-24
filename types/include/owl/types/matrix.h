#pragma once

#include <Eigen/Core>


namespace owl {

using Eigen::Vector;
using Eigen::RowVector;
using Eigen::Matrix;

template <int Size>
using Vectorf = Eigen::Vector<float, Size>;
template <int Size>
using Vectord = Eigen::Vector<double, Size>;
template <int Size>
using Vectori = Eigen::Vector<int, Size>;

template <int Size>
using RowVectorf = Eigen::RowVector<float, Size>;
template <int Size>
using RowVectord = Eigen::RowVector<double, Size>;

template <int Rows, int Cols>
using Matrixf = Eigen::Matrix<float, Rows, Cols>;
template <int Rows, int Cols>
using Matrixd = Eigen::Matrix<double, Rows, Cols>;


template <typename Scalar>
using Vector1 = Eigen::Vector<Scalar, 1>;
using Vector1f = Vector1<float>;
using Vector1d = Vector1<double>;
using Vector1i = Vector1<int>;

using Eigen::Vector2;
using Eigen::Vector2f;
using Eigen::Vector2d;
using Eigen::Vector2i;

using Eigen::Vector3;
using Eigen::Vector3f;
using Eigen::Vector3d;
using Eigen::Vector3i;

using Eigen::Vector4;
using Eigen::Vector4f;
using Eigen::Vector4d;
using Eigen::Vector4i;

using Eigen::VectorX;
using Eigen::VectorXf;
using Eigen::VectorXd;
using Eigen::VectorXi;


template <typename Scalar>
using RowVector1 = Eigen::RowVector<Scalar, 1>;
using RowVector1f = RowVector1<float>;
using RowVector1d = RowVector1<double>;

using Eigen::RowVector2;
using Eigen::RowVector2f;
using Eigen::RowVector2d;

using Eigen::RowVector3;
using Eigen::RowVector3f;
using Eigen::RowVector3d;

using Eigen::RowVector4;
using Eigen::RowVector4f;
using Eigen::RowVector4d;

using Eigen::RowVectorX;
using Eigen::RowVectorXf;
using Eigen::RowVectorXd;


template <typename Scalar>
using Matrix1 = Eigen::Matrix<Scalar, 1, 1>;
using Matrix1f = Matrix1<float>;
using Matrix1d = Matrix1<double>;

using Eigen::Matrix2;
using Eigen::Matrix2f;
using Eigen::Matrix2d;
using Eigen::Matrix2Xf;
using Eigen::Matrix2Xd;
using Eigen::MatrixX2f;
using Eigen::MatrixX2d;

using Eigen::Matrix3;
using Eigen::Matrix3f;
using Eigen::Matrix3d;
using Eigen::Matrix3Xf;
using Eigen::Matrix3Xd;
using Eigen::MatrixX3f;
using Eigen::MatrixX3d;

using Eigen::Matrix4;
using Eigen::Matrix4f;
using Eigen::Matrix4d;
using Eigen::Matrix4Xf;
using Eigen::Matrix4Xd;
using Eigen::MatrixX4f;
using Eigen::MatrixX4d;

using Eigen::MatrixX;
using Eigen::MatrixXf;
using Eigen::MatrixXd;


template <typename Scalar>
using Vector6 = Vector<Scalar, 6>;
typedef Vector6<float> Vector6f;
typedef Vector6<double> Vector6d;

template <typename Scalar>
using Matrix6 = Matrix<Scalar, 6, 6>;
typedef Matrix6<float> Matrix6f;
typedef Matrix6<double> Matrix6d;

template <typename Scalar>
using MatrixX6 = Matrix<Scalar, Eigen::Dynamic, 6>;
typedef MatrixX6<float> MatrixX6f;
typedef MatrixX6<double> MatrixX6d;

template <typename Scalar>
using Matrix6X = Matrix<Scalar, 6, Eigen::Dynamic>;
typedef Matrix6X<float> Matrix6Xf;
typedef Matrix6X<double> Matrix6Xd;

} // namespace owl
