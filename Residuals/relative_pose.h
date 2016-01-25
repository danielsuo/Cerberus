//
// Created by Fisher Yu on 1/24/16.
//

#ifndef CERBERUS_RELATIVE_POSE_H
#define CERBERUS_RELATIVE_POSE_H

#include <Eigen/Dense>

class RelativeDistanceCostFunctor {
 public:
  RelativeDistanceCostFunctor(const Eigen::Vector3d &position0,
                              const Eigen::Vector3d &position1,
                              double weight = 1) {
    distance_ = (position0 - position1).norm();
    weight_ = weight;
  }

  template <typename T>
  bool operator () (const T* const p0, const T* const p1, T* residuals) const {
    residuals[0] = T(0);
    for (int i = 0; i < 3; ++i) {
      residuals[0] += (p0[i] - p1[i]) * (p0[i] - p1[i]);
    }
    residuals[0] = weight_ * (ceres::sqrt(residuals[0]) - T(distance_));
    return true;
  }

 private:
  double distance_;
  double weight_;
};

// Square matrix multiplication. The first matrix is transposed.
template <int Size, typename T>
void mtxm(const T* const m0, const T* const m1, T* p) {
  for (int r = 0; r < Size; ++r) {
    for (int c = 0; c < Size; ++c) {
      *p = T(0);
      for (int i = 0; i < Size; ++i) {
        *p += m0[i * Size + r] * m1[i * Size + c];
      }
      ++p;
    }
  }
}

// Calculate angular distance between relative rotation
class RelativeRotationCostFunctor {
 public:
  RelativeRotationCostFunctor(const Eigen::Vector3d &angle_axis_0,
                              const Eigen::Vector3d &angle_axis_1,
                              double weight = 1) {
    relative_rotation_.resize(9);
    std::vector<double> r0(9), r1(9);
    ceres::AngleAxisToRotationMatrix(angle_axis_0.data(), r0.data());
    ceres::AngleAxisToRotationMatrix(angle_axis_1.data(), r1.data());
    mtxm<3>(r0.data(), r1.data(), relative_rotation_.data());
    weight_ = weight;
  }

  template <typename T>
  bool operator () (const T* const aa0, const T* const aa1,
                    T* residuals) const {
    T r0[9], r1[9], r0tr1[9], rtr[9];
    for (int i = 0; i < 9; ++i) {
      rtr[i] = T(relative_rotation_[i]);
    }
    ceres::AngleAxisToRotationMatrix(aa0, r0);
    ceres::AngleAxisToRotationMatrix(aa1, r1);
    mtxm<3>(r0, r1, r0tr1);
    // Reuse memory
    mtxm<3>(rtr, r0tr1, r0);
    ceres::RotationMatrixToAngleAxis(r0, r1);
    residuals[0] = T(weight_) *
        ceres::sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);

    return true;
  }

 private:
  std::vector<double> relative_rotation_;
  double weight_;
};

class CircleConstraintCostFunctor {
 public:
  CircleConstraintCostFunctor(const Eigen::Vector3d &center,
                              const Eigen::Vector3d &normal,
                              double radius,
                              double weight = 1):
      center_(center), normal_(normal), radius_(radius), weight_(weight) {
    normal_.normalize();
  }

  template <typename T>
  bool operator () (const T* const position_ptr, T* residuals) const {
    typedef Eigen::Matrix<T, 3, 1> VectorT;
    VectorT position, direction, rejection;
    VectorT center_t, normal_t;
    for (int i = 0; i < 3; ++i) {
      position[i] = position_ptr[i];
      center_t[i] = T(center_[i]);
      normal_t[i] = T(normal_[i]);
    }
    direction = position - center_t;
    rejection = direction - direction.dot(normal_t) * normal_t;
    residuals[0] = weight_ * (rejection.norm() - radius_);
    return true;
  }

 private:
  Eigen::Vector3d center_;
  Eigen::Vector3d normal_;
  double radius_;
  double weight_;
};

class LineConstraintCostFunctor {
 public:
  LineConstraintCostFunctor(const Eigen::Vector3d &start,
                            const Eigen::Vector3d &direction,
                            double weight = 1) :
      start_(start), direction_(direction), weight_(weight) {
    direction_.normalize();
  }

  template <typename T>
  bool operator () (const T* const position_ptr, T* residuals) const {
    typedef Eigen::Matrix<T, 3, 1> VectorT;
    VectorT position, direction, rejection;
    VectorT start_t, line_direction_t;
    for (int i = 0; i < 3; ++i) {
      position[i] = position_ptr[i];
      start_t[i] = T(start_[i]);
      line_direction_t[i] = T(direction_[i]);
    }
    direction = position - start_t;
    rejection = direction - direction.dot(line_direction_t) * line_direction_t;
    residuals[0] = weight_ * rejection.norm();
    return true;
  }

 private:
  Eigen::Vector3d start_;
  Eigen::Vector3d direction_;
  double weight_;
};

#endif //CERBERUS_RELATIVE_POSE_H
