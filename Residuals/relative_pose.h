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
                              double weight = 0) {
    distance_ = (position0 - position1).norm();
    weight_ = weight;
  }

  template <typename T>
  bool operator () (const T* const p0, const T* const p1, T* residuals) const {
    residuals[0] = T(0);
    for (int i = 0; i < 3; ++i) {
      residuals[0] += (p0[i] - p1[i]) * (p0[i] - p1[i]);
    }
    residuals[0] = ceres::sqrt(residuals[0]);
    return true;
  }

 private:
  double distance_;
  double weight_;
};

#endif //CERBERUS_RELATIVE_POSE_H
