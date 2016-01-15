#ifndef ALIGNMENT3DRESIDUAL_H
#define ALIGNMENT3DRESIDUAL_H

class Alignment3DResidual {
public:
  Alignment3DResidual(double* m_point_observed) : m_point_observed(m_point_observed) {}

  template <typename T>
  bool operator()(const T* const Rt,
                  const T* const m_point_predicted,
                  T* residuals) const {
                  
    // Rt[0,1,2] are the angle-axis rotation.
    T point_predicted[3];
    
    ceres::AngleAxisRotatePoint(Rt, m_point_predicted, point_predicted);

    // Rt[3,4,5] are the translation.
    point_predicted[0] += Rt[3];
    point_predicted[1] += Rt[4];
    point_predicted[2] += Rt[5];
    
    // The error is the difference between the predicted and m_point_observed position.
    residuals[0] = T(m_point_observed[5]) * (point_predicted[0] - T(m_point_observed[2]));
    residuals[1] = T(m_point_observed[5]) * (point_predicted[1] - T(m_point_observed[3]));
    residuals[2] = T(m_point_observed[5]) * (point_predicted[2] - T(m_point_observed[4]));

    return true;
  }
  double* m_point_observed;
};

#endif