#ifndef BUNDLEADJUSTMENTRESIDUAL_H
#define BUNDLEADJUSTMENTRESIDUAL_H

class BundleAdjustmentResidual {
public:
  BundleAdjustmentResidual(double *m_point_observed, double *m_weight):
    m_point_observed(m_point_observed), m_weight(m_weight) {}

  template <typename T>
  bool operator()(const T* const Rt,
                  const T* const m_point_predicted,
                  T* residuals) const {

    // T-ify the weight and points
    T weight[3] = {T(m_weight[0]), T(m_weight[1]), T(m_weight[2])};
    T point_observed[3] = {T(m_point_observed[0]), T(m_point_observed[1]), T(m_point_observed[2])};
    T point_predicted[3] = {m_point_predicted[0], m_point_predicted[1], m_point_predicted[2]};

    // Rt[0,1,2] are the angle-axis rotation.
    ceres::AngleAxisRotatePoint(Rt, point_observed, point_observed);

    // Rt[3,4,5] are the translation.
    point_observed[0] += Rt[3];
    point_observed[1] += Rt[4];
    point_observed[2] += Rt[5];

    // The error is the difference between the predicted and observed position.
    residuals[0] = weight[0] * (point_predicted[0] - point_observed[0]);
    residuals[1] = weight[1] * (point_predicted[1] - point_observed[1]);
    residuals[2] = weight[2] * (point_predicted[2] - point_observed[2]);

    return true;
  }

  double *m_point_observed;
  double *m_weight;
};

#endif