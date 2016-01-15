#ifndef ALIGNMENTTRIANGULATERESIDUAL_H
#define ALIGNMENTTRIANGULATERESIDUAL_H

class AlignmentTriangulateResidual {
public:
  AlignmentTriangulateResidual(double* m_point_observed, double* Rt, double fx, double fy, double px, double py) :
  m_point_observed(m_point_observed), Rt(Rt), fx(fx), fy(fy), px(px), py(py) {}

  template <typename T>
  bool operator()(const T* const m_point_predicted, T* residuals) const {

    T EPS = T(0.00001);

    // Rt[0,1,2] are the angle-axis rotation.
    T point_predicted[3];
    
    ceres::AngleAxisRotatePoint((T*)(Rt), m_point_predicted, point_predicted);

    // Rt[3,4,5] are the translation.
    point_predicted[0] += T(Rt[3]);
    point_predicted[1] += T(Rt[4]);
    point_predicted[2] += T(Rt[5]);
    
    // let point_predicted[2] ~= 0
    if (T(0.0) <= point_predicted[2]) {
      if (point_predicted[2] < EPS) {
        point_predicted[2] = EPS;
      }
    } else {
      if (point_predicted[2] > -EPS) {
        point_predicted[2] = -EPS;
      }
    }
    
    // project it
    point_predicted[0] = T(fx) * point_predicted[0] / point_predicted[2] + T(px);
    point_predicted[1] = T(fy) * point_predicted[1] / point_predicted[2] + T(py);
    
    // reprojection error
    residuals[0] = (point_predicted[0] - T(m_point_observed[0]));
    residuals[1] = (point_predicted[1] - T(m_point_observed[1]));     

    return true;
  }
  
  double* m_point_observed;
  double* Rt;
  double fx, fy, px, py;
};

#endif