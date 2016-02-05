#ifndef ALIGNMENT3DFIXEDRESIDUAL_H
#define ALIGNMENT3DFIXEDRESIDUAL_H

class Alignment3DFixedResidual {
public:
  Alignment3DFixedResidual(double *observed_in) : observed(observed_in){}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
    const T* const point,
    const T* const calibration,
    T* residuals) const {
    T pref[3];
    T pref2[3];
    T p[3];
    
    // camera_extrinsic[0,1,2] is x,z,angle
    pref[0] = point[0] - camera_extrinsic[0];
    pref[1] = point[1];
    pref[2] = point[2]-camera_extrinsic[1];

    // Rotate the points 
    pref2[0] = pref[0] * cos(camera_extrinsic[2]) - pref[2] * sin(camera_extrinsic[2]);
    pref2[1] = pref[1];
    pref2[2] = pref[0] * sin(camera_extrinsic[2]) + pref[2] * cos(camera_extrinsic[2]);

    // rotate from reference camera
    ceres::AngleAxisRotatePoint(calibration , pref2, p);
    p[0] += calibration[3];
    p[1] += calibration[4];
    p[2] += calibration[5];

    // The error is the difference between the predicted and observed position.
    residuals[0] = T(observed[5]) * (p[0] - T(observed[2]));
    residuals[1] = T(observed[5]) * (p[1] - T(observed[3]));
    residuals[2] = T(observed[5]) * (p[2] - T(observed[4]));
    return true;
  }

  double* observed;
};

#endif