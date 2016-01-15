#ifndef ALIGNMENT2DFOCALLENGTHRESIDUAL_H
#define ALIGNMENT2DFOCALLENGTHRESIDUAL_H

class Alignment2DFocalLengthResidual {
public:
  Alignment2DFocalLengthResidual(double* observed_in, double fx, double fy, double px, double py) : 
  observed(observed_in), fx(fx), fy(fy), px(px), py(py) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
                  const T* const point,
                  const T* const focalLen,
                  T* residuals) const {
    T EPS = T(0.00001);
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint(camera_extrinsic, point, p);
    
    // camera_extrinsic[3,4,5] are the translation.
    p[0] += camera_extrinsic[3];
    p[1] += camera_extrinsic[4];
    p[2] += camera_extrinsic[5];
    
    // let p[2] ~= 0
    if (T(0.0)<=p[2]){
      if(p[2]<EPS){
        p[2] = EPS;
      }
    }else{
      if (p[2]>-EPS){
        p[2] = -EPS;
      }
    }
    
    // project it
    p[0] = T(focalLen[0]) * p[0] / p[2] + T(px);
    p[1] = T(focalLen[1]) * p[1] / p[2] + T(py);
    
    // reprojection error
    residuals[0] = T(observed[5])*(p[0] - T(observed[0]));
    residuals[1] = T(observed[5])*(p[1] - T(observed[1]));     
    
    return true;
  }
  
  double* observed;
  double fx, fy, px, py;
};

#endif