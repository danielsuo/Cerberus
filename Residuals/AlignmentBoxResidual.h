#ifndef ALIGNMENTBOXRESIDUAL_H
#define ALIGNMENTBOXRESIDUAL_H

class AlignmentBoxResidual {
public:
  AlignmentBoxResidual(double* observed_in, double* objectHalfSize, double* objectWeight) : 
  observed(observed_in), objectHalfSize(objectHalfSize), objectWeight(objectWeight) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
    const T* const world2object,
    T* residuals) const {

    T p[3];
    T q[3];  
    T o[3];  
    T r[3];

    p[0] = T(observed[2]) - camera_extrinsic[3];
    p[1] = T(observed[3]) - camera_extrinsic[4];
    p[2] = T(observed[4]) - camera_extrinsic[5];

    r[0] = - camera_extrinsic[0];
    r[1] = - camera_extrinsic[1];
    r[2] = - camera_extrinsic[2];
    
    ceres::AngleAxisRotatePoint(r, p, o);
    
    ceres::AngleAxisRotatePoint(world2object, o, q);

    q[0] += world2object[3];
    q[1] += world2object[4];
    q[2] += world2object[5];   

    // hinge loss
    double* szPtr = objectHalfSize + 3 * int(observed[5]);
    T oW = T(objectWeight[int(observed[5])]);

    if (q[0] < T(-szPtr[0])){
      residuals[0] = oW * (q[0] + T(szPtr[0]));
    }else if (q[0] < T(szPtr[0])){
      residuals[0] = T(0.0);
    }else{
      residuals[0] = oW * (q[0] - T(szPtr[0]));
    }

    if (q[1] < T(-szPtr[1])){
      residuals[1] = oW * (q[1] + T(szPtr[1]));

    }else if (q[1] < T(szPtr[1])){
      residuals[1] = T(0.0);
    }else{
      residuals[1] = oW * (q[1] - T(szPtr[1]));

    }

    if (q[2] < T(-szPtr[2])){
      residuals[2] = oW * (q[2] + T(szPtr[2]));
    }else if (q[2] < T(szPtr[2])){
      residuals[2] = T(0.0);
    }else{
      residuals[2] = oW * (q[2] - T(szPtr[2]));
    }    

    return true;
  }

  double* observed;
  double* objectHalfSize;
  double* objectWeight;
};

#endif