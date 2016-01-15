#include <cmath>
#include <cstdio>
#include <iostream>
#include "ceres/ceres.h"
#include "ceres/rotation.h"

using namespace std;

double* objectHalfSize;
double* objectWeight;
double* focalLen;
double fx, fy, px, py, w3Dv2D;

#define EPS T(0.00001)

int exe_time=0;

struct AlignmentErrorTriangulate {
  AlignmentErrorTriangulate(double* observed_in, double* camera_extrinsic_in): observed(observed_in), camera_extrinsic(camera_extrinsic_in) {}

  template <typename T>
  bool operator()(const T* const point, T* residuals) const {
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint((T*)(camera_extrinsic), point, p);

    // camera_extrinsic[3,4,5] are the translation.
    p[0] += T(camera_extrinsic[3]);
    p[1] += T(camera_extrinsic[4]);
    p[2] += T(camera_extrinsic[5]);
    
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
    p[0] = T(fx) * p[0] / p[2] + T(px);
    p[1] = T(fy) * p[1] / p[2] + T(py);
    
    // reprojection error
    residuals[0] = (p[0] - T(observed[0]));
    residuals[1] = (p[1] - T(observed[1]));     

    /*
      cout<<"p[0]="<<p[0]<<endl;
      cout<<"p[1]="<<p[1]<<endl;
      cout<<"observed[0]="<<observed[0]<<endl;
      cout<<"observed[1]="<<observed[1]<<endl;
      cout<<"residuals[0]="<<residuals[0]<<endl;
      cout<<"residuals[1]="<<residuals[1]<<endl;
      cout<<"--------------------------"<<endl;
    */
    return true;
  }
  
  double* observed;
  double* camera_extrinsic;
};




struct AlignmentError2D {
  AlignmentError2D(double* observed_in): observed(observed_in) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
                  const T* const point,
                  T* residuals) const {
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint(camera_extrinsic, point, p);
    /*
      T x = camera_extrinsic[0];
      T y = camera_extrinsic[1];
      T z = camera_extrinsic[2];
      T x2 = x*x;
      T y2 = y*y;
      T z2 = z*z;
      T w2 = T(1.0) - x2 - y2 - z2;
      T w  = sqrt(w2);
    
      p[0] = point[0]*(w2 + x2 - y2 - z2) - point[1]*(T(2.0)*w*z - T(2.0)*x*y) + point[2]*(T(2.0)*w*y + T(2.0)*x*z);
      p[1] = point[1]*(w2 - x2 + y2 - z2) + point[0]*(T(2.0)*w*z + T(2.0)*x*y) - point[2]*(T(2.0)*w*x - T(2.0)*y*z);
      p[2] = point[2]*(w2 - x2 - y2 + z2) - point[0]*(T(2.0)*w*y - T(2.0)*x*z) + point[1]*(T(2.0)*w*x + T(2.0)*y*z);
    */
    
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
    p[0] = T(fx) * p[0] / p[2] + T(px);
    p[1] = T(fy) * p[1] / p[2] + T(py);
    
    // reprojection error
    residuals[0] = T(observed[5])*(p[0] - T(observed[0]));
    residuals[1] = T(observed[5])*(p[1] - T(observed[1]));     
    
    /*
      if (exe_time<10000){
      exe_time++;
      cout<<"p[0]="<<p[0]<<endl;
      cout<<"p[1]="<<p[1]<<endl;
      cout<<"observed[0]="<<observed[0]<<endl;
      cout<<"observed[1]="<<observed[1]<<endl;
      cout<<"residuals[0]="<<residuals[0]<<endl;
      cout<<"residuals[1]="<<residuals[1]<<endl;
      cout<<"--------------------------"<<endl;
      }
    */
    
    return true;
  }
  
  double* observed;

};


struct AlignmentError3D {
  AlignmentError3D(double* observed_in): observed(observed_in) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
                  const T* const point,
                  T* residuals) const {
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint(camera_extrinsic, point, p);
    /*
      T x = camera_extrinsic[0];
      T y = camera_extrinsic[1];
      T z = camera_extrinsic[2];
      T x2 = x*x;
      T y2 = y*y;
      T z2 = z*z;
      T w2 = T(1.0) - x2 - y2 - z2;
      T w  = sqrt(w2);
    
      p[0] = point[0]*(w2 + x2 - y2 - z2) - point[1]*(T(2.0)*w*z - T(2.0)*x*y) + point[2]*(T(2.0)*w*y + T(2.0)*x*z);
      p[1] = point[1]*(w2 - x2 + y2 - z2) + point[0]*(T(2.0)*w*z + T(2.0)*x*y) - point[2]*(T(2.0)*w*x - T(2.0)*y*z);
      p[2] = point[2]*(w2 - x2 - y2 + z2) - point[0]*(T(2.0)*w*y - T(2.0)*x*z) + point[1]*(T(2.0)*w*x + T(2.0)*y*z);
    */
    
    // camera_extrinsic[3,4,5] are the translation.
    p[0] += camera_extrinsic[3];
    p[1] += camera_extrinsic[4];
    p[2] += camera_extrinsic[5];
    
    // The error is the difference between the predicted and observed position.
    residuals[0] = T(observed[5])*(p[0] - T(observed[2]));
    residuals[1] = T(observed[5])*(p[1] - T(observed[3]));
    residuals[2] = T(observed[5])*(p[2] - T(observed[4]));


    /*    
          if (exe_time<10){
          exe_time ++;
          cout<<"fx="<<fx<<endl;
          cout<<"fy="<<fy<<endl;
          cout<<"px="<<px<<endl;
          cout<<"py="<<py<<endl;
          cout<<"w3Dv2D="<<w3Dv2D<<endl;
          cout<<"p[0]="<<p[0]<<endl;
          cout<<"p[1]="<<p[1]<<endl;
          cout<<"p[2]="<<p[2]<<endl;
          cout<<"observed[0]="<<observed[0]<<endl;
          cout<<"observed[1]="<<observed[1]<<endl;
          cout<<"observed[2]="<<observed[2]<<endl;
          cout<<"observed[3]="<<observed[3]<<endl;
          cout<<"observed[4]="<<observed[4]<<endl;
          cout<<"residuals[0]="<<residuals[0]<<endl;
          cout<<"residuals[1]="<<residuals[1]<<endl;
          cout<<"residuals[2]="<<residuals[2]<<endl;
          cout<<"--------------------------"<<endl;
          }
    */


    return true;
  }
  double* observed;
};


template <class T>
T bandFunction(T x, T halfWin){
  return T(1.0)/(T(1.0)+exp(T(-100.0)*(x-halfWin))) + T(1.0)/(T(1.0)+exp(T(-100.0)*(-x-halfWin)));
}

template <class T>
T hingeLoss(T x, T halfWin){
  //return max(T(0), T(-halfWin)-x) + max(T(0), T(-halfWin)+x);
  if (x< -halfWin){
    //return -halfWin - x;
    return x + halfWin;
  }else if (x<halfWin){
    return 0;
  }else{
    return x-halfWin;
  }
}


struct AlignmentErrorBox {
  AlignmentErrorBox(double* observed_in): observed(observed_in) {}

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

    // quadratic function
    //residuals[0] = q[0];
    //residuals[1] = q[1];
    //residuals[2] = q[2];    

    // hinge loss
    double* szPtr = objectHalfSize + 3 * int(observed[5]);
    T oW = T(objectWeight[int(observed[5])]);

    //residuals[0] = hingeLoss<T>(q[0], T(szPtr[0]));
    //residuals[1] = hingeLoss<T>(q[1], T(szPtr[1]));
    //residuals[2] = hingeLoss<T>(q[2], T(szPtr[2])); 

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

    //double* szPtr = objectHalfSize + 3 * int(observed[5]);
    //residuals[0] = bandFunction<T>(q[0], T(szPtr[0]));
    //residuals[1] = bandFunction<T>(q[1], T(szPtr[1]));
    //residuals[2] = bandFunction<T>(q[2], T(szPtr[2]));

    //cout<<"szPtr[0]="<<szPtr[0]<<endl;
    //cout<<"szPtr[1]="<<szPtr[1]<<endl;
    //cout<<"szPtr[2]="<<szPtr[2]<<endl;
    //cout<<"q[0]="<<q[0]<<endl;
    //cout<<"q[1]="<<q[1]<<endl;
    //cout<<"q[2]="<<q[2]<<endl;

    //cout<<"bandFunction<T>(q[0], T(szPtr[0]))="<<bandFunction<T>(q[0], T(szPtr[0]))<<endl;
    //cout<<"bandFunction<T>(q[0], T(szPtr[0]))="<<bandFunction<T>(q[1], T(szPtr[1]))<<endl;
    //cout<<"bandFunction<T>(q[0], T(szPtr[0]))="<<bandFunction<T>(q[2], T(szPtr[2]))<<endl;
    /*
      if (exe_time<3){
      exe_time ++;

      cout<<"(q[1] < T(-szPtr[1]))"<<endl;

      cout<<"observed[2]="<<observed[2]<<endl;
      cout<<"observed[3]="<<observed[3]<<endl;
      cout<<"observed[4]="<<observed[4]<<endl;

      cout<<"p[0]="<<p[0]<<endl;
      cout<<"p[1]="<<p[1]<<endl;
      cout<<"p[2]="<<p[2]<<endl;

      cout<<"o[0]="<<o[0]<<endl;
      cout<<"o[1]="<<o[1]<<endl;
      cout<<"o[2]="<<o[2]<<endl;

      cout<<"world2object[0]="<<world2object[0]<<endl;
      cout<<"world2object[1]="<<world2object[1]<<endl;
      cout<<"world2object[2]="<<world2object[2]<<endl;
      cout<<"world2object[3]="<<world2object[3]<<endl;
      cout<<"world2object[4]="<<world2object[4]<<endl;
      cout<<"world2object[5]="<<world2object[5]<<endl;

      cout<<"camera_extrinsic[0]="<<camera_extrinsic[0]<<endl;
      cout<<"camera_extrinsic[1]="<<camera_extrinsic[1]<<endl;
      cout<<"camera_extrinsic[2]="<<camera_extrinsic[2]<<endl;
      cout<<"camera_extrinsic[3]="<<camera_extrinsic[3]<<endl;
      cout<<"camera_extrinsic[4]="<<camera_extrinsic[4]<<endl;
      cout<<"camera_extrinsic[5]="<<camera_extrinsic[5]<<endl;

      cout<<"q[0]="<<q[0]<<endl;
      cout<<"q[1]="<<q[1]<<endl;
      cout<<"q[2]="<<q[2]<<endl;

      cout<<endl;
      }
    */

    return true;
  }
  double* observed;
};

struct AlignmentError2D3D {
  AlignmentError2D3D(double* observed_in): observed(observed_in) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
                  const T* const point,
                  T* residuals) const {
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint(camera_extrinsic, point, p);
    /*
      T x = camera_extrinsic[0];
      T y = camera_extrinsic[1];
      T z = camera_extrinsic[2];
      T x2 = x*x;
      T y2 = y*y;
      T z2 = z*z;
      T w2 = T(1.0) - x2 - y2 - z2;
      T w  = sqrt(w2);
    
      p[0] = point[0]*(w2 + x2 - y2 - z2) - point[1]*(T(2.0)*w*z - T(2.0)*x*y) + point[2]*(T(2.0)*w*y + T(2.0)*x*z);
      p[1] = point[1]*(w2 - x2 + y2 - z2) + point[0]*(T(2.0)*w*z + T(2.0)*x*y) - point[2]*(T(2.0)*w*x - T(2.0)*y*z);
      p[2] = point[2]*(w2 - x2 - y2 + z2) - point[0]*(T(2.0)*w*y - T(2.0)*x*z) + point[1]*(T(2.0)*w*x + T(2.0)*y*z);
    */
    
    
    // camera_extrinsic[3,4,5] are the translation.
    p[0] += camera_extrinsic[3];
    p[1] += camera_extrinsic[4];
    p[2] += camera_extrinsic[5];
    
    // The error is the difference between the predicted and observed position.
    residuals[2] = T(observed[5])*(p[0] - T(observed[2]))*w3Dv2D;
    residuals[3] = T(observed[5])*(p[1] - T(observed[3]))*w3Dv2D;
    residuals[4] = T(observed[5])*(p[2] - T(observed[4]))*w3Dv2D;

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
    p[0] = T(fx) * p[0] / p[2] + T(px);
    p[1] = T(fy) * p[1] / p[2] + T(py);
    
    // reprojection error
    residuals[0] = T(observed[5])*(p[0] - T(observed[0]));
    residuals[1] = T(observed[5])*(p[1] - T(observed[1]));     
    
    /*
      if (exe_time<10){
      exe_time ++;
      cout<<"fx="<<fx<<endl;
      cout<<"fy="<<fy<<endl;
      cout<<"px="<<px<<endl;
      cout<<"py="<<py<<endl;
      cout<<"w3Dv2D="<<w3Dv2D<<endl;
      cout<<"p[0]="<<p[0]<<endl;
      cout<<"p[1]="<<p[1]<<endl;
      cout<<"observed[0]="<<observed[0]<<endl;
      cout<<"observed[1]="<<observed[1]<<endl;
      cout<<"residuals[0]="<<residuals[0]<<endl;
      cout<<"residuals[1]="<<residuals[1]<<endl;
      cout<<"residuals[2]="<<residuals[2]<<endl;
      cout<<"residuals[3]="<<residuals[3]<<endl;
      cout<<"residuals[4]="<<residuals[4]<<endl;
      cout<<"--------------------------"<<endl;
      }
    */
    return true;
  }
  double* observed;
};

struct AlignmentError2DfocalLen {
  AlignmentError2DfocalLen(double* observed_in): observed(observed_in) {}

  template <typename T>
  bool operator()(const T* const camera_extrinsic,
                  const T* const point,
                  const T* const focalLen,
                  T* residuals) const {
                  
    // camera_extrinsic[0,1,2] are the angle-axis rotation.
    T p[3];
    
    ceres::AngleAxisRotatePoint(camera_extrinsic, point, p);
    /*
      T x = camera_extrinsic[0];
      T y = camera_extrinsic[1];
      T z = camera_extrinsic[2];
      T x2 = x*x;
      T y2 = y*y;
      T z2 = z*z;
      T w2 = T(1.0) - x2 - y2 - z2;
      T w  = sqrt(w2);
    
      p[0] = point[0]*(w2 + x2 - y2 - z2) - point[1]*(T(2.0)*w*z - T(2.0)*x*y) + point[2]*(T(2.0)*w*y + T(2.0)*x*z);
      p[1] = point[1]*(w2 - x2 + y2 - z2) + point[0]*(T(2.0)*w*z + T(2.0)*x*y) - point[2]*(T(2.0)*w*x - T(2.0)*y*z);
      p[2] = point[2]*(w2 - x2 - y2 + z2) - point[0]*(T(2.0)*w*y - T(2.0)*x*z) + point[1]*(T(2.0)*w*x + T(2.0)*y*z);
    */
    
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
    
    /*
      if (exe_time<10000){
      exe_time++;
      cout<<"p[0]="<<p[0]<<endl;
      cout<<"p[1]="<<p[1]<<endl;
      cout<<"observed[0]="<<observed[0]<<endl;
      cout<<"observed[1]="<<observed[1]<<endl;
      cout<<"residuals[0]="<<residuals[0]<<endl;
      cout<<"residuals[1]="<<residuals[1]<<endl;
      cout<<"--------------------------"<<endl;
      }
    */
    
    return true;
  }
  
  double* observed;

};


int main(int argc, char** argv)
{
  // NOTE: Data from MATLAB is in column-major order

  cout << "Ba2D3D bundle adjuster in 2D and 3D. Writen by Jianxiong Xiao." << endl;
  cout << "Usage: EXE mode(1,2,3,5) w3Dv2D input_file_name output_file_name" << endl;

  // Get bundle adjustment mode
  int mode = atoi(argv[1]);

  // TODO: figure out
  w3Dv2D = atof(argv[2]);

  // Open input file
  FILE* fp = fopen(argv[3],"rb");
  if (fp == NULL) {
    cout << "fail to open file" << endl;
    return false;
  }

  // Get number of camera poses (i.e., frames; assumes that each pose
  // uses the same camera)
  unsigned int nCam;  fread((void*)(&nCam), sizeof(unsigned int), 1, fp);

  // Get number of points in the point cloud
  unsigned int nPts;  fread((void*)(&nPts), sizeof(unsigned int), 1, fp);

  // Get number of SIFT features
  unsigned int nObs;  fread((void*)(&nObs), sizeof(unsigned int), 1, fp);

  // Read camera intrinsic
  fread((void*)(&fx), sizeof(double), 1, fp);
  fread((void*)(&fy), sizeof(double), 1, fp);
  fread((void*)(&px), sizeof(double), 1, fp);
  fread((void*)(&py), sizeof(double), 1, fp);

  focalLen = new double [2];
  focalLen[0] = fx;
  focalLen[1] = fy;

  // Read all of the extrinsic matrices for the nCam camera
  // poses. Converts world coordinates into camera coordinates
  double* cameraRt = new double [12*nCam];
  fread((void*)(cameraRt), sizeof(double), 12*nCam, fp);

  // read initial 3D point position
  double* pointCloud = new double [3*nPts];
  fread((void*)(pointCloud), sizeof(double), 3*nPts, fp);

  // Read all of the SIFT feature transforms (TODO: what is coordinate frame?)
  // observation
  unsigned int* pointObservedIndex = new unsigned int [2*nObs];
  double* pointObservedValue = new double [6*nObs];
  fread((void*)(pointObservedIndex), sizeof(unsigned int), 2*nObs, fp);
  fread((void*)(pointObservedValue), sizeof(double), 6*nObs, fp);

  // finish reading
  fclose(fp);

  // output info
  cout << "Parameters: ";
  cout << "mode=" << mode << " ";
  cout << "w3Dv2D=" << w3Dv2D << endl;

  cout << "Meta Info: ";
  cout << "nCam=" << nCam << " ";
  cout << "nPts=" << nPts << " ";
  cout << "nObs=" << nObs << endl;

  cout << "Camera Intrinsic: ";
  cout << "fx=" << fx << " ";
  cout << "fy=" << fy << " ";
  cout << "px=" << px << " ";
  cout << "py=" << py << endl;

  // Construct camera parameters from camera matrix
  double* cameraParameter = new double [6*nCam];

  // Loop through all poses
  for(int cameraID = 0; cameraID < nCam; ++cameraID) {
    double* cameraPtr = cameraParameter + 6*cameraID;
    double* cameraMat = cameraRt + 12*cameraID;

    if (!(isnan(*cameraPtr))) {
      // Converting column-major order cameraMat into first three
      // elements of cameraPtr
      ceres::RotationMatrixToAngleAxis<double>(cameraMat, cameraPtr);

      // Grabbing translation
      cameraPtr[3] = cameraMat[9];
      cameraPtr[4] = cameraMat[10];
      cameraPtr[5] = cameraMat[11];
    }
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);

  //----------------------------------------------------------------

  for (unsigned int idObs = 0; idObs < nObs; ++idObs) {

    double* cameraPtr = cameraParameter + pointObservedIndex[2*idObs] * 6;
    double* observePtr = pointObservedValue+6*idObs;

    //if (observePtr[5] < 0){

      double* pointPtr  = pointCloud + pointObservedIndex[2*idObs+1] * 3;

      ceres::CostFunction* cost_function;
      switch (mode){
      case 1:
        // 2D triangulation
        cost_function = new ceres::AutoDiffCostFunction<AlignmentErrorTriangulate, 2, 3>(new AlignmentErrorTriangulate(observePtr,cameraPtr));
        problem.AddResidualBlock(cost_function,loss_function,pointPtr);
        break;
      case 2:
        // 2D bundle adjustment
        cost_function = new ceres::AutoDiffCostFunction<AlignmentError2D, 2, 6, 3>(new AlignmentError2D(observePtr));
        problem.AddResidualBlock(cost_function,loss_function,cameraPtr,pointPtr);
        break;
      case 3:
        // 3D bundle adjustment
        cost_function = new ceres::AutoDiffCostFunction<AlignmentError3D, 3, 6, 3>(new AlignmentError3D(observePtr));
        problem.AddResidualBlock(cost_function,loss_function,cameraPtr,pointPtr);
        break;
      case 4:
        // 3D bundle adjustment
        cost_function = new ceres::AutoDiffCostFunction<AlignmentError3D, 3, 6, 3>(new AlignmentError3D(observePtr));
        problem.AddResidualBlock(cost_function,loss_function,cameraPtr,pointPtr);
        break;
      case 5:
        // 5D bundle adjustment
        if (isnan(observePtr[2])){
          cost_function = new ceres::AutoDiffCostFunction<AlignmentError2D,   2, 6, 3>(new AlignmentError2D(observePtr));
        }else{
          cost_function = new ceres::AutoDiffCostFunction<AlignmentError2D3D, 5, 6, 3>(new AlignmentError2D3D(observePtr));
        }
        problem.AddResidualBlock(cost_function,loss_function,cameraPtr,pointPtr);
        break;
      case 6:
        // 2D bundle adjustment with focal length
        cost_function = new ceres::AutoDiffCostFunction<AlignmentError2DfocalLen, 2, 6, 3, 2>(new AlignmentError2DfocalLen(observePtr));
        problem.AddResidualBlock(cost_function,loss_function,cameraPtr,pointPtr,focalLen);
        break;

      }
      /*}else{
      double* objectPtr = objectParameter + int(observePtr[5]) * 6;

      ceres::CostFunction* cost_function;
      cost_function = new ceres::AutoDiffCostFunction<AlignmentErrorBox, 3, 6, 6>(new AlignmentErrorBox(observePtr));
      problem.AddResidualBlock(cost_function,loss_function,cameraPtr,objectPtr);

      }*/



  }

  //----------------------------------------------------------------

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
  // for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.max_num_iterations = 200;  
  options.minimizer_progress_to_stdout = true;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;  //ceres::SPARSE_SCHUR;  //ceres::DENSE_SCHUR;
  //  options.ordering_type = ceres::SCHUR;
  
  /*
    options.linear_solver_type = ceres::DENSE_SCHUR; //ceres::SPARSE_SCHUR; //ceres::DENSE_SCHUR; //ceres::SPARSE_NORMAL_CHOLESKY; //
    options.ordering_type = ceres::SCHUR;
    options.minimizer_progress_to_stdout = true;
    // New options
    //options.preconditioner_type = ceres::JACOBI; // ceres::IDENTITY
    options.num_linear_solver_threads = 12;
    //options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    //options.use_block_amd = true;
    //options.eta=1e-2;
    //options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
    //options.use_nonmonotonic_steps=false;
    */
  
  /*
    options.trust_region_strategy_type =  ceres::LEVENBERG_MARQUARDT; // DEFINE_string(trust_region_strategy, "lm", "Options are: lm, dogleg");
    options.eta = 1e-2; //  DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the accuracy of each linear solve of the truncated newton step. Changing this parameter can affect solve performance ");
    options.linear_solver_type = ceres::SPARSE_SCHUR; //DEFINE_string(solver_type, "sparse_schur", "Options are:  sparse_schur, dense_schur, iterative_schur, sparse_cholesky,  dense_qr, dense_cholesky and conjugate_gradients");
    options.preconditioner_type = ceres::JACOBI; //DEFINE_string(preconditioner_type, "jacobi", "Options are:  identity, jacobi, schur_jacobi, cluster_jacobi,  cluster_tridiagonal");
    options.sparse_linear_algebra_library =  ceres::SUITE_SPARSE; //DEFINE_string(sparse_linear_algebra_library, "suitesparse", "Options are: suitesparse and cxsparse");
    options.ordering_type = ceres::SCHUR; //DEFINE_string(ordering_type, "schur", "Options are: schur, user, natural");
    options.dogleg_type =  ceres::TRADITIONAL_DOGLEG; //DEFINE_string(dogleg_type, "traditional", "Options are: traditional, subspace");
    options.use_block_amd = true; //DEFINE_bool(use_block_amd, true, "Use a block oriented fill reducing ordering.");
    options.num_threads = 1; //DEFINE_int32(num_threads, 1, "Number of threads");
    options.linear_solver_min_num_iterations = 5; //DEFINE_int32(num_iterations, 5, "Number of iterations");
    options.use_nonmonotonic_steps = false; //DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use nonmonotic steps");
    //DEFINE_double(rotation_sigma, 0.0, "Standard deviation of camera rotation perturbation.");
    //DEFINE_double(translation_sigma, 0.0, "Standard deviation of the camera translation perturbation.");
    //DEFINE_double(point_sigma, 0.0, "Standard deviation of the point perturbation");
    //DEFINE_int32(random_seed, 38401, "Random seed used to set the state of the pseudo random number generator used to generate the pertubations.");
    */
  
  //ceres::Solve(options, &problem, NULL);
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  //cout << summary.FullReport() << endl;
  cout << summary.BriefReport() << endl;
  cout<<" fx: "<<focalLen[0]<<" fy: "<<focalLen[1]<<endl;

  // obtain camera matrix from parameters
  for(int cameraID=0; cameraID<nCam; ++cameraID){
    double* cameraPtr = cameraParameter+6*cameraID;
    double* cameraMat = cameraRt+12*cameraID;
    if (!(isnan(*cameraPtr))){
      ceres::AngleAxisToRotationMatrix<double>(cameraPtr, cameraMat);
      cameraMat[9]  = cameraPtr[3];
      cameraMat[10] = cameraPtr[4];
      cameraMat[11] = cameraPtr[5];
      //cout<<"cameraID="<<cameraID<<" : ";
      //cout<<"cameraPtr="<<cameraPtr[0]<<" "<<cameraPtr[1]<<" "<<cameraPtr[2]<<" "<<cameraPtr[3]<<" "<<cameraPtr[4]<<" "<<cameraPtr[5]<<endl;
    }
  }

  // write back result files

  FILE* fpout = fopen(argv[4],"wb");
  fwrite((void*)(&nCam), sizeof(unsigned int), 1, fpout);
  fwrite((void*)(&nPts), sizeof(unsigned int), 1, fpout);

  fwrite((void*)(cameraRt), sizeof(double), 12*nCam, fpout);
  fwrite((void*)(pointCloud), sizeof(double), 3*nPts, fpout);
  fwrite((void*)(focalLen), sizeof(double), 2, fpout);
  
  
  

  fclose (fpout);

  // clean up
  delete [] cameraRt;
  delete [] pointCloud;
  delete [] pointObservedIndex;
  delete [] pointObservedValue;
  delete [] cameraParameter;

  return 0;  
}