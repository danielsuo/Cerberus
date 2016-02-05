#include <cmath>
#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include "ceres/rotation.h"
#include "glog/logging.h"

#include "Cerberus.h"
#include "Residuals/BundleAdjustmentResidual.h"
#include "Residuals/PoseGraphResidual.h"

using namespace std;

int NUM_BA_POINTS;

int main(int argc, char** argv)
{
  google::InitGoogleLogging(argv[0]);

  // NOTE: Data from MATLAB is in column-major order

  cout << "Ba2D3D bundle adjuster in 2D and 3D. Writen by Jianxiong Xiao." << endl;
  cout << "Usage: EXE mode(1,2,3,5) w3Dv2D input_file_name output_file_name" << endl;

  // Get bundle adjustment mode
  int mode = atoi(argv[1]);

  // TODO: figure out
  NUM_BA_POINTS = atof(argv[2]);

  // Open input file
  FILE* fp = fopen(argv[3],"rb");
  if (fp == NULL) {
    cout << "fail to open file" << endl;
    return false;
  }

  // Get number of camera poses (i.e., frames; assumes that each pose
  // uses the same camera)
  uint32_t nCam;  fread((void*)(&nCam), sizeof(uint32_t), 1, fp);

  // Get number of matched pairs
  uint32_t nPairs; fread((void*)(&nPairs), sizeof(uint32_t), 1, fp);

  // Get number of matched key points
  uint32_t nMatchedPoints; fread((void*)(&nMatchedPoints), sizeof(uint32_t), 1, fp);

  // Read all of the extrinsic matrices for the nCam camera
  // poses. Converts camera coordinates into world coordinates
  double* cameraRtC2W = new double [12*nCam];
  cout << "Reading cameraRtC2W" << endl;
  fread((void*)(cameraRtC2W), sizeof(double), 12*nCam, fp);

  // Read all relative poses for each matched pair
  double *cameraRt_ij = new double[12*nPairs];
  cout << "Reading cameraRt_ij" << endl;
  fread((void*)(cameraRt_ij), sizeof(double), 12*nPairs, fp);

  //Read all camera pose indices for matched pairs
  uint32_t *cameraRt_ij_indices = new uint32_t [2*nPairs];
  cout << "Reading cameraRt_ij_indices" << endl;
  fread((void*)(cameraRt_ij_indices), sizeof(uint32_t), 2*nPairs, fp);

  // Read all of the matched key points for a given pair (3D
  // coordinates in i's coordinates)
  double *cameraRt_ij_points_observed_i = new double[3*nMatchedPoints];
  cout << "Reading cameraRt_ij_points_observed_i" << endl;
  fread((void*)(cameraRt_ij_points_observed_i), sizeof(double), 3*nMatchedPoints, fp);

  // Read all of the matched key points for a given pair (3D
  // coordinates in j's coordinates)
  double *cameraRt_ij_points_observed_j = new double[3*nMatchedPoints];
  cout << "Reading cameraRt_ij_points_observed_j" << endl;
  fread((void*)(cameraRt_ij_points_observed_j), sizeof(double), 3*nMatchedPoints, fp);

  // Read all of the matched key points for a given pair (3D
  // coordinates in world frame)
  double *cameraRt_ij_points_predicted = new double[3*nMatchedPoints];
  cout << "Reading cameraRt_ij_points_predicted" << endl;
  fread((void*)(cameraRt_ij_points_predicted), sizeof(double), 3*nMatchedPoints, fp);

  // Read the count of matched key points per pair so we can index
  // correctly
  uint32_t *cameraRt_ij_points_count = new uint32_t [nPairs];
  cout << "Reading cameraRt_ij_points_count" << endl;
  fread((void*)(cameraRt_ij_points_count), sizeof(uint32_t), nPairs, fp);

  // finish reading
  fclose(fp);

  // Construct camera parameters from camera matrix
  double *cameraParameter = new double [6*nCam];

  cout << "Converting camera parameters" << endl;
  // Loop through all poses
  for(int cameraID = 0; cameraID < nCam; cameraID++) {
    double *cameraPtr = cameraParameter + 6 * cameraID;
    double *cameraMat = cameraRtC2W + 12 * cameraID;

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

  double* cameraParameter_ij = new double[6*nPairs];
  for (int cameraPairID = 0; cameraPairID < nPairs; cameraPairID++) {
    double *cameraPtr = cameraParameter_ij + 6 * cameraPairID;
    double *cameraMat = cameraRt_ij + 12 * cameraPairID;

    if (!(isnan(*cameraPtr))) {
      ceres::RotationMatrixToAngleAxis<double>(cameraMat, cameraPtr);
      cameraPtr[3] = cameraMat[9];
      cameraPtr[4] = cameraMat[10];
      cameraPtr[5] = cameraMat[11];
    }
  }

  Cerberus *c1 = new Cerberus();
  Cerberus *c2 = new Cerberus();

  // Which matched point are we on?
  uint32_t points_count = 0;
  uint32_t cameraRt_ij_points_total = 0;
  double *points_observed_i = new double[3 * NUM_BA_POINTS * nPairs];
  double *points_observed_j = new double[3 * NUM_BA_POINTS * nPairs];
  double *points_predicted = new double[3 * NUM_BA_POINTS * nPairs];

  cout << "Building problem" << endl;

  double sum = 0;
  int printCount = 0;

  // TODO: compute inverse outside
  for (int i = 0; i < nPairs; i++) {
    // if (i > 50) continue;
    // if (cameraRt_ij_indices[2*i+1] - cameraRt_ij_indices[2*i] > 1) continue;
    /*
     * Build Pose Graph problem
     */

    // Get initial values for Rt_i and Rt_j. Note that we subtract one
    // from our offset because MATLAB is 1-indexed. Note that these
    // transforms are camera to world coordinates
     double *Rt_i = cameraParameter + 6 * (cameraRt_ij_indices[2 * i] - 1);
     double *Rt_j = cameraParameter + 6 * (cameraRt_ij_indices[2 * i + 1] - 1);

    // Get relative pose between matched frames
     double *Rt_ij = cameraParameter_ij + 6 * i;

    // Overweight time-based alignment
     double w = cameraRt_ij_indices[2 * i + 1] - cameraRt_ij_indices[2 * i] == 1 ? 50.0 : 1.0;
     double weight_PoseGraph[9] = {w, w, w, w, w, w, w, w, w};

     double residual[9];

    // calc_residual(Rt_ij, Rt_i, Rt_j, weight, residual, i < iters);

    // for (int k = 0; k < 9; k++) {
    //   sum += pow(residual[k], 2);
    // }

    // cout << sum << endl;
    // cout << endl;

     PoseGraphResidual::AddResidualBlock(c1, Rt_ij, weight_PoseGraph, Rt_i, Rt_j);

    // if (printCount++ < iters) {
    //   cout << endl;
    // }


    /*
     * Build Bundle Adjustment problem
     */

    // We want to add some number of points to constrain rotations
     int num_points = min(NUM_BA_POINTS, (int)cameraRt_ij_points_count[i]);

    // Overweight loop closures
     w = cameraRt_ij_indices[2 * i + 1] - cameraRt_ij_indices[2 * i] == 1 ? 1 : 1;
     double weight_BundleAdjustment[3] = {w, w, w};

     for (int j = 0; j < num_points; j++) {
      int index = (points_count + j) * 3;
      int t_index = (cameraRt_ij_points_total + j) * 3;

      // Fill up observed points
      points_observed_i[index] = cameraRt_ij_points_observed_i[t_index];
      points_observed_i[index + 1] = cameraRt_ij_points_observed_i[t_index + 1];
      points_observed_i[index + 2] = cameraRt_ij_points_observed_i[t_index + 2];

      points_observed_j[index] = cameraRt_ij_points_observed_j[t_index];
      points_observed_j[index + 1] = cameraRt_ij_points_observed_j[t_index + 1];
      points_observed_j[index + 2] = cameraRt_ij_points_observed_j[t_index + 2];

      // Fill up predicted points
      // points_predicted[index] = cameraRt_ij_points_predicted[t_index];
      // points_predicted[index + 1] = cameraRt_ij_points_predicted[t_index + 1];
      // points_predicted[index + 2] = cameraRt_ij_points_predicted[t_index + 2];

      double *point_observed_i = points_observed_i + index;
      double *point_observed_j = points_observed_j + index;
      double *point_predicted = points_predicted + index;

      AngleAxisRotateAndTranslatePoint(Rt_i, point_observed_i, point_predicted);

      // TODO: replace these with pointer to Cerberus
      BundleAdjustmentResidual::AddResidualBlock(c2, point_observed_i, weight_BundleAdjustment, Rt_i, point_predicted);
      BundleAdjustmentResidual::AddResidualBlock(c2, point_observed_j, weight_BundleAdjustment, Rt_j, point_predicted);
    }

    points_count += num_points;
    cameraRt_ij_points_total += cameraRt_ij_points_count[i];
  }

  //----------------------------------------------------------------

  // cout << "Starting full pose graph solver" << endl;
  c1->solve();

  cout << "Starting pose graph BA solver" << endl;
  c2->solve();

  // obtain camera matrix from parameters
  for(int cameraID = 0; cameraID < nCam; cameraID++){
    double* cameraPtr = cameraParameter + 6 * cameraID;
    double* cameraMat = cameraRtC2W + 12 * cameraID;
    if (!(isnan(*cameraPtr))){
      ceres::AngleAxisToRotationMatrix<double>(cameraPtr, cameraMat);
      cameraMat[9]  = cameraPtr[3];
      cameraMat[10] = cameraPtr[4];
      cameraMat[11] = cameraPtr[5];
    }
  }

  // write back result files
  FILE* fpout = fopen(argv[4],"wb");
  fwrite((void*)(cameraRtC2W), sizeof(double), 12*nCam, fpout);

  fclose (fpout);

  // clean up
  delete [] cameraRtC2W;
  delete [] cameraRt_ij;
  delete [] cameraRt_ij_indices;
  delete [] cameraRt_ij_points_observed_i;
  delete [] cameraRt_ij_points_observed_j;
  delete [] cameraRt_ij_points_predicted;
  delete [] cameraRt_ij_points_count;
  delete [] cameraParameter;
  delete [] cameraParameter_ij;
  delete [] points_observed_i;
  delete [] points_observed_j;
  delete [] points_predicted;

  return 0;
}

 // double t_o[3];
 //    double t_p[3];

 //    ceres::AngleAxisRotatePoint(Rt_i, Rt_ij + 3, t_o);
 //    t_p[0] = Rt_j[3] - Rt_i[3];
 //    t_p[1] = Rt_j[4] - Rt_i[4];
 //    t_p[2] = Rt_j[5] - Rt_i[5];

 //    double IR_i[3] = {-Rt_i[0], -Rt_i[1], -Rt_i[2]};
 //    // Z and y vectors to constrain rotation
 //    double z[3] = {0, 0, 1}; double z_o[3]; double z_p[3];
 //    double y[3] = {0, 1, 0}; double y_o[3]; double y_p[3];

 //    ceres::AngleAxisRotatePoint(Rt_ij, z, z_o);
 //    ceres::AngleAxisRotatePoint(Rt_j, z, z_p);
 //    ceres::AngleAxisRotatePoint(IR_i, z_p, z_p);

 //    ceres::AngleAxisRotatePoint(Rt_ij, y, y_o);
 //    ceres::AngleAxisRotatePoint(Rt_j, y, y_p);
 //    ceres::AngleAxisRotatePoint(IR_i, y_p, y_p);

 //    sum += pow(weight * (z_p[0] - z_o[0]), 2.0);
 //    sum += pow(weight * (z_p[1] - z_o[1]), 2.0);
 //    sum += pow(weight * (z_p[2] - z_o[2]), 2.0);
 //    sum += pow(weight * (y_p[0] - y_o[0]), 2.0);
 //    sum += pow(weight * (y_p[1] - y_o[1]), 2.0);
 //    sum += pow(weight * (y_p[2] - y_o[2]), 2.0);

 //    // Create observed translation vector
 //    double t_ij_observed[3] = {Rt_ij[3], Rt_ij[4], Rt_ij[5]};
 //    ceres::AngleAxisRotatePoint(Rt_i, t_ij_observed, t_ij_observed);

 //    // Create predicted translation vector. Note this should be t_j -
 //    // t_i, but because Rt_j and Rt_i are world-to-camera and Rt_ij is
 //    // camera-to-world
 //    double t_ij_predicted[3] = {Rt_j[3] - Rt_i[3], Rt_j[4] - Rt_i[4], Rt_j[5] - Rt_i[5]};

 //    sum += pow(weight * 2 * (t_ij_predicted[0] - t_ij_observed[0]), 2.0);
 //    sum += pow(weight * 2 * (t_ij_predicted[1] - t_ij_observed[1]), 2.0);
 //    sum += pow(weight * 2 * (t_ij_predicted[2] - t_ij_observed[2]), 2.0);

 //    double p_o[3];
 //    double p_p[3];

 //    // Compute observed view direction vector from world to camera by
 //    // relative rotation between i and j
 //    ceres::AngleAxisRotatePoint(Rt_ij, z, p_o);

 //    // Rotate from according to j's rotation
 //    ceres::AngleAxisRotatePoint(Rt_j, z, p_p);

 //    // Rotate back from according to i's rotation
 //    ceres::AngleAxisRotatePoint(IR_i, p_p, p_p);

 //    sum += t_p[0] - t_o[0] + t_p[1] - t_o[1] + t_p[2] - t_o[2] + p_p[0] - p_o[0] + p_p[1] - p_o[1] + p_p[2] - p_o[2];

 //    cout << sum << " " << weight << endl;
 //    cout << t_o[0] << " " << t_o[1] << " " << t_o[2] << " " << p_o[0] << " " << p_o[1] << " " << p_o[2] << endl;
 //    cout << t_p[0] << " " << t_p[1] << " " << t_p[2] << " " << p_p[0] << " " << p_p[1] << " " << p_p[2] << endl;
 //    cout << endl;
