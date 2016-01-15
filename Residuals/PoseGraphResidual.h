#ifndef POSEGRAPHRESIDUAL_H
#define POSEGRAPHRESIDUAL_H

#include "ceres/ceres.h"
#include "Cerberus.h"

class PoseGraphResidual {
public:
  PoseGraphResidual(double *m_Rt_ij, double *m_weight): m_Rt_ij(m_Rt_ij), m_weight(m_weight) {}
  template <typename T>
  bool operator()(const T* const Rt_i,
    const T* const Rt_j,
    T* residuals) const {

    // TODO: compare results to optimizing angle axis directly rather
    // than converting in residual block

    // Get T-ified weight
    T weight[9];
    for (int i = 0; i < 9; i++) {
      weight[i] = T(m_weight[i]);
    }

    // Create a T-ified Rt_ij for computation. Note that this rotation
    // is camera-to-world
    T Rt_ij[6];
    for (int k = 0; k < 6; k++) {
      Rt_ij[k] = T(m_Rt_ij[k]);
    }

    // Rotation
    T IR_i[3] = {-Rt_i[0], -Rt_i[1], -Rt_i[2]};
    T z[3] = {T(0), T(0), T(1)}; T z_o[3]; T z_p[3];
    T y[3] = {T(0), T(1), T(0)}; T y_o[3]; T y_p[3];

    ceres::AngleAxisRotatePoint(Rt_ij, z, z_o);
    ceres::AngleAxisRotatePoint(Rt_j, z, z_p);
    ceres::AngleAxisRotatePoint(IR_i, z_p, z_p);

    ceres::AngleAxisRotatePoint(Rt_ij, y, y_o);
    ceres::AngleAxisRotatePoint(Rt_j, y, y_p);
    ceres::AngleAxisRotatePoint(IR_i, y_p, y_p);

    // Translation
    T t_o[3] = {Rt_ij[3], Rt_ij[4], Rt_ij[5]};
    T t_p[3] = {Rt_j[3] - Rt_i[3], Rt_j[4] - Rt_i[4], Rt_j[5] - Rt_i[5]};
    ceres::AngleAxisRotatePoint(Rt_i, t_o, t_o);

    // Set residuals
    residuals[0] = weight[0] * (z_p[0] - z_o[0]);
    residuals[1] = weight[1] * (z_p[1] - z_o[1]);
    residuals[2] = weight[2] * (z_p[2] - z_o[2]);
    residuals[3] = weight[3] * (y_p[0] - y_o[0]);
    residuals[4] = weight[4] * (y_p[1] - y_o[1]);
    residuals[5] = weight[5] * (y_p[2] - y_o[2]);
    residuals[6] = weight[6] * (t_p[0] - t_o[0]);
    residuals[7] = weight[7] * (t_p[1] - t_o[1]);
    residuals[8] = weight[8] * (t_p[2] - t_o[2]);

    return true;
  }

  // Convenience functions that are horribly designed.
  static CostFunction *GenerateCostFunction(double *Rt_ij, double *weight) {
    return new AutoDiffCostFunction<PoseGraphResidual, 9, 6, 6>(new PoseGraphResidual(Rt_ij, weight));
  }

  static void AddResidualBlock(Cerberus *solver, double *Rt_ij, double *weight, double *Rt_i, double *Rt_j) {
    CostFunction *cost = PoseGraphResidual::GenerateCostFunction(Rt_ij, weight);
    solver->problem.AddResidualBlock(cost, solver->loss, Rt_i, Rt_j);
  }

  double *m_Rt_ij;
  double *m_weight;
};

#endif