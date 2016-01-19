#ifndef CERBERUS_H
#define CERBERUS_H

#include <iostream>
#include "ceres/ceres.h"
#include "ceres/rotation.h"

using namespace std;
using namespace ceres;

class Cerberus {
public:
  Cerberus();

  Solver::Options options;
  LossFunction *loss;
  Problem problem;
  Solver::Summary summary;

  void solve();
};

// Rt is a 4x3 transformation matrix of the form [R|t]. result is a
// 6-dimensional vector with first three terms representing axis, magnitude
// angle, and the last three terms represent translation. Two typenames
// because Ceres requires doubles.
template <typename T, typename U> void TransformationMatrixToAngleAxisAndTranslation(T *Rt, U *result) {
  U tmp[12];
  for (int i = 0; i < 12; i++) {
    tmp[i] = U(Rt[i]);
  }

  TransformationMatrixToAngleAxisAndTranslation(tmp, result);
}

template <typename T> void TransformationMatrixToAngleAxisAndTranslation(T *Rt, T *result) {
  RotationMatrixToAngleAxis<T>(MatrixAdapter<const T, 4, 1>(Rt), result);
  result[3] = Rt[3];
  result[4] = Rt[7];
  result[5] = Rt[11];
}

// At is a 6-dimensional vector with first three terms representing axis,
// magnitude angle, and the last three terms represent translation. Rt is a
// 4x3 transformation matrix of the form [R|t]. Two typenames because Ceres
// requires doubles.
template <typename T, typename U> void AngleAxisAndTranslationToTransformationMatrix(const U *At, T *result) {
  U tmp[12];

  AngleAxisAndTranslationToTransformationMatrix(At, tmp);

  for (int i = 0; i < 12; i++) {
    result[i] = T(tmp[i]);
  }
}

template <typename T> void AngleAxisAndTranslationToTransformationMatrix(const T *At, T *result) {
  AngleAxisToRotationMatrix<T>(At, MatrixAdapter<T, 4, 1>(result));
  result[3] = At[3];
  result[7] = At[4];
  result[11] = At[5];
}

#endif