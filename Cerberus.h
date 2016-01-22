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
template <typename T, typename U> void TransformationMatrixToAngleAxisAndTranslation(T *Rt, U *result, bool column_major = false) {
  U tmp[12];
  for (int i = 0; i < 12; i++) {
    tmp[i] = U(Rt[i]);
  }

  if (!column_major) {
    RotationMatrixToAngleAxis<U>(MatrixAdapter<const U, 4, 1>(tmp), result);
    result[3] = tmp[3];
    result[4] = tmp[7];
    result[5] = tmp[11];
  } else {
    RotationMatrixToAngleAxis<U>(tmp, result);
    result[3] = tmp[9];
    result[4] = tmp[10];
    result[5] = tmp[11];
  }
}

// At is a 6-dimensional vector with first three terms representing axis,
// magnitude angle, and the last three terms represent translation. Rt is a
// 4x3 transformation matrix of the form [R|t]. Two typenames
// because Ceres requires doubles.
template <typename T, typename U> void AngleAxisAndTranslationToTransformationMatrix(const U *At, T *result, bool column_major = false) {
  U tmp[12];

  if (!column_major) {
    AngleAxisToRotationMatrix<U>(At, MatrixAdapter<U, 4, 1>(tmp));
    tmp[3] = At[3];
    tmp[7] = At[4];
    tmp[11] = At[5];
  } else {
    AngleAxisToRotationMatrix<U>(At, tmp);
    tmp[9] = At[3];
    tmp[10] = At[4];
    tmp[11] = At[5];
  }

  for (int i = 0; i < 12; i++) {
    result[i] = T(tmp[i]);
  }
}

// Rotate a point about axis and translate
template <typename T, typename U> void AngleAxisRotateAndTranslatePoint(const U *Rt, const T *from, U *to, bool inverse = false) {
  U tmp[3] = { U(from[0]), U(from[1]), U(from[2]) };

  if (!inverse) {
    AngleAxisRotatePoint(Rt, tmp, to);

    to[0] += Rt[3];
    to[1] += Rt[4];
    to[2] += Rt[5];
  } else {
    to[0] = tmp[0] - Rt[3];
    to[1] = tmp[1] - Rt[4];
    to[2] = tmp[2] - Rt[5];

    U IRt[3] = { -Rt[0], -Rt[1], -Rt[2] };
    AngleAxisRotatePoint(IRt, to, to);
  }
}

template <typename T> void InverseAngleAxis(T *Rt) {
  Rt[0] = -Rt[0];
  Rt[1] = -Rt[1];
  Rt[2] = -Rt[2];
  Rt[3] = -Rt[3];
  Rt[4] = -Rt[4];
  Rt[5] = -Rt[5];

  AngleAxisRotatePoint(Rt, Rt + 3, Rt + 3);
}

#endif