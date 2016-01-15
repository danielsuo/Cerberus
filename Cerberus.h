#ifndef CERBERUS_H
#define CERBERUS_H

#include <iostream>
#include "ceres/ceres.h"

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

#endif