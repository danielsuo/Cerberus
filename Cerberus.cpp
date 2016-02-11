#include "Cerberus.h"

Cerberus::Cerberus() {
  // Set default solver options
  options.max_num_iterations = 2000;
  options.minimizer_progress_to_stdout = true;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.trust_region_strategy_type = LEVENBERG_MARQUARDT;
  options.num_threads = 8;

  // Set default loss function
  loss = new TrivialLoss();
}

void Cerberus::solve() {
  Solve(options, &problem, &summary);
  cerr << summary.BriefReport() << endl;
}