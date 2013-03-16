/*
 * Solver.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: chaako
 */

#include "Solver.h"

Solver::Solver() {
	// TODO Auto-generated constructor stub

}

Solver::~Solver() {
	// TODO Auto-generated destructor stub
}

void Solver::compute(Eigen::SparseMatrix<double> A) {
	solver.compute(A);
}

Eigen::VectorXd Solver::solve(Eigen::VectorXd b) {
	return solver.solve(b);
}

int Solver::info() {
	return solver.info();
}
