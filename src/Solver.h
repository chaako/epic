/*
 * Solver.h
 *
 *  Created on: Mar 16, 2013
 *      Author: chaako
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "typesAndDefinitions.h"

class Solver {
public:
	Solver();
	virtual ~Solver();
	void compute(Eigen::SparseMatrix<double> A);
	Eigen::VectorXd solve(Eigen::VectorXd b);
	int info();

#ifndef MESHER
//	Eigen::SuperLU<Eigen::SparseMatrix<double> > preconditioner;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double>,Eigen::SuperILU<Eigen::SparseMatrix<double> > > solver;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double>,Eigen::SuperLU<Eigen::SparseMatrix<double> > > solver;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
	Eigen::SuperLU<Eigen::SparseMatrix<double> > solver;
#else
	Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
#endif

};

#endif /* SOLVER_H_ */
