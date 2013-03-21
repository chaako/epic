/*
 * Solver.h
 *
 *  Created on: Mar 16, 2013
 *      Author: chaako
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "typesAndDefinitions.h"

#ifdef HAVE_MPI
#ifndef MESHER
#include <superlu_ddefs.h>
#endif
#endif

class Solver {
public:
	Solver();
	virtual ~Solver();
	void compute(Eigen::SparseMatrix<double> A);
	Eigen::VectorXd solve(Eigen::VectorXd b);
	int info();

#ifndef MESHER
#ifndef HAVE_MPI
//	Eigen::SuperLU<Eigen::SparseMatrix<double> > preconditioner;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double>,Eigen::SuperILU<Eigen::SparseMatrix<double> > > solver;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double>,Eigen::SuperLU<Eigen::SparseMatrix<double> > > solver;
//	Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
	Eigen::SuperLU<Eigen::SparseMatrix<double> > solver;
#endif
#else
	Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
#endif

	int currentStatus;
	double *currentSolution;
	double *currentError;
	int currentSize;
	struct {
		union {int nnz;int lda;};
	    void *values;
		int *innerInd;
		int *outerInd;
	} storage;
	int nnz;
    double *values;
	int *innerInd;
	int *outerInd;

#ifdef HAVE_MPI
#ifndef MESHER

	gridinfo_t grid;
	bool partOfGrid;
	bool factored;
	bool needToFree;

	superlu_options_t options;
	SuperMatrix Aslu;
	ScalePermstruct_t ScalePermstruct;
	int nrhs;
	LUstruct_t LUstruct;
	SOLVEstruct_t SOLVEstruct;
    double *berr;
    SuperLUStat_t stats;
    int infoSlu;

#endif
#endif

};


#endif /* SOLVER_H_ */
