/*
 * Solver.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: chaako
 */

#include "Solver.h"

#ifdef HAVE_MPI
#ifndef MESHER

Solver::Solver() {
	currentSize = 1;
	currentSolution = new double[currentSize];
	currentError = new double[currentSize];
	currentStatus = Eigen::Success;

    int nprow=sqrt(MPI::COMM_WORLD.Get_size());
    int npcol=nprow;
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
    partOfGrid = (grid.iam < nprow*npcol);

    nrhs = 1;
    needToFree = false;
}

Solver::~Solver() {
	delete[] currentSolution;
	delete[] currentError;
	if (needToFree) {
		delete[] values;
		delete[] innerInd;
		delete[] outerInd;
	    ScalePermstructFree(&ScalePermstruct);
	    Destroy_LU(currentSize, &grid, &LUstruct);
	    LUstructFree(&LUstruct);
	}
}

extern "C" {
	extern void pdgssvx(superlu_options_t *options, SuperMatrix *Aslu,
	ScalePermstruct_t *ScalePermstruct,
	double B[], int ldb, int nrhs, gridinfo_t *grid,
	LUstruct_t *LUstruct, SOLVEstruct_t *SOLVEstruct, double *berr,
	SuperLUStat_t *stat, int *info);
}

void Solver::compute(Eigen::SparseMatrix<double> A) {

	currentSize = A.rows();
    // TODO: better way to handle new/delete?
	delete[] currentSolution;
	currentSolution = new double[currentSize];
	delete[] currentError;
	currentError = new double[currentSize];

	if ( partOfGrid ){

		set_default_options_dist(&options);
		// TODO: explore other reorderings for sparsity-preservation (e.g. METIS),
		//       since NATURAL gives almost an order of magnitude more non-zeros
//		options.ColPerm = NATURAL;
		options.ColPerm = MMD_AT_PLUS_A;

		needToFree = true;

//		values = A.valuePtr();
//		innerInd = A.innerIndexPtr();
//		outerInd = A.outerIndexPtr();

		nnz = A.nonZeros();
//		// TODO: debugging
//		cout << "vals: " << values[nnz-1] << " " << values[nnz-2] << endl;
//		cout << "inner: " << innerInd[nnz-1] << " " << innerInd[nnz-2] << endl;
//		cout << "outer: " << outerInd[nnz-1] << " " << outerInd[nnz-2] << endl;

		values = new double[nnz];
		innerInd = new int[nnz];
		outerInd = new int[nnz];
		for (int i=0; i<nnz; i++) {
			values[i] = *((double*)(A.valuePtr()) + i);
			innerInd[i] = *((int*)(A.innerIndexPtr()) + i);
			outerInd[i] = *((int*)(A.outerIndexPtr()) + i);
		}

//		// TODO: debugging
//		cout << "vals: " << values[nnz-1] << " " << values[nnz-2] << endl;
//		cout << "inner: " << innerInd[nnz-1] << " " << innerInd[nnz-2] << endl;
//		cout << "outer: " << outerInd[nnz-1] << " " << outerInd[nnz-2] << endl;

		// Map sparse Eigen matrix to SuperLU matrix
		Aslu.Stype = SLU_NC;
		Aslu.nrow = A.rows();
		Aslu.ncol = A.cols();
		Aslu.Mtype = SLU_GE;
//		storage.nnz = A.nonZeros();
//		storage.values = A.valuePtr();
//		storage.innerInd = A.innerIndexPtr();
//		storage.outerInd = A.outerIndexPtr();
		storage.nnz = nnz;
		storage.values = values;
		storage.innerInd = innerInd;
		storage.outerInd = outerInd;
		Aslu.Store = &storage;
		Aslu.Dtype = SLU_D;

//    // TODO: debugging
//    cout << "storage: " << storage.nnz << " " << storage.lda << " " <<
//    		storage.values << " " << storage.innerInd << " " <<
//    		storage.outerInd << endl;
//    cout << "storage: " << A.nonZeros() << " " << A.rows() << " " <<
//    		A.valuePtr() << " " << A.innerIndexPtr() << " " <<
//    		A.outerIndexPtr() << endl;


//    int m,n,nnz;
//    SuperMatrix Atest;
//    double *a, *rhs;
//    double s, u, p, e, r, l;
//    int *asub, *xa;
//    m = n = 5;
//    nnz = 12;
//    a = new double[nnz];
//    asub = new int[nnz];
//    xa = new int[n+1];
//    s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
//    a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
//    a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
//    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
//    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
//    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
//    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;
//    dCreate_CompCol_Matrix_dist(&Atest, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
//    rhs = new double[m * nrhs];
//    for (int i = 0; i < m; ++i) rhs[i] = 1.0;
//    ScalePermstructInit(m, n, &ScalePermstruct);
//    LUstructInit(m, n, &LUstruct);
//    PStatInit(&stats);
//	pdgssvx_ABglobal(&options, &Atest, &ScalePermstruct, rhs, n, 1, &grid,
//			&LUstruct, berr, &stats, &info);
//
//	if (info!=0)
//		throw;
//    PStatFree(&stats);
//    ScalePermstructFree(&ScalePermstruct);
//    Destroy_LU(A.cols(), &grid, &LUstruct);
//    LUstructFree(&LUstruct);
//    delete[] a;
//    delete[] asub;
//    delete[] xa;
//    delete[] rhs;

		ScalePermstructInit(A.rows(), A.cols(), &ScalePermstruct);
		LUstructInit(A.rows(), A.cols(), &LUstruct);


		PStatInit(&stats);
		infoSlu = 0;
		pdgssvx_ABglobal(&options, &Aslu, &ScalePermstruct, currentSolution, currentSize, nrhs, &grid,
				&LUstruct, currentError, &stats, &infoSlu);
	    PStatFree(&stats);

		if (infoSlu!=0) {
			factored = false;
			throw;
		} else {
			factored = true;
		}
		if ( options.SolveInitialized ) {
			dSolveFinalize(&options, &SOLVEstruct);
		}

		this->currentStatus = (infoSlu == 0) ? Eigen::Success : Eigen::NumericalIssue;
    }
}

Eigen::VectorXd Solver::solve(Eigen::VectorXd b) {
	assert(currentSize==b.rows());
	for (int i=0; i<currentSize; i++) {
		currentSolution[i] = b[i];
	}
	if (!factored)
		throw;
	if (partOfGrid) {
//		set_default_options_dist(&options);
		// TODO: explore other reorderings for sparsity-preservation (e.g. METIS),
		//       since NATURAL gives almost an order of magnitude more non-zeros
//		options.ColPerm = NATURAL;
//		options.ColPerm = MMD_AT_PLUS_A;
		options.Trans = NOTRANS;
		options.Fact = FACTORED;
		options.IterRefine = NOREFINE;
		PStatInit(&stats);
		infoSlu = 0;
		pdgssvx_ABglobal(&options, &Aslu, &ScalePermstruct, currentSolution, currentSize, nrhs, &grid,
				&LUstruct, currentError, &stats, &infoSlu);
	    PStatFree(&stats);
		this->currentStatus = (infoSlu == 0) ? Eigen::Success : Eigen::NumericalIssue;
	}

	// TODO: Safe to assume master is part of SuperLU processor grid?
	int mpiId = 0;
	mpiId = MPI::COMM_WORLD.Get_rank();
	MPI::COMM_WORLD.Bcast(currentSolution, currentSize, MPI::DOUBLE, 0);
	Eigen::VectorXd x(currentSize);
	for (int i=0; i<currentSize; i++) {
		x[i] = currentSolution[i];
//		// TODO: debugging
//		cout << " " << x[i];
	}

	return x;
}

int Solver::info() {
	return currentStatus;
}
#endif
#else

Solver::Solver() {
}

Solver::~Solver() {
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

#endif
