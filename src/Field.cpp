/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Field.h"

ElectricField::ElectricField(Mesh *inputMesh_ptr, string inputName)
		: Field<vect3d>(inputMesh_ptr, inputName, iBase_VERTEX) {
}

ElectricField::ElectricField(Mesh *inputMesh_ptr, string inputName,
		CodeField vertexType)
		: Field<vect3d>(inputMesh_ptr, inputName, iBase_VERTEX) {
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank() == 0)
#endif
		cout << "Computing the Poisson solve LU decomposition..." << endl << endl;
	int nVerts = entities.size();
	int m = nVerts + nVerts*NDIM;
	vector<Eigen::Triplet<double> > coefficients;
	// TODO: don't hard-code estimated number of adjacent tets
	coefficients.reserve(20*m);
	vector<bool> potentialBoundaryVertexSet(nVerts);
	vector<bool> eFieldBoundaryVertexSet(nVerts);
	for (int j=0; j<nVerts; j++) {
		potentialBoundaryVertexSet[j] = false;
		eFieldBoundaryVertexSet[j] = false;
	}
	for (int i=0; i<mesh_ptr->entitiesVectors[iBase_REGION].size(); i++) {
		vector<vect3d> vVs = mesh_ptr->getVertexVectors(i,iBase_REGION);
		double volume = mesh_ptr->getTetVolume(vVs);
		vect3d centroid(0.,0.,0.);
		centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
		vector<int> adjacentVertices =
				mesh_ptr->adjacentEntitiesVectors[iBase_REGION][i][iBase_VERTEX];
		Eigen::Matrix<double,NDIM,NDIM+1> basisDerivatives;
		mesh_ptr->evaluateLinearBasisFunctionDerivatives(centroid, i,
				&basisDerivatives);
		// The order of the basis derivatives corresponds to that of adjacentVertices
		for (int j=0; j<adjacentVertices.size(); j++) {
			int jj = adjacentVertices[j];
			vect3d vertexPosition = mesh_ptr->getCoordinates(entities[jj]);
			double potential = 0.;
			// TODO: don't hard-code boundary codes
			if (vertexType[jj]==4) {
				// Set potential at object surface
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, jj, 1.));
					potentialBoundaryVertexSet[jj] = true;
				}
			} else if (vertexType[jj]==5) {
				// Set potential at object surface
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, jj, 1.));
					potentialBoundaryVertexSet[jj] = true;
				}
			}
			double debyeLength = 1.e0;
			for (int k=0; k<adjacentVertices.size(); k++) {
				int kk = adjacentVertices[k];
				double basisBasisCoefficient = volume/20;
				if (jj==kk)
					basisBasisCoefficient *= 2.;
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, kk,
							basisBasisCoefficient/debyeLength/debyeLength*
							exp(potential)));
				}
				for (int l=0; l<NDIM; l++) {
					// Coefficients multiplying eField
					if (!eFieldBoundaryVertexSet[jj]) {
						// TODO: Make functions for indexing
						coefficients.push_back(Eigen::Triplet<double>(
								nVerts+jj*NDIM+l, nVerts+kk*NDIM+l,
								basisBasisCoefficient));
					}
					double basisDerivativeBasisCoefficient =
							-basisDerivatives(l,k)*volume/4;
					if (!potentialBoundaryVertexSet[jj]) {
						coefficients.push_back(Eigen::Triplet<double>(
								jj, nVerts+kk*NDIM+l,
								basisDerivativeBasisCoefficient));
					}
					// Coefficients multiplying potential
					if (!eFieldBoundaryVertexSet[jj]) {
						coefficients.push_back(Eigen::Triplet<double>(
								nVerts+jj*NDIM+l, kk,
								basisDerivativeBasisCoefficient));
					}
				}
			}
		}
	}

	Eigen::SparseMatrix<double> fixedPartOfA(m,m);
	fixedPartOfA.setFromTriplets(coefficients.begin(), coefficients.end());
//	// TODO: Don't hard-code tolerance
//	solver.setTolerance(1.e-2);
//	solver.setMaxIterations(10);
	solver.compute(fixedPartOfA);
//	solver.analyzePattern(fixedPartOfA);
	if(solver.info()!=Eigen::Success) {
		// decomposition failed
		throw;
	}
}

void ElectricField::calcField(PotentialField potentialField) {
	int m = entities.size()*NDIM;
	// Find components of E by finite element method (solving Ax=b)
	vector<Eigen::Triplet<double> > coefficients;
	Eigen::VectorXd b(m);
	for (int i=0; i<m; i++)
		b[i] = 0.;
	// TODO: don't hard-code estimated number of adjacent tets
	coefficients.reserve(20*m);
	for (int i=0; i<mesh_ptr->entitiesVectors[iBase_REGION].size(); i++) {
		vector<vect3d> vVs = mesh_ptr->getVertexVectors(i,iBase_REGION);
		double volume = mesh_ptr->getTetVolume(vVs);
		vect3d centroid(0.,0.,0.);
		centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
		double potential=0.;
		vect3d potentialDerivative(0.,0.,0.);
		int elementIndex=i;
		potentialField.evalFieldAndDeriv(&potential, &potentialDerivative,
				centroid, &elementIndex, 1);
		assert(elementIndex==i);
		vector<int> adjacentVertices =
				mesh_ptr->adjacentEntitiesVectors[iBase_REGION][i][iBase_VERTEX];
		for (int l=0; l<NDIM; l++) {
			for (int j=0; j<adjacentVertices.size(); j++) {
				int jj = adjacentVertices[j];
				b[jj*NDIM+l] += potentialDerivative[l]*volume/4.;
				for (int k=0; k<adjacentVertices.size(); k++) {
					int kk = adjacentVertices[k];
					double coefficient = volume/20;
					if (jj==kk)
						coefficient *= 2.;
					// TODO: Make functions for indexing
					coefficients.push_back(Eigen::Triplet<double>(
							jj*NDIM+l, kk*NDIM+l, coefficient));
				}
			}
		}
	}

	Eigen::SparseMatrix<double> A(m,m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	// A is SPD, so can use Cholesky factorization (or Biconjugate Gradient)
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > choleskyFactorization(A);
	if(choleskyFactorization.info()!=Eigen::Success) {
		// decomposition failed
		throw;
	}
	Eigen::VectorXd x = choleskyFactorization.solve(b);
	if(choleskyFactorization.info()!=Eigen::Success) {
		// solving failed
		throw;
	}
//	for (int i=0; i<m; i++) {
//		cout << "A[" << i << "," << i << "] = " << A.coeff(i,i) << " ";
//		cout << "x[" << i << "] = " << x[i] << " ";
//		cout << "b[" << i << "] = " << b[i] << endl;
//	}

	for (int i=0; i<entities.size(); i++) {
		vect3d eField(0.,0.,0.);
		// TODO: Make functions for indexing
		// TODO: Check sign of eField
		for (int j=0; j<NDIM; j++)
			eField[j] = -x[i*NDIM+j];
		this->setField(entities[i], eField);
	}

}

void ElectricField::calcField(PotentialField *potentialField_ptr, CodeField vertexType,
		DensityField ionDensity) {
	int nVerts = entities.size();
	int m = nVerts + nVerts*NDIM;
	double potentialChange=1.;
	int numberOfSolves=1;
	// Iterate linearized FEM solve to get Boltzmann electron response
	while (numberOfSolves<20 && potentialChange>2.e-5) {
	potentialChange=0.;
	// Find potential and components of E by finite element method (solving Ax=b)
	vector<Eigen::Triplet<double> > coefficients;
	Eigen::VectorXd b(m);
	for (int i=0; i<m; i++)
		b[i] = 0.;
	// TODO: don't hard-code estimated number of adjacent tets
	coefficients.reserve(20*m);
	vector<bool> potentialBoundaryVertexSet(nVerts);
	vector<bool> eFieldBoundaryVertexSet(nVerts);
	for (int j=0; j<nVerts; j++) {
		potentialBoundaryVertexSet[j] = false;
		eFieldBoundaryVertexSet[j] = false;
	}
	for (int i=0; i<mesh_ptr->entitiesVectors[iBase_REGION].size(); i++) {
		vector<vect3d> vVs = mesh_ptr->getVertexVectors(i,iBase_REGION);
		double volume = mesh_ptr->getTetVolume(vVs);
		vect3d centroid(0.,0.,0.);
		centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
		vector<int> adjacentVertices =
				mesh_ptr->adjacentEntitiesVectors[iBase_REGION][i][iBase_VERTEX];
		Eigen::Matrix<double,NDIM,NDIM+1> basisDerivatives;
		mesh_ptr->evaluateLinearBasisFunctionDerivatives(centroid, i,
				&basisDerivatives);
		// The order of the basis derivatives corresponds to that of adjacentVertices
		for (int j=0; j<adjacentVertices.size(); j++) {
			int jj = adjacentVertices[j];
			vect3d vertexPosition = mesh_ptr->getCoordinates(entities[jj]);
			double potential = potentialField_ptr->operator[](jj);
			// TODO: don't hard-code boundary codes
			if (vertexType[jj]==4) {
				// Set potential at object surface
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, jj, 1.));
					// TODO: hard-coding zero potential of ExB-field as origin here
					double shieldedPotential=vertexPosition.dot(E);
					// TODO: Don't hard-code object potential
					b[jj] = -4. - shieldedPotential;
//					b[jj] = potential;
					potentialBoundaryVertexSet[jj] = true;
				}
			} else if (vertexType[jj]==5) {
				// Set potential at object surface
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, jj, 1.));
					// TODO: Don't hard-code uperturbed potential
					b[jj] = -0.;
//					b[jj] = potential;
					potentialBoundaryVertexSet[jj] = true;
				}
//			} else if (vertexType[jj]==5) {
//				// Set electric field to zero
//				if (!eFieldBoundaryVertexSet[jj]) {
//					for (int l=0; l<NDIM; l++) {
//						coefficients.push_back(Eigen::Triplet<double>(
//								nVerts+jj*NDIM+l, nVerts+jj*NDIM+l, 1.));
//						b[nVerts+jj*NDIM+l] = 0.;
//					}
//					eFieldBoundaryVertexSet[jj] = true;
//				}
			}
			double debyeLength = 1.e0;
			for (int k=0; k<adjacentVertices.size(); k++) {
				int kk = adjacentVertices[k];
				double basisBasisCoefficient = volume/20;
				if (jj==kk)
					basisBasisCoefficient *= 2.;
				if (!potentialBoundaryVertexSet[jj]) {
					coefficients.push_back(Eigen::Triplet<double>(jj, kk,
							basisBasisCoefficient/debyeLength/debyeLength*
							exp(potential)));
				}
				for (int l=0; l<NDIM; l++) {
					// Coefficients multiplying eField
					if (!eFieldBoundaryVertexSet[jj]) {
						// TODO: Make functions for indexing
						coefficients.push_back(Eigen::Triplet<double>(
								nVerts+jj*NDIM+l, nVerts+kk*NDIM+l,
								basisBasisCoefficient));
					}
					double basisDerivativeBasisCoefficient =
							-basisDerivatives(l,k)*volume/4;
					if (!potentialBoundaryVertexSet[jj]) {
						coefficients.push_back(Eigen::Triplet<double>(
								jj, nVerts+kk*NDIM+l,
								basisDerivativeBasisCoefficient));
					}
					// Coefficients multiplying potential
					if (!eFieldBoundaryVertexSet[jj]) {
						coefficients.push_back(Eigen::Triplet<double>(
								nVerts+jj*NDIM+l, kk,
								basisDerivativeBasisCoefficient));
					}
				}
			}
			// TODO: add charge source term
			if (!potentialBoundaryVertexSet[jj]) {
				b[jj] += 1/debyeLength/debyeLength*
//						(ionDensity[jj] - exp(potential)*(1.-potential))*
						// TODO: check that this works for potential~1
						(ionDensity[jj] - (exp(potential)-potential))*
						volume/4.;
			}
		}
	}

	Eigen::SparseMatrix<double> A(m,m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

//	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(A);
#ifndef MESHER
//	Eigen::SuperLU<Eigen::SparseMatrix<double> > solver;
//	solver.compute(A);
//	if(solver.info()!=Eigen::Success) {
//		// decomposition failed
//		throw;
//	}
//	for (int i=0; i<m; i++) {
//		cout << "A[" << i << "," << i << "] = " << A.coeff(i,i) << " ";
////		cout << "A[" << 3 << "," << i << "] = " << A.coeff(3,i) << " ";
////		cout << "x[" << i << "] = " << x[i] << " ";
//		cout << "b[" << i << "] = " << b[i] << endl;
//	}
//	cout << A << endl << endl;
//	cout << b.transpose() << endl << endl;
//	solver.factorize(A);
	Eigen::VectorXd x = solver.solve(b);
//	Eigen::VectorXd x = solver.solveWithGuess(b,currentSolution);
	if(solver.info()!=Eigen::Success) {
		// solving failed
		cout << "Poisson solve failed with error: " << solver.info() << endl;
		throw;
	}
//	cout << x.transpose() << endl << endl;

	for (int i=0; i<nVerts; i++) {
		double dPot = fabs(x[i]-potentialField_ptr->operator[](i));
		potentialChange = (potentialChange>dPot) ?
				potentialChange	: dPot;
	}
	for (int i=0; i<nVerts; i++) {
		vect3d eField(0.,0.,0.);
		vect3d currentEField=this->operator[](i);
		// TODO: Make functions for indexing
		// TODO: Do more advanced under-relaxation?
		// TODO: Fix sign of eField in solve rather than here
		for (int j=0; j<NDIM; j++)
			eField[j] = 0.5*(-x[nVerts+i*NDIM+j]+currentEField[j]);
		this->setField(entities[i], eField);
		potentialField_ptr->setField(entities[i],
				0.5*(x[i]+potentialField_ptr->operator[](i)));
	}
#endif
#ifdef HAVE_MPI
		if (MPI::COMM_WORLD.Get_rank() == 0)
#endif
			cout << "linearized FEM iter " << numberOfSolves <<
				" : potentialChange = " << potentialChange << endl;
	numberOfSolves++;
	}
}

void ElectricField::calcField_Gatsonis(PotentialField potentialField) {
	for (int i=0; i<entities.size(); i++) {
		vector<entHandle> superCellFaces =
				mesh_ptr->getSuperCellFaces(entities[i]);
		vect3d eField(0.,0.,0.);
		vect3d point = mesh_ptr->getCoordinates(entities[i]);
		double volume=0;
		for (int j=0; j<superCellFaces.size(); j++)  {
			vect3d surfaceVector =
					mesh_ptr->getSurfaceVector(superCellFaces[j], point);
			double potential = potentialField.getAverageField(superCellFaces[j]);
			volume += mesh_ptr->getTetVolume(point, superCellFaces[j]);
			eField += potential*surfaceVector;
		}
		eField /= volume;
		this->setField(entities[i], eField);
	}
}


PotentialField::PotentialField(Mesh *inputMesh_ptr, string inputName)
	: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX) {
}

PotentialField::PotentialField(PotentialField potential, string inputName)
	: Field<double>(potential.mesh_ptr, inputName, iBase_VERTEX) {
	for (int i=0; i<entities.size(); i++) {
		this->setField(entities[i], potential.getField(entities[i]));
	}
}


void PotentialField::calcField(CodeField vertexType) {
	for (int i=0; i<entities.size(); i++) {
		vect3d point = mesh_ptr->getCoordinates(entities[i]);
		// TODO: Change this from just a test Coulomb field
//		double potential = -1./point.norm();
		double potential;
		if (vertexType.getField(entities[i])==4) {
			// TODO: don't hard-code sheath potential
			potential = -1./2.;
		} else if (vertexType.getField(entities[i])==5) {
			// TODO: don't hard-code boudnary potential
			potential = 0;
		} else {
			potential = 0;
//			potential = -0.77;
		}
		this->setField(entities[i], potential);
	}
}

void PotentialField::calcField(DensityField ionDensity,
		DensityField electronDensity, CodeField vertexType, FILE *outFile) {
	assert(ionDensity.mesh_ptr == mesh_ptr);
	assert(electronDensity.mesh_ptr == mesh_ptr);
	for (int i=0; i<entities.size(); i++) {
		// TODO: remove this restriction on vertices
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);
		vect3d xzPosition = nodePosition;

		double potential;
		// TODO: don't hard-code boundary type and quasi-neutral operation
//		if (vertexType.getField(entities[i])==4) {
//			// TODO: don't hard-code sheath potential
//			potential = -1./2.;
//		} else
		if (vertexType.getField(entities[i])==5) {
			// TODO: don't hard-code boudnary potential
			potential = 0;
		} else {
			potential = this->getField(entities[i]);
			potential += 1./2.*log(ionDensity.getField(entities[i])/
					electronDensity.getField(entities[i]));
		}
		this->setField(entities[i], potential);

//#ifdef HAVE_MPI
//		if (MPI::COMM_WORLD.Get_rank() == 0)
//#endif
//			cout << nodePosition.norm() << " " << potential << endl;
		if (outFile)
			fprintf(outFile, "%g %g\n", nodePosition.norm(), potential);
//		}
	}
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

void PotentialField::calcField(DensityField ionDensity,
		CodeField vertexType, FILE *outFile) {
	assert(ionDensity.mesh_ptr == mesh_ptr);
	for (int i=0; i<entities.size(); i++) {
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);

		double potential;
		// TODO: don't hard-code boundary type and quasi-neutral operation
//		if (vertexType.getField(entities[i])==4) {
//			// TODO: don't hard-code sheath potential
//			potential = -1./2.;
//		} else
		if (vertexType.getField(entities[i])==5) {
			// TODO: don't hard-code boundary potential
			potential = 0;
		} else {
			double currentPotential = this->getField(entities[i]);
			double boltzmannPotential = log(ionDensity.getField(entities[i]));
			potential = (currentPotential+boltzmannPotential)/2.;
		}
		this->setField(entities[i], potential);

		if (outFile)
			fprintf(outFile, "%g %g\n", nodePosition.norm(), potential);
	}
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

void PotentialField::calcField(DensityField ionDensity,
		DensityField ionDensityPP, DensityField ionDensityNP,
		DensityField electronDensity,
		DensityField electronDensityPP, DensityField electronDensityNP,
		CodeField vertexType, double positivePotentialPerturbation,
		double negativePotentialPerturbation, FILE *outFile) {
	assert(ionDensity.mesh_ptr == mesh_ptr);
	assert(electronDensity.mesh_ptr == mesh_ptr);
	for (int i=0; i<entities.size(); i++) {
		// TODO: remove this restriction on vertices
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);
		vect3d xzPosition = nodePosition;

		double potential;
		// TODO: don't hard-code boundary type and quasi-neutral operation
//		if (vertexType.getField(entities[i])==4) {
//			// TODO: don't hard-code sheath potential
//			potential = -1./2.;
//		} else
		if (vertexType.getField(entities[i])==5) {
			// TODO: don't hard-code boudnary potential
			potential = 0;
		} else {
			double ionDerivativePP = (ionDensityPP[i]-ionDensity[i])/
					positivePotentialPerturbation;
			double ionDerivativeNP = (ionDensityNP[i]-ionDensity[i])/
					negativePotentialPerturbation;
			double electronDerivativePP = (electronDensityPP[i]-electronDensity[i])/
					positivePotentialPerturbation;
			double electronDerivativeNP = (electronDensityNP[i]-electronDensity[i])/
					negativePotentialPerturbation;
			double potentialCorrectionPP = (ionDensity[i]-electronDensity[i])/
					(electronDerivativePP-ionDerivativePP);
			double potentialCorrectionNP = (ionDensity[i]-electronDensity[i])/
					(electronDerivativeNP-ionDerivativeNP);
			// TODO: set the max extrapolation multiplier globally?
			double limitedCorrectionPP = min(potentialCorrectionPP,
					2.*positivePotentialPerturbation);
			double limitedCorrectionNP = max(potentialCorrectionNP,
					2.*negativePotentialPerturbation);
			potential = this->getField(entities[i]);
			if (potentialCorrectionPP>0. && potentialCorrectionNP>0.) {
				// Only one relevant solution
				potential += limitedCorrectionPP;
			} else if (potentialCorrectionPP<0. && potentialCorrectionNP<0.) {
				// Only one relevant solution
					potential += limitedCorrectionNP;
			} else if (potentialCorrectionPP>0. && potentialCorrectionNP<0.) {
				// Choose closest solution
				potential += (potentialCorrectionPP<-potentialCorrectionNP) ?
						limitedCorrectionPP : limitedCorrectionNP;
			}
//			potential += 1./2.*log(ionDensity.getField(entities[i])/
//					electronDensity.getField(entities[i]));
		}
		this->setField(entities[i], potential);

//#ifdef HAVE_MPI
//		if (MPI::COMM_WORLD.Get_rank() == 0)
//#endif
//			cout << nodePosition.norm() << " " << potential << endl;
		if (outFile)
			fprintf(outFile, "%g %g\n", nodePosition.norm(), potential);
//		}
	}
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

DensityField::DensityField(Mesh *inputMesh_ptr, string inputName)
		: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX) {
}

void DensityField::calcField() {}

void DensityField::calcField(CodeField vertexType,
		PotentialField potential, double charge) {
	// TODO: get rid of hard-coding here
	for (int i=0; i<entities.size(); i++) {
		double density=1.;
		if (vertexType.getField(entities[i])==4) {
			density = -1./2.;
		} else if (vertexType.getField(entities[i])==5) {
			density = 1.;
		} else {
			density = exp(charge*potential[i]);
		}
		this->setField(entities[i], density);
	}
}

void DensityField::calcField(DensityField ionDensity,
		DensityField electronDensity) {
	for (int i=0; i<entities.size(); i++) {
		double density = ionDensity.getField(entities[i]) -
				electronDensity.getField(entities[i]);
		this->setField(entities[i], density);
	}
}

void DensityField::calcField(ElectricField& electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdge, double charge,
		double potentialPerturbation, FILE *outFile) {
	int mpiId = 0;
#ifdef HAVE_MPI
	double *density = new double[entities.size()];
	// TODO: Am assuming here that fields are identical on master and slaves
	mpiId = MPI::COMM_WORLD.Get_rank();
	if (mpiId == 0) {
		DensityField::requestDensityFromSlaves(electricField,
				potentialField, faceType, vertexType,
				shortestEdge, charge, potentialPerturbation, outFile);
		for (int node=0; node<entities.size(); node++) {
			density[node] = this->getField(entities[node]);
		}
	} else {
		DensityField::processDensityRequests(electricField,
				potentialField, faceType, vertexType,
				shortestEdge, charge, potentialPerturbation);
	}
	MPI::COMM_WORLD.Bcast(density, entities.size(), MPI::DOUBLE, 0);
	for (int node=0; node<entities.size(); node++) {
		this->setField(entities[node], density[node]);
	}
#else
	extern_findTet=0;
	clock_t startClock = clock(); // timing
	for (int node=0; node<entities.size(); node++) {
		double error;
	    double density = this->calculateDensity(node, electricField,
				potentialField, faceType, vertexType, shortestEdge,
				charge, potentialPerturbation, &error);
		this->setField(entities[node], density);
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
//		cout << nodePosition.norm() << " " << density << " " << error << endl;
		if (outFile)
			fprintf(outFile, "%g %g %g\n", nodePosition.norm(), density, error);
	}
	clock_t endClock = clock(); // timing
	cout << "calcField total (s)= "
			<< (double)(endClock-startClock)/(double)CLOCKS_PER_SEC << endl; // timing
	cout << "findTet total (s)= "
			<< (double)extern_findTet/(double)CLOCKS_PER_SEC << endl; // timing
//	cout << "checkIfInNewTet (s)= "
//			<< (double)extern_checkIfInNewTet/(double)CLOCKS_PER_SEC << endl; // timing
#endif
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

void DensityField::poissonCubeTest() {
	for (int i=0; i<entities.size(); i++) {
		vect3d xyz = mesh_ptr->getCoordinates(entities[i]);
		double chargeDensity = -0.9*sin(xyz[0]*M_PI)*sin(xyz[1]*M_PI)*sin(xyz[2]*M_PI);
		// TODO: don't hard-code debye length
		double debyeLength=1.;
		double density = debyeLength*debyeLength*chargeDensity
				+ exp(1./M_PI/M_PI/3.*chargeDensity);
//		cout << density << " " << xyz.transpose() << endl;
		this->setField(entities[i], density);
	}
}

double DensityField::calculateDensity(int node, ElectricField& electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdgeField, double charge,
		double potentialPerturbation, double *error) {
	// TODO: Need unified way of specifying unperturbed boundary plasma
	if (vertexType[node]==5) {
		return 1.;
//	} else if (vertexType[node]==4) {
//		// TODO: don't need sheath entrance density if specifying potential,
//		//       but might be interested in it or other moments later
//		// TODO: could make below analytic expression for planar electron dens.
//		return exp(-0.5);
	}
//	// TODO: remove this
//	if (node%100!=1) {
//		return 1.;
//	}
	double density=0.;
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
	IntegrandContainer integrandContainer;
	integrandContainer.mesh_ptr = mesh_ptr;
	integrandContainer.node = entities[node];
	integrandContainer.electricField_ptr = &electricField;
	integrandContainer.potentialField_ptr = &potentialField;
	integrandContainer.faceTypeField_ptr = &faceType;
	integrandContainer.vertexTypeField_ptr = &vertexType;
	integrandContainer.shortestEdgeField_ptr = &shortestEdgeField;
	stringstream fileNameStream;
//	fileNameStream << "distFunc/distributionFunction_r" << nodePosition.norm()
//			<< "_vert" << node << ".p3d";
	fileNameStream << "distFunc/distributionFunction_q" << charge << "_r"
			<< nodePosition.norm() << "_vert" << node << ".p3d";
	integrandContainer.outFile = NULL;
	stringstream fileNameStreamOrbit;
	fileNameStreamOrbit << "orbits/orbits_q" << charge << "_r"
			<< nodePosition.norm() << "_vert" << node << ".p3d";
	integrandContainer.orbitOutFile = NULL;
	bool doThisNode = true;
	// TODO: don't hard-code output nodes
//	if (node<5) {
//	if (node==0) {
//	if (node==5 || node==2540) {
//	if (node==105) {
//	if (nodePosition[2]<0.5 && nodePosition[2]>-0. && nodePosition[0]<2.5 &&
//			nodePosition[0]>1.8 && nodePosition[1]<1.4 && nodePosition[1]>0.8
//			&& charge>0.) {
//	if (vertexType[node]==4 && charge>0.) {
//	vect3d desiredNodePosition(0.908757, 0.411679, 0.0684252);
//	if ((desiredNodePosition-nodePosition).norm()<1e-2 && charge>0.) {
//		doThisNode = true;
//		integrandContainer.outFile = fopen(fileNameStream.str().c_str(), "w");
//		fprintf(integrandContainer.outFile, "x y z f\n");
//		integrandContainer.orbitOutFile =
//				fopen(fileNameStreamOrbit.str().c_str(), "w");
//		fprintf(integrandContainer.orbitOutFile, "x y z energy\n");
//	} else {
//		// TODO: decouple do and record
//		doThisNode = false;
//	}
	integrandContainer.charge = charge;
	int vdim=3;
	double xmin[vdim], xmax[vdim];
	for (int j=0; j<vdim; j++) {
		xmin[j] = -1.;
		xmax[j] = 1.;
	}
	// TODO: should make number of orbits adaptive
	int numberOfOrbits=100;
//	if (charge<0.)
//	if (charge>0.)
//		numberOfOrbits*=10;
	// TODO: make potential perturbation more robust, transparent, and flexible
	potentialField[node] += potentialPerturbation;
	int actualNumberOfOrbits=0;
	int failureType=0;
	double probabilityThatTrueError=0.;
	if (doThisNode) {
//		adapt_integrate(1, &distributionFunctionFromBoundary, (void*)&integrandContainer,
//				vdim, xmin, xmax, numberOfOrbits, 1.e-5, 1.e-5, &density, error);
		// TODO: Turn off smoothing flag bit
//		Vegas(NDIM, 1, &distributionFunctionFromBoundaryCuba,
//				(void*)&integrandContainer, 1.e-5, 1.e-5, 0, 0,
//				numberOfOrbits, numberOfOrbits, min(numberOfOrbits,1000), 1000, 1000, 0,
//				NULL, &actualNumberOfOrbits, &failureType,
//				&density, error, &probabilityThatTrueError);
		double moments[5];
		double errors[5];
		double probabilities[5];
		Vegas(NDIM, 1, &distributionFunctionFromBoundaryCuba,
				(void*)&integrandContainer, 1.e-5, 1.e-5, 0, 0,
				numberOfOrbits, numberOfOrbits, min(numberOfOrbits,1000), 1000, 1000, 0,
				NULL, &actualNumberOfOrbits, &failureType,
				moments, errors, probabilities);
		density = moments[0];
		*error = errors[0];
		probabilityThatTrueError = probabilities[0];
	}
	// TODO: make potential perturbation more robust, transparent, and flexible
	potentialField[node] -= potentialPerturbation;
	if (integrandContainer.outFile)
		fclose(integrandContainer.outFile);
	if (integrandContainer.orbitOutFile)
		fclose(integrandContainer.orbitOutFile);
	return density;
}

#ifdef HAVE_MPI
void DensityField::requestDensityFromSlaves(ElectricField& electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdge, double charge,
		double potentialPerturbation, FILE *outFile) {
	int nProcesses = MPI::COMM_WORLD.Get_size();
	MPI::Status status;
	int node=0;

	// Send one node to each process
	for (int rank=1; rank<nProcesses; ++rank) {
		if (node<entities.size()) {
			MPI::COMM_WORLD.Send(&node, 1, MPI::INT, rank, WORKTAG);
			node++;
		}
	}

	// Process incoming density and send new nodes until all done
	while (node<entities.size()) {
		MPI::Status status = this->receiveDensity(outFile);
		MPI::COMM_WORLD.Send(&node, 1, MPI::INT, status.Get_source(),
				WORKTAG);
		node++;
	}

	// Process any outstanding densities
	// TODO: could run into problems here if fewer nodes than processes?
	assert(nProcesses<entities.size());
	for (int rank=1; rank<nProcesses; ++rank) {
		this->receiveDensity(outFile);
	}

	// Send empty message with DIETAG to signal done with nodes
	for (int rank=1; rank<nProcesses; ++rank) {
		MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, DIETAG);
	}
}
#endif

#ifdef HAVE_MPI
MPI::Status DensityField::receiveDensity(FILE *outFile) {
	MPI::Status status;
	double density[2];
	MPI::COMM_WORLD.Recv(&density[0], 2, MPI::DOUBLE, MPI_ANY_SOURCE,
			MPI_ANY_TAG, status);
	int incomingNode = status.Get_tag();
	if (incomingNode>=0 && incomingNode<entities.size())
		this->setField(entities[incomingNode], density[0]);
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[incomingNode]);
//	cout << nodePosition.norm() << " " << density[0] << " " <<
//			density[1] << endl;
	if (outFile)
		fprintf(outFile, "%g %g %g\n", nodePosition.norm(), density[0], density[1]);
	return status;
}
#endif

#ifdef HAVE_MPI
void DensityField::processDensityRequests(ElectricField& electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdge, double charge,
		double potentialPerturbation) {
	MPI::Status status;
	int node;

	while (1) {
		MPI::COMM_WORLD.Recv(&node, 1, MPI::INT, 0, MPI_ANY_TAG,
				status);
		if (status.Get_tag() == DIETAG) {
			return;
		}
		double density[2];
		density[0] = this->calculateDensity(node, electricField,
				potentialField, faceType, vertexType,
				shortestEdge, charge, potentialPerturbation,
				&density[1]);
		MPI::COMM_WORLD.Send(&density[0], 2, MPI::DOUBLE, 0, node);
	}
}
#endif

CodeField::CodeField(Mesh *inputMesh_ptr, string inputName, int elementType)
		: Field<int>(inputMesh_ptr, inputName, elementType) {
}

void CodeField::calcField(Field<int> faceTypeField) {
	for (int i=0; i<entities.size(); i++) {
		// TODO: this does not tag regions with an edge (but not whole face)
		//       on the boundary...could do 2ndAdjacent through vertex
		vector<entHandle> faces =
				mesh_ptr->getAdjacentEntities(entities[i],iBase_FACE);
		int elementType=0;
		for (int j=0; j<faces.size(); j++) {
			int faceType = faceTypeField.getField(faces[j]);
			// TODO: handle case where two adjacent faces are different boundaries?
			if (faceType>0)
				elementType=faceType;
		}
		this->setField(entities[i],elementType);
	}
}

ShortestEdgeField::ShortestEdgeField(Mesh *inputMesh_ptr, string inputName)
	: Field<double>(inputMesh_ptr, inputName, iBase_REGION) {
}

void ShortestEdgeField::calcField() {
	for (int i=0; i<entities.size(); i++) {
		vector<vect3d> vertexVectors = mesh_ptr->getVertexVectors(entities[i]);
		double shortestEdge = 1./DELTA_LENGTH;
		for (int j=0; j<vertexVectors.size()-1; j++) {
			for (int k=j+1; k<vertexVectors.size(); k++) {
				vect3d relativePosition = vertexVectors[j]-vertexVectors[k];
				if (relativePosition.norm() < shortestEdge)
					shortestEdge = relativePosition.norm();
			}
		}
		this->setField(entities[i], shortestEdge);
	}
}
