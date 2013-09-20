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
		CodeField vertexType, double debyeLength, bool doLuDecomposition=true)
		: Field<vect3d>(inputMesh_ptr, inputName, iBase_VERTEX) {
	if (doLuDecomposition && debyeLength>0.) {
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
	} else {
#ifdef HAVE_MPI
	if (MPI::COMM_WORLD.Get_rank() == 0)
#endif
		cout << "Not computing the Poisson solve LU decomposition: doLuDecomposition = " <<
		doLuDecomposition << ", debyeLength = " << debyeLength << endl << endl;
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
		if (elementIndex!=i)
			throw string("problem in ElectricField::calcField");
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
		DensityField ionDensity, double debyeLength) {
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
					double shieldedPotential=vertexPosition.dot(extern_E);
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
			for (int k=0; k<adjacentVertices.size(); k++) {
				int kk = adjacentVertices[k];
				double basisBasisCoefficient = volume/20;
				if (jj==kk)
					basisBasisCoefficient *= 2.;
				// TODO: use referenceElectronDensity in electron response
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
						// TODO: use referenceElectronDensity in electron response
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


PotentialField::PotentialField(Mesh *inputMesh_ptr, string inputName,
		double boundaryPotential, double sheathPotential,
		bool fixSheathPotential, int numberOfComponents)
	: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX, numberOfComponents),
	  boundaryPotential(boundaryPotential),
	  sheathPotential(sheathPotential),
	  fixSheathPotential(fixSheathPotential) {
	referenceElectronDensity_ptr = NULL;
	referenceElectronTemperature_ptr = NULL;
}

PotentialField::PotentialField(PotentialField potential, string inputName)
	: Field<double>(potential.mesh_ptr, inputName, iBase_VERTEX),
	  boundaryPotential(potential.boundaryPotential),
	  sheathPotential(potential.sheathPotential),
	  fixSheathPotential(potential.fixSheathPotential) {
	for (int i=0; i<entities.size(); i++) {
		this->setField(entities[i], potential.getField(entities[i]));
	}
	referenceElectronDensity_ptr = NULL;
	referenceElectronTemperature_ptr = NULL;
}


void PotentialField::calcField(CodeField vertexType, double debyeLength,
		double boundaryPotential, double surfacePotential,
		double sheathPotential) {
	for (int i=0; i<entities.size(); i++) {
		vect3d point = mesh_ptr->getCoordinates(entities[i]);
		double potential;
		if (debyeLength==0.) {
			if (vertexType.getField(entities[i])==4) {
				potential = sheathPotential;
			} else if (vertexType.getField(entities[i])==5) {
				potential = boundaryPotential;
			} else {
				potential = boundaryPotential;
			}
		} else {
			// TODO: Change this from just a test Coulomb field
			cout << "WARNING: Assuming spherical mesh in PotentialField::calcField" << endl;
			double potential = surfacePotential/point.norm();
		}
		this->setField(entities[i], potential);
	}
}

void PotentialField::calcField(DensityField ionDensity,
		DensityField electronDensity, CodeField vertexType, FILE *outFile) {
	if (ionDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
	if (electronDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
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
			// TODO: don't hard-code boundary potential
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

void PotentialField::calcField(DensityField& ionDensity, Field<vect3d>& ionVelocity,
		CodeField& vertexType, FILE *outFile, double boundaryPotential,
		double sheathPotential, bool fixSheathPotential) {
	if (ionDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
	for (int i=0; i<entities.size(); i++) {
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);

		double potential;
		// TODO: don't hard-code boundary type
		if (vertexType.getField(entities[i])==5) {
			potential = boundaryPotential;
		} else {
			double currentPotential = this->getField(entities[i]);
			double referenceDensity = referenceElectronDensity_ptr->getField(entities[i]);
			double referenceElectronTemperature =
					referenceElectronTemperature_ptr->operator()(nodePosition);
			double boltzmannPotential = boundaryPotential +
					referenceElectronTemperature*
					log(ionDensity.getField(entities[i])/referenceDensity);
//			// Beam-like distribution doesn't give correct potential update
//			// TODO: don't hard-code parallel velocity cutoff
//			// TODO: replace this hack
//			if (fabs(ionVelocity.getField(entities[i]).dot(extern_B))<=0.5*extern_B.norm() ||
//					boltzmannPotential<currentPotential) {
				potential = (currentPotential+boltzmannPotential)/2.;
//			} else {
//				potential = currentPotential;
//			}
		}
		if (fixSheathPotential) {
			// TODO: make this function since do multiple places
			vect3d vertexPosition = mesh_ptr->getCoordinates(entities[i]);
			// TODO: hard-coding zero potential of ExB-field as origin here
			// TODO: sign of E appears to be wrong in mover, so flip here
			//       (should be + if one - from shielding and one from -grad(phi) )
			double shieldedPotential=sheathPotential-vertexPosition.dot(extern_E);
			// TODO: don't hard-code boundary type
			if (vertexType.getField(entities[i])==4) {
				potential = shieldedPotential;
			}
		}
		this->setField(entities[i], potential);

		if (outFile)
			fprintf(outFile, "%g %g\n", nodePosition.norm(), potential);
	}
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

void PotentialField::calcFieldAtNode(entHandle entity, double ionDensity, int vertexType,
		double boundaryPotential, double sheathPotential, bool fixSheathPotential) {
	vect3d nodePosition = mesh_ptr->getCoordinates(entity);
	double potential;
	// TODO: don't hard-code boundary type
	if (vertexType==5) {
		potential = boundaryPotential;
	} else {
		// TODO: separate this out to function since here and in calcField?
		double currentPotential = this->getField(entity);
		double referenceDensity = referenceElectronDensity_ptr->getField(entity);
		double referenceElectronTemperature =
				referenceElectronTemperature_ptr->operator()(nodePosition);
		double boltzmannPotential = boundaryPotential +
				referenceElectronTemperature*
				log(ionDensity/referenceDensity);
		potential = (currentPotential+boltzmannPotential)/2.;
	}
	if (fixSheathPotential) {
		// TODO: make this function since do multiple places
		vect3d vertexPosition = mesh_ptr->getCoordinates(entity);
		// TODO: hard-coding zero potential of ExB-field as origin here
		// TODO: sign of E appears to be wrong in mover, so flip here
		//       (should be + if one - from shielding and one from -grad(phi) )
		double shieldedPotential=sheathPotential-vertexPosition.dot(extern_E);
		// TODO: don't hard-code boundary type
		if (vertexType==4) {
			potential = shieldedPotential;
		}
	}
	this->setField(entity, potential);
}

void PotentialField::calcField(DensityField& ionDensity, DerivativeField& ionDensityDerivative,
		CodeField& vertexType, FILE *outFile, double boundaryPotential,
		double sheathPotential, bool fixSheathPotential) {
	if (ionDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
	for (int i=0; i<entities.size(); i++) {
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);

		double potential;
		// TODO: don't hard-code boundary type
		if (vertexType.getField(entities[i])==5) {
			potential = boundaryPotential;
		} else {
			double referenceDensity = referenceElectronDensity_ptr->getField(entities[i]);
			double referenceElectronTemperature =
					referenceElectronTemperature_ptr->operator()(nodePosition);
			double currentPotential = this->getField(entities[i]);
			double chargeDensityDerivative = ionDensityDerivative.getField(entities[i])
					-referenceDensity*exp(currentPotential/referenceElectronTemperature);
			double fractionToApply = 0.5;
			double denominator = chargeDensityDerivative;
			// TODO: bad parameter name
			if (fabs(denominator)<SMALL_DENSITY_DERIVATIVE)
				denominator = copysign(SMALL_DENSITY_DERIVATIVE,chargeDensityDerivative);
			double potentialCorrection = -(ionDensity.getField(entities[i]) -
					referenceDensity*exp(currentPotential/referenceElectronTemperature))/
					denominator;
			if (fabs(potentialCorrection)>LARGE_POTENTIAL_CHANGE)
				potentialCorrection = copysign(LARGE_POTENTIAL_CHANGE,potentialCorrection);
			potential = currentPotential + fractionToApply*potentialCorrection;
		}
		if (fixSheathPotential) {
			// TODO: think about possible influence of background E field
			// TODO: don't hard-code boundary type
			if (vertexType.getField(entities[i])==4) {
				potential = sheathPotential;
			}
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
	if (ionDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
	if (electronDensity.mesh_ptr!=mesh_ptr)
		throw string("mesh pointers not the same in potentialField::calcField");
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

void PotentialField::computePerturbedPotentials(double negativePerturbation,
		double positivePerturbation, double minPotential, double maxPotential) {
	boost::scoped_array<double> potentials_ptr(new double[numberOfComponents]);
	for (int i=0; i<entities.size(); i++) {
		double referencePotential = this->getField(entities[i]);
		potentials_ptr[0] = referencePotential;
		if (numberOfComponents>2) {
			double lowerLimit = min(minPotential,referencePotential+negativePerturbation);
			double upperLimit = max(maxPotential,referencePotential+positivePerturbation);
			double step = (upperLimit-lowerLimit)/(numberOfComponents-2);
			for (int j=1; j<numberOfComponents; j++) {
				potentials_ptr[j] = lowerLimit+(j-1)*step;
			}
		} else if (numberOfComponents==2) {
			// TODO: better behaviour here?
			potentials_ptr[1] = referencePotential+negativePerturbation;
		}
		this->setField(entities[i],potentials_ptr.get());
	}
}

void PotentialField::setReferenceElectronDensity(DensityField& referenceElectronDensity){
	this->referenceElectronDensity_ptr = &referenceElectronDensity;
}

void PotentialField::setReferenceElectronTemperature(SpatialDependence& referenceElectronTemperature){
	this->referenceElectronTemperature_ptr = &referenceElectronTemperature;
}


DensityField::DensityField(Mesh *inputMesh_ptr, string inputName,
		Field<vect3d> *inputAverageVelocity_ptr, Field<double> *inputTemperature_ptr,
		int numberOfComponents)
		: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX, numberOfComponents) {
	averageVelocity_ptr = inputAverageVelocity_ptr;
	temperature_ptr = inputTemperature_ptr;
	distributionFunction_ptr = NULL;
}

void DensityField::calcField() {}

void DensityField::calcField(CodeField& vertexType,
		PotentialField& potential, DensityField& referenceDensity,
		SpatialDependence& referenceTemperature, double charge) {
	// TODO: get rid of hard-coding here
	for (int i=0; i<entities.size(); i++) {
		vect3d position=mesh_ptr->getCoordinates(entities[i]);
		double density=referenceDensity[entities[i]];
		double temperature=referenceTemperature(position);
		// TODO: this isn't very general
		if (vertexType.getField(entities[i])==4) {
			density *= -1./2.;
		} else if (vertexType.getField(entities[i])==5) {
			density *= 1.;
		} else {
			density *= exp(potential[i]/temperature);
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
		PotentialField *potentialField_ptr, DensityField& referenceDensity,
		Field<int>& faceType, CodeField& vertexType,
		ShortestEdgeField& shortestEdge, double charge,
		double potentialPerturbation, bool doAllNodes,
		double unconvergednessThreshold, FILE *outFile) {
	vector<pair<int,double> > unconvergednessPairs;
	vector<double> absDensityDifferences;
	for (int i=0; i<entities.size(); i++) {
		double boltzmannDensity = referenceDensity.getField(entities[i])*
				exp(potentialField_ptr->getField(entities[i]));
		double densityDifference = this->getField(entities[i])-boltzmannDensity;
		absDensityDifferences.push_back(fabs(densityDifference));
//		// TODO: better to use fractional difference?
//		unconvergednessPairs.push_back(make_pair(i,fabs(densityDifference)));
		// TODO: clean up if really doing random ordering
		unconvergednessPairs.push_back(make_pair(i,rand()/double(RAND_MAX)));
	}
	sort(unconvergednessPairs.begin(),unconvergednessPairs.end(),
			boost::bind(&std::pair<int, double>::second, _1) >
			boost::bind(&std::pair<int, double>::second, _2));
	vector<int> sortedNodes;
	for (int i=0; i<entities.size(); i++) {
		double absDensityDifference=absDensityDifferences[unconvergednessPairs[i].first];
		if (doAllNodes || absDensityDifference>unconvergednessThreshold)
			sortedNodes.push_back(unconvergednessPairs[i].first);
	}
	int mpiId = 0;
#ifdef HAVE_MPI
	// TODO: Change this to using smart pointer
	boost::scoped_array<double> density(new double[entities.size()]);
	boost::scoped_array<vect3d> averageVelocity(new vect3d[entities.size()]);
	boost::scoped_array<double> temperature(new double[entities.size()]);
	// TODO: Am assuming here that fields are identical on master and slaves
	mpiId = MPI::COMM_WORLD.Get_rank();
	if (mpiId == 0) {
		DensityField::requestDensityFromSlaves(electricField,
				potentialField_ptr, sortedNodes, faceType, vertexType,
				shortestEdge, charge, potentialPerturbation, outFile);
		for (int node=0; node<entities.size(); node++) {
			density[node] = this->getField(entities[node]);
			averageVelocity[node] = averageVelocity_ptr->getField(entities[node]);
			temperature[node] = temperature_ptr->getField(entities[node]);
		}
	} else {
		DensityField::processDensityRequests(electricField,
				potentialField_ptr, referenceDensity, faceType, vertexType,
				shortestEdge, charge, potentialPerturbation);
	}
	MPI::COMM_WORLD.Bcast(density.get(), entities.size(), MPI::DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(averageVelocity.get(), entities.size()*sizeof(vect3d), MPI::BYTE, 0);
	MPI::COMM_WORLD.Bcast(temperature.get(), entities.size(), MPI::DOUBLE, 0);
	for (int node=0; node<entities.size(); node++) {
		this->setField(entities[node], density[node]);
		averageVelocity_ptr->setField(entities[node], averageVelocity[node]);
		temperature_ptr->setField(entities[node], temperature[node]);
	}
#else
	extern_findTet=0;
	clock_t startClock = clock(); // timing
	for (int i=0; i<sortedNodes.size(); i++) {
		int node=sortedNodes[i];
		boost::scoped_array<double> densities(new double[numberOfComponents]);
		boost::scoped_array<double> densityErrors(new double[numberOfComponents]);
		boost::scoped_array<vect3d> averageVelocities(new vect3d[numberOfComponents]);
		boost::scoped_array<vect3d> averageVelocityErrors(new vect3d[numberOfComponents]);
		boost::scoped_array<double> temperatures(new double[numberOfComponents]);
		boost::scoped_array<double> temperatureErrors(new double[numberOfComponents]);
	    this->calculateDensity(node, electricField,
				*potentialField_ptr, referenceDensity,
				faceType, vertexType, shortestEdge,
				charge, potentialPerturbation, densities.get(), densityErrors.get(),
				averageVelocities.get(), averageVelocityErrors.get(),
				temperatures.get(), temperatureErrors.get());
		this->setField(entities[node], densities.get());
		averageVelocity_ptr->setField(entities[node], averageVelocities.get());
		temperature_ptr->setField(entities[node], temperatures.get());
		// TODO: handle multi-component case or treat in more transparent manner?
//		if (numberOfComponents==1) {
//			potentialField_ptr->calcFieldAtNode(entities[node],densities[0],vertexType[node],
//					potentialField_ptr->boundaryPotential,potentialField_ptr->sheathPotential,
//					potentialField_ptr->fixSheathPotential);
//		}
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
//		cout << nodePosition.norm() << " " << density << " " << error << endl;
		if (outFile)
			fprintf(outFile, "%g %g %g\n", nodePosition.norm(), densities[0], densityErrors[0]);
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

// TODO: figure out why get error if put include in typesAndDefinitions.h
#include <boost/numeric/quadrature/adaptive.hpp> // Not true boost library
//#include <boost/numeric/quadrature/kronrodgauss.hpp> // Not true boost library
// TODO: make reference density its own class?
void DensityField::calcField(CodeField& vertexType,
		DistributionFunction& distributionFunction, double charge) {
	for (int i=0; i<entities.size(); i++) {
		vect3d position = mesh_ptr->getCoordinates(entities[i]);
//		DistributionFunctionVy distFuncVy;
//		distFuncVy.distFunc_ptr = &distributionFunction;
//		distFuncVy.pos = position;
		// TODO: shouldn't really be altering things passed by reference
		distributionFunction.setPresetPosition(position);
		double density, densityError;
//		boost::numeric::quadrature::adaptive()(distFuncVy,-1.,1.,
//				density, densityError);
		boost::numeric::quadrature::adaptive()(distributionFunction,-1.,1.,
				density, densityError);
//		boost::numeric::quadrature::adaptive()(dynamic_cast<const Maxwellian&>(distributionFunction),-1.,1.,
//				density, densityError);
//		// TODO: figure out why 15 works but not 20 (only odd?)
//		//       "N is defined for 15, 21, 31, 41 and 51"
//		boost::numeric::quadrature::kronrod_gauss<15> integrator;
//		integrator(dynamic_cast<const DistributionFunction&>(distribFunc),-1.,1.,
//				density);
		if (densityError/density>0.0001) {
			cout << "Large reference density error: " << densityError <<
					", n= " << density << ", pos= " << position.transpose() << endl;
		}
		this->setField(entities[i], density);
	}

}

void DensityField::poissonCubeTest(double debyeLength) {
	for (int i=0; i<entities.size(); i++) {
		vect3d xyz = mesh_ptr->getCoordinates(entities[i]);
		double chargeDensity = -0.9*sin(xyz[0]*M_PI)*sin(xyz[1]*M_PI)*sin(xyz[2]*M_PI);
		double density = debyeLength*debyeLength*chargeDensity
				+ exp(1./M_PI/M_PI/3.*chargeDensity);
//		cout << density << " " << xyz.transpose() << endl;
		this->setField(entities[i], density);
	}
}

void DensityField::calculateDensity(int node, ElectricField& electricField,
		PotentialField& potentialField, DensityField& referenceDensity,
		Field<int>& faceType, CodeField& vertexType,
		ShortestEdgeField& shortestEdgeField, double charge,
		double potentialPerturbation, double *density, double *error,
		vect3d *averageVelocity, vect3d *averageVelocityError,
		double *temperature, double *temperatureError) {
	for (int k=0; k<numberOfComponents; k++)
		density[k] = 0.;
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
	bool doThisNode = false;
	bool recordThisNode = false;
	if (extern_evalPositions_ptr->size()==0) {
		doThisNode = true;
	} else {
		// TODO: checking every evalPosition for every node is not efficient
		//       (evalPosition list probably short enough not to matter)
		for (int i=0; i<extern_evalPositions_ptr->size(); i++) {
			vect3d desiredNodePosition=extern_evalPositions_ptr->operator[](i);
			// TODO: don't hard-code distance threshold
			doThisNode = doThisNode ||
					(desiredNodePosition-nodePosition).norm()<NODE_DISTANCE_THRESHOLD;
			// TODO: decouple do and record
			recordThisNode = doThisNode;
		}
	}
	// TODO: Need unified way of specifying unperturbed boundary plasma
	if (vertexType[node]==5) {
		for (int k=0; k<numberOfComponents; k++) {
			averageVelocity[k] = -extern_VEXB;
			averageVelocityError[k] = vect3d(0.,0.,0.);
			temperature[k] = 1.;
			temperatureError[k] = 0.;
			density[k] = referenceDensity[node];
		}
//	} else if (vertexType[node]==4) {
//		// TODO: don't need sheath entrance density if specifying potential,
//		//       but might be interested in it or other moments later
//		// TODO: could make below analytic expression for planar electron dens.
//		return exp(-0.5);
	} else {
		IntegrandContainer integrandContainer;
		integrandContainer.mesh_ptr = mesh_ptr;
		integrandContainer.node = entities[node];
		integrandContainer.electricField_ptr = &electricField;
		integrandContainer.potentialField_ptr = &potentialField;
		integrandContainer.faceTypeField_ptr = &faceType;
		integrandContainer.vertexTypeField_ptr = &vertexType;
		integrandContainer.shortestEdgeField_ptr = &shortestEdgeField;
		integrandContainer.distributionFunction_ptr = distributionFunction_ptr;
		integrandContainer.outFile = NULL;
		integrandContainer.orbitOutFile = NULL;
		if (recordThisNode) {
			stringstream fileNameStream;
//			fileNameStream << "distFunc/distributionFunction_r" << nodePosition.norm()
//					<< "_vert" << node << ".p3d";
			fileNameStream << "distFunc/distributionFunction_q" << charge << "_r"
					<< nodePosition.norm() << "_vert" << node << ".p3d";
			integrandContainer.outFile = fopen(fileNameStream.str().c_str(), "w");
			fprintf(integrandContainer.outFile, "x y z f\n");
		}
		if (extern_saveOrbits) {
			stringstream fileNameStreamOrbit;
			fileNameStreamOrbit << "orbits/orbits_q" << charge << "_r"
					<< nodePosition.norm() << "_vert" << node << ".p3d";
			integrandContainer.orbitOutFile =
					fopen(fileNameStreamOrbit.str().c_str(), "w");
			fprintf(integrandContainer.orbitOutFile, "x y z energy\n");
		}
		integrandContainer.charge = charge;
		int vdim=3;
		double xmin[vdim], xmax[vdim];
		for (int j=0; j<vdim; j++) {
			xmin[j] = -1.;
			xmax[j] = 1.;
		}
		// TODO: should make number of orbits adaptive
		int numberOfOrbits=1000;
//		if (charge<0.)
//		if (charge>0.)
//			numberOfOrbits*=10;
		// TODO: make potential perturbation more robust, transparent, and flexible
		potentialField[node] += potentialPerturbation;
		int actualNumberOfOrbits=0;
		int failureType=0;
		double probabilityThatTrueError=0.;
		if (doThisNode) {
//			adapt_integrate(1, &distributionFunctionFromBoundary, (void*)&integrandContainer,
//					vdim, xmin, xmax, numberOfOrbits, 1.e-5, 1.e-5, &density, error);
			// TODO: Turn off smoothing flag bit
//			Vegas(NDIM, 1, &distributionFunctionFromBoundaryCuba,
//					(void*)&integrandContainer, 1.e-5, 1.e-5, 0, 0,
//					numberOfOrbits, numberOfOrbits, min(numberOfOrbits,1000), 1000, 1000, 0,
//					NULL, &actualNumberOfOrbits, &failureType,
//					&density, error, &probabilityThatTrueError);
			if (potentialField.numberOfComponents!=numberOfComponents)
				throw string("number of components differ in density and potential");
			boost::scoped_array<double> potentials_ptr(new double[numberOfComponents]);
			potentialField.getField(entities[node],potentials_ptr.get());
			for (int k=0; k<numberOfComponents; k++) {
				potentialField[node] = potentials_ptr[k];
				double moments[5];
				double errors[5];
				double probabilities[5];
				Vegas(NDIM, 1+NDIM+1, &distributionFunctionFromBoundaryCuba,
						(void*)&integrandContainer, 1.e-5, 1.e-5, 0, 0,
						numberOfOrbits, numberOfOrbits, min(numberOfOrbits,1000), 1000, 1000, 0,
						NULL, &actualNumberOfOrbits, &failureType,
						moments, errors, probabilities);
				// TODO: just set NCOMP higher?
				if (failureType<0)
					cout << "failureType: " << failureType <<
					" (Probably need to change NCOMP in cuba/Makefile.am)" << endl;
				density[k] = moments[0];
				error[k] = errors[0];
				// TODO: assuming NDIM==3
				if (NDIM!=3)
					throw string("can only handle NDIM==3");
				for (int i=0; i<NDIM; i++) {
					averageVelocity[k][i] = moments[i+1]/density[k];
					// TODO: include density error
					averageVelocityError[k][i] = errors[i+1]/density[k];
				}
				// Subtract ordered kinetic energy from second moment to get temperature
				temperature[k] = moments[4]/density[k] - pow(averageVelocity[k].norm(),2.)/3.;
				// TODO: calculate temperature error more carefully
				temperatureError[k] = sqrt(pow(errors[4]/density[k],2.) + pow(errors[0],2.) +
						pow(0.5*averageVelocityError[k].norm(),2.));
				probabilityThatTrueError = probabilities[0];
			}
			potentialField[node] = potentials_ptr[0];
		}
		// TODO: make potential perturbation more robust, transparent, and flexible
		potentialField[node] -= potentialPerturbation;
		if (integrandContainer.outFile)
			fclose(integrandContainer.outFile);
		if (integrandContainer.orbitOutFile)
			fclose(integrandContainer.orbitOutFile);
	}
//	// TODO: debugging
//	if (recordThisNode) {
//		cout << nodePosition.transpose() << " " << node << " " << vertexType[node] <<
//				" " << vertexType[entities[node]] << " " << density << endl;
//	}
}

#ifdef HAVE_MPI
void DensityField::requestDensityFromSlaves(ElectricField& electricField,
		PotentialField *potentialField_ptr, vector<int>& sortedNodes,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdge, double charge,
		double potentialPerturbation, FILE *outFile) {
	int nProcesses = MPI::COMM_WORLD.Get_size();
	MPI::Status status;
	int nodeCounter=0;
	int node=-1;
	double potential;
//	boost::scoped_array<double> potentials(new double[entities.size()]);
//	for (int i=0; i<entities.size(); i++) {
//		potentials[i] = potentialField_ptr->getField(entities[i]);
//	}

	// Send one node to each process (some may not get one)
	for (int rank=1; rank<min(nProcesses,int(sortedNodes.size()+1)); ++rank) {
		if (nodeCounter<sortedNodes.size()) {
			node = sortedNodes[nodeCounter];
			nodeCounter++;
			MPI::COMM_WORLD.Send(&node, 1, MPI::INT, rank, WORKTAG);
//			MPI::COMM_WORLD.Send(potentials.get(), entities.size(), MPI::DOUBLE,
//					rank, WORKTAG);
		} else {
			throw string("nodeCounter out of bounds in requestDensityFromSlaves");
		}
	}

	// Process incoming density and send new nodes until all done
	while (nodeCounter<sortedNodes.size()) {
		node = sortedNodes[nodeCounter];
		nodeCounter++;
		status = this->receiveDensity(&potential,outFile);
//		// TODO: deal with multi-component case or make more transparent
//		if (numberOfComponents==1) {
//			potentialField_ptr->setField(entities[status.Get_tag()],potential);
//			potentials[status.Get_tag()] = potential;
//		}
		MPI::COMM_WORLD.Send(&node, 1, MPI::INT, status.Get_source(),
				WORKTAG);
//		MPI::COMM_WORLD.Send(potentials.get(), entities.size(), MPI::DOUBLE,
//				status.Get_source(), WORKTAG);
	}

	// Process any outstanding densities
	for (int rank=1; rank<min(nProcesses,int(sortedNodes.size()+1)); ++rank) {
		status = this->receiveDensity(&potential,outFile);
//		// TODO: deal with multi-component case or make more transparent
//		if (numberOfComponents==1) {
//			potentialField_ptr->setField(entities[status.Get_tag()],potential);
//		}
	}

	// Send empty message with DIETAG to signal done with nodes
	for (int rank=1; rank<nProcesses; ++rank) {
		MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, DIETAG);
	}
}
#endif

#ifdef HAVE_MPI
MPI::Status DensityField::receiveDensity(double *potential, FILE *outFile) {
	MPI::Status status;
	boost::scoped_array<double> densities(new double[numberOfComponents]);
	boost::scoped_array<double> densityErrors(new double[numberOfComponents]);
	boost::scoped_array<vect3d> averageVelocities(new vect3d[numberOfComponents]);
	boost::scoped_array<vect3d> averageVelocityErrors(new vect3d[numberOfComponents]);
	vect3d incomingNodePosition;
	boost::scoped_array<double> temperatures(new double[numberOfComponents]);
	boost::scoped_array<double> temperatureErrors(new double[numberOfComponents]);
	MPI::COMM_WORLD.Recv(densities.get(), numberOfComponents, MPI::DOUBLE, MPI_ANY_SOURCE,
			MPI_ANY_TAG, status);
	int node = status.Get_tag();
	int source = status.Get_source();
	MPI::COMM_WORLD.Recv(densityErrors.get(), numberOfComponents, MPI::DOUBLE, source, node, status);
	MPI::COMM_WORLD.Recv(averageVelocities.get(), numberOfComponents*sizeof(vect3d), MPI::BYTE,
			source, node, status);
	MPI::COMM_WORLD.Recv(averageVelocityErrors.get(), numberOfComponents*sizeof(vect3d), MPI::BYTE,
			source, node, status);
	MPI::COMM_WORLD.Recv(&incomingNodePosition, sizeof(vect3d), MPI::BYTE,
			source, node, status);
	MPI::COMM_WORLD.Recv(temperatures.get(), numberOfComponents, MPI::DOUBLE, source, node, status);
	MPI::COMM_WORLD.Recv(temperatureErrors.get(), numberOfComponents, MPI::DOUBLE, source, node, status);
	if (node>=0 && node<entities.size()) {
		this->setField(entities[node], densities.get());
		averageVelocity_ptr->setField(entities[node], averageVelocities.get());
		temperature_ptr->setField(entities[node], temperatures.get());
//		// TODO: handle multi-component case or treat in more transparent manner?
//		if (numberOfComponents==1) {
//			MPI::COMM_WORLD.Recv(potential, numberOfComponents, MPI::DOUBLE, source, node, status);
//		}
	}
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
	if ((nodePosition-incomingNodePosition).norm()>NODE_DISTANCE_THRESHOLD) {
		cout << "ERROR: node position on slave different than on master " <<
				nodePosition.transpose() << " vs " <<
				incomingNodePosition.transpose() << endl;
	}
//	cout << nodePosition.norm() << " " << density[0] << " " <<
//			density[1] << endl;
	if (outFile)
		fprintf(outFile, "%g %g %g\n", nodePosition.norm(), densities[0], densityErrors[0]);

	return status;
}
#endif

#ifdef HAVE_MPI
void DensityField::processDensityRequests(ElectricField& electricField,
		PotentialField *potentialField_ptr, DensityField& referenceDensity,
		Field<int> faceType, CodeField vertexType,
		ShortestEdgeField shortestEdge, double charge,
		double potentialPerturbation) {
	MPI::Status status;
	int node;
//	boost::scoped_array<double> potentials(new double[entities.size()]);

	while (1) {
		MPI::COMM_WORLD.Recv(&node, 1, MPI::INT, 0, MPI_ANY_TAG,
				status);
		if (status.Get_tag() == DIETAG) {
			return;
		}
//		MPI::COMM_WORLD.Recv(potentials.get(), entities.size(), MPI::DOUBLE, 0, MPI_ANY_TAG,
//				status);
//		// TODO: handle multi-component case or treat in more transparent manner?
//		if (numberOfComponents==1) {
//			for (int i=0; i<entities.size(); i++) {
//				potentialField_ptr->setField(entities[i],potentials[i]);
//			}
//		}
		boost::scoped_array<double> densities(new double[numberOfComponents]);
		boost::scoped_array<double> densityErrors(new double[numberOfComponents]);
		boost::scoped_array<vect3d> averageVelocities(new vect3d[numberOfComponents]);
		boost::scoped_array<vect3d> averageVelocityErrors(new vect3d[numberOfComponents]);
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
		boost::scoped_array<double> temperatures(new double[numberOfComponents]);
		boost::scoped_array<double> temperatureErrors(new double[numberOfComponents]);
		this->calculateDensity(node, electricField,
				*potentialField_ptr, referenceDensity,
				faceType, vertexType,
				shortestEdge, charge, potentialPerturbation,
				densities.get(), densityErrors.get(),
				averageVelocities.get(), averageVelocityErrors.get(),
				temperatures.get(), temperatureErrors.get());
		// TODO: just send one big data-structure?
		MPI::COMM_WORLD.Send(densities.get(), numberOfComponents, MPI::DOUBLE, 0, node);
		MPI::COMM_WORLD.Send(densityErrors.get(), numberOfComponents, MPI::DOUBLE, 0, node);
		MPI::COMM_WORLD.Send(averageVelocities.get(), numberOfComponents*sizeof(vect3d), MPI::BYTE, 0, node);
		MPI::COMM_WORLD.Send(averageVelocityErrors.get(), numberOfComponents*sizeof(vect3d), MPI::BYTE, 0, node);
		MPI::COMM_WORLD.Send(&nodePosition, sizeof(vect3d), MPI::BYTE, 0, node);
		MPI::COMM_WORLD.Send(temperatures.get(), numberOfComponents, MPI::DOUBLE, 0, node);
		MPI::COMM_WORLD.Send(temperatureErrors.get(), numberOfComponents, MPI::DOUBLE, 0, node);
//		// TODO: handle multi-component case or treat in more transparent manner?
//		if (numberOfComponents==1) {
//			potentialField_ptr->calcFieldAtNode(entities[node],densities[0],vertexType[node],
//					potentialField_ptr->boundaryPotential,potentialField_ptr->sheathPotential,
//					potentialField_ptr->fixSheathPotential);
//			double potential=potentialField_ptr->getField(entities[node]);
//			MPI::COMM_WORLD.Send(&potential, numberOfComponents, MPI::DOUBLE, 0, node);
//		}
	}
}
#endif

void DensityField::setDistributionFunction(DistributionFunction& distributionFunction){
	this->distributionFunction_ptr = &distributionFunction;
}


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

DerivativeField::DerivativeField(Mesh *inputMesh_ptr, string inputName, int elementType)
		: Field<double>(inputMesh_ptr, inputName, elementType) {
}

void DerivativeField::calcField(Field<double>& numeratorValue, Field<double>& numeratorReference,
		Field<double>& denominatorValue, Field<double>& denominatorReference) {
	for (int i=0; i<entities.size(); i++) {
		double numerator = numeratorValue[i]-numeratorReference[i];
		double denominator = denominatorValue[i]-denominatorReference[i];
		if (fabs(denominator)<SMALL_DENOMINATOR)
			denominator = copysign(SMALL_DENOMINATOR,denominator);
		double derivative = numerator/denominator;
		this->setField(entities[i],derivative);
	}
}

DensityDerivativeField::DensityDerivativeField(Mesh *inputMesh_ptr, string inputName, int elementType)
		: DerivativeField(inputMesh_ptr, inputName, elementType) {
}

void DensityDerivativeField::calcField(Field<double>& numeratorValue, Field<double>& numeratorReference,
		Field<double>& denominatorValue, Field<double>& denominatorReference) {
	for (int i=0; i<entities.size(); i++) {
		double numerator = numeratorValue[i]-numeratorReference[i];
		double denominator = denominatorValue[i]-denominatorReference[i];
		// TODO: might be taking derivative wrt to something other than potential?
		double derivative;
		if ((fabs(numerator)<SMALL_DENSITY_CHANGE) || (fabs(denominator)<SMALL_POTENTIAL_CHANGE)) {
			// TODO: this forces small steps (and oscillations if sign wrong) near convergence?
			derivative = -1.;
		} else {
			derivative = numerator/denominator;
		}
		this->setField(entities[i],derivative);
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
