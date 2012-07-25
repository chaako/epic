/*
 * Field.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef FIELD_H_
#define FIELD_H_

#include "typesAndDefinitions.h"
//#include <stdio.h>
////#include <type_traits> // Requires -std=c++0x compiler flag
//#include <boost/type_traits.hpp>
//
//#ifdef HAVE_MPI
//#include "mpi.h"
//#endif
//
//#include "iMesh.h"
//#include "Eigen/Dense"
#include "Mesh.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, string inputName,
			int entityDimension);
	// TODO: Does an empty destructor cause a memory leak?
	virtual ~Field() {}

	T& operator[](entHandle& entity);
	T& operator[](int& entityIndex);

	T getField(vect3d position, entHandle *entity=NULL,
			int interpolationOrder=INTERPOLATIONORDER);
	T getField(entHandle node);
	T getFieldFromMeshDB(entHandle entityHandle);
	T getAverageField(entHandle element);
	void evalFieldAndDeriv(T *fieldValue,
			Eigen::Matrix<T,NDIM,1> *fieldDeriv,
			const vect3d &position, int *entityIndex=-1,
			int interpolationOrder=INTERPOLATIONORDER);
	void setField(entHandle node, T field);
	Eigen::VectorXd getErrorCoefficients(entHandle element,
			int interpolationOrder);
	Eigen::VectorXd getErrorCoefficients(int elementIndex,
			int interpolationOrder);
	void getErrorCoefficients(int elementIndex,
			int interpolationOrder,
			Eigen::VectorXd *errorCoefficients);

	Mesh *mesh_ptr;
	string name;
	iBase_TagHandle tag;
	vector<entHandle> &entities;
	vector<T> values;
	// TODO: if precalculating errorCoefficients,
	//       need to update when field changed
//	vector<Eigen::VectorXd> errorCoefficientsVectors;

	int currentRegionIndex;
	vector<int> currentVerticesIndices;
	int currentInterpolationElementIndex;
	int currentInterpolationOrder;
	vect3d currentPosition;
	T currentFieldValue;
	Eigen::Matrix<T,NDIM,1> currentFieldDerivatives;
//	Eigen::Vector4d currentLinearBasisFunctions;
	Eigen::VectorXd currentErrorBases;

	entHandle currentElement;
	vector<entHandle> currentVertices;
	vector<T> currentFields;
	entHandle currentInterpolationElement;
	Eigen::VectorXd currentErrorCoefficients;

	entHandle previousElement;
	vector<entHandle> previousVertices;
	vector<T> previousFields;
	entHandle previousInterpolationElement;
	Eigen::VectorXd previousErrorCoefficients;
};

class ElectricField : public Field<vect3d> {
public:
	ElectricField(Mesh *inputMesh_ptr, string inputName);
	virtual ~ElectricField() {}

	void calcField(PotentialField potentialField);
};

class CodeField : public Field<int> {
public:
	CodeField(Mesh *inputMesh_ptr, string inputName,
			int elementType);
	virtual ~CodeField() {}

	void calcField(Field<int> faceTypeField);
};

class DensityField : public Field<double> {
	// Charge density
public:
	DensityField(Mesh *inputMesh_ptr, string inputName);
	virtual ~DensityField() {}

	void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity);
	void calcField(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation, FILE *outFile=NULL);
	double calculateDensity(int node, ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation, double *error);
#ifdef HAVE_MPI
	void requestDensityFromSlaves(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation, FILE *outFile);
	MPI::Status receiveDensity(FILE *outFile);
	void processDensityRequests(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation);
#endif

};

class PotentialField : public Field<double> {
public:
	PotentialField(Mesh *inputMesh_ptr, string inputName);
	PotentialField(PotentialField potential, string inputName);
	virtual ~PotentialField() {}

	void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity,
			CodeField vertexType, FILE *outFile);
	void calcField(DensityField ionDensity,
			DensityField ionDensityPP, DensityField ionDensityNP,
			DensityField electronDensity,
			DensityField electronDensityPP, DensityField electronDensityNP,
			CodeField vertexType, FILE *outFile);

};

class ShortestEdgeField : public Field<double> {
public:
	ShortestEdgeField(Mesh *inputMesh_ptr, string inputName);
	virtual ~ShortestEdgeField() {}

	void calcField();
};

// gcc doesn't implement the export keyword, so define template functions here
template <class T>
Field<T>::Field(Mesh *inputMesh_ptr, string inputName,
		int entityDimension) : entities(
				inputMesh_ptr->entitiesVectors[entityDimension]) {
	mesh_ptr = inputMesh_ptr;
	name = inputName;
	int ierr;
	int size, type;
//	if (is_same<T,double>::value) {
	if (boost::is_same<T,double>::value) {
		size = 1;
		type = iBase_DOUBLE;
//	} else if (is_same<T,int>::value) {
	} else if (boost::is_same<T,int>::value) {
		size = 1;
		type = iBase_INTEGER;
	} else {
		size = (int)sizeof(T);
		type = iBase_BYTES;
	}
	tag = mesh_ptr->createTagHandle(name, size, type);
	// TODO: consider more transparent handling of exiting tag here
	tag = mesh_ptr->getTagHandle(name);
//	entities = mesh_ptr->getEntities(entityDimension);
	values = vector<T>(entities.size());
	int numberOfUninitializedFields=0;
	int numberOfInitializedFields=0;
	for (int i=0; i<entities.size(); i++) {
		try {
			values[i] = this->getFieldFromMeshDB(entities[i]);
			numberOfInitializedFields++;
		} catch (int FAILURE_GETTING_FIELD) {
			numberOfUninitializedFields++;
		}
	}
//	cout << "The field '" << name << "' was created with " <<
//			numberOfInitializedFields << " values set from mesh DB, and " <<
//			numberOfUninitializedFields << " uninitialized values" << endl;

	// TODO: nVertices should be standardized across functions
	int nVertices = 4;
	currentRegionIndex=-1;
	currentVerticesIndices = vector<int>(nVertices);
	currentInterpolationElementIndex = -1;
	currentInterpolationOrder = -1;
	currentPosition = vect3d(0.,0.,0.);
	currentFieldValue = T();
	currentFieldDerivatives = Eigen::Matrix<T,NDIM,1>(T(),T(),T());
//	currentLinearBasisFunctions = Eigen::Vector4d(0.,0.,0.,0.);
	currentErrorBases = Eigen::VectorXd(N_BASES_QUADRATIC);
	currentErrorCoefficients = Eigen::VectorXd(N_BASES_QUADRATIC);

	currentElement = NULL;
	currentFields = vector<T>(nVertices);
	currentInterpolationElement = NULL;

	previousElement = NULL;
	previousFields = vector<T>(nVertices);
	previousInterpolationElement = NULL;

//	for (int i=0; i<entities.size(); i++) {
//		errorCoefficientsVectors.push_back(
//				this->getErrorCoefficients(i, INTERPOLATIONORDER));
//	}
}

template <class T>
T& Field<T>::operator[](entHandle& entityHandle) {
	// TODO: The pointers (handles) may be to consectutive memory,
	//       in which case there may be a faster way (linear interpolation?)
	//       to find the index (and then just check entities[index])
	int entityIndex = mesh_ptr->indicesOfEntities[entityHandle];
//	cout << "field[" << entityHandle << "] called" << endl;
	return this->operator[](entityIndex);
}

template <class T>
T& Field<T>::operator[](int& entityIndex) {
//	cout << "field[" << entityIndex << "] called" << endl;
	return this->values[entityIndex];
}

template <class T>
T Field<T>::getField(vect3d position, entHandle *entity,
		int interpolationOrder) {
	// TODO: this doesn't work for non-vertex fields
	T field=T();
	bool inElement=false;
	Eigen::Vector4d linearBasisFunctions;
	Eigen::VectorXd errorBases, errorCoefficients;
//	cout << this << " " << *entity << " " << currentElement << " " <<
//			position.transpose() << endl;
	if (*entity==NULL) {
//		cout << this << " " << *entity << " " << currentElement << " " <<
//				position.transpose() << endl;
		// TODO: handle case with no entity hint
		vector<entHandle> adjacentEntities =
				mesh_ptr->getAdjacentEntities(entities[0],iBase_REGION);
		*entity = adjacentEntities[0];
//		cout << this << " " << *entity << " " << currentElement << " " <<
//				position.transpose() << endl << endl;
	} else if (*entity==currentElement) {
//		int dimension=mesh_ptr->getEntityDimension(*entity);
//		if (dimension==iBase_REGION) {
		// TODO: should probably reset currentElement somehow at beginning
		//       of each orbit to prevent data from an old field being used
		//       (probably very unlikely since end tet of one orbit would have
		//       to be start tet of the next)
			inElement=true;
			// TODO: should replace use of linear coeffs with fast checkIfInTet
			linearBasisFunctions =
					mesh_ptr->evaluateLinearBasisFunctions(position, currentElement);
			for (int i=0; i<linearBasisFunctions.rows(); i++)
				if (linearBasisFunctions[i]<0.)
					inElement=false;
//		} else {
//			cout << currentElement << ": " << dimension << " " <<
//					position.transpose() <<
//					" " << position.norm() << endl;
//		}
	}
	bool isElement;
	if (!inElement) {
//		if (currentElement==NULL) {
//			cout << this << " " << *entity << " " << currentElement << " " <<
//					position.transpose() << endl;
//			cout << " " << mesh_ptr->indexOfVertices[*entity] <<
//					" " << mesh_ptr->indexOfFaces[*entity] <<
//					" " << mesh_ptr->indexOfElements[*entity] << endl;
//		}
		vect3d centroid(0.,0.,0.);
		int dimension=mesh_ptr->getEntityDimension(*entity);
		if (dimension!=iBase_REGION) {
			// TODO: find adjacent element
	//		*entity = entities[0];
			isElement = false;
		} else {
			isElement = true;
			vector<vect3d> vVs =
					mesh_ptr->getVertexVectors(*entity);
			centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
		}
		currentElement = *entity;
		bool foundTet=false;
		currentElement = mesh_ptr->findTet(centroid,
				position, currentElement, &foundTet, isElement);
		// TODO: failure to find tet doesn't really mean outside domain yet
		if (!foundTet)
			throw int(OUTSIDE_DOMAIN);
		// TODO: should ensure order of handles and vectors is same
		currentVertices = mesh_ptr->getVertices(currentElement);
		for (int i=0;i<currentVertices.size();i++) {
			currentFields[i] = this->getField(currentVertices[i]);
		}
	}
	linearBasisFunctions =
			mesh_ptr->evaluateLinearBasisFunctions(position, currentElement);
	if (interpolationOrder==2) {
		errorBases =
				mesh_ptr->evaluateQuadraticErrorBases(linearBasisFunctions);
	} else if (interpolationOrder==3) {
		errorBases =
				mesh_ptr->evaluateCubicErrorBases(linearBasisFunctions);
	}
	field *= 0.;
	assert(currentFields.size()==linearBasisFunctions.rows());
	for (int i=0;i<linearBasisFunctions.rows();i++) {
		field += linearBasisFunctions[i]*currentFields[i];
	}
	if (interpolationOrder>1) {
		errorCoefficients = this->getErrorCoefficients(
				currentElement, interpolationOrder);
//		cout << errorBases.rows() << " " << errorCoefficients.rows() <<
//				" " << errorBases.cols() << " " << errorCoefficients.cols() <<
//				endl;
		assert(errorBases.rows()==errorCoefficients.rows());
//		cout << errorCoefficients.transpose() << endl;
//		cout << errorBases.transpose() << endl << endl;
//		cout << field << endl;
//		cout << errorCoefficients.dot(errorBases) << endl;
		field += errorCoefficients.dot(errorBases);
//		for (int i=0;i<errorBases.rows();i++) {
//			field += errorCoefficients[i]*errorBases[i];
//		}
	}
//	cout << field << endl;
	*entity = currentElement;
	return field;
}

template <class T>
void Field<T>::evalFieldAndDeriv(T *fieldValue,
		Eigen::Matrix<T,NDIM,1> *fieldDeriv,
		const vect3d &position, int *entityIndex,
		int interpolationOrder) {
	if (position==currentPosition &&
			interpolationOrder==currentInterpolationOrder) {
		*fieldValue = currentFieldValue;
		*fieldDeriv = currentFieldDerivatives;
	} else {
		bool inElement=false;
//		Eigen::Vector4d &linearBasisFunctions = currentLinearBasisFunctions;
		Eigen::Vector4d linearBasisFunctions(0.,0.,0.,0.);
		int nErrorBases=N_BASES_QUADRATIC;
		if (interpolationOrder==2) {
			nErrorBases = N_BASES_QUADRATIC;
		} else if (interpolationOrder==3) {
			nErrorBases = N_BASES_CUBIC;
		}
		Eigen::VectorXd &errorBases = currentErrorBases;
		Eigen::VectorXd &errorCoefficients = currentErrorCoefficients;
		if (errorBases.rows()!=nErrorBases) {
			errorBases = Eigen::VectorXd(nErrorBases);
			errorCoefficients = Eigen::VectorXd(nErrorBases);
		}
		if (*entityIndex<0) {
			// TODO: handle case with no entity hint
			vector<entHandle> adjacentEntities =
					mesh_ptr->getAdjacentEntities(entities[0],iBase_REGION);
			*entityIndex = mesh_ptr->indicesOfEntities[adjacentEntities[0]];
		} else if (*entityIndex==currentRegionIndex) {
				inElement=true;
				// TODO: should replace use of linear coeffs with fast checkIfInTet
				mesh_ptr->evaluateLinearBasisFunctions(position, currentRegionIndex,
						&linearBasisFunctions);
				for (int i=0; i<linearBasisFunctions.rows(); i++)
					if (linearBasisFunctions[i]<0.-VOLUME_TOLERANCE)
						inElement=false;
		}
		currentRegionIndex = *entityIndex;
		if (!inElement) {
			vect3d centroid(0.,0.,0.);
			vector<vect3d> vVs(NDIM+1);
			mesh_ptr->getVertexVectors(currentRegionIndex, iBase_REGION, &vVs);
			centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
//			currentElement = mesh_ptr->entitiesVectors[iBase_REGION][currentRegionIndex];
			bool foundTet=false;
			currentRegionIndex = mesh_ptr->findTet(centroid,
					position, currentRegionIndex, &foundTet);
			// TODO: failure to find tet doesn't really mean outside domain yet
			if (!foundTet)
				throw int(OUTSIDE_DOMAIN);
			currentVerticesIndices = mesh_ptr->
					adjacentEntitiesVectors[iBase_REGION][currentRegionIndex][iBase_VERTEX];
			for (int i=0;i<currentVerticesIndices.size();i++) {
				currentFields[i] = this->operator[](currentVerticesIndices[i]);
			}
		}
		mesh_ptr->evaluateLinearBasisFunctions(position, currentRegionIndex,
				&linearBasisFunctions);
		if (interpolationOrder==2) {
			mesh_ptr->evaluateQuadraticErrorBases(linearBasisFunctions, &errorBases);
		} else if (interpolationOrder==3) {
			mesh_ptr->evaluateCubicErrorBases(linearBasisFunctions, &errorBases);
		}
		// TODO: vect3d is not initialized to 0, but double is
		*fieldValue = T();
		// TODO: Problem here if *fieldValue=NaN since NaN*0=NaN
		*fieldValue *= 0.;
		assert(currentFields.size()==linearBasisFunctions.rows());
		for (int i=0;i<linearBasisFunctions.rows();i++) {
			*fieldValue += linearBasisFunctions[i]*currentFields[i];
		}
		if (interpolationOrder>1) {
			this->getErrorCoefficients(currentRegionIndex,
					interpolationOrder, &errorCoefficients);
			assert(errorBases.rows()==errorCoefficients.rows());
			*fieldValue += errorCoefficients.dot(errorBases);
		}
		*entityIndex = currentRegionIndex;

		// TODO: remove this check and cast (just for debugging)
		assert(!isnan((double)*fieldValue));

		for (int j=0; j<NDIM; j++) {
			vect3d perturbedPosition = position +
					vect3d::Unit(j)*DELTA_LENGTH;
			mesh_ptr->evaluateLinearBasisFunctions(perturbedPosition,
					currentRegionIndex, &linearBasisFunctions);
			if (interpolationOrder==2) {
				mesh_ptr->evaluateQuadraticErrorBases(linearBasisFunctions,
						&errorBases);
			} else if (interpolationOrder==3) {
				mesh_ptr->evaluateCubicErrorBases(linearBasisFunctions,
						&errorBases);
			}
			// TODO: Too obscure to temp. use fieldDeriv for perturbed field?
			// TODO: vect3d is not initialized to 0, but double is
			(*fieldDeriv)[j] = T();
			// TODO: Problem here if *fieldValue=NaN since NaN*0=NaN
			(*fieldDeriv)[j] *= 0.;
			for (int i=0;i<linearBasisFunctions.rows();i++) {
				(*fieldDeriv)[j] += linearBasisFunctions[i]*currentFields[i];
			}
			if (interpolationOrder>1) {
				(*fieldDeriv)[j] += errorCoefficients.dot(errorBases);
			}
			(*fieldDeriv)[j] -= *fieldValue;
			(*fieldDeriv)[j] /= DELTA_LENGTH;
		}
		currentInterpolationOrder = interpolationOrder;
		currentPosition = position;
		currentFieldValue = *fieldValue;
		currentFieldDerivatives = *fieldDeriv;
	}
}

template <class T>
T Field<T>::getField(entHandle node) {
	T field;
	T *field_ptr=&field;
	int field_alloc = sizeof(T);
	int field_size = sizeof(T);
	int ierr;
	iMesh_getData(mesh_ptr->meshInstance, node, tag, &field_ptr,
			&field_alloc, &field_size, &ierr);
	CHECK("Failure getting field");
	return field;
//	return this->operator[](node);
}

template <class T>
T Field<T>::getFieldFromMeshDB(entHandle entityHandle) {
	T field;
	T *field_ptr=&field;
	int field_alloc = sizeof(T);
	int field_size = sizeof(T);
	int ierr;
	iMesh_getData(mesh_ptr->meshInstance, entityHandle, tag, &field_ptr,
			&field_alloc, &field_size, &ierr);
//	CHECK("Failure getting field");
	if (ierr != iBase_SUCCESS)
		throw FAILURE_GETTING_FIELD;
	return field;
}

template <class T>
T Field<T>::getAverageField(entHandle element) {
	// TODO: this doesn't work for non-vertex fields
	vector<entHandle> vertices =
			mesh_ptr->getAdjacentEntities(element, iBase_VERTEX);
	T average=0;

	for (int i=0; i<vertices.size(); i++) {
		average += this->getField(vertices[i]);
	}
	average /= (double)vertices.size();

	return average;
}

template <class T>
void Field<T>::setField(entHandle node, T field) {
	this->operator[](node) = field;
	int ierr;
	int field_size = sizeof(T);
	iMesh_setData(mesh_ptr->meshInstance, node, tag, &field,
			field_size, &ierr);
}

template <class T>
Eigen::VectorXd Field<T>::getErrorCoefficients(
		entHandle element, int interpolationOrder) {
	Eigen::VectorXd errorCoefficients;
	if (element==currentInterpolationElement && element!=NULL) {
		errorCoefficients = currentErrorCoefficients;
	} else {
		// TODO: do this more elegantly?
		int errorBasesSize;
		if (interpolationOrder==2) {
			errorBasesSize = 6;
		} else if (interpolationOrder==3) {
			errorBasesSize = 16;
//			errorBasesSize = 22;
		} else {
			// TODO: specify exception
			throw;
		}
		// TODO: only dealing with doubles for now
		if (boost::is_same<T,double>::value) {
			vector<entHandle> surroundingVertices =
					mesh_ptr->surroundingVertsMap[element];
			Eigen::MatrixXd evaluatedErrorBases(surroundingVertices.size(),
					errorBasesSize);
			Eigen::VectorXd interpolationError(surroundingVertices.size());
			for (int i=0; i<surroundingVertices.size(); i++) {
				vect3d position =
						mesh_ptr->getCoordinates(surroundingVertices[i]);
				Eigen::Vector4d linearBasisFunctions =
						mesh_ptr->evaluateLinearBasisFunctions(
								position, element);
				Eigen::VectorXd errorBases;
				// TODO: make wrapper to take care of order selection
				if (interpolationOrder==2) {
					errorBases = mesh_ptr->evaluateQuadraticErrorBases(
							linearBasisFunctions);
				} else if (interpolationOrder==3) {
					errorBases = mesh_ptr->evaluateCubicErrorBases(
							linearBasisFunctions);
				} else {
					// TODO: specify exception
					throw;
				}
				// TODO: could probably use block operations here
				assert(errorBases.rows()==evaluatedErrorBases.cols());
				for (int j=0;j<errorBases.rows();j++) {
					evaluatedErrorBases(i,j) = errorBases[j];
				}
				T fieldValue = this->getField(surroundingVertices[i]);
				T interpolatedField;
				interpolatedField *= 0.;
				// TODO: could get fields even if element!=currentElement
				assert(element==currentElement);
				assert(currentFields.size()==linearBasisFunctions.rows());
				// TODO: could make this loop a function
				for (int j=0;j<linearBasisFunctions.rows();j++) {
					interpolatedField +=
							linearBasisFunctions[j]*currentFields[j];
				}
				interpolationError[i] = fieldValue - interpolatedField;
			}
			Eigen::MatrixXd leastSquaresMatrix =
					evaluatedErrorBases.transpose()*evaluatedErrorBases;
			Eigen::VectorXd leastSquaresVector =
					evaluatedErrorBases.transpose()*interpolationError;
			// TODO: could perhapse store intermediate matrix if solve is expensive
			errorCoefficients =
					leastSquaresMatrix.colPivHouseholderQr().solve(
							leastSquaresVector);
		} else {
			// TODO: specify exception
			throw;
		}
		currentErrorCoefficients = errorCoefficients;
		currentInterpolationElement = element;
	}
	return errorCoefficients;
}

template <class T>
Eigen::VectorXd Field<T>::getErrorCoefficients(
		int elementIndex, int interpolationOrder) {
	// TODO: do this more elegantly?
	int errorBasesSize;
	if (interpolationOrder==2) {
		errorBasesSize = 6;
	} else if (interpolationOrder==3) {
		errorBasesSize = 16;
//			errorBasesSize = 22;
	} else {
		// TODO: specify exception
		throw;
	}
	Eigen::VectorXd errorCoefficients(errorBasesSize);
	this->getErrorCoefficients(elementIndex, interpolationOrder,
			&errorCoefficients);
	return errorCoefficients;
}

template <class T>
void Field<T>::getErrorCoefficients(
		int elementIndex, int interpolationOrder,
		Eigen::VectorXd *errorCoefficients) {
	if (elementIndex==currentInterpolationElementIndex &&
			elementIndex>=0) {
		*errorCoefficients = currentErrorCoefficients;
	} else {
		// TODO: do this more elegantly?
		int errorBasesSize;
		if (interpolationOrder==2) {
			errorBasesSize = 6;
		} else if (interpolationOrder==3) {
			errorBasesSize = 16;
//			errorBasesSize = 22;
		} else {
			// TODO: specify exception
			throw;
		}
		assert(errorCoefficients->rows()==errorBasesSize);
		// TODO: only dealing with doubles for now
		if (boost::is_same<T,double>::value) {
			vector<int> &surroundingVertices =
					mesh_ptr->verticesSurroundingRegions[elementIndex];
			Eigen::MatrixXd evaluatedErrorBases(surroundingVertices.size(),
					errorBasesSize);
			Eigen::VectorXd interpolationError(surroundingVertices.size());
			for (int i=0; i<surroundingVertices.size(); i++) {
				vect3d &position = mesh_ptr->
						verticesPositions[surroundingVertices[i]];
				Eigen::Vector4d linearBasisFunctions =
						mesh_ptr->evaluateLinearBasisFunctions(
								position, elementIndex);
				Eigen::VectorXd errorBases;
				// TODO: make wrapper to take care of order selection
				if (interpolationOrder==2) {
					errorBases = mesh_ptr->evaluateQuadraticErrorBases(
							linearBasisFunctions);
				} else if (interpolationOrder==3) {
					errorBases = mesh_ptr->evaluateCubicErrorBases(
							linearBasisFunctions);
				} else {
					// TODO: specify exception
					throw;
				}
				// TODO: could probably use block operations here
				assert(errorBases.rows()==evaluatedErrorBases.cols());
				for (int j=0;j<errorBases.rows();j++) {
					evaluatedErrorBases(i,j) = errorBases[j];
				}
				T fieldValue = this->operator[](surroundingVertices[i]);
				T interpolatedField;
				interpolatedField *= 0.;
				// TODO: could get fields even if element!=currentElement
				assert(elementIndex==currentRegionIndex);
				assert(currentFields.size()==linearBasisFunctions.rows());
				// TODO: could make this loop a function
				for (int j=0;j<linearBasisFunctions.rows();j++) {
					interpolatedField +=
							linearBasisFunctions[j]*currentFields[j];
				}
				interpolationError[i] = fieldValue - interpolatedField;
			}
			Eigen::MatrixXd leastSquaresMatrix =
					evaluatedErrorBases.transpose()*evaluatedErrorBases;
			Eigen::VectorXd leastSquaresVector =
					evaluatedErrorBases.transpose()*interpolationError;
			// TODO: could perhapse store intermediate matrix if solve is expensive
			*errorCoefficients =
					leastSquaresMatrix.colPivHouseholderQr().solve(
							leastSquaresVector);
		} else {
			// TODO: specify exception
			throw;
		}
		currentErrorCoefficients = *errorCoefficients;
		currentInterpolationElementIndex = elementIndex;
	}
}

#endif /* FIELD_H_ */


