/*
 * Field.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef FIELD_H_
#define FIELD_H_

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "typesAndDefinitions.h"
#include "variables.h"
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
#include "SpatialDependence.h"
#include "DistributionFunction.h"
#include "Mesh.h"
#include "Solver.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

class DensityField;
class PotentialField;
class ShortestEdgeField;

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, string inputName,
			int entityDimension, int numberOfComponents=1);
	// TODO: Does an empty destructor cause a memory leak?
	virtual ~Field() {}

	void copyValues(Field<T>& fieldToCopy, int targetSlot=0);

	T& operator[](entHandle& entity);
	T& operator[](int& entityIndex);
	// TODO: operator+= should probably return Field<T>&
	void operator+=(T& valueToAdd);

	T getField(vect3d position, entHandle *entity=NULL,
			int interpolationOrder=INTERPOLATIONORDER);
	T getField(entHandle node);
	void getField(entHandle entityHandle, T *field_ptr);
	vector<T> getFieldVector(entHandle entityHandle);
	void getFieldVector(entHandle entityHandle, vector<T> *fieldVector_ptr);
	T getFieldFromMeshDB(entHandle entityHandle);
	void getFieldFromMeshDB(entHandle entityHandle, T *field_ptr,
			int requestedNumberOfComponents=1);
	T getAverageField(entHandle element);
	void evalField(T *fieldValue,
			const vect3d &position, int *entityIndex=-1,
			int interpolationOrder=INTERPOLATIONORDER);
	void evalFieldAndDeriv(T *fieldValue,
			Eigen::Matrix<T,NDIM,1> *fieldDeriv,
			const vect3d &position, int *entityIndex=-1,
			int interpolationOrder=INTERPOLATIONORDER);
	void setField(entHandle node, T field, int targetSlot=0);
	void setField(entHandle node, T *field_ptr);
	void setFieldVector(entHandle entityHandle, vector<T>& fieldVector);
	void setVtkField(entHandle node, T *field_ptr);
	Eigen::VectorXd getErrorCoefficients(entHandle element,
			int interpolationOrder);
	Eigen::VectorXd getErrorCoefficients(int elementIndex,
			int interpolationOrder);
	void getErrorCoefficients(int elementIndex,
			int interpolationOrder,
			Eigen::VectorXd *errorCoefficients);
	void updateTagHandle();

	void computeWeightedAverage(Field<T> *averageField_ptr,
			vector<double> weights);

	void print();

	Mesh *mesh_ptr;
//	// TODO: can this be done better?
//	typedef boost::conditional<boost::is_same<T,int>::value,
//								vtkIntArray,vtkDoubleArray>::type vtkFieldType;
//	vtkSmartPointer<vtkFieldType> vtkField_ptr;
	vtkSmartPointer<vtkDataArray> vtkField_ptr;
	bool storeVtkField;
	int numberOfComponents;
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
	ElectricField(Mesh *inputMesh_ptr, string inputName,
			CodeField vertexType, double debyeLength, bool doLuDecomposition);
	virtual ~ElectricField() {}

	void calcField(PotentialField potentialField);
	void calcField(PotentialField *potentialField_ptr, CodeField vertexType,
			DensityField ionDensity, double debyeLength);
	void calcField_Gatsonis(PotentialField potentialField);

	Solver solver;
};

class CodeField : public Field<int> {
public:
	CodeField(Mesh *inputMesh_ptr, string inputName,
			int elementType);
	virtual ~CodeField() {}

	void calcField(Field<int> faceTypeField);
};

class DerivativeField : public Field<double> {
public:
	DerivativeField(Mesh *inputMesh_ptr, string inputName,
			int elementType);
	virtual ~DerivativeField() {}

	void calcField(Field<double>& numeratorValue, Field<double>& numeratorReference,
			Field<double>& denominatorValue, Field<double>& denominatorReference);
};

class DensityDerivativeField : public DerivativeField {
public:
	DensityDerivativeField(Mesh *inputMesh_ptr, string inputName,
			int elementType);
	virtual ~DensityDerivativeField() {}

	void calcField(Field<double>& numeratorValue, Field<double>& numeratorReference,
			Field<double>& denominatorValue, Field<double>& denominatorReference);
};

class DensityField : public Field<double> {
	// Charge density
public:
	DensityField(Mesh *inputMesh_ptr, string inputName,
			Field<vect3d> *inputAverageVelocity_ptr=NULL,
			Field<double> *inputTemperature_ptr=NULL,
			int numberOfComponents=1);
	virtual ~DensityField() {}

	void calcField();
	void calcField(CodeField& vertexType,
			PotentialField& potential, DensityField& referenceDensity,
			SpatialDependence& referenceTemperature, double charge);
	void calcField(DensityField ionDensity, DensityField electronDensity);
	void calcField(ElectricField& electricField,
			PotentialField *potentialField_ptr, DensityField& referenceDensity,
			Field<int>& faceType, CodeField& vertexType,
			ShortestEdgeField& shortestEdge, double charge,
			double potentialPerturbation, bool doAllNodes=false,
			double unconvergednessThreshold=0.03, FILE *outFile=NULL);
	// TODO: make reference density its own class?
	void calcField(CodeField& vertexType,
			DistributionFunction& distributionFunction, double charge);
	void poissonCubeTest(double debyeLength);
	void calculateDensity(int node, ElectricField& electricField,
			PotentialField& potentialField, DensityField& referenceDensity,
			Field<int>& faceType, CodeField& vertexType,
			ShortestEdgeField& shortestEdge, double charge,
			double potentialPerturbation, double *density, double *error,
			vect3d *averageVelocity, vect3d *averageVelocityError,
			double *temperature, double *temperatureError);
#ifdef HAVE_MPI
	void requestDensityFromSlaves(ElectricField& electricField,
			PotentialField *potentialField_ptr, vector<int>& sortedNodes,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation, FILE *outFile);
	MPI::Status receiveDensity(double *potential,FILE *outFile);
	void processDensityRequests(ElectricField& electricField,
			PotentialField *potentialField_ptr, DensityField& referenceDensity,
			Field<int> faceType, CodeField vertexType,
			ShortestEdgeField shortestEdge, double charge,
			double potentialPerturbation);
#endif

	void setDistributionFunction(DistributionFunction& distributionFunction);

	Field<vect3d> *averageVelocity_ptr;
	Field<double> *temperature_ptr;
	DistributionFunction *distributionFunction_ptr;
};

class PotentialField : public Field<double> {
public:
	PotentialField(Mesh *inputMesh_ptr, string inputName, double boundaryPotential=0.,
			double sheathPotential=-0.5, bool fixSheathPotential=false,
			int numberOfComponents=1);
	PotentialField(PotentialField potential, string inputName);
	virtual ~PotentialField() {}

	void calcField(CodeField vertexType, double debyeLength,
			double boundaryPotential, double surfacePotential,
			double sheathPotential);
	void calcField(DensityField ionDensity, DensityField electronDensity,
			CodeField vertexType, FILE *outFile);
	void calcField(DensityField& ionDensity, Field<vect3d>& ionVelocity,
			CodeField& vertexType, FILE *outFile, double boundaryPotential,
			double sheathPotential, bool fixSheathPotential);
	void calcFieldAtNode(entHandle entity, double ionDensity, int vertexType,
			double boundaryPotential, double sheathPotential, bool fixSheathPotential);
	void calcField(DensityField& ionDensity, DerivativeField& ionDensityDerivative,
			CodeField& vertexType, FILE *outFile, double boundaryPotential,
			double sheathPotential, bool fixSheathPotential);
	void calcField(DensityField ionDensity,
			DensityField ionDensityPP, DensityField ionDensityNP,
			DensityField electronDensity,
			DensityField electronDensityPP, DensityField electronDensityNP,
			CodeField vertexType,
			double positivePotentialPerturbation,
			double negativePotentialPerturbation,
			FILE *outFile);

	// TODO: not very transparent to use inverted defaults to disable min/max
	void computePerturbedPotentials(double negativePerturbation,
			double positivePertubation, double minPotential=10., double maxPotential=-10.);

	void setReferenceElectronDensity(DensityField& referenceElectronDensity);
	void setReferenceElectronTemperature(SpatialDependence& referenceElectronTemperature);

	double boundaryPotential;
	double sheathPotential;
	bool fixSheathPotential;

	DensityField *referenceElectronDensity_ptr;
	SpatialDependence *referenceElectronTemperature_ptr;
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
		int entityDimension, int numberOfComponents) : entities(
				inputMesh_ptr->entitiesVectors[entityDimension]),
				numberOfComponents(numberOfComponents) {
	mesh_ptr = inputMesh_ptr;
	name = inputName;
	int ierr;
	int size, type;
//	if (is_same<T,double>::value) {
	if (boost::is_same<T,double>::value) {
		size = numberOfComponents;
		type = iBase_DOUBLE;
		vtkField_ptr = vtkSmartPointer<vtkDoubleArray>::New();
//	} else if (is_same<T,int>::value) {
	} else if (boost::is_same<T,int>::value) {
		size = numberOfComponents;
		type = iBase_INTEGER;
		vtkField_ptr = vtkSmartPointer<vtkIntArray>::New();
	} else if (boost::is_same<T,vect3d>::value) {
		size = NDIM*numberOfComponents;
		type = iBase_DOUBLE;
		vtkField_ptr = vtkSmartPointer<vtkDoubleArray>::New();
	} else {
		size = numberOfComponents*(int)sizeof(T);
		type = iBase_BYTES;
//		// TODO: figure out why this doesn't work
//		vtkField_ptr = vtkSmartPointer<vtkDataArray>::New();
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

	// TODO: don't assume entityDimension==0
//	throw string("nonzero entity dimension in SurfaceField not implemented yet");
	if (entityDimension==0) {
		storeVtkField = true;
		if (entities.size() != mesh_ptr->vtkMesh_ptr->GetNumberOfPoints())
			throw string("vtkMesh and mesh don't have the same number of elements");
		// TODO: add try/catch?
		vtkField_ptr->SetName(name.c_str());
		vtkField_ptr->SetNumberOfComponents(size);
		vtkField_ptr->SetNumberOfTuples(entities.size());
		// TODO: verify that this doens't make a copy, i.e. can use field_ptr to modify
		mesh_ptr->vtkMesh_ptr->GetPointData()->AddArray(vtkField_ptr);
		for (int i=0; i<entities.size(); i++) {
			this->setVtkField(entities[i],&(values[i]));
		}
	} else {
		storeVtkField = false;
	}

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
void Field<T>::copyValues(Field<T>& fieldToCopy, int targetSlot) {
	for (int i=0; i<entities.size(); i++) {
		this->setField(entities[i],fieldToCopy.getField(entities[i]),targetSlot);
	}
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
void Field<T>::operator+=(T& valueToAdd) {
	for (int i=0; i<entities.size(); i++) {
		this->setField(entities[i],this->getField(entities[i])+valueToAdd);
	}
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
	entHandle entityIfNull=NULL;
	// TODO: figure out why *entity=NULL initialization can give entity=NULL...opt?
	if (entity==NULL)
		entity = &entityIfNull;
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
	if (currentFields.size()!=linearBasisFunctions.rows())
		throw string("problem in getField");
	for (int i=0;i<linearBasisFunctions.rows();i++) {
		field += linearBasisFunctions[i]*currentFields[i];
	}
	if (interpolationOrder>1) {
		errorCoefficients = this->getErrorCoefficients(
				currentElement, interpolationOrder);
//		cout << errorBases.rows() << " " << errorCoefficients.rows() <<
//				" " << errorBases.cols() << " " << errorCoefficients.cols() <<
//				endl;
		if (errorBases.rows()!=errorCoefficients.rows())
			throw string("problem in getField()");
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
void Field<T>::evalField(T *fieldValue,
		const vect3d &position, int *entityIndex,
		int interpolationOrder) {
	// TODO: should evalField just be getField?
	// TODO: too much code-duplication with evalFieldAndDeriv
	// Can only do 1st order interpolation of non-doubles
	if (interpolationOrder!=1)
		throw string("can only do first order interpolation of non-doubles");
	if (position==currentPosition &&
			interpolationOrder==currentInterpolationOrder) {
		*fieldValue = currentFieldValue;
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
		if (currentFields.size()!=linearBasisFunctions.rows())
			throw string("problem in evalField");
		for (int i=0;i<linearBasisFunctions.rows();i++) {
			*fieldValue += linearBasisFunctions[i]*currentFields[i];
		}
		// TODO: should be able to call evalField for doubles with higher order
//		if (interpolationOrder>1) {
//			this->getErrorCoefficients(currentRegionIndex,
//					interpolationOrder, &errorCoefficients);
//			assert(errorBases.rows()==errorCoefficients.rows());
//			*fieldValue += errorCoefficients.dot(errorBases);
//		}
		*entityIndex = currentRegionIndex;

		currentInterpolationOrder = interpolationOrder;
		currentPosition = position;
		currentFieldValue = *fieldValue;
	}
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
		if (currentFields.size()!=linearBasisFunctions.rows())
			throw string("problem in evalFieldAndDeriv");
		for (int i=0;i<linearBasisFunctions.rows();i++) {
			*fieldValue += linearBasisFunctions[i]*currentFields[i];
		}
		if (interpolationOrder>1) {
			this->getErrorCoefficients(currentRegionIndex,
					interpolationOrder, &errorCoefficients);
			if (errorBases.rows()!=errorCoefficients.rows())
				throw string("problem in evelFieldAndDeriv()");
			*fieldValue += errorCoefficients.dot(errorBases);
		}
		*entityIndex = currentRegionIndex;

		// TODO: remove this check and cast (just for debugging)
		if (isnan((double)*fieldValue))
			throw string("NaN in evalFieldAndDeriv");

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
	try {
		return this->getFieldFromMeshDB(node);
	} catch (int FAILURE_GETTING_FIELD) {
		throw string("Failure getting field.");
	}
//	return this->operator[](node);
}

template <class T>
void Field<T>::getField(entHandle entityHandle, T *field_ptr) {
	this->getFieldFromMeshDB(entityHandle, field_ptr, numberOfComponents);
}

template <class T>
T Field<T>::getFieldFromMeshDB(entHandle entityHandle) {
	T field;
	this->getFieldFromMeshDB(entityHandle, &field, 1);
	return field;
}

template <class T>
vector<T> Field<T>::getFieldVector(entHandle entityHandle) {
	vector<T> fieldVector(numberOfComponents);
	this->getFieldVector(entityHandle, &fieldVector);
	return fieldVector;
}

template <class T>
void Field<T>::getFieldVector(entHandle entityHandle,
		vector<T> *fieldVector_ptr) {
	boost::scoped_array<T> buffer(new T[numberOfComponents]);
	T *field_ptr = buffer.get();
	this->getFieldFromMeshDB(entityHandle, field_ptr, numberOfComponents);
	if (fieldVector_ptr->size()!=numberOfComponents)
		throw string("fieldVector_ptr->size() != numberOfComponents in setFieldVector");
	for (int i=0; i<fieldVector_ptr->size(); i++) {
		fieldVector_ptr->operator[](i) = field_ptr[i];
	}
}

template <class T>
void Field<T>::getFieldFromMeshDB(entHandle entityHandle, T *fieldOutput_ptr,
		int requestedNumberOfComponents) {
	boost::scoped_array<T> buffer(new T[numberOfComponents]);
	T *field_ptr = buffer.get();
	int field_alloc;
	int field_size;
	if (boost::is_same<T,vect3d>::value) {
		field_alloc = numberOfComponents*3*sizeof(double);
	} else {
		field_alloc = numberOfComponents*sizeof(T);
	}
	field_size = field_alloc;
	int ierr;
	iMesh_getData(mesh_ptr->meshInstance, entityHandle, tag, &field_ptr,
			&field_alloc, &field_size, &ierr);
	if (ierr != iBase_SUCCESS) {
		throw int(FAILURE_GETTING_FIELD);
	}
	for (int i=0; i<min(requestedNumberOfComponents,numberOfComponents); i++) {
		fieldOutput_ptr[i] = field_ptr[i];
	}
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
void Field<T>::setField(entHandle node, T field, int targetSlot) {
	if (numberOfComponents==1) {
		this->setField(node,&field);
	} else {
		boost::scoped_array<T> buffer(new T[numberOfComponents]);
		T *field_ptr = buffer.get();
		this->getField(node,field_ptr);
		if (targetSlot>=numberOfComponents)
			throw string("targetSlot >= numberOfComponents in setField");
		field_ptr[targetSlot] = field;
		this->setField(node,field_ptr);
	}
}

template <class T>
void Field<T>::setField(entHandle node, T *field_ptr) {
	this->setVtkField(node,field_ptr);
	// TODO: cache more than first component?
	this->operator[](node) = *field_ptr;
	int ierr;
	void *components_ptr;
	int field_size;
	boost::scoped_array<vect3d> vect3dBuffer(new vect3d[NDIM*numberOfComponents]);
	if (boost::is_same<T,vect3d>::value) {
		for (int k=0; k<numberOfComponents; k++) {
			for (int i=0; i<NDIM; i++)
				vect3dBuffer[k][i] = (((vect3d*)field_ptr)[k])[i];
		}
		components_ptr = (T*)(vect3dBuffer.get());
//		// TODO: figure out why using .data() breaks in optimization
//		//       (adding print .transpose() prevents optimization and works)
//		vect3d vector(field);
//		field_ptr = vector.data();
		field_size = NDIM*numberOfComponents*sizeof(double);
	} else {
		components_ptr = field_ptr;
		field_size = numberOfComponents*sizeof(T);
	}
	iMesh_setData(mesh_ptr->meshInstance, node, tag, components_ptr,
			field_size, &ierr);
	if (ierr != iBase_SUCCESS)
		throw string("Failure setting field.");
}

template <class T>
void Field<T>::setFieldVector(entHandle entityHandle, vector<T>& fieldVector) {
	boost::scoped_array<T> buffer(new T[numberOfComponents]);
	T *field_ptr = buffer.get();
	if (fieldVector.size()!=numberOfComponents)
		throw string("fieldVector.size() != numberOfComponents in setFieldVector");
	for (int i=0; i<fieldVector.size(); i++) {
		field_ptr[i] = fieldVector[i];
	}
	this->setField(entityHandle,field_ptr);
}

template <class T>
void Field<T>::setVtkField(entHandle node, T *field_ptr) {
	void *components_ptr;
	boost::scoped_array<vect3d> vect3dBuffer(new vect3d[NDIM*numberOfComponents]);
	if (boost::is_same<T,vect3d>::value) {
		for (int k=0; k<numberOfComponents; k++) {
			for (int i=0; i<NDIM; i++)
				vect3dBuffer[k][i] = (((vect3d*)field_ptr)[k])[i];
		}
		components_ptr = (T*)(vect3dBuffer.get());
	} else {
		components_ptr = field_ptr;
	}
	if (storeVtkField) {
		vtkIdType vtkNode = mesh_ptr->iMeshToVtk[node];
		// TODO: other cases than int or double?
		if (boost::is_same<T,int>::value) {
			vtkIntArray::SafeDownCast(vtkField_ptr)->SetTupleValue(vtkNode, (int*)components_ptr);
		} else {
			vtkDoubleArray::SafeDownCast(vtkField_ptr)->SetTupleValue(vtkNode, (double*)components_ptr);
		}
	}
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
					throw string("problem in getErrorCoefficients");
				}
				// TODO: could probably use block operations here
				if (errorBases.rows()!=evaluatedErrorBases.cols())
					throw string("problem in getErrorCoefficients");
				for (int j=0;j<errorBases.rows();j++) {
					evaluatedErrorBases(i,j) = errorBases[j];
				}
				T fieldValue = this->getField(surroundingVertices[i]);
				T interpolatedField;
				interpolatedField *= 0.;
				// TODO: could get fields even if element!=currentElement
				if (element!=currentElement)
					throw string("problem in getErrorCoefficients");
				if (currentFields.size()!=linearBasisFunctions.rows())
					throw string("problem in getErrorCoefficients");
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
		if (errorCoefficients->rows()!=errorBasesSize)
			throw string("problem in getErrorCoefficients");
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
				if (errorBases.rows()!=evaluatedErrorBases.cols())
					throw string("problem in getErrorCoefficients");
				for (int j=0;j<errorBases.rows();j++) {
					evaluatedErrorBases(i,j) = errorBases[j];
				}
				T fieldValue = this->operator[](surroundingVertices[i]);
				T interpolatedField;
				interpolatedField *= 0.;
				// TODO: could get fields even if element!=currentElement
				if (elementIndex!=currentRegionIndex)
					throw string("problem in getErrorCoefficients");
				if (currentFields.size()!=linearBasisFunctions.rows())
					throw string("problem in getErrorCoefficients");
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

template <class T>
void Field<T>::updateTagHandle() {
	tag = mesh_ptr->getTagHandle(name);
}

template <class T>
void Field<T>::computeWeightedAverage(Field<T> *averageField_ptr,
		vector<double> weights) {
	vector<T> fields(numberOfComponents);
	for (int i=0; i<entities.size(); i++) {
		fields = this->getFieldVector(entities[i]);
		// TODO: improve initialization?
		T weightedAverage=0;
		weightedAverage *= 0.;
		if (weights.size()!=numberOfComponents)
			throw string("weights.size() != numberOfComponents in computeWeightedAverage");
		for (int j=0; j<weights.size(); j++) {
			weightedAverage += weights[j]*fields[j];
		}
		averageField_ptr->setField(entities[i],weightedAverage);
	}
}

template <class T>
void Field<T>::print() {
	for (int i=0; i<entities.size(); i++)
		cout << this->getField(entities[i]) << "..." << endl;
}

#endif /* FIELD_H_ */


