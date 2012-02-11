/*
 * Field.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef FIELD_H_
#define FIELD_H_

class Mesh;

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	~Field() {}

//	virtual T getField(Eigen::Vector3d position);
//	virtual T getField(iBase_EntityHandle node);

//	virtual void calcField();

	Mesh *mesh_ptr;
	std::string name;
	iBase_TagHandle tag;

};

class ScalarField : public Field<double> {
public:
	ScalarField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	~ScalarField() {}

//	double getField(Eigen::Vector3d position);
//	double getField(iBase_EntityHandle node);

//	virtual void calcField();

};

class VectorField : public Field<Eigen::Vector3d> {
public:
	VectorField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	~VectorField() {}

//	Eigen::Vector3d getField(Eigen::Vector3d position);
//	Eigen::Vector3d getField(iBase_EntityHandle node);

//	virtual void calcField();

};

class ElectricField : public VectorField {
public:
	ElectricField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	~ElectricField() {}

//	virtual void calcField();

};

class DensityField : public ScalarField {
public:
	DensityField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	~DensityField() {}

//	virtual void calcField();

};

class PotentialField : public ScalarField {
public:
	PotentialField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
//	virtual ~Field();

//	virtual void calcField();

};

#endif /* FIELD_H_ */


