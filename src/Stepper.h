/*
 * Stepper.h
 *
 *  Created on: Jun 4, 2012
 *      Author: chaako
 */

#ifndef STEPPER_H_
#define STEPPER_H_

#include "typesAndDefinitions.h"
//#include <boost/array.hpp>
//#include "Eigen/Dense"

template<int N> class VelocityVerletStepper{
public:
	typedef boost::array<Eigen::Matrix<double,N,1>, 2> state_type;

    static unsigned short order(){return 2;} // Assuming a(x,v)=a(x)

    template<class OdeSystem>
//    void do_step(OdeSystem odeSystem, state_type &x, double t,
//    		state_type &xOut, double dt) {
	void do_step(OdeSystem odeSystem, state_type &x, double t, double dt) {
//	void do_step(O deSystem &odeSystem, state_type &x, double t, double dt) {
    	state_type dxdt;
//    	// Store input x so that can reset state if exception thrown
//    	state_type xIn;
//    	xIn[0] = x[0];
//    	xIn[1] = x[1];
//		try {
//			// Advance velocity half timestep
//			odeSystem(x, dxdt);
//			x[1] += dt/2. * dxdt[1];
//			// Advance position
//			odeSystem(x, dxdt);
//			x[0] += dt * dxdt[0];
//			// Advance velocity half timestep
//			odeSystem(x, dxdt);
//			x[1] += dt/2. * dxdt[1];

			// Advance position half timestep
//			boost::unwrap_ref(odeSystem)(x, dxdt);
//    		odeSystem(x, dxdt);
//			x[0] += dt/2. * dxdt[0];
			x[0] += dt/2. * x[1];
			// Advance velocity
			boost::unwrap_ref(odeSystem)(x, dxdt);
//    		odeSystem(x, dxdt);
			x[1] += dt * dxdt[1];
			// Advance position half timestep
//			boost::unwrap_ref(odeSystem)(x, dxdt);
//    		odeSystem(x, dxdt);
//			x[0] += dt/2. * dxdt[0];
			x[0] += dt/2. * x[1];

//			xOut[0] = x[0];
//			xOut[1] = x[1];
//			// Advance position half timestep
//    		boost::unwrap_ref(odeSystem)(xOut, dxdt);
//			xOut[0] += dt/2. * dxdt[0];
//			// Advance velocity
//    		boost::unwrap_ref(odeSystem)(xOut, dxdt);
//			xOut[1] += dt * dxdt[1];
//			// Advance position half timestep
//    		boost::unwrap_ref(odeSystem)(xOut, dxdt);
//			xOut[0] += dt/2. * dxdt[0];

//			// Advance position
//			odeSystem(x, dxdt);
//			for (int i=0; i<x.size()/2; ++i )
//				x[i] += dt * dxdt[i];
//			// Advance velocity half timestep
//			odeSystem(x, dxdt);
//			for (int i=x.size()/2; i<x.size(); ++i )
//				x[i] += dt/2. * dxdt[i];
//		} catch (int signal) {
//			x[0] = xIn[0];
//			x[1] = xIn[1];
//			throw signal;
//		}
    }
};

class CyclotronicStepper{
public:
	typedef boost::array<vect3d, 2> state_type;

    static unsigned short order(){return 2;}

    template<class OdeSystem>
	void do_step(OdeSystem odeSystem, state_type &x, double t, double dt) {
		state_type dxdt;
		// TODO: should get B from odeSystem, though Cyclotronic integrator
		//       is meant for uniform B so external is fine in that sense
		vect3d unitB = B/B.norm();
		double charge = boost::unwrap_ref(odeSystem).charge;
		// TODO: figure out how to distinguish electron and ion masses
		double mass = 1.;
		// TODO: replace this opaque hack to account for mass difference
		if (charge<0)
			charge *= sqrt(1836.); // sqrt to also account for v_th difference
		// TODO: think through rotation signs since backwards in time (appears correct)
		double rotationAngle = -dt*charge*B.norm()/mass;
		// TODO: determine how Eigen rotation handles large angles
//		int wholePeriods = rotationAngle/(2.*M_PI);
//		rotationAngle -= 2.*M_PI*wholePeriods;
		vect3d larmorVector, rotatedLarmorVector;
//		cout << B.transpose() << ", m " << mass << ", c " << charge << endl;
//		cout << x[0].transpose() << ", v " << x[1].transpose() << endl;

		// Advance position half timestep (drift)
		x[0] += dt/2. * x[1].dot(unitB) * unitB;
		// TODO: handle EXB differently?
		x[0] += VEXB*dt/2.;
		larmorVector = -mass*x[1].cross(unitB)/charge/B.norm();
		rotatedLarmorVector =
				Eigen::AngleAxisd(rotationAngle/2., unitB) * larmorVector;
		x[0] += rotatedLarmorVector - larmorVector;
		x[1] = Eigen::AngleAxisd(-rotationAngle/2., unitB) * x[1];

		// Advance velocity (kick)
		boost::unwrap_ref(odeSystem)(x, dxdt);
		x[1] += dt * dxdt[1];

		// Advance position half timestep (drift)
		x[0] += dt/2. * x[1].dot(unitB) * unitB;
		// TODO: handle EXB differently?
		x[0] += VEXB*dt/2.;
		larmorVector = -mass*x[1].cross(unitB)/charge/B.norm();
		rotatedLarmorVector =
				Eigen::AngleAxisd(rotationAngle/2., unitB) * larmorVector;
		x[0] += rotatedLarmorVector - larmorVector;
		x[1] = Eigen::AngleAxisd(-rotationAngle/2., unitB) * x[1];
    }
};

class TaylorStepper{
public:
	typedef boost::array<vect3d, 2> state_type;

    static unsigned short order(){return 2;}

    template<class OdeSystem>
	void do_step(OdeSystem odeSystem, state_type &x, double t, double dt) {
//    	cout << x[0].transpose() << endl;
//    	cout << x[1].transpose() << endl;
//    	x[0] += vect3d(1.,2.,3.);
//    	x[1][0] += 1.;
//    	x[1][1] += 2.;
//    	x[1][2] += 3.;
//    	cout << x[0].transpose() << endl;
//    	cout << x[1].transpose() << endl;
//    	exit(0);

		state_type dxdt;
		// TODO: should get B from odeSystem, though Taylor integrator
		//       is meant for uniform B so external is fine in that sense
		vect3d unitB = B/B.norm();
		double charge = boost::unwrap_ref(odeSystem).charge;
		// TODO: figure out how to distinguish electron and ion masses
		double mass = 1.;
		// TODO: replace this opaque hack to account for mass difference
		if (charge<0)
			charge *= sqrt(1836.); // sqrt to also account for v_th difference
		double omega = -charge*B.norm()/mass; // minus sign since integrating backwards

//		// TODO: remove this
//		double rotationAngle = dt*charge*B.norm()/mass;
//		vect3d larmorVector, rotatedLarmorVector;
////		x[0] += dt * x[1].dot(unitB) * unitB;
//		larmorVector = -mass*x[1].cross(unitB)/charge/B.norm();
//		rotatedLarmorVector =
//				Eigen::AngleAxisd(rotationAngle, unitB) * larmorVector;
////		x[0] += rotatedLarmorVector - larmorVector;
////		x[1] = Eigen::AngleAxisd(rotationAngle, unitB) * x[1];

//		// TODO: remove this
//		double cc = cos(omega*dt) - 1.;
//		x[0][0] += (x[1][0]*sin(omega*dt)-x[1][1]*cc)/omega;
//		x[0][1] += (x[1][1]*sin(omega*dt)+x[1][0]*cc)/omega;
//		x[0][2] += x[1][2]*dt;
//		double vx = x[1][0];
//		double vy = x[1][1];
//		x[1][0] = vx*cos(omega*dt)+vy*sin(omega*dt);
//		x[1][1] = vy*cos(omega*dt)-vx*sin(omega*dt);
////		double vx = x[1][0]*cos(omega*dt)+x[1][1]*sin(omega*dt);
////		double vy = x[1][1]*cos(omega*dt)-x[1][0]*sin(omega*dt);
////		x[1] = vect3d(vx,vy,x[1][2]);

		// TODO: generalize to any unitB
		assert(unitB[0]==0. && unitB[1]==0);
		vect3d pos = x[0];
		vect3d vel = x[1];
		boost::unwrap_ref(odeSystem)(x, dxdt);
		vect3d accel = dxdt[1];
		vect3d velExB = VEXB;

		double s = sin(omega*dt) - omega*dt;
		double sm = sin(-omega*dt) - (-omega)*dt;
		double c = cos(omega*dt) - 1.;
		double cm = cos(-omega*dt) - 1.;

		// Advance position
		pos[0] += 1./omega*(vel[0]*sin(omega*dt)-vel[1]*c)
				+1./omega/omega*(-accel[0]*c-accel[1]*s);
		pos[1] += 1./(-omega)*(vel[1]*sin(-omega*dt)-vel[0]*cm)
				+1./omega/omega*(-accel[1]*cm-accel[0]*sm);
		pos[2] += dt * vel[2] + 1./2.*dt*dt*accel[2];
		// TODO: handle EXB differently?
		pos += velExB*dt;

		// Advance velocity
		// TODO: generalize to any unitB
		x[0] = pos;
		x[1] = vel;
		boost::unwrap_ref(odeSystem)(x, dxdt);
		pos = x[0];
		vel = x[1];
		vect3d origAccel = accel;
		accel = dxdt[1];
		vect3d dadt = (accel-origAccel)/dt;
		vect3d updatedVel = vel;
		updatedVel[0] = vel[0]*cos(omega*dt)+vel[1]*sin(omega*dt)
				+1./omega*(origAccel[0]*sin(omega*dt)-origAccel[1]*c)
				+1./omega/omega*(-dadt[0]*c-dadt[1]*s);
		updatedVel[1] = vel[1]*cos(-omega*dt)+vel[0]*sin(-omega*dt)
				+1./(-omega)*(origAccel[1]*sin(-omega*dt)-origAccel[0]*cm)
				+1./omega/omega*(-dadt[1]*cm-dadt[0]*sm);
		updatedVel[2] += 1./2.*dt*(origAccel[2]+accel[2]);

		// TODO: generalize to any unitB
		x[0] = pos;
		x[1] = updatedVel;
    }
};

#endif /* STEPPER_H_ */
