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
		double rotationAngle = dt*charge*B.norm()/mass;
		// TODO: determine how Eigen rotation handles large angles
//		int wholePeriods = rotationAngle/(2.*M_PI);
//		rotationAngle -= 2.*M_PI*wholePeriods;
		vect3d larmorVector, rotatedLarmorVector;
//		cout << B.transpose() << ", m " << mass << ", c " << charge << endl;
//		cout << x[0].transpose() << ", v " << x[1].transpose() << endl;

		// Advance position half timestep (drift)
		x[0] += dt/2. * x[1].dot(unitB) * unitB;
		larmorVector = -mass*x[1].cross(unitB)/charge/B.norm();
		rotatedLarmorVector =
				Eigen::AngleAxisd(rotationAngle/2., unitB) * larmorVector;
		x[0] += rotatedLarmorVector - larmorVector;
		x[1] = Eigen::AngleAxisd(rotationAngle/2., unitB) * x[1];

		// Advance velocity (kick)
		boost::unwrap_ref(odeSystem)(x, dxdt);
		x[1] += dt * dxdt[1];

		// Advance position half timestep (drift)
		x[0] += dt/2. * x[1].dot(unitB) * unitB;
		larmorVector = -mass*x[1].cross(unitB)/charge/B.norm();
		rotatedLarmorVector =
				Eigen::AngleAxisd(rotationAngle/2., unitB) * larmorVector;
		x[0] += rotatedLarmorVector - larmorVector;
		x[1] = Eigen::AngleAxisd(rotationAngle/2., unitB) * x[1];
    }
};

#endif /* STEPPER_H_ */
