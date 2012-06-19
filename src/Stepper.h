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
//    void do_step(OdeSystem odeSystem, state_type &x, double t, double dt) const {
	void do_step(OdeSystem &odeSystem, state_type &x, double t, double dt) {
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
//			odeSystem(x, dxdt);
//			x[0] += dt/2. * dxdt[0];
			x[0] += dt/2. * x[1];
			// Advance velocity
			odeSystem(x, dxdt);
			x[1] += dt * dxdt[1];
			// Advance position half timestep
//			odeSystem(x, dxdt);
//			x[0] += dt/2. * dxdt[0];
			x[0] += dt/2. * x[1];

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

#endif /* STEPPER_H_ */
