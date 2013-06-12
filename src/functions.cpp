#include "functions.h"

bool vect3dLessThan(vect3d a, vect3d b) {
	bool aLessThanB = false;
	// TODO: sort differently than by components?
	// TODO: do components with recursive function?
	if (a[0]<b[0]) {
		aLessThanB = true;
	} else if (a[0]==b[0]) {
		if (a[1]<b[1]) {
			aLessThanB = true;
		} else if (a[1]==b[1]) {
			if (a[2]<b[2]) {
				aLessThanB = true;
			}
		}
	}
	return aLessThanB;
}
