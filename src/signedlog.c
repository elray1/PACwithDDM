#include <R.h>
#include <Rmath.h>

#include "utility.h"
#include "signedlog.h"

signedlog* new_signedlog(double slval, int slsign) {
	signedlog *result = Calloc(1, signedlog);
	result->slval = slval;
	result->slsign = slsign;
	return result;
}

void free_signedlog(signedlog *sl) {
	Free(sl);
}

double signedlog_to_dbl(signedlog *sl) {
	return ((double)sl->slsign) * exp(sl->slval);
}

void logspace_add_dbl_to_signedlog(signedlog *sl, double dbl) {
	logspace_add_signedlogs(sl->slval, sl->slsign, dbl, 1, &(sl->slval), &(sl->slsign));
}

void logspace_sub_dbl_from_signedlog(signedlog *sl, double dbl) {
	logspace_add_signedlogs(sl->slval, sl->slsign, dbl, -1, &(sl->slval), &(sl->slsign));
}

void logspace_add_signedlogs(double slval1, int slsign1, double slval2, int slsign2, double *slvaltarget, int *slsigntarget) {
	if(slval1 == -INFINITY && slval2 == -INFINITY) {
		*slvaltarget = -INFINITY;
		*slsigntarget = 1;
	} else if(slsign1 == slsign2) {
		*slvaltarget = logspace_add_safe(slval1, slval2);
		if(*slvaltarget == -INFINITY) {
			*slsigntarget = 1;
		} else {
			*slsigntarget = slsign1;
		}
	} else if(slval1 > slval2) {
		*slvaltarget = logspace_sub(slval1, slval2);
		if(*slvaltarget == -INFINITY) {
			*slsigntarget = 1;
		} else {
			*slsigntarget = slsign1;
		}
	} else {
		*slvaltarget = logspace_sub(slval2, slval1);
		if(*slvaltarget == -INFINITY) {
			*slsigntarget = 1;
		} else {
			*slsigntarget = slsign2;
		}
	}
}

void logspace_mult_signedlogs(double slval1, int slsign1, double slval2, int slsign2, double *slvaltarget, int *slsigntarget) {
	*slvaltarget = slval1 + slval2;
	if(*slvaltarget == -INFINITY) {
		*slsigntarget = 1;
	} else {
		*slsigntarget = slsign1 * slsign2;
	}
}

