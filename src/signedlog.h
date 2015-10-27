typedef struct signedlog {
	double slval;
	int slsign;
} signedlog;

signedlog* new_signedlog(double val, int sign);

void free_signedlog(signedlog *sl);

double signedlog_to_dbl(signedlog *sl);

void logspace_add_dbl_to_signedlog(signedlog *sl, double dbl);

void logspace_sub_dbl_from_signedlog(signedlog *sl, double dbl);

void logspace_add_signedlogs(double val1, int sign1, double val2, int sign2, double *valtarget, int *signtarget);

void logspace_mult_signedlogs(double val1, int sign1, double val2, int sign2, double *valtarget, int *signtarget);
