/* File:   curvefit.hpp
 * Role:   auxiliary classes to fit analytic function to raw data
 * Input:  MMSP data file from tricrystal simulation
 * Output: vertex position, dihedral angle, grain size
 */

#ifndef _CURVEFIT_H_

struct data {
    int n;
    double* x;
    double* y;
};

class curvefitter
{
public:
    // constructor
    curvefitter(int N, double* X, double* Y);
    // destructor
    ~curvefitter();
    // accessor
    void fit(double& a, double& theta);

private:
    const int n;
    const int p;
    const int maxiter;
    gsl_matrix* J;
    gsl_matrix* covar;
    gsl_vector* k;
    gsl_vector* res_f;
    double chi;
    double chi0;

    struct data d;

	gsl_multifit_function_fdf f;

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	int info;
};

#endif
