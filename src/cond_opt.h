
#ifndef cond_opt_H
#define cond_opt_H

#include "green.h"


/*! \file cond_opt.h
Definition of class cond_opt
 (C) Dominic Bergeron
*/

extern "C++" {


/*! \class cond_opt
calculate the optical conductivity in the Two-particle self-consistent approach*/
class cond_opt: public green
{
public:


//@{
//! Constructor
	cond_opt(int =8, int =256, double =1.0, double =1.0, double =0.1, double*  =NULL);
//! Destructor
    ~cond_opt();
//@}

//! calculate the bubble term only optical conductivity with FFT
	void calc_chijj_bulle();
	void calc_chijj_bulle_1();
	
//! calculate the contributions to the bubble term optical conductivity with FFT for different angles in the Brillouin zone
	void calc_chijj_bulle(int N_theta_k, int *k_chijj, int nk_chijj);	
	
//! calculate all the terms in the optical conductivity with FFT and cubic spline. That version uses less memory than the other one.
	void calc_chijj_vertex_corr_all_optim(bool chi3, int mmin, int N_theta_k, int *k_chijj, int nk_chijj);
	
//! calculate all the terms in the current-current correlation function with FFT and cubic spline
	void calc_chijj_vertex_corr_all(bool chi3, int mmin, int N_theta_k, int *k_chijj, int nk_chijj);		



protected:

	double eta0;

	dcomplex *sigma_array;

	double *chijj_bulle_qn;
	
	double *chijj1_qn, *chijj1_qn_corr, *chijj2_qn, *chijj2_qn_corr, *chijj3_qn, *chijj3_qn_corr;
	double *chijj1_qn4, *chijj2_qn6, *chijj3_qn4;
	double chijj1_inf2, chijj1_inf4, chijj1_inf6, chijj1_w0, chijj2_inf2, chijj2_inf4, chijj2_inf6, chijj2_inf8; 

	
	double *V_array;

	double *ImChijj0, *ReSigma0;

	double chijj_w0, kx0;


	int nbr_Akw_neg;
	double  sum_Akw_neg;
	
	char paramsFileName[200];

};

} /* extern "C++" */

#endif /* !defined cond_opt_H */
