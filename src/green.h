
#ifndef GREEN_H
#define GREEN_H

#include "chi.h"



/*! \file green.h
Definition of class green
 (C) Dominic Bergeron
*/

extern "C++" {


/*! \class green
interacting Green's function */
class green: public chi
{
public:


//@{
//! Constructor
	green(int =8, int =32, double =1.0, double =1.0, double =0.1, double*  =NULL);
//! Destructor
    ~green();
//@}


//@{
//! Get the tight-binding parameters, density, chemical potential
	void get_params() const;
	double get_mu(){return mu;}
//@}


//@{
//! calculate the self-energy with fast Fourier transform (imaginary time discretized)
	void calc_Self_FFT();
	
//! calculate the self-energy using FFTs and cubic spline
	void calc_Self_FFT_spline();
	
//! Save the first Nwn Matsubara frequencies of the self-energy for a given set of wave vector
	void save_self_ikn(int *k_chijj, int nk_chijj, int Nwn);
	
//! Save the self-energy for a given set of wave vectors and a non-uniform Matsubara frequency grid defined by integer r, at r=0 all the Matsubara frequencies are saved and the frequency grid becomes more sparse as r increases. The sparsity of the grid increases with frequency.
	void save_self_non_uniform_freq(int *k_save, int nk_save, int r);

//! calculate the self-energy with the "brute force" method and return the value
	dcomplex Self_kw(long int [], long int);
	dcomplex Self_kw(double [], long int);
	
//@}

//@{
	//! save the Matsubara Green's function in ascii format 
	//!format: freq. index		kx index		ky index		real(G)		imag(G)
	void save_Green(int *k);
	
//! calculate the spectral weight for a given k=(k_x, k_y) with 0<k_x<nbk-1 for the energy vector w of size nw0
//! using the N points Pade approximant, the result is stored in the vector spectr_w
	void calc_spectral_weight(int *k, double *w, int nw0, int N, double *spectr_w);
	void calc_spectral_weight(double *x, double *w, int nw0, double *spectr_w)
		{ int k[2];double eps=0.5/nbk;double f=(x[0]+eps)*nbk;k[0]=(int)f;f=(x[1]+eps)*nbk;k[1]=(int)f;
		  calc_spectral_weight(k, w, nw0, nbw/2, spectr_w); }
//! using N points Pade approximants, calculate the spectral weight for a given k=(k_x, k_y) with 0<k_x<nbk-1 and 0<k_y<=k_x for the energy vector w of size nw0
//! the result is stored in the vector spectr_w with the index of eta as the leading one
//! z0 contains the Matsubara frequencies used to compute de Pade parameters, coef and func are size N dcomplex vectors
	void calc_spectral_weight(int *k, double *w, int nw0, int N, double *eta_p, int N_eta, double *spectr_w);

//! using Pade approximants, obtain energy distribution curves for all k between k0=(kx[0],ky[0]) and kf=(kx[Nk-1],kyf[Nk-1]) defined as integers
//! if k0=NULL or kf=NULL, the spectral weight is calculated for all k in the reduced Brillouin zone
//! the default value of eta is used
	void EDC(int *kx, int *ky, int Nk, char *fname);
	
//! create a non-uniform frequency index grid, sparsity of the grid increases with m
	void create_non_uniform_ind_freq(int m);
	
//! using Pade approximants, obtain energy distribution curves for all k between k0=(kx[0],ky[0]) and kf=(kx[Nk-1],kyf[Nk-1]) defined as integers
//! if k0=NULL or kf=NULL, the spectral weight is calculated for all k in the reduced Brillouin zone
//! w is the energy vector of size nw0, z0 is the Matsubara frequency vector of size NP, eta_p constains N_eta values of the imaginary part of energy
//! fname contains a prefix for the output file name
	void EDC(int *kx, int *ky, int Nk, double *w, int nw0, int m_freq, int NP_max, double *eta_p, int N_eta, char *fname);
	
//! obtain momentum distribution curve for the energy w, n is a vector containing the indices of the Matsubara frequencies,
//! Nn is the number of Matsubara frequencies to be used (3 or 4),
//! eta is a vector containing different values of the imaginary part of the energy, Neta is the number of those values
	void MDC(int *n, int Nn, double w, double *eta, int Neta, bool use_Pade=false);
	
//! integrate the Green's function over k
	void Gwn_local();
	
//! find the Fermi surface with interactions
	void find_FS_inter();
	
//! Calculate the real frequency self energy with Padé approximants for a given k and the vector w taking into account the Matsubara frequencies with indices in the vector n with the vector eta containing the imaginary part of w
	void calc_Self_Re_w(int *k, int *n, int Nn, double *w, int Nw, double *eta, int Neta);	
	
//! Call the function calc_Self_Re_w()  for all k vectors in the reduced Brillouin zone
	void calc_Self_Re_w_BZ(int *n, int Nn, double *w, int Nw, double *eta, int Neta);
	
	void calc_Self_Re_w_file(int N0, int mp, double *eta2, int Neta);


//@}

//! trace Sigma*G1
	void traceSelfG();

//! trace Sigma*G2
	void traceSelfG2();

//! Find chemical potential
	void find_mu();
		
//! calculate the total density using the spectral weight
	void calc_dens_Re_w(int N0, int mp, double *eta2, int Neta);	
	
//! calculate the total density using the spectral weight obtained with the Pade applied on the total self-energy
	void calc_dens_Re_w_Self_tot(int N0, int mp, double *eta2, int Neta);
	
//! save the Matsubara self-energy in ascii format
//!format: freq. index		kx index		ky index		real(Sigma)		imag(Sigma)
	void save_Self(int Nfreq);

protected:
	
//! return a value of the spectral weight in real frequency using Pade approximant for the third floor of a continued fraction form for the self-energy 
//! the Pade coefficients must have been calculated before with pade_cont_frac_coef, 
//! params[0]: epsilon(k)-mu, params[1] : number of points used for the Pade
//! coefficients calculation, params[2]: small imaginary part in the energy, params[3]...params[6]: coefficients in the continued fraction
	double spectral_weight(double w, double params[]);	
	
//! return a value of the spectral weight in real frequency using Pade approximant for the self-energy minus the 
//! asymptotic part (the terms in powers of (1/(iwn))^j for j=1...4), the Pade coefficients must have been
//! calculated before with pade_cont_frac_coef, params[0]: epsilon(k)-mu, params[1] : number of points used for the Pade
//! coefficients calculation, params[2]: small imaginary part in the energy
	double spectral_weight_Self(double w, double params[]);	
	
//! integrant dans le calcul de n(k) en frequences reelles
	double integrand_density(double w, double params[]);
	
//!load the Matsubara self-energy saved in ascii format
	void load_Self();	
	
//! save the Matsubara self-energy in binary file
	void save_Self_bin();
	
//! load the Matsubara self-energy frome the binary file
	void load_Self_bin();

	//! trace of G-G_0
	double traceDiffG(double pmu, double par[]);

	int *ind_freq;
	int N_freq, m_freq;
	
	double mu;
	bool mu_calcule;
	
	dcomplex *pade_coef;
	dcomplex *pade_z0;
	dcomplex *pade_func;
	
	int nbr_Akw_neg;
	double  sum_Akw_neg;	

	double eta;
	
	dcomplex *Self_kw_array, *Self_array_tot;
	double Sigma_inf, *Sigma_inf2, *Sigma_inf3;
	
	int *FS_inter;
	double *Ak_FS_inter;
	
	double  *Sigma_inf4;

};

} /* extern "C++" */

#endif /* !defined GREEN_H */

