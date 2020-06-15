
#ifndef CHI0_H
#define CHI0_H

#include "green0.h"


/*! \file chi0.h
Definition of class chi0 
(C) Dominic Bergeron
*/

extern "C++" {


/*! \class chi0
Lindhard function (non interacting response function) */
class chi0 : public green0
{
public:


//! Constructor
	chi0(int nk0=8, int Nw0=32, double n=1.0, double T=1.0, double params[]=NULL);
//! Destructor
    ~chi0();
	
//!	set keep_chirt
	void set_keep_chirt(bool v) {keep_chirt=v;}
//!	set keep_chiqw
	void set_keep_chiqw(bool v) {keep_chiqw=v;}
//! set save_chi0_q	
	void set_save_chi0_q(bool v) {save_chi0_q=v;}
//! set save_chi0_to_file	
	void set_save_chi0_to_file(bool v) {save_chi0_to_file=v;}	
	
//! In Fourier space
//! Calculate chi(q,iw) for one value of (q,w), par=(qy, w) on an infinite lattice
	double chiqw(double qx, double par[]);
	
//! Calculate chi(q,iwn) on the finite lattice
	double chiqw_discret(int *q, int n);	

//! calculate an array of chi0 on a grid using fast Fourier transform. imaginary time is discretized
//	void calc_chiqwFFT(int nq0, int nw0);
 
//! calculate an array of chi0 on a grid using fast Fourier transform and cubic splines for the time transform
	void calc_chiqw_FFT_spline(int nq0, int nw0);	
	
//! calculer chi0 par interpolation par splines cubiques et par transformee de Fourier rapide en enlevant les nw0/(2^r) dernières fréquences
	void calc_chiqw_FFT_spline_optim(int nq0, int nw0, int r);
	
//! calculate chi0(q,w) in the region delimited by q1=(qx1,qy1), q2=(qx2,qy2) at the frequency wn=2*n*pi*tem. qlims={qx1, qy1, qx2, qy2}. Nq={Nqx,Nqy} is the number of points in qx and qy directions
	void calc_chi0_qlims(int n, double qlims[], int Nq[]);
	
//! save chi0(q, w) for a given real frequency w
	void save_chi0q(double w);
	
//! save chi0(q,iqn)
	void save_chi();
	
//! trace of chi0(q,iwn)
	void traceChi0();	


	double fabs(dcomplex a){ return sqrt( real(a)*real(a) + imag(a)*imag(a) ); }
	
protected:

	double ekmu_diag(double kx, double params[]);
	
//! load chi0 from a file if it exists
	void loadChi0();

//! free chiqwRe_array
	void free_chi_array(){ if (chiqw_array) {delete [] chiqw_array; chiqw_array=NULL; }}

//! free chirt
	void free_chirt(){ if (chirt) {delete [] chirt; chirt=NULL; }}
	
//! free dtau_chirt0
	void free_dtau_chirt0(){if (dtau_chirt0) {delete [] dtau_chirt0; dtau_chirt0=NULL; }}
	
	void free_dtau_Grt0p(){if (dtau_Grt0p) {delete [] dtau_Grt0p; dtau_Grt0p=NULL; }}
	
//! integrand in the calculation of chi0
	double chi0IntegRe(double kx, double params[]);


	fstream chi0File;

	double tolChi0Integ;

	bool keep_chirt;
	bool keep_chiqw;
	bool save_chi0_q;
	bool save_chi0_to_file;
	
	double *chiqw_array, *dtau_chirt0, *dtau_Grt0p, *chirt;
	double *chi0_inf, *chi0_inf2;
	
	double *chiqw_array_c, *grid;

	double chi0max;
	int qmax[2];

	int nbq;
	int nqc;
	int nwr;
	
//	double **qpeaks_coord;

};

} /* extern "C++" */

#endif /* !defined CHI0_H */

