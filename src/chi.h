
#ifndef CHI_H
#define CHI_H


#include "chi0.h"


/*! \file chi.h
Definition of class chi 
(C) Dominic Bergeron
*/

extern "C++" {


/*! \class chi
spin and charge response functions in Two-Particle Self-Consistent approach */
class chi : public chi0
{
public:

//@{
//! Constructor
	chi(int =8, int =32, double =0.0, double =1.0, double* =NULL);
//! Destructor
    ~chi();
//@}


//! Get the Hubbard model parameters, density, chemical potential
	void get_params() const;
    
//! find Usp and Uch
	void find_Usp_Uch();
//! find Usp and Uch for chi0 defined on a given grid, Usc=(Usp,Uch) after the call
	void find_Usp_Uch(int nq0, int nw0, double Usc[]);

//! calculate the correlation length for a set of densities and temperature, nD is the number of densities in the vector dens and nT, the number of temperatures in the vector Temp
	void calc_corr_length(double dens[], int nD, double Temp[], int nT);
	
	//! obtain the spin correlation length
	void sp_corr_length();
	
protected:	
	
//@{

//! Evaluate chi
	void chi_Re_w(int nk, int Nwn, double dens, double T, double *w, int Nwr, int N0, double *eta, int Neta, int m, char *rep_in, char *rep_out, char *nom, int type_fic);

//! return a value of chisp(q,iw), params=(qy,w)
	double chisp(double qx, double params[]);

//! return a value of chich(q,iw),  params=(qy,w)
	double chich(double qx, double params[]);

//! calculate V(r,t)=3*Usp*chisp(r,tau)+Uch*chich(r,tau)
	void calc_chirt();

//! trace of chisp(q,iw)-chi0(q,iw) as a function of Usp=Us
	double traceDiffChispChi0(double Us);

//! trace of chich(q,iw)-chi0(q,iw) as a function of Uch=Uc
	double traceDiffChichChi0(double Uc);

//! function of which the root is Usp
	double sum_rule_sp(double Us, double par[]);

//! function of which the root is Uch
	double sum_rule_ch(double Uc, double par[]);

	double fabs(dcomplex a){ return sqrt( real(a)*real(a) + imag(a)*imag(a) ); }

//! find kx_Fermi for a given ky, return kx+ky-PI and put kx in par[0]
	double kx_plus_ky_Fermi(double y, double par[]);

 
//! Calculate the Fermi velocity at the hot spots
	void vF_hot_spots();


//protected:

	void free_vertex_array(){ if (vertex_rt_array) delete [] vertex_rt_array; vertex_rt_array=NULL;}
	
	double Uch, errUc, Usp, errUs;
	double dblOcc;
	double spCorrLgthS_x, spCorrLgthS_y, spCorrLgth_x, spCorrLgth_y, spCorrLgthAs_x, spCorrLgthAs_y, xi0_x, xi0_y;
	double vF;

	double nsr;

	double *vertex_rt_array;

	char paramsFileName[200];

};

} /* extern "C++" */

#endif /* !defined CHI_H */

