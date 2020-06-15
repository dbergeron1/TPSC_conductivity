
#ifndef GREEN0_H
#define GREEN0_H

#include "hamiltonien.h"


/*! \file green0.h
Definition of class green0 
(C) Dominic Bergeron 2007	
*/

extern "C++" {


/*! \class green0		
non interacting Green's function */
class green0: public hamiltonien
{
public:

//@{
//! Constructor
	green0(int =8, int =32, double =0.0, double =1.0, double[] =NULL);
//! Destructor
    ~green0();
//@}


//@{

//! Set the tight-binding parameter
	void set_TB_params(double *params) {hamiltonien::set_TB_params(params); find_mu();}
//! Set temperature
	void set_tem(double T){hamiltonien::set_tem(T); find_mu();}
	
//! Set the density
	void set_density (double n){hamiltonien::set_density (n); find_mu();}
	
//! set the chemical potential mu0 by hand
	void set_mu0(double mu0p);
	
//! set the system size
	void set_grid(int nk0, int Nw0) {set_system_size(nk0); nbw=Nw0; nwmax=Nw0/2; find_mu();}

//@}

//@{
//! Get the tight-binding parameters, density, chemical potential
	void get_params() const;
	double get_mu() const {return mu0;}
//@}


//@{

//! Evaluate the Green's function
//! In Fourier space
	dcomplex Gkw(double *k, dcomplex w)
	{return ((double)1.0)/(w - dispk(k[0],k[1]) + mu0);}
//	{return 1.0/(w + 2*(cos(k[0]) + cos(k[1])) + 4*tp*cos(k[0])*cos(k[1]) + 2*tpp*(cos(2*k[0]) + cos(2*k[1])) + mu0);}
//! Integrand for traceG()
	dcomplex GkwInteg(double kx, double params[])
	{ double k[]={kx,params[0]}; dcomplex w(params[1],params[2]); return Gkw(k,w);}
//! Real part
	double ReGkw(double kx, double params[])
	{ double k[]={kx,params[0]}; dcomplex w(params[1],params[2]); double reG=real(Gkw(k,w));
		return reG;}
//! Imaginary part	
	double ImGkw(double kx, double params[])
	{ double k[]={kx,params[0]}; dcomplex w(params[1],params[2]); double imG=imag(Gkw(k,w));
		return imG;}

//! integrate the Green's function over k
	void G0wn_local();

//! Calulate Green's function on a grid in different spaces
	void calc_Grt(int nr, int nt);
	void calc_Gkt(int nk, int nt);
	void calc_Grw(int nr, int nw);
	void calc_Gkw(int nk, int nw);
	
//! In real space
	dcomplex Grt(int *r, int t);
	double GrtRe(int *r, int t){return real(Grt(r,t));}
	double GrtIm(int *r, int t){return imag(Grt(r,t));}
//! Mixed spaces
	dcomplex Gkt(int *k, int t);
	
	dcomplex Grw(int *r, int nw);
	double GrwRe(int *r, int nw){return real(Grw(r,nw));}
	double GrwIm(int *r, int nw){return imag(Grw(r,nw));}
//@}

//! Trace of Green's functions
	void traceG();

//! Fermi-Dirac distribution
	double Fermi(double kx, double params[]);
	
//! Determine the Fermi surface
//	void find_fermi_surf(int nk);
	void find_FS(int nk);
	
//! Find the non interacting Fermi surface for the finite system
	void find_FS_PBC();	

 protected:
	
//! 
	double ekmu_diag(double kx, double params[]);
	
	
//! Find chemical potential with Fermi function	
	void find_mu();
//find the chemical potential for the infinite system
	void find_mu_inf_syst();
	
//! difference between the actual density and the density calculated for a given chemical potential for the infinite system
	double ddensity_mu(double, double[]);
	
//! difference between the actual density and the density calculated for a given chemical potential for the finite system (periodic boundary conditions) 
	double ddensity_mu_PBC(double m, double params[]);

	void allocate_Gkw(int, int);
	void allocate_Grt(int, int);
	void allocate_Grw(int, int);
	void free_Gkw(){delete [] Gkw_array[0];	delete [] Gkw_array; Gkw_array=NULL; }
	void free_Grt(){delete [] Grt_array[0];	delete [] Grt_array; Grt_array=NULL; }
	void free_Grw(){delete [] Grw_array[0];	delete [] Grw_array; Grw_array=NULL; }
	void free_Gkt(){delete [] Gkt_array[0];	delete [] Gkt_array; Gkt_array=NULL; }

//! 

	double mu0;

	dcomplex **Gkw_array, **Grw_array;
	double **Gkt_array, **Grt_array;
	
	int nbr, nbt, nbw, nwmax;
	
	int nkFS;
	double *FS;
	
//	double *fermi_surf, kFmin, **FS_coord, FSlength;
//	int nksurf, nkFS;

};

} /* extern "C++" */

#endif /* !defined GREEN0_H */

