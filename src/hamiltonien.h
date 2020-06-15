
#ifndef HAMILTONIEN_H
#define HAMILTONIEN_H

#include "generique.h"


/*! \file hamiltonien.h
Definition of class hamiltonien 
 (C) Dominic Bergeron 2011	
*/

extern "C++" {


/*! \class hamiltonien */
class hamiltonien: public generique
{
public:

//@{
//! Constructor
	hamiltonien(int =8, double =1.0, double =1.0, double[] =NULL);
//! Destructor
    ~hamiltonien();
//@}


//@{
//! Set the tight-binding parameter
	void set_TB_params(double *params) {tp=params[0], tpp=params[1], t3=params[2]; er[2]=-tp, er[3]=-tpp, er[4]=-t3; calc_ek();}
//! Set U
	void set_U(double Up) {U=Up;}	
//! Set temperature
	void set_tem(double T){tem=T;}
//! Set the density
	void set_density (double n){density=n;}
//! set the system size
	void set_system_size(int nk0);
//@}

//@{
//! Get the tight-binding parameters, density and temperature
	void get_params() const;
//@}


protected:
	
//! Calculate the values of the dispersion relation on the discrete grid of the finite system
	void calc_ek();
	
//! give the value of the dispersion relation at the integer wave vector (kx, ky)
	double ek_val(int kx, int ky);
//! give the value of the dispersion relation at an arbitrary wave vector (kx,ky)
	double dispk(double kx, double ky);
//! give the gradient of ek. k=(k[0],k[1])
	void grad_ek(double k[], double gr_ek[]);
//! give an element of the gradient of ek. k=(k[0],k[1]), dir=0 =>x, dir!=0 =>y
	double grad_ek(double k[], int dir);	
//! give an element of the Hessian of ek, k=(k[0],k[1]), ind=(row, column)
	double Hessian_ek(double k[], int ind[]);
	
//! 
	double density;
	
	double U, tp, tpp, t3;

	double tem;

	int nbk;
	
	double *ek, *er;
	
};

} /* extern "C++" */

#endif /* !defined HAMILTONIEN_H */

