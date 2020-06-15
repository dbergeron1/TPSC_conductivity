/*
 *  hamiltonien.cpp
 *  TPSC
 *
 *  Created by Dominic Bergeron on 11-06-02.
 *
 */


#include "includeDef.h"
#include "hamiltonien.h"


hamiltonien::hamiltonien(int nk0, double n, double T, double  params[])
{
	nbk=nk0/2+1;
	density=n;
	tem=T;
	if (params) { tp=params[0], tpp=params[1], t3=params[2], U=params[3];}
	else { tp=0, tpp=0, t3=0, U=1;}
	
	int nk8th=(nbk*(nbk+1))/2;
	er=new double[nk8th];
	for (int j=0; j<nk8th; j++) er[j]=0;
	
	er[1]=-1;
	er[2]=-tp;
	er[3]=-tpp;
	er[4]=-t3;
	
	ek=NULL;
	calc_ek();

}


hamiltonien::~hamiltonien()
{
	if (ek) delete [] ek;
	if (er) delete [] er;
}

void hamiltonien::get_params() const
{
	cout<<setiosflags(ios::left)<<"U: "<<setw(10)<<U<<"tp: "<<setw(10)<<tp<<"tpp: "<<setw(10)<<tpp<<"t3: "<<setw(10)<<t3<<"density: "<<setw(10)<<density<<"T: "<<tem<<endl;
}

//! Calculate the values of the dispersion relation on the discrete grid of the finite system
void hamiltonien::calc_ek()
{
	if (ek) delete [] ek;
	
	ek=new double[nbk*nbk];
	
	int l,m;
	double kx, ky;
	
//	cout<<setiosflags(ios::left);
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		for(m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			ek[m+l*nbk]=dispk(kx, ky);
//			cout<<setw(20)<<l<<setw(20)<<m<<ek[m+l*nbk]<<endl;
		}
	}
	
	int nk0=2*(nbk-1);
	
	const char *nameForm="./e_k_%5.3f_%5.3f_%5.3f_%1d.dat";
	char name[200];
	
	fstream file;
	sprintf(name, nameForm , (double)tp, (double)tpp, (double)t3, nk0);
	file.open(name, ios::out);
	
	file<<setiosflags(ios::left)<<setprecision(15);
	for (l=0; l<nbk; l++)
		for (m=0; m<nbk; m++)
			file<<setw(10)<<l<<setw(10)<<m<<ek[m+l*nbk]<<endl;
	
	file.close();
}

void hamiltonien::set_system_size(int nk0) 
{
	if (ek) delete [] ek;
	ek=NULL;
	if (er) delete [] er;
	
	nbk=nk0/2+1;
	int nk8th=(nbk*(nbk+1))/2;
	er=new double[nk8th];
	for (int j=0; j<nk8th; j++) er[j]=0;
	
	er[1]=-1;
	er[2]=-tp;
	er[3]=-tpp;
	er[4]=-t3;
	
	calc_ek();
}

//! give the value of the dispersion relation at the integer wave vector (kx, ky)
double hamiltonien::ek_val(int kx, int ky)
{
	int Nk=2*(nbk-1);
	while (kx<0)	kx+=Nk;
	while (ky<0)	ky+=Nk;
	if (kx>nbk-1) kx=Nk-kx;
	if (ky>nbk-1) ky=Nk-ky;
	
	return ek[ky+kx*nbk];
}

//! give the value of the dispersion relation at an arbitrary wave vector (kx,ky)
double hamiltonien::dispk(double kx, double ky)
{
	return -2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - 4*t3*(cos(2*kx)*cos(ky)+cos(2*ky)*cos(kx));
}

//! give the gradient of ek. k=(k[0],k[1])
void hamiltonien::grad_ek(double k[], double gr_ek[])
{
	double kx=k[0], ky=k[1];
	
	gr_ek[0]=2.0*sin(kx)+4.0*tp*sin(kx)*cos(ky)+4.0*tpp*sin(2*kx)+4*t3*(2.0*sin(2*kx)*cos(ky)+cos(2*ky)*sin(kx));
	gr_ek[1]=2.0*sin(ky)+4.0*tp*cos(kx)*sin(ky)+4.0*tpp*sin(2*ky)+4*t3*(cos(2*kx)*sin(ky)+2.0*sin(2*ky)*cos(kx));
}

//! give an element of the gradient of ek. k=(k[0],k[1]), dir=0 =>x, dir!=0 =>y
double hamiltonien::grad_ek(double k[], int dir)
{
	double kx=k[0], ky=k[1];
	
	if (!dir)
		return 2.0*sin(kx)+4.0*tp*sin(kx)*cos(ky)+4.0*tpp*sin(2*kx)+4*t3*(2.0*sin(2*kx)*cos(ky)+cos(2*ky)*sin(kx));
	else
		return 2.0*sin(ky)+4.0*tp*cos(kx)*sin(ky)+4.0*tpp*sin(2*ky)+4*t3*(cos(2*kx)*sin(ky)+2.0*sin(2*ky)*cos(kx));
}

//! give an element of the Hessian of ek, k=(k[0],k[1]), ind=(row, column)
double hamiltonien::Hessian_ek(double k[], int ind[])
{
	double kx=k[0], ky=k[1];
	
	if (!ind[0] && !ind[1])
		return 2.0*cos(kx)+4.0*tp*cos(kx)*cos(ky)+8.0*tpp*cos(2*kx)+4*t3*(4.0*cos(2*kx)*cos(ky)+cos(2*ky)*cos(kx));
	else if (ind[0] && ind[1])
		return 2.0*cos(ky)+4.0*tp*cos(kx)*cos(ky)+8.0*tpp*cos(2*ky)+4*t3*(cos(2*kx)*cos(ky)+4.0*cos(2*ky)*cos(kx));
	else
		return -4.0*tp*sin(kx)*sin(ky)-8.0*t3*(sin(2*kx)*sin(ky)+sin(2*ky)*sin(kx));
}



