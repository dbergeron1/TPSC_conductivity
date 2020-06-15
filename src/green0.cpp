
/*
 *  green0.cpp
 *  TPSC
 *
 *  Created by Dominic Bergeron
 *
 */

#include "includeDef.h"
#include "green0.h"

green0::green0(int nk0, int Nw0, double n, double T, double  params[]): hamiltonien(nk0, n, T, params) 
{
	nbw=Nw0;
	nwmax=(nbw-2)/2;
	
	if (n!=1.0 || tp || tpp || t3) 
	{
		find_mu_inf_syst();
		find_mu();
	}
	else mu0=0;
	
	Gkw_array=NULL;
	Gkt_array=NULL;
	Grt_array=NULL;
	Grw_array=NULL;
	nbr=0;
	nbt=0;
	
	nkFS=0;
	FS=NULL;
	
	find_FS(10001);
	find_FS_PBC();
}


green0::~green0()
{
	if (Gkt_array) free_Gkt();
	if (Grt_array) free_Grt();
	if (Grw_array) free_Grw();
	if (Gkw_array) free_Gkw();
	
	if (FS) delete [] FS;
}


//! find the non interacting Fermi surface for the finite system
void green0::find_FS_PBC()
{
	cout<<"calcul des coordonnees de la surface de Fermi sans interaction sur le reseau discret\n";
	
	long int l, m, x, y;
	
	dcomplex z;
	double Akw, ekmu;
	int FS;
	double Ak_FS;
	double eta1=0.00001;
	
	fstream file, file2;
	const char *nameForm="./FS_discret_AkF_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d.dat";
	const char *nameForm2="./FS_discret_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d.dat";
	char name[150];
	
	sprintf(name, nameForm, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1));
	
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(8);
	
	sprintf(name, nameForm2, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1));
	
	file2.open(name,ios::out);
	file2<<setiosflags(ios::left)<<setprecision(8);
		
	for (l=0; l<nbk; l++)
	{
		Ak_FS=0;
		FS=0;
		for (m=0; m<nbk-l; m++)
		{	
			x=l+m;
			y=m;
			
			ekmu=ek[x+y*nbk]-mu0;
						
			z=dcomplex(0,eta1);
			
			Akw=-imag(((double)1.0)/(z-ekmu))/PI;
			
			if (Akw>Ak_FS)
			{
				Ak_FS=Akw;
				FS=m;
			}	
			
		}
		file<<setw(10)<<l+FS<<setw(10)<<FS<<Ak_FS<<endl;
		file2<<setw(10)<<l+FS<<setw(10)<<FS<<endl;
	}
	
	file.close();
	file2.close();
	
}



//! integrate the Green's function over k
void green0::G0wn_local()
{
	
	int nk=nbk;
	int nk02=4*(nk-1)*(nk-1);
	int nk8th=(nk*(nk+1))/2;
	
	dcomplex *Gloc=new dcomplex[nbw/2];
	for (int i=0; i<nbw/2; i++) Gloc[i]=0;
	
	fstream file;
	
	const char *nameForm="./G0wn_loc_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	sprintf(name, nameForm, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file.open(name,ios::out);
	file<<setprecision(10)<<setiosflags(ios::left);
	
#pragma omp parallel
	{
		dcomplex z;
		double ektmp;
		int j, l, m, w;
#pragma omp for
		for (j=0; j<nbw/2; j++)
		{
			z=dcomplex(0,(2*j+1)*PI*tem);
			for (l=0; l<nk; l++)
				for(m=0; m<=l; m++)
				{
					ektmp=ek[m+l*nk];
					w=8;
					if (l==0 || l==nk-1) w=w/2;
					if (m==0 || m==nk-1) w=w/2;
					if (m==l) w=w/2;
					
					Gloc[j]+=((double)(1.0*w))/(((double)(1.0*nk02))*(z-ektmp+mu0));
				}
		}
	}
	
	for (int j=0; j<nbw/2; j++)		file<<setw(10)<<j<<setw(25)<<Gloc[j].real()<<Gloc[j].imag()<<endl;
	
	file.close();
	
	delete [] Gloc;
}

void green0::get_params() const
{
	hamiltonien::get_params();
	cout<<"mu0: "<<mu0<<'\n';
}


double green0::Fermi(double kx, double params[])
{
	double ky=params[0], mupar=params[1];

	double ektmp=dispk(kx, ky);

//	cout<<"green0: tem:  "<<tem<<endl;
	
	return 1.0/( exp((ektmp-mupar)/tem) + 1.0 );
}

//! difference between the actual density and the density calculated for a given chemical potential for the finite system (periodic boundary conditions)
double green0::ddensity_mu_PBC(double mup, double params[])
{
	double dens=0;	
	double w;
	int Nk2=4*(nbk-1)*(nbk-1);

	int l,m;
 	for (l=0; l<nbk; l++)
 	{
 		for (m=0; m<=l; m++)
 		{
			w=8;
			if (l==0 || l==nbk-1) w=w/2;
			if (m==0 || m==nbk-1) w=w/2;
			if (m==l) w=w/2;

			dens+=2*w/((exp((ek[m+l*nbk]-mup)/tem)+1.0)*Nk2);
		 }
 	}

	return dens-density;
}

//! set the chemical potential mu0 by hand
void green0::set_mu0(double mu0p)
{
	mu0=mu0p;
	
	double params[]={0, 0};
	
	double dens=ddensity_mu(mu0p, params)+density;
	
	density=dens;
	
	cout<<setiosflags(ios::left)<<setprecision(15);
	cout<<"density for mu0= "<<mu0<<"  :  "<<density<<endl;
}

//! difference between the actual density and the density calculated for a given chemical potential
double green0::ddensity_mu(double mup, double params[])
{
	double tol=1.0e-10;
	
	double limx[]={0.0, PI}, limy[]={0.0, PI};
	int nbEval[]={0};
	double dens;
	
	params[1]=mup;
	
	fctPtr FermiPtr=static_cast<fctPtr>(&green0::Fermi);
	
	dens=1.0/(PI*PI)*2*quadInteg2D(FermiPtr, limx, limy, tol, nbEval, params);
	
	return dens-density;
}


//find the chemical potential for the finite system
void green0::find_mu()
{		
	if (density==1.0 && tp==0 && tpp==0)
	{
		mu0=0.0;
		return;
	}

	fctPtr densPtr=static_cast<fctPtr>(&green0::ddensity_mu_PBC);

	double root[]={0.0};
	double init[]={0.0, -8.0, 8.0};
	double *params=NULL;

	if ( find_zero(densPtr, init, params, root) )
		mu0=root[0];
	else
		cout<<"mu0 non trouve\n";

	cout<<"mu0 du systeme fini: "<<mu0<<endl;
//	cout<<"densite calculee avec mu0: "<<ddensity_mu_PBC(mu0, params)+density<<endl;
	
}

//find the chemical potential for the infinite system
void green0::find_mu_inf_syst()
{
	if (density==1.0 && tp==0 && tpp==0)
	{
		mu0=0.0;
		return;
	}

	fctPtr densPtr=static_cast<fctPtr>(&green0::ddensity_mu);

	double root[]={0.0};
	double init[]={0.0, -8.0, 8.0};
	double params[]={0.0, 0.0};

	if ( find_zero(densPtr, init, params, root) )
		mu0=root[0];
	else
		cout<<"mu0 non trouve\n";
	
	cout<<"mu0 du systeme infini: "<<mu0<<endl;
//	cout<<"densite calculee avec mu0: "<<ddensity_mu(mu0, params)+density<<endl;

}

/*
// calculate the trace of G(k,ikn)
void green0::traceG()
{
 double params[]={0.0,0.0,0.0}, limx[]={0.0, PI}, limy[]={0.0, PI}, tol=1e-6;
 int nbEvalRe[]={0}, nbEvalIm[]={0};
 int nbEval[]={0};
 dcomplex retrG=0.0, imtrG=0.0;
 dcomplex trG=0.0, convFact=0.0;

// fctPtr rePtr=static_cast<fctPtr>(&green0::ReGkw);
// fctPtr imPtr=static_cast<fctPtr>(&green0::ImGkw);

 complexFctPtr Ptr=static_cast<complexFctPtr>(&green0::GkwInteg);

 int nmax=1000000;

 double eta=1e-3;


 for (int n=-nmax-1; n<=nmax; n++)
 {
	 params[2]=(2*n+1)*PI*tem;
	 convFact=exp(dcomplex(0.0,(2*n+1)*PI*tem*eta));

	 nbEval[0]=0;
	 nbEvalRe[0]=0;
	 nbEvalIm[0]=0;

	 trG+=tem*convFact/(PI*PI)*2.0*quadInteg2D(Ptr, limx, limy, tol, nbEval, params);
 }
 cout<<"traceG: "<<trG<<'\n';
 double n=real(retrG)-imag(imtrG);

// cout<<"retrG: "<<retrG<<'\n';
// cout<<"imtrG: "<<imtrG<<'\n';
 cout<<"density: "<<density<<'\n';
}
*/

void green0::allocate_Gkw(int nk, int nw)
{
	set_grid(nk, nw);
	
// nbk=nk;
// nbw=nw;
// nwmax=(nw-2)/2;

	int j;
	int nk2=nbk*nbk;


	Gkw_array=new dcomplex*[nbw];

	long int array_size=(long int) nk2*nbw;

	Gkw_array[0] = new dcomplex[array_size];	
//	Gkw_array[0] = new dcomplex[nk2*nbw];

    for(j=1;j<nbw;j++)
        Gkw_array[j]=Gkw_array[j-1]+nk2;

	for (j=0; j<nbw; j++)
		for (int l=0; l<nk2; l++)
			Gkw_array[j][l]=0.0;
}

void green0::allocate_Grt(int nr, int nt)
{
	nbr=nr;
	nbt=nt;

	int j;
	int nr2=nbr*nbr;


	Grt_array=new double*[nbt];

	long int array_size=(long int) nr2*nbt;

	Grt_array[0] = new double[array_size];
//	Grt_array[0] = new double[nr2*nbt];

    for(j=1;j<nbt;j++)
        Grt_array[j]=Grt_array[j-1]+nr2;

	for (j=0; j<nbt; j++)
		for (int l=0; l<nr2; l++)
			Grt_array[j][l]=0.0;
}


void green0::allocate_Grw(int nr, int nw)
{
	nbr=nr;
	nbw=nw;
	nwmax=(nw-2)/2;
	
	int j;
	int nr2=nbr*nbr;

	Grw_array=new dcomplex*[nbw];
	
	long int array_size=(long int) nr2*nbw;

	Grw_array[0] = new dcomplex[array_size];
//	Grw_array[0] = new dcomplex[nr2*nbw];

    for(j=1;j<nbw;j++)
        Grw_array[j]=Grw_array[j-1]+nr2;

	for (j=0; j<nbw; j++)
		for (int l=0; l<nr2; l++)
			Grw_array[j][l]=0.0;
}

/*
void green0::calc_Gkt(int nk, int nt)
{
	if (Gkt_array)
		free_Gkt();

	nbk=nk;
	nbt=nt;
	int nk2=nk*nk;

	Gkt_array[0]=new double[nk2*nt];

	int j,l,m;
	for (j=1; j<nt; j++)
		Gkt_array[j]=Gkt_array[j-1]+nk2;

	for (j=0; j<nt; j++)
		for (l=0; l<nk2; l++)
			Gkt_array[j][l]=0.0;

	double kx, ky, t, ekmu;

	for (l=0; l<nk; l++)
	{
		kx=l*PI/(nk-1);
		for(m=0; m<nk; m++)
		{
			ky=m*PI/(nk-1);
//			ekmu=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu0;
			ekmu=ek[m+l*nbk] - mu0;
			Gkt_array[0][m + l*nk]=-1.0/(exp(-ekmu/tem)+1);

		}
	}

	for (j=1; j<nt; j++)
	{
		t=j/(2.0*(nt-1)*tem);
		for (l=0; l<nk; l++)
		{
			kx=l*PI/(nk-1);
			for(m=0; m<nk; m++)
			{
				ky=m*PI/(nk-1);
//				ekmu=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu0;
				ekmu=ek[m+l*nbk] - mu0;
				Gkt_array[j][m + l*nk]=exp(-t*ekmu)*Gkt_array[0][m + l*nk];
			}
		}
	}

}
*/

void green0::calc_Gkw(int nk, int nw)
{
	if (fmod(nw,2.0)!=0)
	{
		 cerr<<"le nombre de valeurs de wn doit etre un multiple de deux\n";
		 return;
	}

	if (!Gkw_array)
		allocate_Gkw(nk, nw);
	else if (nbk!=nk || nbw!=nw)
	{
		free_Gkw();
		allocate_Gkw(nk, nw);
	}

	 int n,lx,ly;

	dcomplex w;
	double k[]={0.0, 0.0};

	for (n=-nwmax-1; n<=nwmax; n++)
	{
		w=dcomplex(0.0,(2.0*n+1)*PI*tem);
		for (ly=0; ly<nbk; ly++)
		{
			k[1]=ly*2*PI/nbk;
			for (lx=0; lx<nbk; lx++)
			{
				k[0]=lx*2*PI/nbk;

			 	Gkw_array[n+nwmax+1][ly*nbk+lx]=Gkw(k,w);
			}
		}
	}

}


void green0::calc_Grt(int nr, int nt)
{
	fftw_plan fftplan;

	if (fmod(nt,2.0)!=0)
	{
		cerr<<"le nombre de valeurs de t doit etre un multiple de deux\n";
		return;
	}

	int nwtmp=nbw;

	if (!Gkw_array)
		allocate_Gkw(nr, nt);
	else if (nbk!=nr || nbw!=nt)
	{
		free_Gkw();
		allocate_Gkw(nr, nt);
	}

	if (!Grw_array)
		allocate_Grw(nr, nt);
	else if (nbr!=nr || nwtmp!=nt)
	{
		free_Grw();
		allocate_Grw(nr, nt);
	}

	int nr2=nr*nr;
	int j,l;

	calc_Gkw(nr, nt);

	dcomplex *Gkw_tmp, *Grw_tmp;

	Gkw_tmp=new dcomplex[nr2];
	Grw_tmp=new dcomplex[nr2];

	fftplan=fftw_plan_dft_2d(nr, nr, reinterpret_cast<fftw_complex*>(Gkw_tmp), reinterpret_cast<fftw_complex*>(Grw_tmp),
						 FFTW_BACKWARD, FFTW_ESTIMATE);

	for (j=0; j<nt;j++)
	{
		for (l=0; l<nr2; l++)
			Gkw_tmp[l]=Gkw_array[j][l];

		fftw_execute(fftplan);

		for (l=0; l<nr2; l++)
			Grw_array[j][l]=((double)1.0/nr)*Grw_tmp[l];
	}

	fftw_destroy_plan(fftplan);

	delete [] Grw_tmp;
	delete [] Gkw_tmp;

	if (!Grt_array)
		allocate_Grt(nr, nt);
	else if (nbr!=nr || nbt!=nt)
	{
		free_Grt();
		allocate_Grt(nr, nt);
	}

	dcomplex *Grt_tmp;
	Grt_tmp=new dcomplex[nt];
	Grw_tmp=new dcomplex[nt];

	fftplan=fftw_plan_dft_1d(nt, reinterpret_cast<fftw_complex*>(Grw_tmp), reinterpret_cast<fftw_complex*>(Grt_tmp),
						 FFTW_FORWARD, FFTW_ESTIMATE);

	for (j=0; j<nr2; j++)
	{
		for (l=0; l<nt; l++)
			Grw_tmp[l]=Grw_array[l][j];

		fftw_execute(fftplan);

		for (l=0; l<nt; l++)
			Grt_array[l][j]=1.0/sqrt((double)nt)*real(exp(dcomplex(0.0,((nt-1)*l*PI)/nt))*Grt_tmp[l]);
	}

	fftw_destroy_plan(fftplan);
	delete [] Grw_tmp;
	delete [] Grt_tmp;

}

/*
dcomplex green0::Grw(int *r, int nw)
{
	if (!Grw_array)
	{
		cerr<<"Grw_array n'existe pas";
		return 0.0;
	}
	if (r[0]>=0 && r[1]>=0 && r[0]<nbr && r[1]<nbr && nw>=0 && nw<nbw)
		return Grw_array[nw][r[1]*nbr+r[0]];
	else
	{
		cerr<<"indices hors limites";
		return 0.0;
	}
}
*/

dcomplex green0::Grt(int *r, int t)
{
	if (!Grt_array)
	{
		cerr<<"Grt_array n'existe pas";
		return 0.0;
	}
	if (r[0]>=0 && r[1]>=0 && r[0]<nbr && r[1]<nbr && t>=0 && t<nbt)
		return Grt_array[t][r[1]*nbr+r[0]];
	else
	{
		cerr<<"indices hors limites";
		return 0.0;
	}
}

double green0::ekmu_diag(double kx, double kx0[])
{	
	return dispk(kx,kx-kx0[0])-mu0;
}

void green0::find_FS(int nk)
{
	if (FS) delete [] FS;
	
	cout<<"calcul des coordonnees de la surface de Fermi sans interaction\n";
	
	nkFS=nk;
	FS= new double[2*nk];
	int j;
	for (j=0; j<nk; j++) 
	{
		FS[2*j]=-1;
		FS[2*j+1]=-1;
	}
	
	fctPtr ptr=static_cast<fctPtr> (&green0::ekmu_diag);
	
	double init[]={0,0,PI}, lims[]={0,PI};
	double root[1];
	double kx0[1];
	
	fstream file;
	const char *nameForm="./FS_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d.dat";
	char name[150];
	
	sprintf(name, nameForm, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk);
	
	file.open(name,ios::out);
	file<<setiosflags(ios::left);
	file<<setprecision(10);

	cout<<setiosflags(ios::left);
	for (j=0; j<nk; j++)
	{
		kx0[0]=j*PI/(nk-1);
		init[1]=kx0[0];
		init[0]=0.5*(init[1]+init[2]);
		lims[0]=kx0[0];
		if (find_zero(ptr, init, kx0, root, lims))
		{
			FS[2*j]=root[0];
			FS[2*j+1]=root[0]-kx0[0];
		}
		else
			break;
		
		file<<setw(25)<<FS[2*j]<<FS[2*j+1]<<endl;
//		cout<<setw(20)<<FS[2*j]<<setw(20)<<FS[2*j+1]<<FS[2*j]+FS[2*j+1]<<endl;
	}
	file.close();
}


