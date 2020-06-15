/*
 *  green.cpp
 *  TPSC
 *
 *  Created by Dominic Bergeron
 *
 */


#include "includeDef.h"
#include "green.h"

green::green(int nk0, int Nw0, double n, double T, double eta0, double  params[]):chi(nk0, Nw0, n, T, params)
{
		
	eta=eta0;
	
	mu=0;
	mu_calcule=false;
	
	Self_kw_array=NULL;
	Self_array_tot=NULL;
	Usp=0.0;
	Uch=0.0;
	Sigma_inf2=NULL;
	Sigma_inf3=NULL;
	FS_inter=NULL;
	Ak_FS_inter=NULL;
	
	Sigma_inf4=NULL;
	
	pade_coef=NULL;
	pade_z0=NULL;
	pade_func=NULL;

	ind_freq=NULL;
}

green::~green()
{
	 
	if (Self_kw_array) delete [] Self_kw_array;
	if (Self_array_tot) delete [] Self_array_tot;
	if (Sigma_inf2) delete [] Sigma_inf2;
	if (Sigma_inf3) delete [] Sigma_inf3;
	if (FS_inter) delete [] FS_inter;	
	if (Ak_FS_inter) delete [] Ak_FS_inter;
	
	if (Sigma_inf4) delete [] Sigma_inf4;
	 
 	if (pade_coef) delete [] pade_coef;
 	if (pade_z0) delete [] pade_z0;
	if (pade_func) delete[] pade_func;
	
	if (ind_freq) delete [] ind_freq;
}

//! Save the self_energy for a given set of wave vector
void green::save_self_non_uniform_freq(int *k_save, int nk_save, int r)
{
	int i, j, l, p;
	int nw=nbw/2;
	
	if (!Self_kw_array)
	{
		set_keep_chiqw(false);
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}
	
	int N1=nw/((int)pow(2,r));
	int N2=N1/2;
	int Nf=N1+r*N2+1;
	int *n=new int[Nf];
	
	for (j=0; j<N1; j++) n[j]=j;
	for (j=0; j<Nf-N1; j++)
	{
		l=j/N2;
		n[i+N1]=(j % N2)*((int)pow(2,l+1)) + N1*((int)pow(2,l));
	}
	
	fstream file;
	char name[200];
	const char *nameForm_self_k="./self_k_ikn_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	int nk8th=(nbk*(nbk+1))/2;
	dcomplex selftmp;
	long int x, y, tmp;
	double C2, C3, kn;
	for (j=0; j<=Nf; j++)
	{
		kn=(2.0*n[j]+1)*PI*tem;
		for (p=0; p<nk_save; p++)
		{
			l=0;
			//			for (l=-1; l<=1; l++)
			{
				x=k_save[2*p]+l;
				y=k_save[2*p+1]+l;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				if (x>=0 && x<nbk && y>=0)
				{
					C2=Sigma_inf2[y+(x*(x+1))/2];
					C3=Sigma_inf3[y+(x*(x+1))/2];
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+n[j]*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+n[j]*nk8th])/(kn*kn*kn*kn));
					
					file<<setw(10)<<n[j]<<setw(10)<<x<<setw(10)<<y<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
				}
			}
		}
	}
	file.close();
}

//! Save the self_energy for a given set of wave vector
void green::save_self_ikn(int *k_save, int nk_save, int Nwn)
{	
	int nw=nbw/2;
	
	if (!Self_kw_array)
	{
		set_keep_chiqw(false);
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();		
	}
	
	fstream file;
	char name[200];
	const char *nameForm_self_k="./self_k_ikn_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	int nk8th=(nbk*(nbk+1))/2;
	int j, l, p;
	dcomplex selftmp;
	long int x, y, tmp;
	double C2, C3, kn;
	for (j=0; j<=Nwn; j++)
	{
		kn=(2.0*j+1)*PI*tem;
		for (p=0; p<nk_save; p++)
		{
			l=0;
//			for (l=-1; l<=1; l++)
			{
				x=k_save[2*p]+l;
				y=k_save[2*p+1]+l;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				if (x>=0 && x<nbk && y>=0)
				{
					C2=Sigma_inf2[y+(x*(x+1))/2];
					C3=Sigma_inf3[y+(x*(x+1))/2];
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
					
					file<<setw(10)<<j<<setw(10)<<x<<setw(10)<<y<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
				}
			}
		}
	}
	file.close();
	
/*	
	for (j=0; j<=nbw/8; j++)
	{
		kn=(2.0*j+1)*PI*tem;
		for (p=0; p<nk_save; p++)
		{
			x=k_save[2*p];
			y=k_save[2*p+1];
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			C2=Sigma_inf2[y+(x*(x+1))/2];
			C3=Sigma_inf3[y+(x*(x+1))/2];
			
			selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
			
			file<<setw(10)<<j<<setw(10)<<k_save[2*p]<<setw(10)<<k_save[2*p+1]<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
		}
	}
	file.close();
 */
}

//! integrate the Green's function over k
void green::Gwn_local()
{
	if (!mu_calcule || !Self_kw_array)
	{
		cout<<"G_local(): mu non calcule ou Self_kw_array inexistant\n";
		return;
	}
	
	int nk=nbk;
	int nk02=4*(nk-1)*(nk-1);
	int nk8th=(nk*(nk+1))/2;

	dcomplex *Gloc=new dcomplex[nbw/2];
	for (int i=0; i<nbw/2; i++) Gloc[i]=0;
	
	fstream file;
	
	const char *nameForm="./Gwn_loc_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";
	char name[100];
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file.open(name,ios::out);
	file<<setprecision(10)<<setiosflags(ios::left);
	
#pragma omp parallel
	{
		dcomplex z, selftmp;
		double ektmp;
		long int j, l, m;
		int w;
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
					
					selftmp=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
					
					Gloc[j]+=((double)w/nk02)/(z-ektmp+mu-selftmp);
				}
		}
	}
	
	for (int j=0; j<nbw/2; j++)		file<<setw(10)<<j<<setw(25)<<Gloc[j].real()<<Gloc[j].imag()<<endl;
	
	file.close();
	
	delete [] Gloc;
	
}

void green::calc_Self_Re_w_file(int N0, int mp, double *eta2, int Neta)
{

	fstream file;
//	const char *forme_nom_in="/home/result_Bourbonnais/kperp_%1d_neq.dat";
	const char *forme_nom_in="/home/result_Bourbonnais/Selfw_n12_05_09.dat";
	const char *forme_nom_out="/home/result_Bourbonnais/Self_wRe_kperp_%1d_%1d_neq.dat";
	const char *forme_nom_eta="/home/result_Bourbonnais/imag_freq_%1d_%1.0f_%1.0f_neq.dat";
	char nom[100];

	double T=100;
	
	dcomplex *coef=new dcomplex[N0];
	dcomplex *z0=new dcomplex[N0];
	dcomplex *func=new dcomplex[N0];
	int *n_pade=new int[N0];

	pade_coef=coef;
	pade_z0=z0;
	
	int j, l, m, n;

	int Nwn=40, Np, Np0=37;
	int nkperp=1;
	int k;
	double wn;

//	for (j=0; j<Nwn; j++) z0[j]=dcomplex(0,(2*j+1)*PI*T);

	double wf=10000;
	int Nwr=10001;
	double dw=2*wf/(Nwr-1);

	double *w=new double[Nwr];

	w[0]=-wf;
	for (j=1; j<Nwr; j++) w[j]=w[j-1]+dw;

	dcomplex self;
	

	sprintf(nom,forme_nom_eta,Neta,(double)eta2[0],(double)eta2[Neta-1]);
	file.open(nom,ios::out);
	file<<setiosflags(ios::left);
	for (j=0; j<Neta; j++) file<<setw(15)<<eta2[j];
	file.close();


	double n0, nceil, nfloor, nint, wnx;
	cout<<setiosflags(ios::left)<<setprecision(16);
	for (k=0; k<nkperp; k++)
	{
//		sprintf(nom, forme_nom_in , k);
		strcpy(nom,forme_nom_in);
		file.open(nom, ios::in);
		if (!file)
		{
			cout<<"calc_Self_Re_w(): fichier  "<<nom<<"  inexistant\n";
			return;
		}
        double tmpRe, tmpIm;
		for (j=0; j<Nwn; j++)
		{
/*
			file>>wn>>func[j].real()>>func[j].imag();
			n0=(wn/(PI*T)-1.0)/2.0;
			nceil=ceil(n0);
			nfloor=floor(n0);
			if (fabs(nceil-n0)<fabs(nfloor-n0)) nint=nceil;
			else nint=nfloor;
			wnx=(2*nint+1.0)*PI*T;
			z0[j]=dcomplex(0,wnx);
*/
            
//			file>>n>>func[j].real()>>func[j].imag();
            file>>n>>tmpRe>>tmpIm;
            func[j]=dcomplex(tmpRe,tmpIm);
			z0[j]=dcomplex(0,(2*n+1)*PI*T);

//			if (k==0) 
//				cout<<setw(25)<<wn<<wnx<<endl;
//				cout<<setw(25)<<n0<<setw(25)<<nint<<setw(25)<<wn<<wnx<<endl;
		}
		file.close();
		
		for (Np=Np0; Np<=Nwn; Np++ )
		{
//			pade_cont_frac_coef(func, z0, Np, coef);
			pade_cont_frac_coef_rec(func, z0, Np, coef);			
			sprintf(nom, forme_nom_out, k, Np);
			file.open(nom, ios::out);
			file<<setiosflags(ios::left)<<setprecision(16);
			for (j=0; j<Nwr; j++)
			{
				file<<setw(30)<<w[j];
				for (l=0; l<Neta; l++)
				{
					self=pade(dcomplex(w[j],eta2[l]), Np, z0, coef);
					file<<setw(25)<<self.real()<<setw(35)<<self.imag();
				}
				file<<endl;
			}
			file.close();
		}		
	}

	delete [] w;
	
	delete [] n_pade;
	delete [] coef;
	delete [] z0;
	delete [] func;
	
	pade_coef=NULL;
	pade_z0=NULL;	
	
}

//! Call the function calc_Self_Re_w()  for all k vectors in the reduced Brillouin zone
void green::calc_Self_Re_w_BZ(int *n, int Nn, double *w, int Nw, double *eta, int Neta)
{

	int nk2=(nbk-1)/2+1;
	
#pragma omp parallel
	{
		int l,m, k[2];
		
#pragma omp for		
		for (l=0; l<nk2; l++)
			for (m=0; l+m<=nbk-m-1; m++)
			{
				k[0]=nbk-m-1;
				k[1]=l+m;
				calc_Self_Re_w(k, n, Nn, w, Nw, eta, Neta);
			}
	}
	
}

//! Calculate the real frequency self energy with Padé approximants for a given k and the vector w taking into account the Matsubara frequencies with indices in the 
// vector n with the vector eta containing the imaginary part of w
void green::calc_Self_Re_w(int *k, int *n, int Nn, double *w, int Nw, double *eta_p, int Neta)
{
	dcomplex *coef=new dcomplex[Nn];
	dcomplex *z0=new dcomplex[Nn];
	dcomplex *func=new dcomplex[Nn];
	
	int nk8th=(nbk*(nbk+1))/2;
	
	long int x=k[0];
	long int y=k[1];
	if (y>x)
	{
		y=x;
		x=k[1];
	}
	
	long int j;
	dcomplex z;
	for (j=0; j<Nn; j++)
	{
		z0[j]=dcomplex(0.0, (2*n[j]+1)*PI*tem);
		z=z0[j];
		func[j]=Sigma_inf/z+Sigma_inf2[y + (x*(x+1))/2]/(z*z)+Sigma_inf3[y + (x*(x+1))/2]/(z*z*z)+Self_kw_array[y + (x*(x+1))/2+j*nk8th]/(z*z*z*z);
	}
	
	pade_cont_frac_coef_rec(func, z0, Nn, coef);
	
	fstream fileSelf, fileAkw;
	const char nameFormSelf[]="Self_Re_w_U%4.2f_tp%5.3f_tpp%5.3f_tt%5.3f_dens%7.5f_tem%6.4f_Nk%1d_Nwn%1d_%1d_%1d.dat";
	const char nameFormAkw[]="spectr_wgt_U%4.2f_tp%5.3f_tpp%5.3f_tt%5.3f_dens%7.5f_tem%6.4f_Nk%1d_Nwn%1d_%1d_%1d.dat";
	char name[200];
	
	sprintf(name, nameFormSelf , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw, k[0], k[1]);
	fileSelf.open(name,ios::out);
	fileSelf<<setiosflags(ios::left)<<setprecision(10);
	sprintf(name, nameFormAkw , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw, k[0], k[1]);
	fileAkw.open(name,ios::out);
	fileAkw<<setiosflags(ios::left)<<setprecision(10);
	
	double kx=k[0]*PI/(nbk-1),  ky=k[1]*PI/(nbk-1);
	
//	cout<<"kx:  "<<kx<<"    ky:  "<<ky<<endl;
//	double ekmu=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu;
	double ekmu=ek_val(k[0],k[1]) - mu;	

	int l;
	dcomplex self_w;
	double Akw;
	
	dcomplex *selftmp1=new dcomplex[Nw];

	for (j=0; j<Nw; j++)
	{
		z=dcomplex(w[j],eta_p[0]);
		
		self_w=pade(z, Nn, z0, coef);
		selftmp1[j]=self_w;
		
		Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w) );
		
		fileSelf<<setw(25)<<w[j]<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
		fileAkw<<setw(25)<<w[j]<<setw(20)<<Akw;
		
		for (l=1; l<Neta; l++)
		{
			z=dcomplex(w[j],eta_p[l]);
			
			self_w=pade(z, Nn, z0, coef);
			
			Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w) );
			
			fileSelf<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
			fileAkw<<setw(20)<<Akw;
		}
		fileSelf<<endl;
		fileAkw<<endl;
	}
	fileSelf.close();
	fileAkw.close();
	
	double c1=Sigma_inf, c2=Sigma_inf2[y + (x*(x+1))/2], c3=Sigma_inf3[y + (x*(x+1))/2];
	double a1=c1, b1=-c2/c1;
	double a2=-c3/c1+b1*b1;

	dcomplex fk;
	int Nn2=Nn/8;
	if (Nn2<4) Nn2=4;
	for (j=0; j<Nn2; j++)
	{
		z=z0[j];
		fk=Self_kw_array[y + (x*(x+1))/2+j*nk8th];
		func[j]=(b1*fk + (fk+b1*c3)*z + (c3-a1*b1*b1)*z*z)/(-fk - c3*z + a1*b1*z*z - a1*z*z*z);
	}
	
	pade_cont_frac_coef_rec(func, z0, Nn2, coef);

	const char nameFormSelf_Jfrac[]="Self_Re_w_Jfrac_U%4.2f_tp%5.3f_tpp%5.3f_tt%5.3f_dens%7.5f_tem%6.4f_Nk%1d_Nwn%1d_%1d_%1d.dat";
	sprintf(name, nameFormSelf_Jfrac, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw, k[0], k[1]);
	fileSelf.open(name,ios::out);
	fileSelf<<setiosflags(ios::left)<<setprecision(10);
	
	const char nameFormAkw_Jfrac[]="spectr_wgt_Jfrac_U%4.2f_tp%5.3f_tpp%5.3f_tt%5.3f_dens%7.5f_tem%6.4f_Nk%1d_Nwn%1d_%1d_%1d.dat";
	sprintf(name, nameFormAkw_Jfrac , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw, k[0], k[1]);
	fileAkw.open(name,ios::out);
	fileAkw<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex sub_self;
	for (j=0; j<Nw; j++)
	{
		z=dcomplex(w[j],eta_p[0]);
		
		sub_self=pade(z, Nn2, z0, coef);
		self_w=a1/(z + b1 + sub_self);
		
		Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w) );
		
		fileSelf<<setw(25)<<w[j]<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
		fileAkw<<setw(25)<<w[j]<<setw(20)<<Akw;
		
		for (l=1; l<Neta; l++)
		{
			z=dcomplex(w[j],eta_p[l]);
			
			sub_self=pade(z, Nn2, z0, coef);
			self_w=a1/(z + b1 + sub_self);		
			
			Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w) );
			
			fileSelf<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
			fileAkw<<setw(20)<<Akw;
		}
		fileSelf<<endl;
		fileAkw<<endl;
	}
	fileSelf.close();
	fileAkw.close();
	
	delete [] coef;
	delete [] z0;
	delete [] func;		
	
}

//! create a non-uniform frequency index grid, sparsity of the grid increases with m
void green::create_non_uniform_ind_freq(int m_freq_p)
{
	m_freq=m_freq_p;
	int N1=nbw/((int)pow(2,m_freq+1));
	int N2=N1/2;
	
	N_freq=N1+m_freq*N2;
	
	cout<<"non uniform frequency grid parameters:\n";
	cout<<"N1: "<<N1<<endl;
	cout<<"N_freq: "<<N_freq<<endl;
	
	ind_freq=new int[N_freq];
	
	int j, lj, p2lj;
	
	for (j=0; j<N1; j++)
	{
		ind_freq[j]=j;
	//	cout<<setw(10)<<j<<ind_freq[j]<<'\n';
	}
	for (j=N1; j<N_freq; j++)
	{
		lj=(j-N1)/N2;
		p2lj=(int)pow(2,lj);
		ind_freq[j]=p2lj*N1+p2lj*2*((j-N1)%N2);
	//	cout<<setw(10)<<j<<ind_freq[j]<<'\n';
	}
	
}

//! using Pade approximants, obtain energy distribution curves for all k between k0=(kx[0],ky[0]) and kf=(kx[Nk-1],kyf[Nk-1]) defined as integers
//! if k0=NULL or kf=NULL, the spectral weight is calculated for all k in the reduced Brillouin zone
//! w is the energy vector of size nw0, z0 is the Matsubara frequency vector of size NP, eta_p constains N_eta values of the imaginary part of energy
//! fname contains a prefix for the output file name
void green::EDC(int *kx, int *ky, int Nk, double *w, int nw0, int m_freq_p, int NP_max, double *eta_p, int N_eta, char *fname)
{
	if (!Self_kw_array)
	{
		cout<<"EDC(): la self-energie n'a pas ete calculee\n";
		return;
	}
	
	int j, l, m, n;
	
	create_non_uniform_ind_freq(m_freq_p);
	
	pade_z0=new dcomplex[N_freq];
	
	for (j=0; j<N_freq; j++)
	{
		pade_z0[j]=dcomplex(0,(2*ind_freq[j]+1)*PI*tem);
	}
	
	int NP=N_freq;
	if (NP>NP_max) NP=NP_max;
	
	cout<<setprecision(5);
	cout<<"frequence maximale incluse dans le prolongement analytique: "<<(2*ind_freq[NP-1]+1)*PI*tem<<endl;

	char nameFormAkw[200];
	if (fname) strcpy(nameFormAkw,fname);
	else strcpy(nameFormAkw,"spectral_weight");
	strcat(nameFormAkw,"_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d_%1d_%3.1e_%3.1e.dat");
	
	if (!kx || !ky)
	{
		cout<<"calcul du poids spectral dans la zone reduite\n";
		
		int lmin=0;
		int lmax=nbk;
		
#pragma omp parallel private(l,m,j)
		{
			double *spectr_w=new double[nw0*N_eta];
			fstream Akw_file;
			char name[200];
			int k[2];
			pade_coef=new dcomplex[N_freq];
			pade_func=new dcomplex[N_freq];
			
#pragma omp for
			for (l=lmin; l<lmax; l+=2)
			{
				k[0]=l;
				for (m=0; m<=l; m+=2)
				{
					k[1]=m;
					calc_spectral_weight(k, w, nw0, NP, eta_p, N_eta, spectr_w);
					
					sprintf(name, nameFormAkw, U, tp, tpp, t3, density, tem, 2*(nbk-1), nbw, k[0], k[1], w[0], w[nw0-1]);
					Akw_file.open(name,ios::out);
					Akw_file<<setprecision(10)<<setiosflags(ios::left);
					for (j=0; j<nw0; j++)
					{
						Akw_file<<setw(20)<<w[j];
						for (n=0; n<N_eta; n++)
						{
							Akw_file<<setw(20)<<spectr_w[n+j*N_eta];
						}
						Akw_file<<endl;
					}
					Akw_file.close();
				}
			}
			delete [] pade_coef;
			pade_coef=NULL;
			delete [] pade_func;
			pade_func=NULL;
			delete [] spectr_w;
		}
	}
	else
	{
		cout<<"calcul du poids spectral pour les vecteurs d'onde entre k=("<<(PI*kx[0])/(nbk-1)<<","<<(PI*ky[0])/(nbk-1)<<") et ("<<(PI*kx[Nk-1])/(nbk-1)<<","<<(PI*ky[Nk-1])/(nbk-1)<<")\n";
#pragma omp parallel private(l,j)
		{
			double *spectr_w=new double[nw0*N_eta];
			fstream Akw_file;
			char name[200];
			int k[2];
			pade_coef=new dcomplex[N_freq];
			pade_func=new dcomplex[N_freq];
#pragma omp for
			for (l=0; l<Nk; l++)
			{
				k[0]=kx[l];
				k[1]=ky[l];
				
				calc_spectral_weight(k, w, nw0, NP, eta_p, N_eta, spectr_w);
				
				sprintf(name, nameFormAkw, U, tp, tpp, t3, density, tem, 2*(nbk-1), nbw, k[0], k[1], w[0], w[nw0-1]);
				Akw_file.open(name,ios::out);
				Akw_file<<setprecision(10)<<setiosflags(ios::left);
				for (j=0; j<nw0; j++)
				{
					Akw_file<<setw(20)<<w[j];
					for (n=0; n<N_eta; n++)
					{
						Akw_file<<setw(20)<<spectr_w[n+j*N_eta];
					}
					Akw_file<<endl;
				}
				Akw_file.close();
			}
			delete [] pade_coef;
			pade_coef=NULL;
			delete [] pade_func;
			pade_func=NULL;
			delete [] spectr_w;
		}
	}
	
	delete [] pade_z0;
	pade_z0=NULL;
	
}

//! obtain energy distribution curves for all k between k0=(kx[0],ky[0]) and kf=(kx[Nk-1],kyf[Nk-1]) defined as integers
//! if k0=NULL or kf=NULL, the spectral weight is calculated for all k in the reduced Brillouin zone
void green::EDC(int *kx, int *ky, int Nk, char *fname)
{
	if (!Self_kw_array)
	{
		cout<<"EDC(): la self-energie n'a pas ete calculee\n";
		return;
	}

	int j, l, m;

	double Dw=10.0;
	double dw=0.01;
	int nw0=(int)(2*Dw/dw)+1;

	double *w=new double[nw0];

	w[0]=-Dw;
	for (j=1; j<nw0; j++)
	{
		w[j]=w[0]+j*dw;
//		cout<<w[j]<<'\n';
	}

	int N=nbw/2;

	int lmin=0;
	int lmax=nbk;


	if (!kx || !ky)
	{
		cout<<"calcul du poids spectral dans la zone reduite\n";

#pragma omp parallel private(l,m,j)
		{
			double spectr_w[nw0];
			fstream Akw_file;
			char nameFormAkw[100], name[100];
			strcpy(nameFormAkw,fname);
			strcat(nameFormAkw,"_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%3.1e_%1d_%1d_%1d_%1d.dat");
			int k[2];
#pragma omp for
			for (l=lmin; l<lmax; l+=2)
			{
				k[0]=l;
				for (m=0; m<=l; m+=2)
				{
					k[1]=m;
					calc_spectral_weight(k, w, nw0, N, spectr_w);

					sprintf(name, nameFormAkw, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)eta, 2*(nbk-1), nbw, k[0], k[1]);
					Akw_file.open(name,ios::out);
					Akw_file<<setprecision(15)<<setiosflags(ios::left);

					for (j=0; j<nw0; j++)
					{
						Akw_file<<setw(30)<<w[j]<<spectr_w[j]<<'\n';
					}
					Akw_file.close();
				}
			}
		}
	}
	else
	{
		cout<<"calcul du poids spectral pour les vecteurs d'onde donnes\n";
#pragma omp parallel private(l,j)
		{
			double spectr_w[nw0];
			fstream Akw_file;
			char nameFormAkw[100], name[100];
			strcpy(nameFormAkw,fname);
			strcat(nameFormAkw,"_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%3.1e_%1d_%1d_%1d_%1d.dat");
			int k[2];
#pragma omp for
			for (l=0; l<Nk; l++)
			{
				k[0]=kx[l];
				k[1]=ky[l];

				calc_spectral_weight(k, w, nw0, N, spectr_w);

				sprintf(name, nameFormAkw, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)eta, 2*(nbk-1), nbw, k[0], k[1]);
				Akw_file.open(name,ios::out);
				Akw_file<<setprecision(15)<<setiosflags(ios::left);

				for (j=0; j<nw0; j++)
				{
					Akw_file<<setw(30)<<w[j]<<spectr_w[j]<<'\n';
				}
				Akw_file.close();
			}
		}
	}

}

//! obtain momentum distribution curve for the energy w, n is a vector containing the indices of the Matsubara frequencies, 
//! Nn is the number of Matsubara frequencies to be used (3 or 4), 
//! eta is a vector containing different values of the imaginary part of the energy, Neta is the number of those values
void green::MDC(int *n, int Nn, double w, double *eta_p, int Neta, bool use_Pade)
{
	if (!Self_kw_array)
	{
		cout<<"MDC(): la self-energie n'a pas ete calculee\n";
		return;
	}
	
	if (Nn<2) Nn=2;
	
	char name[100];
	
	//	int N=nbw/2;
	
	dcomplex *z0=new dcomplex[Nn];
	dcomplex *self_ikn=new dcomplex[Nn];
	
	long int j, l, m, x, y;
	int nk8th=(nbk*(nbk+1))/2;
	
	double kx, ky, ekmu, Akw;

	double c1=Sigma_inf, c2, c3, a1, b1, a2;
	
	dcomplex z;
	
	clock_t     begin;
	clock_t     total;
	double temps_sec;
	
	if (!use_Pade)
	{
		dcomplex self_w_poly;
		//	dcomplex self_w;
		
		fstream Akw_poly_file, Self_poly_file;
		
		const char *nameFormAkw_poly="Akw_poly_k_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%5.3e_%1d_%1d.dat";
		
		sprintf(name, nameFormAkw_poly, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw);
		Akw_poly_file.open(name,ios::out);
		Akw_poly_file<<setprecision(12)<<setiosflags(ios::left);
		
		const char *nameFormSelf_poly="Self_poly_k_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%5.3e_%1d_%1d.dat";
		
		sprintf(name, nameFormSelf_poly, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw);
		Self_poly_file.open(name,ios::out);
		Self_poly_file<<setprecision(12)<<setiosflags(ios::left);
		
		for (l=0; l<nbk; l++)
		{
			for (m=0; m<=l; m++)
			{
				ekmu=ek[m+l*nbk] - mu;
				
				if (l==0 && m==0) begin=clock();
				
				for (j=0; j<Nn; j++)
				{
					z0[j]=dcomplex(0.0, (2*n[j]+1)*PI*tem);
					z=z0[j];
					self_ikn[j]=Sigma_inf/z+Sigma_inf2[m+(l*(l+1))/2]/(z*z)+Sigma_inf3[m+(l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
				}
				
				z=dcomplex(w,eta_p[0]);
				
				
				if (Nn>=4)
					self_w_poly=self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))*((z-z0[3])/(z0[0]-z0[3]))
					+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))*((z-z0[3])/(z0[1]-z0[3]))
					+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]))*((z-z0[3])/(z0[2]-z0[3]))
					+self_ikn[3]*((z-z0[0])/(z0[3]-z0[0]))*((z-z0[1])/(z0[3]-z0[1]))*((z-z0[2])/(z0[3]-z0[2]));
				else
					self_w_poly=self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))
					+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))
					+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]));
				
				Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w_poly) );
				Self_poly_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<real(self_w_poly)<<setw(25)<<imag(self_w_poly);
				Akw_poly_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<Akw;
				
				for (j=1; j<Neta; j++)
				{
					
					z=dcomplex(w,eta_p[j]);
					if (Nn>=4)
						self_w_poly=self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))*((z-z0[3])/(z0[0]-z0[3]))
						+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))*((z-z0[3])/(z0[1]-z0[3]))
						+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]))*((z-z0[3])/(z0[2]-z0[3]))
						+self_ikn[3]*((z-z0[0])/(z0[3]-z0[0]))*((z-z0[1])/(z0[3]-z0[1]))*((z-z0[2])/(z0[3]-z0[2]));
					else
						self_w_poly=self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))
						+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))
						+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]));
					
					Akw=-imag( ((double)1.0/PI)/(z-ekmu-self_w_poly) );
					Self_poly_file<<setw(20)<<real(self_w_poly)<<setw(25)<<imag(self_w_poly);
					Akw_poly_file<<setw(20)<<Akw;
				}
				Self_poly_file<<endl;
				Akw_poly_file<<endl;
				
				if (l==0 && m==0)
				{
					total = clock() - begin;
					temps_sec=(double) total / CLOCKS_PER_SEC;
					
					//				cout<<"MDC():"<<'\n'<<"temps pour un k (sec):  "<<temps_sec<<'\n';
					//				cout<<"temps pour un k (min):  "<<temps_sec/60<<'\n';
					//				cout<<"temps estime pour tous les k (min):  "<<(nbk*(nbk+1))/2*temps_sec/60<<'\n';
				}
			}
		}
		
		Self_poly_file.close();
		Akw_poly_file.close();
		
		delete [] self_ikn;
		delete [] z0;
	}
	else
	{
		dcomplex self_w_Pade;
		dcomplex *coef_Pade=new dcomplex[Nn];
	//	dcomplex *coef=new dcomplex[Nn];
	//	dcomplex *func=new dcomplex[Nn];
		
		fstream Akw_file, Self_file;
		
	/*
		const char *nameFormAkw="Akw_k_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d_%1d.dat";
		const char *nameFormSelf="Self_k_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameFormAkw, (double)tp, (double)tpp, (double)U, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw, Nn);
		Akw_file.open(name,ios::out);
		Akw_file<<setprecision(12)<<setiosflags(ios::left);
		
		sprintf(name, nameFormSelf, (double)tp, (double)tpp, (double)U, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw, Nn);
		Self_file.open(name,ios::out);
		Self_file<<setprecision(12)<<setiosflags(ios::left);
	*/
		
		fstream Akw_Pade_file, Self_Pade_file;
		
		const char *nameFormAkw_Pade="Akw_Pade_k_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%5.3e_%1d_%1d.dat";
		const char *nameFormSelf_Pade="Self_Pade_k_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%5.3e_%1d_%1d.dat";
		
		sprintf(name, nameFormAkw_Pade, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw);
		Akw_Pade_file.open(name,ios::out);
		Akw_Pade_file<<setprecision(12)<<setiosflags(ios::left);
		
		sprintf(name, nameFormSelf_Pade, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)w, 2*(nbk-1), nbw);
		Self_Pade_file.open(name,ios::out);
		Self_Pade_file<<setprecision(12)<<setiosflags(ios::left);
		
		//	dcomplex fk,fztmp;
		//	dcomplex sub_self;
		
		for (l=0; l<nbk; l++)
		{
			for (m=0; m<=l; m++)
			{
				ekmu=ek[m+l*nbk] - mu;
				
				if (l==0 && m==0) begin=clock();
				/*
				 x=l;
				 y=m;
				 c2=Sigma_inf2[y + (x*(x+1))/2];
				 c3=Sigma_inf3[y + (x*(x+1))/2];
				 a1=c1;
				 b1=-c2/c1;
				 a2=-c3/c1+b1*b1;
				 */
				for (j=0; j<Nn; j++)
				{
					z0[j]=dcomplex(0.0, (2*n[j]+1)*PI*tem);
					z=z0[j];
					self_ikn[j]=Sigma_inf/z+Sigma_inf2[m+(l*(l+1))/2]/(z*z)+Sigma_inf3[m+(l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
					
					//				fk=Self_kw_array[y + (x*(x+1))/2+j*nk8th];
					//				func[j]=(b1*fk + (fk+b1*c3)*z + (c3-a1*b1*b1)*z*z)/(-fk - c3*z + a1*b1*z*z - a1*z*z*z);
				}
				
				
			//	pade_cont_frac_coef_rec(func, z0, Nn, coef);
				pade_cont_frac_coef_rec(self_ikn, z0, Nn, coef_Pade);
				
				z=dcomplex(w,eta_p[0]);
				
				/*
				 sub_self=pade(z, Nn, z0, coef);
				 self_w=a1/(z + b1 + sub_self);
				 Akw=-imag(1.0/(z-ekmu-self_w))/PI;
				 Self_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
				 Akw_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<Akw;
				*/
				
				 self_w_Pade=pade(z, Nn, z0, coef_Pade);
				 Akw=-imag(1.0/(z-ekmu-self_w_Pade))/PI;
				 Self_Pade_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<real(self_w_Pade)<<setw(25)<<imag(self_w_Pade);
				 Akw_Pade_file<<setw(10)<<l<<setw(20)<<m<<setw(20)<<Akw;
				
				for (j=1; j<Neta; j++)
				{
					//				self_w=pade(dcomplex(w,eta_p[j]), Nn, z0, coef);
					
					z=dcomplex(w,eta_p[j]);
					/*
					 sub_self=pade(z, Nn, z0, coef);
					 self_w=a1/(z + b1 + sub_self);
					 Akw=-imag(1.0/(z-ekmu-self_w))/PI;
					 Self_file<<setw(20)<<real(self_w)<<setw(25)<<imag(self_w);
					 Akw_file<<setw(20)<<Akw;
					*/
					
					 self_w_Pade=pade(z, Nn, z0, coef_Pade);
					 Akw=-imag(1.0/(z-ekmu-self_w_Pade))/PI;
					 Self_Pade_file<<setw(20)<<real(self_w_Pade)<<setw(25)<<imag(self_w_Pade);
					 Akw_Pade_file<<setw(20)<<Akw;
					
				}
				//			Self_file<<endl;
				//			Akw_file<<endl;
				
				Self_Pade_file<<endl;
				Akw_Pade_file<<endl;
				
				if (l==0 && m==0)
				{
					total = clock() - begin;
					temps_sec=(double) total / CLOCKS_PER_SEC;
					
					//				cout<<"MDC():"<<'\n'<<"temps pour un k (sec):  "<<temps_sec<<'\n';
					//				cout<<"temps pour un k (min):  "<<temps_sec/60<<'\n';
					//				cout<<"temps estime pour tous les k (min):  "<<(nbk*(nbk+1))/2*temps_sec/60<<'\n';
				}
			}
		}
		
		//	Self_file.close();
		//	Akw_file.close();
		
		Self_Pade_file.close();
		Akw_Pade_file.close();
		
		delete [] self_ikn;
		delete [] coef_Pade;
		delete [] z0;
		//	delete [] coef;
		//	delete [] func;
	}
		
	
}

//! using N points Pade approximants, calculate the spectral weight for a given k=(k_x, k_y) with 0<k_x<nbk-1 and 0<k_y<=k_x for the energy vector w of size nw0
//! the result is stored in the vector spectr_w with the index of eta as the leading one
//! z0 contains the Matsubara frequencies used to compute de Pade parameters, coef and func are size N dcomplex vectors
void green::calc_spectral_weight(int *k, double *w, int nw0, int N, double *eta_p, int N_eta, double *spectr_w)
{
	double wi=w[0], wf=w[nw0-1];
	
	if (wi>wf)
	{
		cout<<"calc_spectral_weight: wi>wf!\n";
		return;
	}
	
	if (N>nbw/2)
	{
		cout<<"calc_spectral_weight: Value N larger than nw="<<nbw/2<<'\n';
		return;
	}
	
	if (N<=0)
	{
		cout<<"calc_spectral_weight: Value N too small: "<<N<<'\n';
		return;
	}
	
	int Nk=2*(nbk-1);
	while (k[0]<0)	k[0]+=Nk;
	while (k[1]<0)	k[1]+=Nk;
	if (k[0]>nbk-1) k[0]=Nk-k[0];
	if (k[1]>nbk-1) k[1]=Nk-k[1];
	
//	if (k[0]<0) k[0]*=-1;
//	if (k[1]<0) k[1]*=-1;
//	if (k[0]>=nbk || k[1]>=nbk)
//	{
//		cout<<"calc_spectral_weight: invalid wave vector k=("<<k[0]<<","<<k[1]<<")"<<'\n';
//		return;
//	}
	
	if (!Self_kw_array)
	{
		cout<<"calc_spectral_weight: Self-energy not computed"<<'\n';
		return;
	}
	
	int nk8th=(nbk*(nbk+1))/2;
	
	int x=k[0];
	int y=k[1];
	if (y>x)
	{
		y=x;
		x=k[1];
	}

	int j, n;
	dcomplex z;
	for (n=0; n<N; n++)
	{
		z=pade_z0[n];
		pade_func[n]=Sigma_inf/z+Sigma_inf2[y + (x*(x+1))/2]/(z*z)+Sigma_inf3[y + (x*(x+1))/2]/(z*z*z)+Self_kw_array[y + (x*(x+1))/2+n*nk8th]/(z*z*z*z);
	}

	dcomplex self_w;
	
	double ekmu=ek_val(k[0], k[1]) - mu;
	pade_cont_frac_coef_rec(pade_func, pade_z0, N, pade_coef);
	
	for (n=0; n<nw0; n++)
	{
		for (j=0; j<N_eta; j++)
		{
			self_w=pade(dcomplex(w[n], eta_p[j]), N, pade_z0, pade_coef);
		
			spectr_w[j+n*N_eta]=-imag(((double)1.0)/(dcomplex(w[n]-ekmu,eta_p[j])-self_w))/PI;
		}
	}
	
}

//! calculate the spectral weight for a given k=(k_x, k_y) with 0<k_x<nbk-1 and 0<k_y<=k_x for the energy vector w of size nw0
//using the N points Pade approximant, the result is stored in the vector spectr_w
void green::calc_spectral_weight(int *k, double *w, int nw0, int N, double *spectr_w)
{
	double wi=w[0], wf=w[nw0-1];

	if (wi>wf)
	{
		cout<<"calc_spectral_weight: wi>wf!\n";
		return;
	}

	if (N>nbw/2)
	{
		cout<<"calc_spectral_weight: Value N larger than nw="<<nbw/2<<'\n';
		return;
	}

	if (N<=0)
	{
		cout<<"calc_spectral_weight: Value N too small: "<<N<<'\n';
		return;
	}

	if (k[0]<0) k[0]*=-1;
	if (k[1]<0) k[1]*=-1;
	if (k[0]>=nbk || k[1]>=nbk)
	{
		cout<<"calc_spectral_weight: invalid wave vector k=("<<k[0]<<","<<k[1]<<")"<<'\n';
		return;
	}

	if (!Self_kw_array)
	{
		cout<<"calc_spectral_weight: Self-energy not computed"<<'\n';
		return;
	}

//	if (!k[0] && !k[1]) cout<<"NP: "<<N<<endl;
	
	ptrdiff_t w_local = N;
	ptrdiff_t w_start = 0;
	int nk8th=(nbk*(nbk+1))/2;

	dcomplex *coef=new dcomplex[N];
	dcomplex *z0=new dcomplex[N];
	dcomplex *func;
	if (w_local>0) func=new dcomplex[w_local];

	long int n;
	long int x=k[0];
	long int y=k[1];
	if (y>x)
	{
		y=x;
		x=k[1];
	}

//#pragma omp parallel for private(n)
	for (n=0; n<N; n++)
	{
		z0[n]=dcomplex(0.0, (2*n+1)*PI*tem);
	}
//#pragma omp parallel private(n)
	{
		dcomplex z;
//#pragma omp for
		for (n=0; n<w_local; n++)
		{
			z=z0[n];
	//		func[n]=Self_kw_array[y + (x*(x+1))/2 + n*nk8th];
			func[n]=Sigma_inf/z+Sigma_inf2[y + (x*(x+1))/2]/(z*z)+Sigma_inf3[y + (x*(x+1))/2]/(z*z*z)+Self_kw_array[y + (x*(x+1))/2+n*nk8th]/(z*z*z*z);
		}
	}

	double kx=k[0]*PI/(nbk-1),  ky=k[1]*PI/(nbk-1);
//	double ekmu=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu;
	double ekmu=ek_val(k[0], k[1]) - mu;	

	pade_cont_frac_coef_rec(func, z0, (int)(w_local+w_start), coef);
	if (w_local>0) delete [] func;

	w_local = nw0;
	dcomplex self_w;
	double *spec_ptr=spectr_w;

//	fstream file;
//	const char nameForm[]="spectr_wgt_U%4.2f_tp%5.3f_tpp%5.3f_tt%5.3f_dens%7.5f_tem%6.4f_Nk%1d_Nwn%1d_%1d_%1d.dat";
//	char name[200];
	
//	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw, k[0], k[1]);
//	file.open(name,ios::out);
//	file<<setiosflags(ios::left)<<setprecision(10);
	dcomplex A, B;
//#pragma omp parallel for private(n,self_w)
	for (n=0; n<w_local; n++)
	{
		self_w=pade(dcomplex(w[n],eta), N, z0, coef);

		spec_ptr[n]=-imag(((double)1.0)/(dcomplex(w[n]-ekmu,eta)-self_w))/PI;
		
//		file<<setw(25)<<w[n]<<spec_ptr[n]<<endl;
	}
//	file.close();
	
	delete [] coef;
	delete [] z0;
	
/*	
	cout<<" Spectre pour le vecteur ("<<k[0]<<","<<k[1]<<")\n";
	n = ( nw0/2-8<0 ? 0 : nw0/2-8 );
	for (; n<nw0/2+8&&n<nw0; n++)
	{
		cout<<" "<<w[n]<<" "<<spectr_w[n]<<'\n';
	}
*/
}


//! calculate the self-energy with the "brute force" method for a given k=(kx,ky) and frequency index n (kn=(2n+1)*PI*tem)
dcomplex green::Self_kw(double k[], long int n)
{
	if (Usp==0.0 || Uch==0.0)
	{
		cout<<"calc_Self(): Usp et Uch doivent etre calcules d'abord\n";
		return 0.0;
	}
	if (k[0]>nbk-1 || k[1]>nbk-1 || k[0]<0 || k[1]<0)
	{
		cout<<"Self_kw():  k[0] (kx) et k[1] (ky) doivent etre entre 0 et nk0/2\n";
		return 0.0;
	}

	int nk0=2*(nbk-1);
	int Nk=nk0*nk0;
	int nk8th=(nbk*(nbk+1))/2;

	long int j, l, m, p, q, tmp;
	long int  x, y, wn;
	
	double kx=k[0], ky=k[1];
	double qn, kn=(2*n+1)*PI*tem;
	double qx, qy, ekq;
	double chi0tmp, chich, chisp;
	dcomplex G;
	
	dcomplex selfkw=0.0;

	double wt;
	
	for (j=0; j<2*nbw-1; j++)
	{
		if (j<nbw)
		{
			qn=2.0*j*PI*tem;
			wn=j;
		}
		else
		{
			qn=2.0*(j-2*nbw+1)*PI*tem;
			wn=2*nbw-j-1;
		}

		for (l=0; l<nk0; l++)
		{
			qx=l*2.0*PI/nk0;			
			for (m=0; m<nk0; m++)
			{
				qy=m*2.0*PI/nk0;

//				ekq=-2.0*(cos(kx+qx) + cos(ky+qy)) - 4.0*tp*cos(kx+qx)*cos(ky+qy) - 2.0*tpp*(cos(2*(kx+qx)) + cos(2*(ky+qy)));
				ekq=dispk(kx+qx, ky+qy);
				
				x=l;
				if (x>nbk-1) x=nk0-x;
				y=m;
				if (y>nbk-1) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}

				if (wn<=nbw/2)
					chi0tmp=chiqw_array[y+(x*(x+1))/2+wn*nk8th];
				else
					chi0tmp=chi0_inf[y+(x*(x+1))/2]/(qn*qn)+chi0_inf2[y+(x*(x+1))/2]/(qn*qn*qn*qn);
//					chi0tmp=chi_0->chi0_inf[y+(x*(x+1))/2]/(qn*qn)+chi_0->chi0_inf2[y+(x*(x+1))/2]/(qn*qn*qn*qn);
				chisp=chi0tmp/(1.0-Usp*chi0tmp/2.0);
				chich=chi0tmp/(1.0+Uch*chi0tmp/2.0);

				G=((double)1.0)/(dcomplex(0.0, kn+qn)-ekq+mu0);

				selfkw+=(U/8.0)*(tem/Nk)*(3.0*Usp*chisp+Uch*chich)*G;
			}
		}
	}
	

/*	
	for (j=0; j<nbw; j++)
	{
		if (j<nbw/2)
		{
			qn=2.0*j*PI*tem;
			wn=j;
		}
		else
		{
			qn=2.0*(j-nbw)*PI*tem;
			wn=nbw-j;
		}

		for (l=0; l<nk0; l++)
		{
			qx=l*2.0*PI/nk0;

			for (m=0; m<nk0; m++)
			{
				qy=m*2.0*PI/nk0;

				ekq=-2.0*(cos(kx+qx) + cos(ky+qy)) - 4.0*tp*cos(kx+qx)*cos(ky+qy) - 2.0*tpp*(cos(2*(kx+qx)) + cos(2*(ky+qy)));

				x=l;
				if (x>nbk-1) x=nk0-x;
				y=m;
				if (y>nbk-1) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}

				chi0tmp=chiqw_array[y+(x*(x+1))/2+wn*nk8th];
				chisp=chi0tmp/(1.0-Usp*chi0tmp/2.0);
				chich=chi0tmp/(1.0+Uch*chi0tmp/2.0);

				G=((double)1.0)/(dcomplex(0.0, kn+qn)-ekq+mu0);

				selfkw+=(U/8.0)*(tem/Nk)*(3.0*Usp*chisp+Uch*chich)*G;
			}
		}
	}	
*/
	
	return selfkw;
}

//! calculate the self-energy with the "brute force" method for a given k=(kx,ky) and frequency index n (kn=(2n+1)*PI*tem)
dcomplex green::Self_kw(long int k[], long int n)
{
	if (Usp==0.0 || Uch==0.0)
	{
		cout<<"calc_Self(): Usp et Uch doivent etre calcules d'abord\n";
		return 0.0;
	}
	if (k[0]>nbk-1 || k[1]>nbk-1 || k[0]<0 || k[1]<0)
	{
		cout<<"Self_kw():  k[0] (kx) et k[1] (ky) doivent etre entre 0 et nk0/2\n";
		return 0.0;
	}

	int nk0=2*(nbk-1);
	int Nk=nk0*nk0;
	int nk8th=(nbk*(nbk+1))/2;

	long int j, l, m, p, q, tmp; 
	long int x, y, wn;
	long int kx=k[0], ky=k[1];
	
	double qn, kn=(2*n+1)*PI*tem;
	double qx, qy, ekq;
	double chi0tmp, chich, chisp;
	dcomplex G;
	
	dcomplex selfkw=0.0;
	
	double wt;
	
	for (j=0; j<2*nbw-1; j++)
	{
		if (j<nbw)
		{
			qn=2.0*j*PI*tem;
			wn=j;
		}
		else
		{
			qn=2.0*(j-2*nbw+1)*PI*tem;
			wn=2*nbw-j-1;
		}

		for (l=0; l<nk0; l++)
		{
//			qx=l*2.0*PI/nk0;

			p=l+kx;
			if (p>nk0) p=p-nk0;
			if (p>nbk-1) p=nk0-p;
			
			for (m=0; m<nk0; m++)
			{
//				qy=m*2.0*PI/nk0;

//				ekq0=-2.0*(cos(kx2+qx) + cos(ky2+qy)) - 4.0*tp*cos(kx2+qx)*cos(ky2+qy) - 2.0*tpp*(cos(2*(kx2+qx)) + cos(2*(ky2+qy)));
				
				q=m+ky;
				if (q>nk0) q=q-nk0;
				if (q>nbk-1) q=nk0-q;
				
				ekq=ek[q+p*nbk];
				
//				if (j==0) cout<<ekq-ekq0<<'\n';
//				ekq=ekq0;
				
				x=l;
				if (x>nbk-1) x=nk0-x;
				y=m;
				if (y>nbk-1) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}

				if (wn<=nbw/2)
					chi0tmp=chiqw_array[y+(x*(x+1))/2+wn*nk8th];
				else
					chi0tmp=chi0_inf[y+(x*(x+1))/2]/(qn*qn)+chi0_inf2[y+(x*(x+1))/2]/(qn*qn*qn*qn);
//					chi0tmp=chi_0->chi0_inf[y+(x*(x+1))/2]/(qn*qn)+chi_0->chi0_inf2[y+(x*(x+1))/2]/(qn*qn*qn*qn);
				chisp=chi0tmp/(1.0-Usp*chi0tmp/2.0);
				chich=chi0tmp/(1.0+Uch*chi0tmp/2.0);

				G=((double)1.0)/(dcomplex(0.0, kn+qn)-ekq+mu0);

				selfkw+=(U/8.0)*(tem/Nk)*(3.0*Usp*chisp+Uch*chich)*G;
			}
		}
	}	
	
	return selfkw;
}

void green::save_Self_bin()
{
	if (!Self_kw_array)
	{
		cout<<"Self_kw_array n'existe pas\n";
		return;
	}
	
	int nk0=2*(nbk-1);
	long int nk8th=(long int)(nbk*(nbk+1))/2;
	
	const char *nameForm="./self_bin_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	fstream file;
	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&Self_kw_array[0],  (nbw/2+1)*nk8th*sizeof(dcomplex));
	file.close();
}

void green::load_Self_bin()
{
	int nk0=2*(nbk-1);
	long int nk8th=(long int)(nbk*(nbk+1))/2;
	
	if (Self_kw_array) delete [] Self_kw_array;
	
	
	Self_kw_array=new dcomplex[ (nbw/2+1)*nk8th];
	
	const char *nameForm="./self_bin_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	fstream file;
	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::in | ios::binary);
	if (file)
	{
		file.read((char*)&Self_kw_array[0], (nbw/2+1)*nk8th*sizeof(dcomplex));
		file.close();
	}
	else
		cout<<"fichier  "<<name<<"  inexistant\n";
}

//! save the self-energy in ascii format 
//!format: freq. index		kx index		ky index		real(Sigma)		imag(Sigma)
void green::save_Self(int Nfreq)
{
	if (!Self_kw_array)
	{
		cout<<"save_Self(): Self_kw_array n'existe pas\n";
		return;
	}
	
	if (Nfreq>nbw/2 || Nfreq==0) Nfreq=nbw/2;
	
	int nk0=2*(nbk-1);
	int nk8th=(nbk*(nbk+1))/2;
	
	const char *nameForm="./self_tpsc_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[200];
	
	fstream file;
	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out);
	
	long int j, l, m;
	file<<setiosflags(ios::left)<<setiosflags(ios::scientific)<<setprecision(16);
	for (j=0; j<Nfreq; j++)
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
				file<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(25)<<real(Self_kw_array[m + (l*(l+1))/2 + j*nk8th])
						<<imag(Self_kw_array[m + (l*(l+1))/2 + j*nk8th])<<'\n';
	
	file.close();
}

//! save the Matsubara Green's function in ascii format 
//!format: freq. index		real(G)		imag(G)
void green::save_Green(int *k)
{
	if (!Self_kw_array)
	{
		cout<<"save_Green(): Self_kw_array n'existe pas\n";
		return;
	}
	
	int nk0=2*(nbk-1);
	int nk8th=(nbk*(nbk+1))/2;
	
	long int j, l, m, tmp;
	l=k[0];
	m=k[1];
	if (l<0) l=-l;
	if (m<0) m=-m;
	if (l>nk0 || m>nk0)
	{
		cout<<"save_Green(): k[0] et k[1] doivent etre entre 0 et nk0\n";
		return;
	}
	
	const char *nameForm="./Green_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d_%1d.dat";
	char name[200];
	
	dcomplex z, selftmp, Gtmp;
	
	
	if (l>nbk-1) l=nk0-l;
	if (m>nbk-1) m=nk0-m;
	if (m>l)
	{
		tmp=m;
		m=l;
		l=tmp;
	}
	
	double kx, ky;
	kx=l*PI/(nbk-1);
	ky=m*PI/(nbk-1);
//	double ektmp=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu;
	double ektmp=ek_val(k[0], k[1]) - mu;	
	
	fstream file;
	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, l, m);
	file.open(name, ios::out);
	
	file<<setiosflags(ios::left)<<setprecision(16);
	for (j=0; j<nbw/2; j++)
	{
		z=dcomplex(0, (2*j+1)*PI*tem);
		selftmp=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
		Gtmp=((double)1.0)/(z-ektmp-selftmp);
		file<<setw(10)<<j<<setw(25)<<real(Gtmp)<<imag(Gtmp)<<'\n';
	}
	
	file.close();
}

//!load the Matsubara self-energy saved in ascii format
void green::load_Self()
{
	int nw=nbw/2+1;
	int nk0=2*(nbk-1);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (Self_kw_array) delete [] Self_kw_array;
	Self_kw_array=new dcomplex[nw*nk8th];
	
	const char *nameForm="./self_tpsc_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	fstream file;
	sprintf(name, nameForm , (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::in);
	if (!file)
	{
		cout<<"load_Self():  fichier  "<<name<<"  inexistant\n";
		return;
	}
	
    double tmpRe, tmpIm;
	long int j, l, m;
	for (j=0; j<nw-1; j++)
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
            {
                file>>j>>l>>m>>tmpRe>>tmpIm;
                Self_kw_array[m + (l*(l+1))/2 + j*nk8th]=dcomplex(tmpRe,tmpIm);
//                file>>j>>l>>m>>Self_kw_array[m + (l*(l+1))/2 + j*nk8th].real()>>Self_kw_array[m + (l*(l+1))/2 + j*nk8th].imag();
            }
	
	file.close();
}

//calculate the self-energy using FFTs and cubic spline
void green::calc_Self_FFT_spline()
{
	bool compare_force_brute=false;
	
	keep_chirt=true;
	find_Usp_Uch();
//	Usp=U;
//	Uch=U;
	sp_corr_length();
	calc_chirt();	

	cout<<"calcul de la self-energie par calc_Self_FFT_spline()\n";

	double *ekmu, *Gkt_tmp;
	
	double *Grt_tmp;

	int nk=nbq;
//	int nk=chi1->nbq;
//	nbk=nk;
	int nk2=nk*nk;
	int nk0=2*(nk-1);
	int nk02=nk0*nk0;
	int nk8th=(nk*(nk+1))/2;
	
	int nw=nwmax+1;
//	int nw=chi1->nwmax+1;
	nbw=2*(nw-1);

	Grt_tmp=new double[nk2];
	ekmu=new double[nk2];

	long int array_size=(long int)nw*nk8th;
	
	cout<<"nombre d'elements dans la matrice de la self-energie (complexes):  "<<array_size<<'\n';
	cout<<"taille de la matrice de la self-energie (MO):  "<<array_size*32/pow(2.0,20)<<'\n';
	cout<<"taille de la matrice de la self-energie (GO):  "<<array_size*32/pow(2.0,30)<<'\n';
	
	if (Self_kw_array)	delete [] Self_kw_array;
//	Self_kw_array=new dcomplex[nw*nk8th];
	Self_kw_array=new dcomplex[array_size];
	
	unsigned flag = FFTW_MEASURE;
	
	fftw_plan fftplan_r, fftplan_t;
	
	fftplan_r=fftw_plan_r2r_2d(nk, nk, Grt_tmp, Grt_tmp, FFTW_REDFT00, FFTW_REDFT00, flag);

	long int j,l,m;

	double kx, ky, t;

//#pragma omp parallel for private(l,m,kx,ky)
	for (l=0; l<nk; l++)
	{
//		kx=l*PI/(nk-1);
		for(m=0; m<nk; m++)
		{
//			ky=m*PI/(nk-1);
//			ekmu[m + l*nk]=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky)) - mu0;
			ekmu[m + l*nk]=ek[m+l*nk]-mu0;
		}
	}

	double *vertex_array = vertex_rt_array;

	double TPSC_fact=U*(3.0*Usp + Uch)/8.0;
	double Gtmp, chitmp, dG, dchi;
	double *dtau_Sigma_qt=new double[2*nk8th];

// calcul de d/dtau Sigma(q,tau) a tau=0 et tau=beta
//tau=0	
		for (l=0; l<nk; l++)
			for(m=0; m<nk; m++)
					Grt_tmp[m + l*nk]=-1.0/(exp(-ekmu[m + l*nk]/tem)+1.0)/nk02;

		fftw_execute(fftplan_r);	
	
//#pragma omp parallel for private(l,m)
		for (l=0; l<nk; l++)
		{
			for(m=0; m<l; m++)
			{
				Gtmp=Grt_tmp[m + l*nk];
				dchi=TPSC_fact*dtau_chirt0[m + (l*(l+1))/2];
				dG=dtau_Grt0p[m + (l*(l+1))/2];				
//				dchi=TPSC_fact*chi_0->dtau_chirt0[m + (l*(l+1))/2];
//				dG=chi_0->dtau_Grt0p[m + (l*(l+1))/2];
				Grt_tmp[m + l*nk]=dchi*Gtmp+vertex_array[m + (l*(l+1))/2]*dG;
				Grt_tmp[l + m*nk]=Grt_tmp[m + l*nk];
			}
			m=l;
			Gtmp=Grt_tmp[m + l*nk];
			dchi=TPSC_fact*dtau_chirt0[m + (l*(l+1))/2];
			dG=dtau_Grt0p[m + (l*(l+1))/2];			
//			dchi=TPSC_fact*chi_0->dtau_chirt0[m + (l*(l+1))/2];
//			dG=chi_0->dtau_Grt0p[m + (l*(l+1))/2];
			Grt_tmp[m + l*nk]=dchi*Gtmp+vertex_array[m + (l*(l+1))/2]*dG;
		}

		fftw_execute(fftplan_r);
	
		for (l=0; l<nk; l++)
			for(m=0; m<=l; m++)
				dtau_Sigma_qt[m + (l*(l+1))/2]=Grt_tmp[m + l*nk];
		
//tau=beta		
		for (l=0; l<nk; l++)
			for(m=0; m<nk; m++)
					Grt_tmp[m + l*nk]=-1.0/(exp(ekmu[m + l*nk]/tem)+1.0)/nk02;

		fftw_execute(fftplan_r);		
		
		double *xi_r=new double[nk8th];
		xi_r[0]=-mu0;
		for (l=1; l<nk8th; l++) xi_r[l]=er[l];

 
//#pragma omp parallel for private(l,m)		
		for (l=0; l<nk; l++)
		{
			for(m=0; m<l; m++)
			{
				Gtmp=Grt_tmp[m + l*nk];
				dchi=-TPSC_fact*dtau_chirt0[m + (l*(l+1))/2];
				dG=-dtau_Grt0p[m + (l*(l+1))/2]+xi_r[m + (l*(l+1))/2];				
//				dchi=-TPSC_fact*chi_0->dtau_chirt0[m + (l*(l+1))/2];
//				dG=-chi_0->dtau_Grt0p[m + (l*(l+1))/2]+xi_r[m + (l*(l+1))/2];
				Grt_tmp[m + l*nk]=dchi*Gtmp+vertex_array[m + (l*(l+1))/2]*dG;
				Grt_tmp[l + m*nk]=Grt_tmp[m + l*nk];
			}
			m=l;
			Gtmp=Grt_tmp[m + l*nk];
			dchi=-TPSC_fact*dtau_chirt0[m + (l*(l+1))/2];
			dG=-dtau_Grt0p[m + (l*(l+1))/2]+xi_r[m + (l*(l+1))/2];			
//			dchi=-TPSC_fact*chi_0->dtau_chirt0[m + (l*(l+1))/2];
//			dG=-chi_0->dtau_Grt0p[m + (l*(l+1))/2]+xi_r[m + (l*(l+1))/2];
			Grt_tmp[m + l*nk]=dchi*Gtmp+vertex_array[m + (l*(l+1))/2]*dG;
		}

		fftw_execute(fftplan_r);
	
		for (l=0; l<nk; l++)
			for(m=0; m<=l; m++)
				dtau_Sigma_qt[m + (l*(l+1))/2+nk8th]=Grt_tmp[m + l*nk];
	
	delete [] xi_r;	
	delete [] Grt_tmp;
	
#pragma omp parallel private(j,t,l,m,Grt_tmp)	
	{
		Grt_tmp=new double[nk2];
#pragma omp for		
		for (j=0; j<=nbw; j++)
		{				
			t=j/(nbw*tem);
			
			for (l=0; l<nk; l++)
			{
				for(m=0; m<nk; m++)
				{
					if (ekmu[m + l*nk]>0)
						Grt_tmp[m + l*nk]=-exp(-t*ekmu[m + l*nk])/(exp(-ekmu[m + l*nk]/tem)+1.0)/nk02;
					else
						Grt_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + l*nk])/(exp(ekmu[m + l*nk]/tem)+1.0)/nk02;
				}
			}

			fftw_execute_r2r(fftplan_r,Grt_tmp,Grt_tmp);
			
			if (j<nw)
			{
				for (l=0; l<nk; l++)
				{
					for(m=0; m<l; m++)
					{
						Grt_tmp[m + l*nk]*=vertex_array[m + (l*(l+1))/2 + j*nk8th];
						Grt_tmp[l + m*nk]=Grt_tmp[m + l*nk];
					}
					m=l;
					Grt_tmp[m + l*nk]=vertex_array[m + (l*(l+1))/2 + j*nk8th]*Grt_tmp[m + l*nk];
				}
			}
			else
			{
				for (l=0; l<nk; l++)
				{
					for(m=0; m<l; m++)
					{
						Grt_tmp[m + l*nk] *= vertex_array[m + (l*(l+1))/2 + (nbw-j)*nk8th];
						Grt_tmp[l + m*nk]=Grt_tmp[m + l*nk];
					}
					m=l;
					Grt_tmp[m + l*nk] *= vertex_array[m + (l*(l+1))/2 + (nbw-j)*nk8th];
				}
			}

			fftw_execute_r2r(fftplan_r,Grt_tmp,Grt_tmp);
			
			if (j<nw)
			{
				for (l=0; l<nk; l++)
					for(m=0; m<=l; m++)
                        Self_kw_array[m + (l*(l+1))/2 + j*nk8th]=dcomplex(Grt_tmp[m + l*nk], Self_kw_array[m + (l*(l+1))/2 + j*nk8th].imag());
//							Self_kw_array[m + (l*(l+1))/2 + j*nk8th].real()=Grt_tmp[m + l*nk];
			}
			else
			{
				for (l=0; l<nk; l++)
					for(m=0; m<=l; m++)
                        Self_kw_array[m + (l*(l+1))/2 + (j-nw)*nk8th]=dcomplex(Self_kw_array[m + (l*(l+1))/2 + (j-nw)*nk8th].real(),Grt_tmp[m + l*nk]);
//							Self_kw_array[m + (l*(l+1))/2 + (j-nw)*nk8th].imag()=Grt_tmp[m + l*nk];
			}		
//			for (l=0; l<nk; l++)
//				for(m=0; m<=l; m++)
//						Self_kw_array[m + (l*(l+1))/2 + j*nk8th]=Grt_tmp[m + l*nk];
				
		}
		delete [] Grt_tmp;
	}
		
	delete [] ekmu;

	fftw_destroy_plan(fftplan_r);
	
	
// TF avec spline sur tau	
	dcomplex *selftmp=new dcomplex[nbw];
	
//	cout<<"nbw: "<<nbw<<endl;

	flag = FFTW_MEASURE;
	fftplan_t=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (selftmp), reinterpret_cast<fftw_complex*> (selftmp), FFTW_BACKWARD, flag);
	
	delete [] selftmp;
	
	double *coeffs;
	double *tau;
	double *selfqt_tmp;
	double wn;
	double d2S0, d2Sbeta;
	double coeff1, coeff2, coeff3;
	Sigma_inf2=new double[nk8th];
	Sigma_inf3=new double[nk8th];
	Sigma_inf4=new double[nk8th];
	
	int NS0=nbw+1;
	
	if (Self_array_tot)	delete [] Self_array_tot;
//	Self_array_tot=new dcomplex[nw*nk8th];
	
	cout<<setiosflags(ios::left)<<setprecision(15);
#pragma omp parallel private(l,j,tau,selfqt_tmp,coeffs,selftmp,d2S0,d2Sbeta,coeff1,coeff2,coeff3,wn)
	{
		tau=new double[NS0];
		selfqt_tmp=new double[NS0];
		coeffs=new double[4*nbw];
		selftmp=new dcomplex[nbw];
		dcomplex tmp;
		double *FP=new double[2];

#pragma omp for
		for (l=0; l<nk8th; l++)
			{							
				for (j=0; j<nw; j++)
				{
					tau[j]=j/(nbw*tem);					
					selfqt_tmp[j]=real(Self_kw_array[l + j*nk8th]);
				}
				for (j=nw; j<=nbw; j++)
				{
					tau[j]=j/(nbw*tem);
					selfqt_tmp[j]=imag(Self_kw_array[l + (j-nw)*nk8th]);
				}
				
				coeffs[0]=dtau_Sigma_qt[l];
				coeffs[1]=dtau_Sigma_qt[l+nk8th];				
				spline_coeffs_rel(tau, selfqt_tmp, NS0, coeffs);
				
				d2S0=2*coeffs[1];
				d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];
				
				for (j=0; j<nbw; j++)
				{
					selftmp[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs[4*j];
				}
				
				fftw_execute_dft(fftplan_t,reinterpret_cast<fftw_complex*>(selftmp),reinterpret_cast<fftw_complex*>(selftmp));
/*
				FP[0]=dtau_Sigma_qt[l];
				FP[1]=dtau_Sigma_qt[l+nk8th];
				spline_coefficients(coeffs,tau,selfqt_tmp, FP, NS0);
				
				d2S0=2*coeffs[2];
				d2Sbeta=6.0*coeffs[4*nbw-1]/(nbw*tem)+2.0*coeffs[4*nbw-2];
				
				for (j=0; j<nbw; j++)
				{
					selftmp[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs[4*j+3];
				}
				
				fftw_execute_dft(fftplan_t,reinterpret_cast<fftw_complex*>(selftmp),reinterpret_cast<fftw_complex*>(selftmp));
*/
/*			    
				coeffs[0]=dtau_Sigma_qt[l];
				coeffs[1]=dtau_Sigma_qt[l+nk8th];				
				spline_coeffs(tau, selfqt_tmp, NS0, coeffs);
				
				for (j=0; j<nbw; j++)
				{
					tmp=exp(dcomplex(0,j*PI/nbw))*6.0*coeffs[4*j];
					selftmp[j]=dcomplex(tmp.real(), tmp.imag());
				}
				
				fftw_execute_dft(fftplan_t,reinterpret_cast<fftw_complex*>(selftmp),reinterpret_cast<fftw_complex*>(selftmp));
								
				d2S0=2*coeffs[1];
				d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
*/						
				coeff1=selfqt_tmp[0]+selfqt_tmp[nbw];
				if (l==0)	Sigma_inf=-coeff1;
				coeff2=dtau_Sigma_qt[l]+dtau_Sigma_qt[l+nk8th];
				Sigma_inf2[l]=coeff2;
				coeff3=d2S0+d2Sbeta;
				Sigma_inf3[l]=-coeff3;

				for (j=0; j<nbw/2; j++)				
				{
					wn=(2*j+1)*PI*tem;
					
//					Self_array_tot[l + j*nk8th]=I*coeff1/wn - coeff2/(wn*wn) - I*coeff3/(wn*wn*wn)
//										+(1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*selftmp[j]/(wn*wn*wn*wn);
					
					Self_kw_array[l + j*nk8th]=((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*selftmp[j];
					
//					tmp.real()=real(selftmp[j]);
//					tmp.imag()=imag(selftmp[j]);
//					Self_kw_array[l + j*nk8th]=(1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*tmp;					
				}
				
				j=nbw/2-1;
//				tmp.real()=real(selftmp[j]);
//				tmp.imag()=imag(selftmp[j]);
//				Sigma_inf4[l]=real((1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*tmp);
				Sigma_inf4[l]=real( ((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*selftmp[j]);
			}
		delete [] tau;
		delete [] selfqt_tmp;
		delete [] coeffs;
		delete [] selftmp;
		delete [] FP;
//		delete [] Selfw4;
	}


	delete [] dtau_Sigma_qt;
	fftw_destroy_plan(fftplan_t);

	cout<<"self-energie calculee\n";
//	cout<<"Sigma_inf:  "<<Sigma_inf<<'\n';

	if (chiqw_array && compare_force_brute)	
	{
		dcomplex selftmp2, selftmp3;
		dcomplex z(0, PI*tem);
		
		j=0;
		l=(nk-1)/2;
		m=(nk-1)/2;
		long int ktmp[]={l,m};
	
		selftmp2=Self_kw(ktmp, j);
		selftmp3=(Sigma_inf+(Sigma_inf2[m + (l*(l+1))/2] + (Sigma_inf3[m + (l*(l+1))/2]+Self_kw_array[m + (l*(l+1))/2 + j*nk8th]/z)/z)/z)/z;
//		selftmp3=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z) + Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m + (l*(l+1))/2 + j*nk8th]/(z*z*z*z);
//		selftmp3=Self_array_tot[m + (l*(l+1))/2 + j*nk8th];
		
		cout<<"valeur a q=("<<l*PI/(nk-1)<<", "<<m*PI/(nk-1)<<") et wn=pi*tem:  "<<selftmp3<<'\n';
		cout<<"valeur ""force brute"" a q=("<<l*PI/(nk-1)<<", "<<m*PI/(nk-1)<<") et wn=pi*tem:  "<<selftmp2<<'\n';
		cout<<"difference:  "<<selftmp3-selftmp2<<'\n';
	}


/*	
	char name[300];
	fstream file, file3, file4, file5, file6, file7;

	dcomplex *Self2=new dcomplex[nbw*nk8th];
	for (l=0; l<nw*nk8th; l++)
		Self2[l]=Self_kw_array[l];	
	
//	save_Self_bin();	
	load_Self_bin();
	
	cout<<setiosflags(ios::left)<<setprecision(16);
	double diff=0;
	dcomplex diffl;
	for (l=0; l<nbw*nk8th; l++)
	{
		diffl=Self2[l]-Self_kw_array[l];
		diff+=sqrt(norm(diffl));
		if (norm(diffl) && l<20)
		{
			cout<<setw(10)<<l<<setw(45)<<Self2[l]<<setw(45)<<Self_kw_array[l]<<sqrt(norm(diffl))/sqrt(norm(Self2[l]))<<'\n';
		}
	}

	cout<<"diff/(nbw*nk8th):  "<<diff/(nbw*nk8th)<<'\n';
	
	delete [] Self2;
*/	
	
//	double nsr=chi1->nsr;
//	double Sigma_inf0=0.5*U*(0.25*(3*Usp+Uch)*nsr-0.25*Uch*nsr*nsr)+Usp*(Usp-Uch)*nsr*nsr/16.0;
//	cout<<'\n'<<"Sigma_inf:   "<<Sigma_inf0<<'\n';	
//	j=nbw-1;
//	l=0;
//	for (l=0; l<nk8th; l+=10)
//			{
//				wn=(2*j+1)*PI*tem;
//				cout<<setw(20)<<Sigma_inf<<I*wn*Self_kw_array[l + j*nk8th]<<'\n';
//			}


/*
// comparer avec le calcul force brute (Self_kw(k[],n))

	dcomplex selftmp2, selftmp3, z;
	int ktmp[2];
	
	l=0;
	m=0;	
	
//	for (l=0; l<nk; l++)
	{
		ktmp[0]=l;
		
//		m=l;
//		for(m=0; m<nk; m++)
		{
			ktmp[1]=m;
		
			for (j=0; j<10; j++)
			{
				z=dcomplex(0, (2*j+1)*PI*tem);
				selftmp2=Self_kw(ktmp, j);
//				selftmp3=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
				selftmp3=Self_array_tot[m + (l*(l+1))/2+j*nk8th];
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(50)
							<<selftmp3<<setw(50)<<selftmp2
							<<setw(25)<<real(selftmp3-selftmp2)/real(selftmp3)<<imag(selftmp3-selftmp2)/imag(selftmp3)<<'\n';
			}
			for (j=nbw/2-10; j<nbw/2; j++)
			{
				z=dcomplex(0, (2*j+1)*PI*tem);
				selftmp2=Self_kw(ktmp, j);
//				selftmp3=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
				selftmp3=Self_array_tot[m + (l*(l+1))/2+j*nk8th];
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(50)
							<<selftmp3<<setw(50)<<selftmp2
							<<setw(25)<<real(selftmp3-selftmp2)/real(selftmp3)<<imag(selftmp3-selftmp2)/imag(selftmp3)<<'\n';
			}
			
			for (j=nbw/2; j<nbw/2+10; j++)
			{
				z=dcomplex(0, (2*j+1)*PI*tem);
				selftmp2=Self_kw(ktmp, j);
				selftmp3=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Sigma_inf4[m+(l*(l+1))/2]/(z*z*z*z)+Sigma_inf5[m+(l*(l+1))/2]/(z*z*z*z*z);
//				wn=(2*j+1)*PI*tem;
//				selftmp3=dcomplex(0,-Sigma_inf/wn+Sigma_inf3[m + (l*(l+1))/2]/(wn*wn*wn))
//						-Sigma_inf2[m + (l*(l+1))/2]/(wn*wn)+Sigma_inf4[m + (l*(l+1))/2]/(wn*wn*wn*wn);
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(50)
							<<selftmp3<<setw(50)<<selftmp2
							<<setw(25)<<real(selftmp3-selftmp2)/real(selftmp3)<<imag(selftmp3-selftmp2)/imag(selftmp3)<<'\n';
			}
			
			for (j=nbw-10; j<nbw; j++)
			{
				selftmp2=Self_kw(ktmp, j);
				z=dcomplex(0, (2*j+1)*PI*tem);
//				wn=(2*j+1)*PI*tem;
				selftmp3=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Sigma_inf4[m+(l*(l+1))/2]/(z*z*z*z)+Sigma_inf5[m+(l*(l+1))/2]/(z*z*z*z*z);
//				selftmp3=dcomplex(0,-Sigma_inf/wn+Sigma_inf3[m + (l*(l+1))/2]/(wn*wn*wn))
//						-Sigma_inf2[m + (l*(l+1))/2]/(wn*wn)+Sigma_inf4[m + (l*(l+1))/2]/(wn*wn*wn*wn);
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(50)
							<<selftmp3<<setw(50)<<selftmp2
							<<setw(25)<<real(selftmp3-selftmp2)/real(selftmp3)<<imag(selftmp3-selftmp2)/imag(selftmp3)<<'\n';
			}		
						
		}	
	}	
*/
/*	
	cout<<'\n';
	
	j=0;
	for (l=0; l<nk; l+=1)
	{
		ktmp[0]=l;
		for(m=0; m<=l; m+=1)
		{		
			ktmp[1]=m;

//			for (j=nbw/2-20; j<nbw/2; j++)
			{
				selftmp2=Self_kw(ktmp, j);
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(45)
							<<Self_kw_array[m + (l*(l+1))/2 + j*nk8th]<<setw(45)<<selftmp2
							<<2*fabs(Self_kw_array[m + (l*(l+1))/2 + j*nk8th]-selftmp2)/fabs(Self_kw_array[m + (l*(l+1))/2 + j*nk8th]+selftmp2)<<'\n';
			}
		}
	}
	
	j=nbw-1;
	for (l=0; l<nk; l+=1)
	{
		ktmp[0]=l;
		for(m=0; m<=l; m+=1)
		{		
			ktmp[1]=m;

//			for (j=nbw/2-20; j<nbw/2; j++)
			{
				selftmp2=Self_kw(ktmp, j);
				wn=(2*j+1)*PI*tem;
				selftmp3=dcomplex(0,-Sigma_inf/wn+Sigma_inf3[m + (l*(l+1))/2]/(wn*wn*wn))
						-Sigma_inf2[m + (l*(l+1))/2]/(wn*wn)+Sigma_inf4[m + (l*(l+1))/2]/(wn*wn*wn*wn);
				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(45)
							<<selftmp3<<setw(45)<<selftmp2
							<<2*fabs(selftmp3-selftmp2)/fabs(selftmp3+selftmp2)<<'\n';	
			}
		}
	}
*/
	
}

//! find the Fermi surface with interactions
void green::find_FS_inter()
{
	cout<<"calcul des coordonnees de la surface de Fermi avec interaction\n";
	
	long int j, l, m, x, y, tmp;
	
	if (!Self_kw_array)
	{
		cout<<"find_FS_inter(): self-energy non calculee\n";
		return;
	}
	
	if (!mu_calcule)
	{
		cout<<"find_FS_inter(): mu n'a pas ete calcule\n";
		return;
	}
	
	int Nn=4, nk8th=(nbk*(nbk+1))/2;
	dcomplex z, z0[4], self_ikn[4], self_w_poly, self_w_Pade;
	double Akw, ekmu;
	double eta1=0.000001;
	
	FS_inter=new int[nbk];
	Ak_FS_inter=new double[nbk];
	
	fstream file, file2, file3, file4;
	const char *nameForm="./FS_inter_AkF_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm2="./FS_inter_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm3="./FS_inter_AkF_Pade_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm4="./FS_inter_Pade_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[150];
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(8);
	
	sprintf(name, nameForm2, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file2.open(name,ios::out);
	file2<<setiosflags(ios::left)<<setprecision(8);
	
	sprintf(name, nameForm3, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file3.open(name,ios::out);
	file3<<setiosflags(ios::left)<<setprecision(8);
	
	sprintf(name, nameForm4, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file4.open(name,ios::out);
	file4<<setiosflags(ios::left)<<setprecision(8);
	
	dcomplex *coeffs=new dcomplex[Nn];
	
	int FS_inter_Pade;
	double Ak_FS_inter_Pade;
	
	for (l=0; l<nbk; l++)
	{
		Ak_FS_inter[l]=0;
		FS_inter[l]=0;
		
		Ak_FS_inter_Pade=0;
		FS_inter_Pade=0;
		for (m=0; m<nbk-l; m++)
		{	
			x=l+m;
			y=m;
			
			ekmu=ek[x+y*nbk]-mu;
			
			for (j=0; j<Nn; j++)
			{				
				z0[j]=dcomplex(0.0, (2*j+1)*PI*tem);
				z=z0[j];
				self_ikn[j]=Sigma_inf/z+Sigma_inf2[y+(x*(x+1))/2]/(z*z)+Sigma_inf3[y+(x*(x+1))/2]/(z*z*z)+Self_kw_array[y+(x*(x+1))/2+j*nk8th]/(z*z*z*z);
			}
			
			pade_cont_frac_coef_rec(self_ikn, z0, Nn, coeffs);
			
			z=dcomplex(0,eta1);
			
			self_w_Pade=pade(z, Nn, z0, coeffs);
			
			Akw=-imag(((double)1.0)/(z-ekmu-self_w_Pade))/PI;
			
			if (Akw>Ak_FS_inter_Pade)
			{
				Ak_FS_inter_Pade=Akw;
				FS_inter_Pade=m;
			}	
			
			//			self_w_poly= self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))*((z-z0[3])/(z0[0]-z0[3]))
			//						+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))*((z-z0[3])/(z0[1]-z0[3]))
			//						+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]))*((z-z0[3])/(z0[2]-z0[3]))
			//						+self_ikn[3]*((z-z0[0])/(z0[3]-z0[0]))*((z-z0[1])/(z0[3]-z0[1]))*((z-z0[2])/(z0[3]-z0[2]));
			
			self_w_poly= self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))
			+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))
			+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]));
			
			Akw=-imag(((double)1.0)/(z-ekmu-self_w_poly))/PI;
			//			cout<<setw(10)<<l<<setw(10)<<m<<Akw<<endl;
			
			if (Akw>Ak_FS_inter[l])
			{
				Ak_FS_inter[l]=Akw;
				FS_inter[l]=m;
			}	
		}
		file<<setw(10)<<l+FS_inter[l]<<setw(10)<<FS_inter[l]<<Ak_FS_inter[l]<<endl;
		file2<<setw(10)<<l+FS_inter[l]<<setw(10)<<FS_inter[l]<<endl;
		file3<<setw(10)<<l+FS_inter_Pade<<setw(10)<<FS_inter_Pade<<Ak_FS_inter_Pade<<endl;
		file4<<setw(10)<<l+FS_inter_Pade<<setw(10)<<FS_inter_Pade<<endl;
	}
	
	file.close();
	file2.close();
	file3.close();
	file4.close();
	
	delete [] coeffs;
	
	//	delete [] FS_inter;
	//	delete [] Ak_FS_inter;
	
	//	cout<<"FS_inter:\n";
	//	for (l=nbk-1; l>=0; l--)
	//		cout<<setw(10)<<l<<setw(10)<<FS_inter[l]<<Ak_FS_inter[l]<<endl;
	
}

/*
//! find the Fermi surface with interactions
void green::find_FS_inter()
{
	cout<<"calcul des coordonnees de la surface de Fermi avec interaction\n";
	
	long int j, l, m, x, y, tmp;
	
	if (!Self_kw_array)
	{
		cout<<"find_FS_inter(): self-energy non calculee\n";
		return;
	}
	
	if (!mu_calcule)
	{
		cout<<"find_FS_inter(): mu n'a pas ete calcule\n";
		return;
	}
	
	int Nn=4, nk8th=(nbk*(nbk+1))/2;
	dcomplex z, z0[4], self_ikn[4], self_w_poly;
	double Akw, ekmu;
	double eta1=0.001;
	
	FS_inter=new int[nbk];
	Ak_FS_inter=new double[nbk];
	
	fstream file, file2;
	const char *nameForm="./FS_inter_AkF_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm2="./FS_inter_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[150];
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file.open(name,ios::out);
	file<<setiosflags(ios::left);
	
	sprintf(name, nameForm2, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	file2.open(name,ios::out);
	file2<<setiosflags(ios::left);
	
	for (l=0; l<nbk; l++)
	{
		Ak_FS_inter[l]=0;
		FS_inter[l]=0;
		for (m=0; m<nbk; m++)
		{	
			ekmu=ek[m+l*nbk]-mu;
			
			x=l;
			y=m;
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			for (j=0; j<Nn; j++)
			{				
				z0[j]=dcomplex(0.0, (2*j+1)*PI*tem);
				z=z0[j];
				self_ikn[j]=Sigma_inf/z+Sigma_inf2[y+(x*(x+1))/2]/(z*z)+Sigma_inf3[y+(x*(x+1))/2]/(z*z*z)+Self_kw_array[y+(x*(x+1))/2+j*nk8th]/(z*z*z*z);
			}
			
						
			z=dcomplex(0,eta1);
						
//			self_w_poly= self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))*((z-z0[3])/(z0[0]-z0[3]))
//						+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))*((z-z0[3])/(z0[1]-z0[3]))
//						+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]))*((z-z0[3])/(z0[2]-z0[3]))
//						+self_ikn[3]*((z-z0[0])/(z0[3]-z0[0]))*((z-z0[1])/(z0[3]-z0[1]))*((z-z0[2])/(z0[3]-z0[2]));
			
			self_w_poly= self_ikn[0]*((z-z0[1])/(z0[0]-z0[1]))*((z-z0[2])/(z0[0]-z0[2]))
						+self_ikn[1]*((z-z0[0])/(z0[1]-z0[0]))*((z-z0[2])/(z0[1]-z0[2]))
						+self_ikn[2]*((z-z0[0])/(z0[2]-z0[0]))*((z-z0[1])/(z0[2]-z0[1]));
			
			Akw=-imag(((double)1.0)/(z-ekmu-self_w_poly))/PI;
//			cout<<setw(10)<<l<<setw(10)<<m<<Akw<<endl;
			
			if (Akw>Ak_FS_inter[l])
			{
				Ak_FS_inter[l]=Akw;
				FS_inter[l]=m;
			}	
		}
		file<<setw(10)<<l<<setw(10)<<FS_inter[l]<<Ak_FS_inter[l]<<endl;
		file2<<setw(10)<<l<<setw(10)<<FS_inter[l]<<endl;
	}
	
	file.close();
	file2.close();
	
//	delete [] FS_inter;
//	delete [] Ak_FS_inter;
	
//	cout<<"FS_inter:\n";
//	for (l=nbk-1; l>=0; l--)
//		cout<<setw(10)<<l<<setw(10)<<FS_inter[l]<<Ak_FS_inter[l]<<endl;

}
*/
/*
//! calculate the total density using the spectral weight
void green::calc_dens_Re_w(int N0, int mp, double *eta2, int Neta)
{
	if (!Self_kw_array)
	{
		cout<<"calc_dens_Re_w(): la self-energie n<a pas ete calculee\n";
		return;
	}
	
	dcomplex *coef=new dcomplex[N0];
	dcomplex *coef2=new dcomplex[N0];
	dcomplex *z0=new dcomplex[N0];
	dcomplex *func=new dcomplex[N0];
	long int *n_pade=new long int[N0];

	pade_coef=coef;
	pade_z0=z0;
	
	long int j, l, m, n;
	
	int Np=mp*((N0-1)/((int)pow(2.0,mp+1)))+(N0-1)/((int)pow(2.0,mp))+1;

	int j0, jf, d;

	jf=(N0-1)/((int)pow(2.0,mp));
	for (j=0; j<=jf; j++)
	{
		z0[j]=dcomplex(0.0, (2.0*j+1)*PI*tem);
		n_pade[j]=j;
	}

	int Nj=(N0-1)/((int)pow(2.0,mp+1));

	n=jf;
	d=1;
	for (l=0; l<mp; l++)
	{
		j0=jf+1;
		jf=j0+Nj-1;
		d=2*d;
		for (j=j0; j<=jf; j++)
		{
			n=n+d;
			n_pade[j]=n;
			z0[j]=dcomplex(0.0, (2.0*n+1)*PI*tem);	
		}
	}
	
	cout<<"calc_dens_Re_w():\n";
	cout<<"N0:  "<<N0<<endl;
	
	if (N0>nbw/2)
	{
		if (mp>0)
		{
			n_pade[Np-1]--;
			z0[Np-1]=dcomplex(0.0, (2.0*n_pade[Np-1]+1)*PI*tem);
		}
		else
		{
			Np-=2;
		}
	}
	
	cout<<"Np:  "<<Np<<endl;
	cout<<"n_pade[Np-1]:  "<<n_pade[Np-1]<<endl;


	int nk0=2*(nbk-1);
	int nk8th=nbk*(nbk+1)/2;
	double normFact=1.0/(nk0*nk0);


	
	fctPtr ptrInteg_dens=static_cast<fctPtr>(&green::integrand_density);
	fctPtr ptrAw=static_cast<fctPtr>(&green::spectral_weight);
	fctPtr ptrAwSelf=static_cast<fctPtr>(&green::spectral_weight_Self);
	
//	cout<<"Np: "<<Np<<endl;
//	Np=nbw/2-1;
//	for (n=0; n<Np; n++)	n_pade[n]=n;
//	for (n=0; n<Np; n++)	z0[n]=dcomplex(0.0, (2*n_pade[n]+1)*PI*tem);

	double params[7];
	params[1]=(double)Np;
	
	int nbEval[]={0};
	double tol=1.0e-6;
	double wmax=200;
	double lims[]={-wmax,wmax};
	double wn, C1, C2, C3, C4, a1, a2, e1, e2;
	
	l=0;
	m=0;	
	C1=Sigma_inf;
	C2=Sigma_inf2[m + (l*(l+1))/2];
	C3=Sigma_inf3[m + (l*(l+1))/2];
	C4=Sigma_inf4[m + (l*(l+1))/2];
	a1=C1;
	e1=C2/C1;
	a2=(C1*C3-C2*C2)/(C1*C1);
	e2=(C1*C1*C4+C2*C2*C2-2*C1*C2*C3)/(C1*(C1*C3-C2*C2));
	dcomplex SigmaHF, SigmaLF, z;
	for (j=0; j<Np; j++) 	
	{
		wn=(2*n_pade[j]+1)*PI*tem;
		z=dcomplex(0,wn);
//		SigmaHF=-I*Sigma_inf/wn-Sigma_inf2[m + (l*(l+1))/2]/(wn*wn)+I*Sigma_inf3[m + (l*(l+1))/2]/(wn*wn*wn)
//			    +Sigma_inf4[m + (l*(l+1))/2]/(wn*wn*wn*wn);
		SigmaLF=Self_kw_array[m + (l*(l+1))/2 + n_pade[j]*nk8th]-C4;
		func[j]=(SigmaLF*(a2-(z-e1)*(z-e2)) + a1*(e1*e1*e1*e1*(z-e2)+a2*a2*(2*e1+e2+z)+a2*(e1*e1*e1-e1*e2*(e2-2.0*z)+e2*e2*z+e1*e1*(3.0*z-2*e2))))/(a1*(e1*e1*e1*e1-a2*(z-e1)*(2*e1+e2+z))-SigmaLF*(z-e1));
	}
//	hx (a2 + (e1 - x) (-e2 + x)) + a1 (e1^4 (-e2 + x) + a2^2 (2 e1 + e2 + x) + a2 (e1^3 - e1 e2 (e2 - 2 x) + e2^2 x + e1^2 (-2 e2 + 3 x)))/(hx (e1 - x) + a1 (e1^4 + a2 (e1 - x) (2 e1 + e2 + x)))

//	pade_cont_frac_coef(func, z0, Np, coef);
	pade_cont_frac_coef_rec(func, z0, Np, coef);
	params[0]=ek[m+l*nbk] - mu;


	
	params[3]=a1, params[4]=e1, params[5]=a2, params[6]=e2;
	nbr_Akw_neg=0;
	sum_Akw_neg=0;
	cout<<"SigmaLF:\n";
	cout<<setprecision(12)<<setiosflags(ios::left);
	cout<<"eta                 "<<"integ(A(0,w))\n";       
	for (j=0; j<Neta; j++)
	{
		params[2]=eta2[j];
		nbEval[0]=0;
		cout<<setw(20)<<eta2[j]<<quadInteg(ptrAw, lims, tol, nbEval, params, NULL, NULL)/(2*PI)<<endl;
	}
	if (nbr_Akw_neg)
	{
		cout<<"nombre de valeurs de poids spectral negatif:  "<<nbr_Akw_neg<<'\n';
		cout<<"valeur moyenne des poids spectraux negatifs:  "<<sum_Akw_neg/nbr_Akw_neg<<'\n';			
	}	
	
	for (j=0; j<Np; j++) 	
	{
		z=dcomplex(0,(2*n_pade[j]+1)*PI*tem);
		func[j]=Sigma_inf/z+Sigma_inf2[m+(l*(l+1))/2]/(z*z)+Sigma_inf3[m+(l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
	}
	
	pade_cont_frac_coef_rec(func, z0, Np, coef2);
//	pade_cont_frac_coef(func, z0, Np, coef2);
	params[0]=ek[m+l*nbk] - mu;
	pade_coef=coef2;
	
	nbEval[0]=0;
	nbr_Akw_neg=0;
	sum_Akw_neg=0;
	cout<<"Sigma total:\n";
	cout<<setprecision(12)<<setiosflags(ios::left);
	cout<<"eta                 "<<"integ(A(0,w))\n";       
	for (j=0; j<Neta; j++)
	{
		params[2]=eta2[j];
		nbEval[0]=0;
		cout<<setw(20)<<eta2[j]<<quadInteg(ptrAwSelf, lims, tol, nbEval, params, NULL, NULL)/(2*PI)<<endl;
	}
	if (nbr_Akw_neg)
	{
		cout<<"nombre de valeurs de poids spectral negatif:  "<<nbr_Akw_neg<<'\n';
		cout<<"valeur moyenne des poids spectraux negatifs:  "<<sum_Akw_neg/nbr_Akw_neg<<'\n';			
	}

	pade_coef=coef;
	double *dens=new double[Neta];
	for (j=0; j<Neta; j++) dens[j]=0;
	
	nbr_Akw_neg=0;
	sum_Akw_neg=0;
	int wgt;
	for (l=0; l<nbk; l++)
	{
		for (m=0; m<=l; m++)
		{
			params[0]=ek[m+l*nbk] - mu;
			
			C1=Sigma_inf;
			C2=Sigma_inf2[m + (l*(l+1))/2];
			C3=Sigma_inf3[m + (l*(l+1))/2];
			C4=Sigma_inf4[m + (l*(l+1))/2];
			a1=C1;
			e1=C2/C1;
			a2=(C1*C3-C2*C2)/(C1*C1);
			e2=(C1*C1*C4+C2*C2*C2-2*C1*C2*C3)/(C1*(C1*C3-C2*C2));
			dcomplex SigmaHF, SigmaLF, z;
			for (j=0; j<Np; j++) 	
			{
				wn=(2*n_pade[j]+1)*PI*tem;
				z=dcomplex(0,wn);
				SigmaLF=Self_kw_array[m + (l*(l+1))/2 + n_pade[j]*nk8th]-C4;
				func[j]=(SigmaLF*(a2 - (z-e1)*(z-e2)) + a1*(e1*e1*e1*e1*(z-e2)+a2*a2*(2*e1+e2+z)+a2*(e1*e1*e1-e1*e2*(e2-2.0*z)+e2*e2*z+e1*e1*(3.0*z-2*e2))))/(a1*(e1*e1*e1*e1-a2*(z-e1)*(2*e1+e2+z))-SigmaLF*(z-e1));
			}

			pade_cont_frac_coef_rec(func, z0, Np, coef);			
//			pade_cont_frac_coef(func, z0, Np, coef);
					
			wgt=8;
			if (m==0 || m==nbk-1) wgt=wgt/2;
			if (l==0 || l==nbk-1) wgt=wgt/2;
			if (l==m) wgt=wgt/2;
			
			params[3]=a1, params[4]=e1, params[5]=a2, params[6]=e2;
			for (j=0; j<Neta; j++)
			{
				params[2]=eta2[j];
				nbEval[0]=0;
				dens[j]+=wgt*normFact*quadInteg(ptrInteg_dens, lims, tol, nbEval, params, NULL, NULL)/PI;
			}
			
		}
	}
	
	cout<<"densite calculee par calc_dens_Re_w():\n";
	cout<<"eta                 "<<"dens\n";       
	for (j=0; j<Neta; j++)
	{
		cout<<setw(20)<<eta2[j]<<dens[j]<<endl;
	}
	
	if (nbr_Akw_neg)
	{
		cout<<"nombre de valeurs de poids spectral negatif:  "<<nbr_Akw_neg<<'\n';
		cout<<"valeur moyenne des poids spectraux negatifs:  "<<sum_Akw_neg/nbr_Akw_neg<<'\n';			
	}
	
	delete [] dens;
	
	delete [] n_pade;
	delete [] coef;
	delete [] z0;
	delete [] func;
	
	pade_coef=NULL;
	pade_z0=NULL;	
	
}
*/

/*
//! calculate the total density using the spectral weight obtained with the Pade applied on the total self-energy
void green::calc_dens_Re_w_Self_tot(int N0, int mp, double *eta2, int Neta)
{
	if (!Self_kw_array)
	{
		cout<<"calc_dens_Re_w(): la self-energie n<a pas ete calculee\n";
		return;
	}
	
	dcomplex *coef=new dcomplex[N0];
	dcomplex *z0=new dcomplex[N0];
	dcomplex *func=new dcomplex[N0];
	int *n_pade=new int[N0];

	pade_coef=coef;
	pade_z0=z0;
	
	int j, l, m, n;
	
	int Np=mp*((N0-1)/((int)pow(2.0,mp+1)))+(N0-1)/((int)pow(2.0,mp))+1;

	int j0, jf, d;

	jf=(N0-1)/((int)pow(2.0,mp));
	for (j=0; j<=jf; j++)
	{
		z0[j]=dcomplex(0.0, (2.0*j+1)*PI*tem);
		n_pade[j]=j;
	}

	int Nj=(N0-1)/((int)pow(2.0,mp+1));

	n=jf;
	d=1;
	for (l=0; l<mp; l++)
	{
		j0=jf+1;
		jf=j0+Nj-1;
		d=2*d;
		for (j=j0; j<=jf; j++)
		{
			n=n+d;
			n_pade[j]=n;
			z0[j]=dcomplex(0.0, (2.0*n+1)*PI*tem);	
		}
	}
	
	cout<<"calc_dens_Re_w():\n";
	cout<<"N0:  "<<N0<<endl;
	
	if (N0>nbw/2)
	{
		if (mp>0)
		{
			n_pade[Np-1]--;
			z0[Np-1]=dcomplex(0.0, (2.0*n_pade[Np-1]+1)*PI*tem);
		}
		else
		{
			Np-=2;
		}
	}
	
	cout<<"Np:  "<<Np<<endl;
	cout<<"n_pade[Np-1]:  "<<n_pade[Np-1]<<endl;


	int nk0=2*(nbk-1);
	int nk8th=nbk*(nbk+1)/2;
	double normFact=1.0/(nk0*nk0);

	double *ek=new double[nbk*nbk];
	double kx, ky;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			ek[m+l*nbk]=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2.0*kx) + cos(2.0*ky));
		}
	}
	
	fctPtr ptrInteg_dens=&green::integrand_density;
	fctPtr ptrAw=&green::spectral_weight_Self;
	
//	cout<<"Np: "<<Np<<endl;
//	Np=nbw/2-1;
//	for (n=0; n<Np; n++)	n_pade[n]=n;
//	for (n=0; n<Np; n++)	z0[n]=dcomplex(0.0, (2*n_pade[n]+1)*PI*tem);

	double params[3];
	params[1]=(double)Np;
	
	int nbEval[]={0};
	double tol=1.0e-6;
	double wmax=200;
	double lims[]={-wmax,wmax};
	
	l=0;
	m=0;
	for (j=0; j<Np; j++) 	func[j]=Self_kw_array[m + (l*(l+1))/2 + n_pade[j]*nk8th];
	pade_cont_frac_coef(func, z0, Np, coef);
	params[0]=ek[m+l*nbk] - mu;
	
	cout<<setprecision(12)<<setiosflags(ios::left);
	cout<<"eta                 "<<"integ(A(0,w))\n";       
	for (j=0; j<Neta; j++)
	{
		params[2]=eta2[j];
		nbEval[0]=0;
		cout<<setw(20)<<eta2[j]<<quadInteg(ptrAw, lims, tol, nbEval, params, NULL, NULL)/(2*PI)<<endl;
	}
	
	double *dens=new double[Neta];
	for (j=0; j<Neta; j++) dens[j]=0;
	
	int wgt;
	for (l=0; l<nbk; l++)
	{
		for (m=0; m<=l; m++)
		{
			params[0]=ek[m+l*nbk] - mu;
			
			for (j=0; j<Np; j++) 	func[j]=Self_kw_array[m + (l*(l+1))/2 + n_pade[j]*nk8th];
		
			pade_cont_frac_coef(func, z0, Np, coef);
					
			wgt=8;
			if (m==0 || m==nbk-1) wgt=wgt/2;
			if (l==0 || l==nbk-1) wgt=wgt/2;
			if (l==m) wgt=wgt/2;
			
			for (j=0; j<Neta; j++)
			{
				params[2]=eta2[j];
				nbEval[0]=0;
				dens[j]+=wgt*normFact*quadInteg(ptrInteg_dens, lims, tol, nbEval, params, NULL, NULL)/PI;
			}
			
		}
	}
	
	cout<<"densite calculee par calc_dens_Re_w():\n";
	cout<<"eta                 "<<"dens\n";       
	for (j=0; j<Neta; j++)
	{
		cout<<setw(20)<<eta2[j]<<dens[j]<<endl;
	}
	
	if (nbr_Akw_neg)
	{
		cout<<"nombre de valeurs de poids spectral negatif:  "<<nbr_Akw_neg<<'\n';
		cout<<"valeur moyenne des poids spectraux negatifs:  "<<sum_Akw_neg/nbr_Akw_neg<<'\n';			
	}	
	
	delete [] n_pade;
	delete [] coef;
	delete [] z0;
	delete [] func;
	
	pade_coef=NULL;
	pade_z0=NULL;	
	
}
*/

//! integrant dans le calcul de n(k) en frequences reelles
double green::integrand_density(double w, double params[])
{
	double f;
	
	if (w<=0)
		f=1.0/(exp(w/tem)+1.0);
	else
	{
		double expw=exp(-w/tem);
		f=expw/(expw+1.0);
	}
	
	return f*spectral_weight(w, params);
}

//! return a value of the spectral weight in real frequency using Pade approximant for the self-energy, 
//! the Pade coefficients must have been calculated before with pade_cont_frac_coef, 
//!params[0]: epsilon(k)-mu, params[1] : number of points used for the Pade
//! coefficients calculation, params[2]: small imaginary part in the energy
double green::spectral_weight_Self(double w, double params[])
{
	double ekmu=params[0];
	int N=(int)params[1];
	double eta2=params[2];
	
	dcomplex selfkw=pade(dcomplex(w,eta), N, pade_z0, pade_coef);
	
	double Akw=-2.0*imag(((double)1.0)/(dcomplex(w,eta2)-ekmu-selfkw));

	if (Akw<0)
	{
		nbr_Akw_neg++;
		sum_Akw_neg+=Akw;
//		if (nbr_Akw_neg<50) cout<<"poids spectral negatif!  ekmu:  "<<ekmu <<"   w:  "<<w<<"   A(w):  "<<Akw<<'\n';
	}

	return Akw;
}

//! return a value of the spectral weight in real frequency using Pade approximant for the third floor of a continued
//! fraction form for the self-energy 
//! the Pade coefficients must have been calculated before with pade_cont_frac_coef, 
//! params[0]: epsilon(k)-mu, 
//! params[1] : number of points used for the Pade coefficients calculation,
//! params[2]: small imaginary part in the energy, 
//! params[3]...params[6]: coefficients in the continued fraction
double green::spectral_weight(double w, double params[])
{
	double ekmu=params[0];
	int N=(int)params[1];
	double eta2=params[2];
	double a1=params[3], e1=params[4], a2=params[5], e2=params[6];
	
	dcomplex fw=pade(dcomplex(w,eta), N, pade_z0, pade_coef);
	
	double Akw=-2.0*imag(((double)1.0)/(dcomplex(w,eta2)-ekmu-a1/(dcomplex(w,eta)-e1-a2/(dcomplex(w,eta)-e2-fw))));

	if (Akw<0)
	{
		nbr_Akw_neg++;
		sum_Akw_neg+=Akw;
//		if (nbr_Akw_neg<50) cout<<"poids spectral negatif!  ekmu:  "<<ekmu <<"   w:  "<<w<<"   A(w):  "<<Akw<<'\n';
	}

	return Akw;
}


//! trace Sigma*G2
void green::traceSelfG2()
{
	double sum=0.0, wn, rSig, iSig, tmp;
	double ektmp, kx, ky;
	
	long int j, l, m;
	int w;
	int nk=nbk;
	int nk02=4*(nk-1)*(nk-1);
	int nk8th=(nk*(nk+1))/2;
	ptrdiff_t w_local = nbw/2+1;
	ptrdiff_t w_start = 0;
	
	if (w_start+w_local==nbw/2+1) w_local--;
	
#pragma omp parallel private(l,m,kx,ky,ektmp,w,j,wn,rSig,iSig,tmp) reduction(+:sum)
	{
		dcomplex z, selftmp;	
#pragma omp for
		for (l=0; l<nk; l++)
		{
			kx=l*PI/(nk-1);
			for(m=0; m<=l; m++)
			{
				ky=m*PI/(nk-1);
//				ek=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky));
				ektmp=ek[m+l*nk];
				w=8;
				if (l==0 || l==nk-1) w=w/2;
				if (m==0 || m==nk-1) w=w/2;
				if (m==l) w=w/2;
				for (j=0; j<w_local; j++)
				{
					wn=(2*(j+w_start)+1)*PI*tem;
					z=dcomplex(0,wn);
					selftmp=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
					rSig=selftmp.real();
					iSig=selftmp.imag();
					
					//				rSig=real( Self_kw_array[m + (l*(l+1))/2 + j*nk8th] );
					//				iSig=imag( Self_kw_array[m + (l*(l+1))/2 + j*nk8th] );
					
					tmp=w*tem*2/( nk02*(  (wn-iSig)*(wn-iSig)+(ektmp-mu+rSig)*(ektmp-mu+rSig) ) );
					sum+=tmp*( iSig*wn - iSig*iSig -(ektmp-mu)*rSig - rSig*rSig );
				}
			}
		}
	}
	
//	cout<<"trace Self*G2: "<<sum<<"  U*dblOcc-U*n^2/4:  "<<U*(chi1->dblOcc)-U*density*density/4.0<<'\n';
	cout<<"trace Self*G2: "<<sum<<"  U*dblOcc-U*n^2/4:  "<<U*dblOcc-U*density*density/4.0<<'\n';	
	cout<<"double occupation (Self*G2): "<<(sum+U*density*density/4.0)/U<<'\n';
	
	fstream paramsFile;
	
	const char *nameForm="./TPSC_infos_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	paramsFile.open(name,ios::out | ios::app);
	paramsFile<<setprecision(16);
	paramsFile<<"trace Self*G2:  "<<sum<<'\n';
	paramsFile<<"double occupation (Self*G2): "<<(sum+U*density*density/4.0)/U<<'\n';
	paramsFile<<"mu0:  "<<mu0<<'\n';
	paramsFile<<"mu:  "<<mu<<'\n';
	paramsFile.close();

}


//! trace Sigma*G1
void green::traceSelfG()
{
	/* Version parallele */
	double sum=0.0, wn, rSig, iSig, tmp;
	double ektmp, kx, ky;
	
	long int j, l, m;
	int w;
	int nk=nbk;
	int nk02=4*(nk-1)*(nk-1);
	int nk8th=(nk*(nk+1))/2;
	ptrdiff_t w_local = nbw/2+1;
	ptrdiff_t w_start = 0;
	
	if (w_start+w_local==nbw/2+1) w_local--;
	
	
#pragma omp parallel private(l,m,kx,ky,ektmp,w,j,wn,rSig,iSig,tmp) reduction(+:sum)
	{
		dcomplex z, selftmp;
#pragma omp for	
		for (l=0; l<nk; l++)
		{
			kx=l*PI/(nk-1);
			for(m=0; m<=l; m++)
			{
				ky=m*PI/(nk-1);
				ektmp=ek[m+l*nk];
//				ek=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky));
				w=8;
				if (l==0 || l==nk-1) w=w/2;
				if (m==0 || m==nk-1) w=w/2;
				if (m==l) w=w/2;
				for (j=0; j<w_local; j++)
				{
					wn=(2*(j+w_start)+1)*PI*tem;
					z=dcomplex(0,wn);
					
					selftmp=Sigma_inf/z+Sigma_inf2[m + (l*(l+1))/2]/(z*z)+Sigma_inf3[m + (l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
					rSig=selftmp.real();
					iSig=selftmp.imag();
					
					//				rSig=real( Self_kw_array[m + (l*(l+1))/2 + j*nk8th] );
					//				iSig=imag( Self_kw_array[m + (l*(l+1))/2 + j*nk8th] );					
					
					tmp=w*tem*2/(nk02*(wn*wn+(ektmp-mu0)*(ektmp-mu0)));
					sum+=tmp*( rSig*(mu0-ektmp) + iSig*wn );
				}
			}
		}
	}
	
	double eps=1.0e-12;
	double expk, sum2=0;
#pragma omp parallel for private(l,m,kx,ky,ektmp,w,j,wn,expk) reduction(+:sum2)
	for (l=0; l<nk; l++)
	{
		kx=l*PI/(nk-1);
		for(m=0; m<=l; m++)
		{
			ky=m*PI/(nk-1);
			ektmp=ek[m+l*nk];
//			ek=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky));
			w=8;
			if (l==0 || l==nk-1) w=w/2;
			if (m==0 || m==nk-1) w=w/2;
			if (m==l) w=w/2;
			for (j=0; j<w_local; j++)
			{
				wn=(2*(j+w_start)+1)*PI*tem;
				sum2+=w*Sigma_inf*tem*2/(nk02*(wn*wn+(ektmp-mu0)*(ektmp-mu0)));
			}
			if (fabs(ektmp-mu0)>eps)
				sum2+=w*Sigma_inf*(1.0/(exp((ektmp-mu0)/tem)+1.0) - 0.5)/(ektmp-mu0)/nk02;
			else
			{
				if ((ektmp-mu0)<0)
					expk=exp((ektmp-mu0)/tem);
				else
					expk=exp(-(ektmp-mu0)/tem);
				sum2-=w*Sigma_inf*expk/( ( expk + 1.0 )*( expk + 1.0 )*nk02*tem);
			}	
		}
	}
	
//	cout<<"U*dblOcc-U*n^2/4:  "<<U*(chi1->dblOcc)-U*density*density/4.0<<'\n';
	cout<<"U*dblOcc-U*n^2/4:  "<<U*dblOcc-U*density*density/4.0<<'\n';
	cout<<"trace Self*G: "<<sum<<'\n';
	cout<<"trace Self*G (avec correction asymptotique): "<<sum+sum2<<'\n';
	
	fstream paramsFile;
	
	const char *nameForm="./TPSC_infos_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
	paramsFile.open(name,ios::out | ios::app);
	paramsFile<<setprecision(16);
	paramsFile<<"U*dblOcc-U*n^2/4:   "<<U*dblOcc-U*density*density/4.0<<'\n';
	paramsFile<<"trace Self*G:  "<<sum<<'\n';
	paramsFile.close();
		
}


//! trace of G-G_0
double green::traceDiffG(double pmu, double par[])
{
	long int j, l, m;
	int nk=nbk;
	//	int nk0=2*(nbk-1);
	int nk8th=(nk*(nk+1))/2;
	int w;
	
	double sum=0, sum1=0, wn, rSig, iSig, tmp;
	double ektmp;
	double kx, ky;
	
	
#pragma omp parallel for private(l,m,kx,ky,ektmp,w,j,wn,rSig,iSig,tmp) reduction(+:sum1)
	for (l=0; l<nk; l++)
	{
		kx=l*PI/(nk-1);
		for(m=0; m<=l; m++)
		{
			ky=m*PI/(nk-1);
			ektmp=ek[m+l*nk];
//			ektmp=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky));
			w=16;
			if (l==0 || l==nk-1) w=w/2;
			if (m==0 || m==nk-1) w=w/2;
			if (m==l) w=w/2;
			for (j=nbw; j>=nbw/2; j--)
			{
				wn=(2*j+1)*PI*tem;
				rSig=-Sigma_inf2[m + (l*(l+1))/2]/(wn*wn)+Sigma_inf4[m + (l*(l+1))/2]/(wn*wn*wn*wn); //-Sigma_inf6[m + (l*(l+1))/2]/(wn*wn*wn*wn*wn*wn);
				iSig=-Sigma_inf/wn+Sigma_inf3[m + (l*(l+1))/2]/(wn*wn*wn); //-Sigma_inf5[m + (l*(l+1))/2]/(wn*wn*wn*wn*wn)+Sigma_inf7[m + (l*(l+1))/2]/(wn*wn*wn*wn*wn*wn*wn);
				tmp=((double) w)/(((wn-iSig)*(wn-iSig)+(ektmp-pmu+rSig)*(ektmp-pmu+rSig))*(wn*wn+(ektmp-mu0)*(ektmp-mu0)));
				sum1+=tmp*((mu0-pmu+rSig)*(ektmp-pmu+rSig)*(ektmp-mu0)+iSig*(iSig-2*wn)*(ektmp-mu0)-wn*wn*(mu0-pmu+rSig));
			}
		}
	}
	
#pragma omp parallel private(l,m,kx,ky,ektmp,w,j,wn,rSig,iSig,tmp) reduction(+:sum)
	{
		dcomplex z, selftmp;
#pragma omp for		
		for (l=0; l<nk; l++)
		{
			kx=l*PI/(nk-1);
			for(m=0; m<=l; m++)
			{
				ky=m*PI/(nk-1);
				ektmp=ek[m+l*nk];
//				ektmp=-2.0*(cos(kx) + cos(ky)) - 4.0*tp*cos(kx)*cos(ky) - 2.0*tpp*(cos(2*kx) + cos(2*ky));
				w=16;
				if (l==0 || l==nk-1) w=w/2;
				if (m==0 || m==nk-1) w=w/2;
				if (m==l) w=w/2;
				for (j=nbw/2-1; j>=0; j--)
				{
					wn=(2*j+1)*PI*tem;
					z=dcomplex(0,wn);
					selftmp=Sigma_inf/z+Sigma_inf2[m+(l*(l+1))/2]/(z*z)+Sigma_inf3[m+(l*(l+1))/2]/(z*z*z)+Self_kw_array[m+(l*(l+1))/2+j*nk8th]/(z*z*z*z);
					rSig=real(selftmp);
					iSig=imag(selftmp);
					
					//				rSig=real(Self_kw_array[m + (l*(l+1))/2 + j*nk8th]);
					//				iSig=imag(Self_kw_array[m + (l*(l+1))/2 + j*nk8th]);
					
					tmp=((double) w)/(((wn-iSig)*(wn-iSig)+(ektmp-pmu+rSig)*(ektmp-pmu+rSig))*(wn*wn+(ektmp-mu0)*(ektmp-mu0)));
					sum+=tmp*((mu0-pmu+rSig)*(ektmp-pmu+rSig)*(ektmp-mu0)+iSig*(iSig-2*wn)*(ektmp-mu0)-wn*wn*(mu0-pmu+rSig));
				}
			}
		}
	}
	
	cout<<"mu: "<<pmu<<"  trace(G-G0): "<<sum+sum1<<'\n';
	
	return sum+sum1;
}


void green::find_mu()
{
	fstream paramsFile;
	fstream muFile;
	
	const char *nameForm="./TPSCparams_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *muFile_nameForm="./mu0_mu_D_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	char paramName[100], paramsFileName[100];
	double param;
	double inexistant = -100000.;
	double uncomputed = 100000.;
	
	
	sprintf(paramsFileName, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	sprintf(paramName, muFile_nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbk-1), nbw);
	
/*
	mu = uncomputed;
 
	paramsFile.open(paramsFileName,ios::in);
 
	if (paramsFile)
	{
		for (int j=0; j<12; j++)
		{
			paramsFile>>paramName>>param;
		}
		
		if (strcmp("mu",paramName)==0)
		{
			mu=param;
			cout<<"mu:  "<<mu<<'\n';
			
			traceSelfG2();
			
			paramsFile.close();
			
			return;
		}
		
		paramsFile.close();
		paramsFile.clear();
	}
	else
	{
		cout<<"green::find_mu():  fichier de parametres inexistant\n";
		
		return;
	}
*/	
	
	fctPtr ptr=static_cast<fctPtr> (&green::traceDiffG);
	
	double init[]={0.0, -10.0, 10.0}, root[]={0.0}, lims[]={-20, 20};

/*	
	if (density==1.0 && tp==0 && tpp==0) 
	{
		mu=0;
		traceDiffG(0,NULL);
		return;
	}
*/	
	
	if (find_zero(ptr, init, NULL, root, lims))
	{
		mu=root[0];
		mu_calcule=true;		
		cout<<"mu:  "<<mu<<'\n';
	}
	else
	{
		cout<<"Attention! mu non calcule!\n";
	}
	
	
	paramsFile.open(paramsFileName, ios::out | ios::app);
	if (paramsFile)
	{
		paramsFile<<setprecision(14);
		paramsFile<<"mu  "<<mu<<'\n';
		
		paramsFile.close();
	}
	else
		cout<<"echec d'ouverture du fichier de parametres\n";
	
	muFile.open(paramName, ios::out);
	if (muFile)
	{
		muFile<<setprecision(15)<<setiosflags(ios::left);
		muFile<<setw(25)<<mu0<<setw(25)<<mu<<dblOcc<<'\n';
		
		muFile.close();
	}
	
	traceSelfG2();
}


//print the parameters on screen
void green::get_params() const
{
	chi::get_params();
	cout<<"mu: "<<mu<<'\n';
}
