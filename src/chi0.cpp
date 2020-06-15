
/*
 *  chi0.cpp
 *  TPSC
 *
 *  Created by Dominic Bergeron
 *
 */

#include "includeDef.h"
#include "chi0.h"
#include <cstdlib>

//! Constructor
chi0::chi0(int nk0, int Nw0, double n, double T, double  params[]):green0(nk0, Nw0, n, T, params)
{
	nbq=nk0/2+1;
	nwmax=Nw0/2;
	
	keep_chirt=true;
	keep_chiqw=true;
	save_chi0_q=true;
	save_chi0_to_file=true;
	
	chi0_inf=NULL;
	chi0_inf2=NULL;
	dtau_Grt0p=NULL;
	dtau_chirt0=NULL;
	chiqw_array=NULL;
	chirt=NULL;
	
	chi0max=0.0;
	
	qmax[0]=0;
	qmax[1]=0;
	
	tolChi0Integ=1.0e-8;
	
	nwr=0;
	
}

//! Destructor
chi0::~chi0()
{
	if (chi0_inf) delete [] chi0_inf;
	
	if (chi0_inf2) delete [] chi0_inf2;
	
	if (chiqw_array) delete [] chiqw_array;

	if (chirt) delete [] chirt;
	
	if (dtau_chirt0) delete [] dtau_chirt0;
	
	if (dtau_Grt0p) delete [] dtau_Grt0p;
	
}

//! save chi0(q,iqn)
void chi0::save_chi()
{
	char chi0FileName[200];
	const char nameForm[]="./chi0_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";

	int nq0=2*(nbq-1);
	int nw0=2*nwmax;
	int nk8th=(nbq*(nbq+1))/2;

	if (!chiqw_array) return;
	sprintf(chi0FileName, nameForm, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nq0, nw0);

	long offset = 0;
	long int j, jmax = nwmax + 1;

	chi0File.open(chi0FileName, ios::out|ios::binary);
	chi0File.seekp(offset);

	if (chi0File.good())
	{
		if (offset==0)
		{
			chi0File.write((char*)&nbq,sizeof(nbq));
			chi0File.write((char*)&nwmax,sizeof(nwmax));
			offset = sizeof(nbq) + sizeof(nwmax);
		}

		for (j=0; j<jmax; j++)
			chi0File.write((char*)&chiqw_array[j*nk8th],nbq*(nbq+1)*sizeof(double)/2);

        	offset += jmax * nbq * (nbq + 1)  * sizeof(double) / 2;
	}
	else std::cout << "Error in chi0::save_chi" << std::endl;
	
	chi0File.close();
}


void chi0::loadChi0()
{
	char chi0FileName[200];
	const char *nameForm="./chi0_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";

	int nq0=2*(nbq-1);
	int nw0=2*nwmax;
	sprintf(chi0FileName, nameForm, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nq0, nw0);
	
	long offset;	
	int nk8th;
	ptrdiff_t local_n;
	ptrdiff_t local_start;
	long int alloc_local;

	fstream file;
	
	file.open(chi0FileName, ios::in|ios::binary);
	if (!file)
	{
		cout << "chi0::loadChi0():  fichier  "<<chi0FileName <<"   inexistant \n";
		return;
	}
	file.seekg(0);
	file.read((char*)&nbq, sizeof(nbq));
	file.read((char*)&nwmax, sizeof(nwmax));
	offset = sizeof(nbq) + sizeof(nwmax);
	nk8th=(nbq*(nbq+1))/2;
	local_n = nwmax+1;
	local_start = 0;
	alloc_local = (long int) (nwmax+1)*nk8th;
	

	if (chiqw_array) delete [] chiqw_array;
	chiqw_array=new double[alloc_local];

	long int j;
	
	for (j=0; j<local_n; j++)
		file.read((char*)&chiqw_array[j*nk8th],(nbq+1)*nbq/2*sizeof(double));
	file.close();

	chi0max=0.0;
	int l,m;
	for (l=0; l<nbq; l++)
		for (m=0; m<=l; m++)
		{
			if (chiqw_array[m+(l*(l+1))/2]>chi0max)
			{
				chi0max=chiqw_array[m+(l*(l+1))/2];
				qmax[0]=l;
				qmax[1]=m;
			}
		}
	cout<<"chi0 charge,  nbq: "<<nbq<<"  nwmax: "<<nwmax<<'\n';
	cout<<"qmax: "<<qmax[0]<<"  "<<qmax[1]<<"  chi0max: "<<chi0max<<'\n';

}

//calculate chi0(q,w) in the region delimited by q1=(qx1,qy1), q2=(qx2,qy2) at the frequency wn=2*n*pi*tem. qlims={qx1, qy1, qx2, qy2}. Nq={Nqx,Nqy} is the number of points in qx and qy directions
void chi0::calc_chi0_qlims(int n, double qlims[], int Nq[])
{
	find_mu_inf_syst();
	
	double qx, par[]={0, 2*n*PI*tem};
	
	int Nqx=Nq[0], Nqy=Nq[1];
	double dqx=0, dqy=0;
	
	if (Nqx>1)	dqx=(qlims[2]-qlims[0])/(Nqx-1);
	if (Nqy>1)	dqy=(qlims[3]-qlims[1])/(Nqy-1);
	
	double chi0tmp;
	
	const char *nameFormq="./chi0_q_%5.3f_%5.3f_%5.3f_%7.5f_%7.5f_%7.5f_%7.5f_%7.5f_%7.5f_%d_%d_%d.dat";
	char name[200];
	
	fstream file;
	
	sprintf(name, nameFormq, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, (double)(qlims[0]/PI), (double)(qlims[1]/PI), (double)(qlims[2]/PI), (double)(qlims[3]/PI), Nqx, Nqy, n);
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(16);
	
	cout<<"calcul de chi0(q,wn="<<par[1]<<")"<<'\n';
	
	qx=qlims[0];
	int jx=1, jy, j=0;
	while (qx<=(qlims[2]+EPSILON))
	{
		par[0]=qlims[1];
		jy=1;
		while (par[0]<=(qlims[3]+EPSILON) && par[0]<=(qx+EPSILON))
		{
			chi0tmp=chiqw(qx, par);
			file<<setw(25)<<qx<<setw(25)<<par[0]<<chi0tmp<<'\n';
//			cout<<setw(25)<<qx/PI<<setw(25)<<par[0]/PI<<chi0tmp<<'\n';
			par[0]=qlims[1]+jy*dqy;
			jy++;
			j++;
		}
		qx=qlims[0]+jx*dqx;
		jx++;
	}

	file.close();	
	
}


// Calculate chi0(q,iw) for one value of (q,w), par=(qy, w)
double chi0::chiqw(double qx, double par[])
{
	double tol=1e-4;

	int nbEval[]={0};
	double params[]={0.0, 0.0, 0.0, 0.0}, f=1.0;
	double limx[]={0.0, 2.0*PI}, limy[]={0.0, 2.0*PI};

	fctPtr Ptr=static_cast<fctPtr>(&chi0::chi0IntegRe);
//params dans chi0IntegRe: ky=params[0], qx=params[1], qy=params[2], wn=params[3];	
	
	params[1]=qx;
	params[2]=par[0];
	params[3]=par[1];

	double chitmp=0.0;

	if (params[1]==0.0  || params[1]==PI)
	{
		limx[1]=PI;
		f=2.0*f;
	}
	if (params[2]==0.0  || params[2]==PI)
	{
		limy[1]=PI;
		f=2.0*f;
	}

	int Ninterv=4;
	
//#pragma omp parallel reduction(+:chitmp)
	{
		int j;
		double limytmp[2];
		
//#pragma omp for	
		for (j=0; j<Ninterv; j++)
		{
			limytmp[0]=j*limy[1]/Ninterv;
			limytmp[1]=(j+1)*limy[1]/Ninterv;
			chitmp-=f/(2.0*PI*PI)*quadInteg2D(Ptr, limx, limytmp, tol, nbEval, params);
		}
	}
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	tol=1e-10*chitmp;
	
	if ( tol<1e-14 ) tol=1e-14;
		
	chitmp=0;
	
//#pragma omp parallel reduction(+:chitmp)
	{
		int j;
		double limytmp[2];
		
//#pragma omp for	
		for (j=0; j<Ninterv; j++)
		{
			limytmp[0]=j*limy[1]/Ninterv;
			limytmp[1]=(j+1)*limy[1]/Ninterv;
			chitmp-=f/(2.0*PI*PI)*quadInteg2D(Ptr, limx, limytmp, tol, nbEval, params);
		}
	}
//	cout<<setw(25)<<params[1]<<setw(25)<<params[2]<<chitmp<<'\n';

	return chitmp;
}

// calculer chi0 par interpolation par splines cubiques et par transformee de Fourier rapide 
void chi0::calc_chiqw_FFT_spline(int nq0, int nw0)
{
	free_chi_array();
	
	//	nwr=0;
	
	int nk=nq0/2+1;
	int nk02=nq0*nq0;
	int nk2=nk*nk;
	int nk8th=(nk*(nk+1))/2;
	int nw=nw0/2+1;
	
	nbq=nk;
	
	cout<<"calcul de chi0 par calc_chiqw_FFT_spline()\n";
	
	//	cout<<"nwmax avant set_grid: "<<nwmax<<endl;
	
	set_grid(nq0, nw0);
	
	//	cout<<"nwmax apres set_grid: "<<nwmax<<endl;
	
	long int array_size=(long int)nw*nk8th;
	
	//	cout<<"array_size:  "<<array_size<<endl;
	
	chiqw_array=new double[array_size];	
	//	chiqw_array=new double[nw*nk8th];
	
	cout<<"nq0: "<<nq0<<"  nw0: "<<nw0<<'\n';
	
	double *Grt_tmp=new double[nk8th];
	double *Gkt_tmp=new double[nk8th];
	double *ekmu=new double[nk8th];
	double *dtau_chiqt0=new double[nk8th];
	double *dtau_Grt0m=new double[nk8th];
	double *chir_tmp=new double[nk2];
	
	long int j,l,m;
	double kx, ky, t;
	
	free_dtau_Grt0p();
	dtau_Grt0p=new double[nk8th];
	
	dtau_Grt0m[0]=-mu0;
	for (j=1; j<nk8th; j++) dtau_Grt0m[j]=er[j];
	
	fftw_plan fftplan_k;
	unsigned flag = FFTW_MEASURE;
	
	// Frequency transform 
	fftw_r2r_kind kind = FFTW_REDFT00;
	
	fftplan_k=fftw_plan_r2r_2d(nk, nk, chir_tmp, chir_tmp, kind, kind, flag);
	
	cout<<setiosflags(ios::left);
	
	//#pragma omp parallel for private(l,m,kx,ky)
	for (l=0; l<nk; l++)
	{
		for(m=0; m<=l; m++)
		{
			ekmu[m + (l*(l+1))/2]=ek[m+l*nbk]-mu0;
			Gkt_tmp[m + (l*(l+1))/2]=-1.0/(exp(-ekmu[m + (l*(l+1))/2]/tem)+1.0);
			//			cout<<setw(20)<<l<<setw(20)<<m<<setw(20)<<ekmu[m + (l*(l+1))/2]<<Gkt_tmp[m + (l*(l+1))/2]<<endl;
		}
	}
	
	//calcul de d/dtau G(r,tau)(tau=0) (pour le calcul de d/dtau chi(q,tau)(tau=0) necessaire aux conditions aux frontieres des splines de chi(q,tau))
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<l; m++)
		{
			//			cout<<setw(50)<<l<<setw(10)<<m<<"m + l*nk: "<<m + l*nk<<endl;
			chir_tmp[m + l*nk]=-ekmu[m + (l*(l+1))/2]*Gkt_tmp[m + (l*(l+1))/2]/nk02;
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
		chir_tmp[l + l*nk]=-ekmu[l + (l*(l+1))/2]*Gkt_tmp[l + (l*(l+1))/2]/nk02;
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
		for(m=0; m<=l; m++)
		{
			dtau_Grt0m[m + (l*(l+1))/2]-=chir_tmp[m + l*nk];
			dtau_Grt0p[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
	
	if (keep_chirt)		chirt=new double[array_size];	
	//	if (keep_chirt)		chirt=new double[nw*nk8th];
	
	// calcul de chi0(q,tau=0) et d/dtau chi(q,tau)(tau=0)
	// calcul de G(r,tau=0)	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=Gkt_tmp[m + (l*(l+1))/2]/nk02;
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
		chir_tmp[l + l*nk]=Gkt_tmp[l + (l*(l+1))/2]/nk02;
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
		for(m=0; m<=l; m++)
		{
			Grt_tmp[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
	
	//calcul de chi(q,tau=0)	
	chir_tmp[0]=-2.0*chir_tmp[0]*(chir_tmp[0]+1.0);
	if (keep_chirt)	chirt[0]=chir_tmp[0];
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=-2.0*chir_tmp[m + l*nk]*chir_tmp[m + l*nk];
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
			if (keep_chirt)	chirt[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
		m=l;
		chir_tmp[l + l*nk]=-2.0*chir_tmp[l + l*nk]*chir_tmp[l + l*nk];
		if (keep_chirt)	chirt[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	//	fftw_execute(fftplan_k);
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			chiqw_array[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	free_dtau_chirt0();
	dtau_chirt0=new double[nk8th];
	
	//calcul de d/dtau chi(q,tau)(tau=0)
	chir_tmp[0]=-2.0*(dtau_Grt0p[0]*(1.0+Grt_tmp[0])+dtau_Grt0m[0]*Grt_tmp[0]);
	dtau_chirt0[0]=chir_tmp[0];
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=-2.0*(dtau_Grt0p[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]+dtau_Grt0m[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]);
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
			dtau_chirt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
		m=l;
		chir_tmp[l + l*nk]=-2.0*(dtau_Grt0p[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]+dtau_Grt0m[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]);
		dtau_chirt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}	
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			dtau_chiqt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	delete [] dtau_Grt0m;
	delete [] chir_tmp;
	delete [] Grt_tmp;
	
	//calcul de chi0(q,tau) pour 0<tau<beta/2	
#pragma omp parallel private(j,t,l,m,chir_tmp,Grt_tmp)	
	{
		chir_tmp=new double[nk2];
		Grt_tmp=new double[nk8th];
		
#pragma omp for		
		for (j=1; j<nw-1;j++)
		{
			t=j/(nw0*tem);
			// calcul de G(r,t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=exp(-t*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + (l*(l+1))/2])/(exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)/nk02;
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// calcul de G(r,-t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					// conserver G(r,t)
					Grt_tmp[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
					
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=-exp((t-1.0/tem)*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=exp(t*ekmu[m + (l*(l+1))/2])/((exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)*nk02);
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			//calcul de chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<l; m++)
				{
					chir_tmp[m + l*nk]=-2.0*chir_tmp[m + l*nk]*Grt_tmp[m + (l*(l+1))/2];
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
					if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
				}
				m=l;
				chir_tmp[l + l*nk]=-2.0*chir_tmp[l + l*nk]*Grt_tmp[l + (l*(l+1))/2];
				if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// garder en memoire chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
					chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chir_tmp[m + l*nk];
			}
		}
		
		delete [] chir_tmp;
		delete [] Grt_tmp;
	}
	
	chir_tmp=new double[nk2];
	
	j=nw-1;
	// calcul de chi0(q,tau=beta/2)
	// calcul de G(r,beta)
	t=1.0/(2.0*tem);
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
		{
			if (ekmu[m + (l*(l+1))/2]>=0.0)
				chir_tmp[m + l*nk]=exp(-t*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
			else
				chir_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + (l*(l+1))/2])/((exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)*nk02);
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
	{
		for(m=0; m<nk; m++)
			chir_tmp[m + l*nk]=2.0*chir_tmp[m + l*nk]*chir_tmp[m + l*nk];
		for(m=0; m<=l; m++)	
			if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chir_tmp[m + l*nk];
	}
	
	delete [] Gkt_tmp;
	delete [] ekmu;
	delete [] chir_tmp;
	fftw_destroy_plan(fftplan_k);
	
	double *coeffs;
	dcomplex *chi_qt=new dcomplex[nw0];
	fftw_plan fftplan_w=fftw_plan_dft_1d(nw0, reinterpret_cast<fftw_complex *>(chi_qt), reinterpret_cast<fftw_complex *>(chi_qt), FFTW_BACKWARD, flag);	
	
	delete [] chi_qt;
	
	double *tau;
	double *chiqt_tmp;
	
	int NS0=nw;
	chi0_inf=new double[nk8th];
	chi0_inf2=new double[nk8th];
	//	chi0_inf3=new double[nk8th];
	
	char name[300];
	
	cout<<setprecision(14)<<setiosflags(ios::left);
	
	double wn, a, b, c, d, x1, x2, integ;	
	
#pragma omp parallel private(l,m,j,tau,chiqt_tmp,coeffs,chi_qt,integ,a,b,c,d,x1,x2,wn)
	{
		tau=new double[NS0];
		chiqt_tmp=new double[NS0];
		coeffs=new double[4*(NS0-1)];
		chi_qt=new dcomplex[nw0];
		double *FP=new double[2];
		double x;
		
		double *coeffs2=new double[4*(NS0-1)];

#pragma omp for
		for (l=0; l<nk; l++)
			for (m=0; m<=l; m++)
			{
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nw0*tem);
					chiqt_tmp[j]=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
				}
					
				coeffs[0]=dtau_chiqt0[m + (l*(l+1))/2];
				coeffs[1]=0;
				spline_coeffs_rel(tau, chiqt_tmp, NS0, coeffs);
				for (j=0; j<nw0/2; j++)
					chi_qt[j]=6.0*coeffs[4*j];
				for (j=nw0/2; j<nw0; j++)
					chi_qt[j]=-chi_qt[nw0-j-1];
				
				fftw_execute_dft(fftplan_w,reinterpret_cast<fftw_complex*>(chi_qt),reinterpret_cast<fftw_complex*>(chi_qt));

				integ=0;
				for (j=0; j<nw0/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x=tau[j+1]-tau[j];
					integ+=a*x*x*x*x/4.0+b*x*x*x/3.0+c*x*x/2.0+d*x;
				}
/*				
				FP[0]=dtau_chiqt0[m + (l*(l+1))/2];
				FP[1]=0;
				spline_coefficients(coeffs2,tau,chiqt_tmp, FP, NS0);
				
				if (l==nk-1 && m==nk-1) 
					for (j=0; j<NS0-1; j++) cout<<setw(30)<<(coeffs2[4*j+3]-coeffs[4*j])/coeffs2[4*j+3]<<setw(30)<<(coeffs2[4*j+2]-coeffs[4*j+1])/coeffs2[4*j+2]<<setw(30)<<(coeffs2[4*j+1]-coeffs[4*j+2])/coeffs2[4*j+1]<<setw(30)<<(coeffs2[4*j]-coeffs[4*j+3])/coeffs2[4*j]<<endl;
								
				for (j=0; j<nw0/2; j++)
					chi_qt[j]=6.0*coeffs2[4*j+3];
				for (j=nw0/2; j<nw0; j++)
					chi_qt[j]=-chi_qt[nw0-j-1];
				
				fftw_execute_dft(fftplan_w,reinterpret_cast<fftw_complex*>(chi_qt),reinterpret_cast<fftw_complex*>(chi_qt));
				
				integ=oned_spline_integral(coeffs2,tau,chiqt_tmp, NS0);
*/
 
/*				
				coeffs[0]=dtau_chiqt0[m + (l*(l+1))/2];
				coeffs[1]=0;				
				spline_coeffs(tau, chiqt_tmp, NS0, coeffs);

				for (j=0; j<nw0/2; j++)
					chi_qt[j]=6.0*coeffs[4*j];
				for (j=nw0/2; j<nw0; j++)
					chi_qt[j]=-chi_qt[nw0-j-1];

				fftw_execute_dft(fftplan_w,reinterpret_cast<fftw_complex*>(chi_qt),reinterpret_cast<fftw_complex*>(chi_qt));
				
				integ=0;
				for (j=0; j<nw0/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j];
					x2=tau[j+1];
					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
*/				
				chiqw_array[m + (l*(l+1))/2]=2*integ;
				
				if ( l!=0 || m!=0)
				{
					chi0_inf[m + (l*(l+1))/2]=-2.0*dtau_chiqt0[m + (l*(l+1))/2];
					for (j=1; j<nw; j++)
					{
						wn=2*j*PI*tem;
						chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chi0_inf[m + (l*(l+1))/2]/(wn*wn)
						+((1.0-cos(2*j*PI/nw0))*real(chi_qt[j])+sin(2*j*PI/nw0)*imag(chi_qt[j]))/(wn*wn*wn*wn);
						//						+real((1.0-exp(dcomplex(0,2*j*PI/nw0)))*chi_qt[j])/(wn*wn*wn*wn);
					}
					
					chi0_inf2[m + (l*(l+1))/2]=2.0*real(chi_qt[nw-1]);
				}
				else
				{
					chi0_inf[0]=0;
					chi0_inf2[0]=0;
					for (j=1; j<nw; j++)
						chiqw_array[ j*nk8th]=0;
				}
			}
		delete [] tau;
		delete [] chiqt_tmp;
		delete [] coeffs;
		delete [] chi_qt;
		delete [] FP;
		delete [] coeffs2;
	}
	
	delete [] dtau_chiqt0;
	fftw_destroy_plan(fftplan_w);
	
	chi0max=0.0;
	
	//	for (l=0; l<nk; l++)
	for (l=nk-1; l>=0; l--)
	{
		for (m=l; m>=0; m--)
		{
			if (chiqw_array[m + (l*(l+1))/2]>chi0max)
			{
				chi0max=chiqw_array[m + (l*(l+1))/2];
				qmax[0]=l;
				qmax[1]=m;
			}
		}
	}
	cout<<"qxmax, qymax: "<<qmax[0]<<",  "<<qmax[1]<<'\n';
	
	if (save_chi0_q)
	{
		fstream file;
		const char *nameFormw="./chi0_qmax_wn_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d_%d_%d.dat";
		const char *nameFormq="./chi0_wn0_q_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";
		
		l=qmax[0];
		m=qmax[1];
		sprintf(name, nameFormw, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nq0, nw0, l, m);
		file.open(name,ios::out);
		file<<setiosflags(ios::left)<<setprecision(16);
		for (j=0; j<nw; j++)
			file<<setw(10)<<j<<chiqw_array[m + (l*(l+1))/2 + j*nk8th]<<'\n';
		file.close();
		
		sprintf(name, nameFormq, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nq0, nw0);
		file.open(name,ios::out);
		file<<setiosflags(ios::left)<<setprecision(16);
		for (l=0; l<nk; l++)
			for (m=0; m<=l; m++)
				file<<setw(10)<<l<<setw(10)<<m<<chiqw_array[m + (l*(l+1))/2]<<'\n';
		file.close();
	}
	
	int q[]={qmax[0],qmax[1]};
	double chitmp=chiqw_discret(q, 0);
	cout<<"erreur relative a (qmax,iwn=0) par rapport au calcul standard:   "<<fabs((chitmp-chi0max)/chitmp)<<'\n';
	
	l=qmax[0];
	m=qmax[1];
	j=1;
	double chitmp2=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	chitmp=chiqw_discret(q, 1);
	cout<<"erreur relative a (qmax,iwn=2*PI*tem) par rapport au calcul standard:   "<<fabs((chitmp-chitmp2)/chitmp)<<'\n';
	//	q[0]=0;
	//	q[1]=0;
	//	chitmp=chiqw_discret(q, 0);
	//	cout<<"erreur relative a (q=0,iwn=0):   "<<fabs((chitmp-chiqw_array[0])/chitmp)<<'\n';	
	traceChi0();
	
	/*	
	 double *chiqw2=new double[nw*nk8th];
	 for (l=0; l<nw*nk8th; l++) chiqw2[l]=chiqw_array[l];
	 
	 //	save_chi();
	 loadChi0();
	 
	 double diff=0;
	 double diffl;
	 for (l=0; l<nw*nk8th; l++)
	 {
	 diffl=chiqw2[l]-chiqw_array[l];
	 diff+=abs(diffl);
	 if (diffl && l<20)
	 {
	 cout<<"l:  "<<l<<"   diffl:  "<<diffl<<'\n';
	 }
	 }
	 
	 cout<<"diff/(nw*nk8th):  "<<diff/(nw*nk8th)<<'\n';
	 
	 delete [] chiqw2;
	 */	
	
	/*	
	 double chitmp1, chitmp2, qx, par[2];
	 
	 l=nk-1;
	 m=nk-1;
	 q[0]=l;
	 q[1]=m;
	 qx=q[0]*PI/(nk-1);
	 par[0]=q[1]*PI/(nk-1);
     */
	 /*	
	 cout<<"chi(pi,pi,iwn):\n"<<setprecision(14)<<setiosflags(ios::left);
	 for (j=0; j<10; j++)
	 {
	 chitmp=chiqw_discret(q, j);
	 chitmp1=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	 par[1]=2*j*PI*tem;
	 chitmp2=chiqw(qx, par);
	 cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //		cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';
	 }	
	 for (j=nw-10; j<nw; j++)
	 {
	 chitmp=chiqw_discret(q, j);
	 chitmp1=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	 par[1]=2*j*PI*tem;
	 chitmp2=chiqw(qx, par);
	 cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //		cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';
	 }
	 for (j=nw; j<nw+10; j++)
	 {
	 wn=2*j*PI*tem;
	 chitmp=chiqw_discret(q, j);
	 par[1]=2*j*PI*tem;
	 chitmp2=chiqw(qx, par);
	 chitmp1=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
	 cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //		cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';
	 }
	 for (j=nw0-10; j<nw0; j++)
	 {
	 wn=2*j*PI*tem;
	 chitmp=chiqw_discret(q, j);
	 par[1]=2*j*PI*tem;
	 chitmp2=chiqw(qx, par);
	 chitmp1=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
	 cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //		cout<<setw(15)<<j<<setw(30)<<chitmp1<<setw(30)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';
	 }	
	 */
	
	/*	
	 int dk=(nk-1)/8;
	 
	 cout<<setprecision(14)<<setiosflags(ios::left);
	 for (j=0; j<1; j++)
	 {
	 cout<<"j:  "<<j<<'\n';
	 par[1]=2*j*PI*tem;
	 for (l=0; l<nk; l+=dk)
	 {
	 q[0]=l;
	 qx=l*PI/(nk-1);
	 for (m=0; m<=l; m+=dk)
	 {
	 q[1]=m;
	 chitmp=chiqw_discret(q, j);
	 chitmp1=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	 par[0]=m*PI/(nk-1);
	 chitmp2=chiqw(qx, par);
	 cout<<setw(10)<<l<<setw(10)<<m<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<chitmp1<<setw(25)<<chitmp<<(chiqw_array[m + (l*(l+1))/2 + j*nk8th]-chitmp)/chitmp<<'\n';				
	 }
	 }
	 cout<<'\n';
	 }
	 
	 
	 j=nw-1;
	 cout<<"j:  "<<j<<'\n';
	 par[1]=2*j*PI*tem;
	 for (l=0; l<nk; l+=dk)
	 {
	 q[0]=l;
	 qx=l*PI/(nk-1);
	 for (m=0; m<=l; m+=dk)
	 {
	 q[1]=m;
	 wn=2*j*PI*tem;
	 chitmp1=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	 chitmp=chiqw_discret(q, j);
	 par[0]=m*PI/(nk-1);
	 chitmp2=chiqw(qx, par);
	 cout<<setw(10)<<l<<setw(10)<<m<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';			
	 //			cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<chitmp1<<setw(25)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';				
	 }
	 }
	 cout<<'\n';
	 
	 j=nw;
	 cout<<"j:  "<<j<<'\n';
	 par[1]=2*j*PI*tem;
	 for (l=0; l<nk; l+=dk)
	 {
	 q[0]=l;
	 qx=l*PI/(nk-1);
	 for (m=0; m<=l; m+=dk)
	 {
	 q[1]=m;
	 wn=2*j*PI*tem;
	 chitmp1=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
	 chitmp=chiqw_discret(q, j);
	 par[0]=m*PI/(nk-1);
	 chitmp2=chiqw(qx, par);
	 cout<<setw(10)<<l<<setw(10)<<m<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';
	 //			cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<chitmp1<<setw(25)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';				
	 }
	 }
	 cout<<'\n';
	 
	 j=nw0-1;
	 cout<<"j:  "<<j<<'\n';
	 par[1]=2*j*PI*tem;
	 for (l=0; l<nk; l+=dk)
	 {
	 q[0]=l;
	 qx=l*PI/(nk-1);
	 for (m=0; m<=l; m+=dk)
	 {
	 q[1]=m;
	 wn=2*j*PI*tem;
	 chitmp1=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
	 chitmp=chiqw_discret(q, j);
	 par[0]=m*PI/(nk-1);
	 chitmp2=chiqw(qx, par);
	 cout<<setw(10)<<l<<setw(10)<<m<<setw(30)<<chitmp1<<setw(30)<<chitmp<<setw(30)<<(chitmp1-chitmp)/chitmp<<setw(30)<<chitmp2<<(chitmp2-chitmp)/chitmp<<'\n';			
	 //			cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<chitmp1<<setw(25)<<chitmp<<(chitmp1-chitmp)/chitmp<<'\n';				
	 }
	 }
	 */
	 /*		
	 double dGtmp, dchi;
	 cout<<setprecision(14)<<setiosflags(ios::left);
	 //	dGtmp=(-Grt[(nw0-1)*nk8th]-Grt[0]-1.0)*nw0*tem;
	 //	cout<<setw(25)<<dGtmp<<setw(25)<<dtau_Grt0m[0]<<fabs((dtau_Grt0m[0]-dGtmp)/dtau_Grt0m[0])<<'\n';	
	 for (l=0; l<nk; l++)
	 for (m=0; m<=l; m++)
	 {
	 dchi=(chirt[m + (l*(l+1))/2+ nk8th]-chirt[m + (l*(l+1))/2])*nw0*tem;
	 cout<<setw(25)<<dchi<<setw(25)<<dtau_chirt0[m + (l*(l+1))/2]<<dtau_chirt0[m + (l*(l+1))/2]-dchi<<'\n';
	 //			dGtmp=(-Grt[m + (l*(l+1))/2+ (nw0-1)*nk8th]-Grt[m + (l*(l+1))/2])*nw0*tem;
	 //			cout<<setw(25)<<dGtmp<<setw(25)<<dtau_Grt0m[m + (l*(l+1))/2]<<fabs((dtau_Grt0m[m + (l*(l+1))/2]-dGtmp)/dtau_Grt0m[m + (l*(l+1))/2])<<'\n';
	 //			dGtmp=(Grt[m + (l*(l+1))/2+ nk8th]-Grt[m + (l*(l+1))/2])*nw0*tem;
	 //			cout<<setw(25)<<dGtmp<<setw(25)<<dtau_Grt0p[m + (l*(l+1))/2]<<dtau_Grt0p[m + (l*(l+1))/2]-dGtmp<<'\n';
	 }
	 */
	
}


/*
// calculer chi0 par interpolation par splines cubiques et par transformee de Fourier rapide en ne gardant pas les nw0/(2^r) derniere frequences pour sauver de la memoire
void chi0::calc_chiqw_FFT_spline_optim(int nq0, int nw0, int r)
{
	free_chi_array();
	
	int nk=nq0/2+1;
	int nk02=nq0*nq0;
	int nk2=nk*nk;
	
	int nk8th=(nk*(nk+1))/2;
	
	int nw=nw0/2+1;
	
	int Dnw=nw0/((int)pow(2.0,r));
	int nwc=nw-Dnw;
	
	nwr=nwc+Dnw/2;
	
	cout<<"nwc: "<<nwc<<'\n';
	cout<<"nwr: "<<nwr<<'\n';	
	
	nbq=nk;
	//	nwmax=nw-1;
	
	cout<<"calcul de chi0 par calc_chiqw_FFT_spline_optim()\n";
	
	//	cout<<"nwmax avant set_grid: "<<nwmax<<endl;
	
	set_grid(nq0, nw0);
	
	//	cout<<"nwmax apres set_grid: "<<nwmax<<endl;

	long int array_size=(long int) nwr*nk8th;
	
	chiqw_array=new double[array_size];
//	chiqw_array=new double[nwr*nk8th];
	
	cout<<"nq0: "<<nq0<<"  nw0: "<<nw0<<'\n';
	
	double *Grt_tmp=new double[nk8th];
	double *Gkt_tmp=new double[nk8th];
	double *ekmu=new double[nk8th];
	double *dtau_chiqt0=new double[nk8th];
	double *dtau_Grt0m=new double[nk8th];
	
	double *chir_tmp;
	chir_tmp=new double[nk2];
	
	long int j,l,m;
	double kx, ky, t;
	
	free_dtau_Grt0p();
	dtau_Grt0p=new double[nk8th];
	dtau_Grt0m[0]=-mu0;
	for (j=1; j<nk8th; j++) dtau_Grt0m[j]=er[j];
	
	fftw_plan fftplan_k;
	unsigned flag = FFTW_MEASURE;
	
	// Frequency transform 
	fftw_r2r_kind kind = FFTW_REDFT00;
	
	fftplan_k=fftw_plan_r2r_2d(nk, nk, chir_tmp, chir_tmp, kind, kind, flag);
	
	cout<<setiosflags(ios::left);
	
	//#pragma omp parallel for private(l,m,kx,ky)
	for (l=0; l<nk; l++)
	{
		for(m=0; m<=l; m++)
		{
			ekmu[m + (l*(l+1))/2]=ek[m+l*nbk]-mu0;
			Gkt_tmp[m + (l*(l+1))/2]=-1.0/(exp(-ekmu[m + (l*(l+1))/2]/tem)+1.0);
			//			cout<<setw(20)<<l<<setw(20)<<m<<setw(20)<<ekmu[m + (l*(l+1))/2]<<Gkt_tmp[m + (l*(l+1))/2]<<endl;
		}
	}
	
	//calcul de d/dtau G(r,tau)(tau=0) (pour le calcul de d/dtau chi(q,tau)(tau=0) necessaire aux conditions aux 
	//frontieres des splines de chi(q,tau))	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=-ekmu[m + (l*(l+1))/2]*Gkt_tmp[m + (l*(l+1))/2]/nk02;
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
		chir_tmp[l + l*nk]=-ekmu[l + (l*(l+1))/2]*Gkt_tmp[l + (l*(l+1))/2]/nk02;
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
		for(m=0; m<=l; m++)
		{
			dtau_Grt0m[m + (l*(l+1))/2]-=chir_tmp[m + l*nk];
			dtau_Grt0p[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
	
//	if (keep_chirt)		chirt=new double[nwr*nk8th];
	if (keep_chirt)		chirt=new double[array_size];
	
	// calcul de chi0(q,tau=0) et d/dtau chi(q,tau)(tau=0)
	// calcul de G(r,tau=0)	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=Gkt_tmp[m + (l*(l+1))/2]/nk02;
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
		chir_tmp[l + l*nk]=Gkt_tmp[l + (l*(l+1))/2]/nk02;
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
		for(m=0; m<=l; m++)
		{
			Grt_tmp[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
	
	//calcul de chi(q,tau=0)	
	chir_tmp[0]=-2.0*chir_tmp[0]*(chir_tmp[0]+1.0);
	if (keep_chirt)	chirt[0]=chir_tmp[0];
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=-2.0*chir_tmp[m + l*nk]*chir_tmp[m + l*nk];
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
			if (keep_chirt)	chirt[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
		m=l;
		chir_tmp[l + l*nk]=-2.0*chir_tmp[l + l*nk]*chir_tmp[l + l*nk];
		if (keep_chirt)	chirt[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	//	fftw_execute(fftplan_k);
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			chiqw_array[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	free_dtau_chirt0();
	dtau_chirt0=new double[nk8th];
	
	//calcul de d/dtau chi(q,tau)(tau=0)
	chir_tmp[0]=-2.0*(dtau_Grt0p[0]*(1.0+Grt_tmp[0])+dtau_Grt0m[0]*Grt_tmp[0]);
	dtau_chirt0[0]=chir_tmp[0];
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>0; l--)
	{
		for(m=0; m<l; m++)
		{
			chir_tmp[m + l*nk]=-2.0*(dtau_Grt0p[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]+dtau_Grt0m[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]);
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
			dtau_chirt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
		}
		m=l;
		chir_tmp[l + l*nk]=-2.0*(dtau_Grt0p[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]+dtau_Grt0m[m + (l*(l+1))/2]*Grt_tmp[m + (l*(l+1))/2]);
		dtau_chirt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}	
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	//#pragma omp parallel for private(l,m)
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			dtau_chiqt0[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
	}
	
	delete [] dtau_Grt0m;
	delete [] chir_tmp;
	delete [] Grt_tmp;
	
	
	
	//calcul de chi0(q,tau) pour 0<tau<beta/2	
#pragma omp parallel private(j,t,l,m,chir_tmp,Grt_tmp)	
	{
		chir_tmp=new double[nk2];
		Grt_tmp=new double[nk8th];
		
#pragma omp for		
		for (j=1; j<nwc;j++)
		{
			t=j/(nw0*tem);
			// calcul de G(r,t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=exp(-t*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + (l*(l+1))/2])/(exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)/nk02;
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// calcul de G(r,-t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					// conserver G(r,t)
					Grt_tmp[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
					
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=-exp((t-1.0/tem)*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=exp(t*ekmu[m + (l*(l+1))/2])/((exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)*nk02);
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			//calcul de chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<l; m++)
				{
					chir_tmp[m + l*nk]=-2.0*chir_tmp[m + l*nk]*Grt_tmp[m + (l*(l+1))/2];
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
					if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
				}
				m=l;
				chir_tmp[l + l*nk]=-2.0*chir_tmp[l + l*nk]*Grt_tmp[l + (l*(l+1))/2];
				if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// garder en memoire chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
					chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chir_tmp[m + l*nk];
			}
		}
		
#pragma omp for		
//		for (j=1; j<(Dnw/2);j++)
		for (j=nwc; j<nwr; j++)
		{
//			t=(nwc-1+2*j)/(nw0*tem);
			t=(2*j-nwc+1)/(nw0*tem);
			// calcul de G(r,t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=exp(-t*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + (l*(l+1))/2])/(exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)/nk02;
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// calcul de G(r,-t)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
				{
					// conserver G(r,t)
					Grt_tmp[m + (l*(l+1))/2]=chir_tmp[m + l*nk];
					
					if (ekmu[m + (l*(l+1))/2]>=0.0)
						chir_tmp[m + l*nk]=-exp((t-1.0/tem)*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
					else
						chir_tmp[m + l*nk]=exp(t*ekmu[m + (l*(l+1))/2])/((exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)*nk02);
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
				}
			}
			
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			//calcul de chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<l; m++)
				{
					chir_tmp[m + l*nk]=-2.0*chir_tmp[m + l*nk]*Grt_tmp[m + (l*(l+1))/2];
					chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
//					if (keep_chirt)	chirt[m + (l*(l+1))/2+(nwc-1+j)*nk8th]=chir_tmp[l + m*nk];
					if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
				}
				m=l;
				chir_tmp[l + l*nk]=-2.0*chir_tmp[l + l*nk]*Grt_tmp[l + (l*(l+1))/2];
//				if (keep_chirt)	chirt[m + (l*(l+1))/2+(nwc-1+j)*nk8th]=chir_tmp[l + m*nk];
				if (keep_chirt)	chirt[m + (l*(l+1))/2+j*nk8th]=chir_tmp[l + m*nk];
			}
			
			fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
			
			// garder en memoire chi0(q,tau)
			for (l=nk-1; l>=0; l--)
			{
				for(m=0; m<=l; m++)
					chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chir_tmp[m + l*nk];
//					chiqw_array[m + (l*(l+1))/2 + (nwc-1+j)*nk8th]=chir_tmp[m + l*nk];
			}
		}		
		
		delete [] chir_tmp;
		delete [] Grt_tmp;
	}
	
	chir_tmp=new double[nk2];
	

	// calcul de chi0(q,tau=beta/2)
	// calcul de G(r,beta)
	t=1.0/(2.0*tem);
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
		{
			if (ekmu[m + (l*(l+1))/2]>=0.0)
				chir_tmp[m + l*nk]=exp(-t*ekmu[m + (l*(l+1))/2])*Gkt_tmp[m + (l*(l+1))/2]/nk02;
			else
				chir_tmp[m + l*nk]=-exp((1.0/tem-t)*ekmu[m + (l*(l+1))/2])/((exp(ekmu[m + (l*(l+1))/2]/tem)+1.0)*nk02);
			chir_tmp[l + m*nk]=chir_tmp[m + l*nk];
		}
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=0; l<nk; l++)
	{
		for(m=0; m<nk; m++)
			chir_tmp[m + l*nk]=2.0*chir_tmp[m + l*nk]*chir_tmp[m + l*nk];
		for(m=0; m<=l; m++)	
			if (keep_chirt)	chirt[m + (l*(l+1))/2+(nwr-1)*nk8th]=chir_tmp[l + m*nk];
//			if (keep_chirt)	chirt[m + (l*(l+1))/2+(nwc-1+Dnw/2)*nk8th]=chir_tmp[l + m*nk];
	}
	
	fftw_execute_r2r(fftplan_k,chir_tmp,chir_tmp);
	
	for (l=nk-1; l>=0; l--)
	{
		for(m=0; m<=l; m++)
			chiqw_array[m + (l*(l+1))/2 + (nwr-1)*nk8th]=chir_tmp[m + l*nk];
//			chiqw_array[m + (l*(l+1))/2 + (nwc-1+Dnw/2)*nk8th]=chir_tmp[m + l*nk];
	}
	
	delete [] Gkt_tmp;
	delete [] ekmu;
	delete [] chir_tmp;
	fftw_destroy_plan(fftplan_k);
	
	double *coeffs;
	dcomplex *chi_qt=new dcomplex[nw0];
	fftw_plan fftplan_w=fftw_plan_dft_1d(nw0, reinterpret_cast<fftw_complex *>(chi_qt), reinterpret_cast<fftw_complex *>(chi_qt), FFTW_BACKWARD, flag);	
	
	delete [] chi_qt;
	
	double *tau;
	double *chiqt_tmp;
	
	int NS0=nwr;
	chi0_inf=new double[nk8th];
	chi0_inf2=new double[nk8th];
	//	chi0_inf3=new double[nk8th];
	
	double chiqtmp1, chiqtmp2;
	
	char name[300];
	
	cout<<setprecision(14)<<setiosflags(ios::left);
	
	double wn, a, b, c, d, x1, x2, integ;	
	
#pragma omp parallel private(l,m,j,tau,chiqt_tmp,coeffs,chi_qt,integ,a,b,c,d,x1,x2,wn)
	{
		tau=new double[NS0];
		chiqt_tmp=new double[NS0];
		coeffs=new double[4*(NS0-1)];
		chi_qt=new dcomplex[nw0];
		
#pragma omp for
		for (l=0; l<nk; l++)
			for (m=0; m<=l; m++)
			{
				for (j=0; j<nwc; j++)
				{
					tau[j]=j/(nw0*tem);
					chiqt_tmp[j]=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
				}
				for (j=nwc; j<NS0; j++)
				{
					tau[j]=(2*j-nwc+1)/(nw0*tem);
					chiqt_tmp[j]=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
				}
				
				coeffs[0]=dtau_chiqt0[m + (l*(l+1))/2];
				coeffs[1]=0;				
				spline_coeffs(tau, chiqt_tmp, NS0, coeffs);
				
				for (j=0; j<nwc-1; j++)
				{
					chi_qt[j]=6.0*coeffs[4*j];
//					if (l==nk-1 && m==nk-1) cout<<setw(10)<<j<<chi_qt[j]<<'\n';
				}
				j=nwc-1;
				chi_qt[2*j-nwc+1]=6.0*(2*coeffs[4*j]+coeffs[4*(j-1)])/3;
				for (j=nwc; j<NS0-1; j++)
				{
					chi_qt[2*j-nwc]=6.0*(3*coeffs[4*(j-1)]+coeffs[4*j])/4;
					chi_qt[2*j-nwc+1]=6.0*(coeffs[4*(j-1)]+3*coeffs[4*j])/4;
				}
				j=NS0-2;
				chi_qt[2*j-nwc+2]=6.0*coeffs[4*j];

//				for (j=nwc-1; j<NS0-1; j++)
//				{
//					chi_qt[2*j-nwc+1]=6.0*coeffs[4*j];
//					chi_qt[2*j-nwc+2]=6.0*coeffs[4*j];
//					if (l==nk-1 && m==nk-1) 
//					{
//						cout<<setw(10)<<2*j-nwc+1<<chi_qt[2*j-nwc+1]<<'\n';
//						cout<<setw(10)<<2*j-nwc+2<<chi_qt[2*j-nwc+2]<<'\n';
//					}
//				}
				
				for (j=nw0/2; j<nw0; j++)
					chi_qt[j]=-chi_qt[nw0-j-1];
				
				fftw_execute_dft(fftplan_w,reinterpret_cast<fftw_complex*>(chi_qt),reinterpret_cast<fftw_complex*>(chi_qt));
				
				integ=0;
				for (j=0; j<NS0-1; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j];
					x2=tau[j+1];
					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
				
				chiqw_array[m + (l*(l+1))/2]=2*integ;
				
				if ( l!=0 || m!=0)
				{
					chi0_inf[m + (l*(l+1))/2]=-2.0*dtau_chiqt0[m + (l*(l+1))/2];
					for (j=1; j<NS0; j++)
					{
						wn=2*j*PI*tem;
						chiqw_array[m + (l*(l+1))/2 + j*nk8th]=chi0_inf[m + (l*(l+1))/2]/(wn*wn) +real((1.0-exp(dcomplex(0,2*j*PI/nw0)))*chi_qt[j]/(wn*wn*wn*wn));
					}
					
					chi0_inf2[m + (l*(l+1))/2]=2.0*real(chi_qt[nw-1]);
					
//					if (l==nk-1 && m==nk-1)
//						for (j=NS0; j<nw; j++)
//						{
//							wn=2*j*PI*tem;
//							chiqtmp1=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+real((1.0-exp(dcomplex(0,2*j*PI/nw0)))*chi_qt[j]/(wn*wn*wn*wn));
//							chiqtmp2=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
//							cout<<setw(10)<<j<<setw(30)<<chiqtmp1<<setw(30)<<chiqtmp2<<(chiqtmp1-chiqtmp2)/chiqtmp1<<'\n';
//						}
				}
				else
				{
					chi0_inf[0]=0;
					chi0_inf2[0]=0;
					for (j=1; j<NS0; j++)
						chiqw_array[ j*nk8th]=0;
				}
			}
		delete [] tau;
		delete [] chiqt_tmp;
		delete [] coeffs;
		delete [] chi_qt;
	}
	
	delete [] dtau_chiqt0;
	fftw_destroy_plan(fftplan_w);
	
	chi0max=0.0;
	
	//	for (l=0; l<nk; l++)
	for (l=nk-1; l>=0; l--)
	{
		for (m=l; m>=0; m--)
		{
			if (chiqw_array[m + (l*(l+1))/2]>chi0max)
			{
				chi0max=chiqw_array[m + (l*(l+1))/2];
				qmax[0]=l;
				qmax[1]=m;
			}
		}
	}
	cout<<"qxmax, qymax: "<<qmax[0]<<",  "<<qmax[1]<<'\n';
	
	fstream file;
	const char *nameFormw="./chi0_qmax_wn_%6.4f_%6.4f_%6.4f_%6.4f_%2d_%2d_%1d_%1d.dat";
	const char *nameFormq="./chi0_wn0_q_%6.4f_%6.4f_%6.4f_%6.4f_%2d_%2d.dat";
	
	l=qmax[0];
	m=qmax[1];
	sprintf(name, nameFormw, (double)tp, (double)tpp, (double)density, (double)tem, nq0, nw0, l, m);
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(16);
	for (j=0; j<NS0; j++)
		file<<setw(10)<<j<<chiqw_array[m + (l*(l+1))/2 + j*nk8th]<<'\n';
	for (j=NS0; j<nw; j++)
	{
		wn=2*j*PI*tem;
		file<<setw(10)<<j<<chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn)<<'\n';
	}
	file.close();
	
	sprintf(name, nameFormq, (double)tp, (double)tpp, (double)density, (double)tem, nq0, nw0);
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(16);
	for (l=0; l<nk; l++)
		for (m=0; m<=l; m++)
			file<<setw(10)<<l<<setw(10)<<m<<chiqw_array[m + (l*(l+1))/2]<<'\n';
	file.close();
	
	int q[]={qmax[0],qmax[1]};
	double chitmp=chiqw_discret(q, 0);
	cout<<"erreur relative a (qmax,iwn=0) par rapport au calcul standard:   "<<fabs((chitmp-chi0max)/chitmp)<<'\n';
	//	q[0]=0;
	//	q[1]=0;
	//	chitmp=chiqw_discret(q, 0);
	//	cout<<"erreur relative a (q=0,iwn=0):   "<<fabs((chitmp-chiqw_array[0])/chitmp)<<'\n';	
	traceChi0();
}
*/
 
//! trace of chi0(q,iwn)
void chi0::traceChi0()
{
	double trChitmp=0.0;
	// double trChitmp2=0.0;
	
	long int j, l, m;
	int w;
	int Nq=4*(nbq-1)*(nbq-1);
	int nk8th=(nbq*(nbq+1))/2;
	double tr1;
	int N1=nwmax+1;
	if (nwr) N1=nwr+1;

#pragma omp parallel for private(j,l,m,w) reduction(+:trChitmp)
		for (j=0; j<N1; j++)
			for (l=0; l<nbq; l++)
				for (m=0; m<=l; m++)
				{
					w=16;
					if (j==0 || j==nwmax) w=w/2;
					if (l==0 || l==nbq-1) w=w/2;
					if (m==0 || m==nbq-1) w=w/2;
					if (m==l) w=w/2;
					
					trChitmp+=w*tem*chiqw_array[m + (l*(l+1))/2 + j*nk8th]/Nq;
					
					//			 q[0]=l;
					//			 q[1]=m;			 
					//			 trChitmp2+=w*tem*chiqw_discret(q,j)/Nq;
					
				}
		tr1=trChitmp;
	
		double chi0tmp, wn;
		if (nwr)
		{
			for (j=N1; j<=nwmax; j++)
				for (l=0; l<nbq; l++)
					for (m=0; m<=l; m++)
					{
						w=16;
						if (j==0 || j==nwmax) w=w/2;
						if (l==0 || l==nbq-1) w=w/2;
						if (m==0 || m==nbq-1) w=w/2;
						if (m==l) w=w/2;
						
						wn=2*PI*j*tem;
						chi0tmp=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
						trChitmp+=w*tem*chi0tmp/Nq;
						
						//			 q[0]=l;
						//			 q[1]=m;			 
						//			 trChitmp2+=w*tem*chiqw_discret(q,j)/Nq;
						
					}
		}
		
		for (j=nwmax+1; j<2*nwmax; j++)
			for (l=0; l<nbq; l++)
				for (m=0; m<=l; m++)
				{
					w=16;
					if (j==0 || j==nwmax) w=w/2;
					if (l==0 || l==nbq-1) w=w/2;
					if (m==0 || m==nbq-1) w=w/2;
					if (m==l) w=w/2;
					
					wn=2*PI*j*tem;
					chi0tmp=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
					trChitmp+=w*tem*chi0tmp/Nq;
					
					//			 q[0]=l;
					//			 q[1]=m;			 
					//			 trChitmp2+=w*tem*chiqw_discret(q,j)/Nq;
					
				}

 
	double tr_exact=density-density*density/2.0;	
	cout<<setprecision(14);
	cout<<"trace(chi0(q,iwn)):   "<<trChitmp<<"	n-n^2/2:   "<<tr_exact<<"	difference relative:  "<<(tr_exact - trChitmp)/tr_exact<<'\n';
// cout<<"trace(chi0(q,iwn))(+correction asymptotique):   "<<tr1+2*(trChitmp-tr1)<<"	n-n^2/2:   "<<density-density*density/2.0<<"	difference:  "<<density-density*density/2.0 - tr1-2*(trChitmp-tr1)<<'\n';
	double tr_corr=trChitmp+N1*(tr1-trChitmp)/(N1-2*nwmax);
	cout<<"trace(chi0(q,iwn))(+correction asymptotique):   "<<tr_corr<<"	n-n^2/2:   "<<tr_exact<<"	difference relative:  "<<(tr_exact - tr_corr)/tr_exact<<'\n';	
// cout<<"trace(chi0_exact(q,iwn)):   "<<trChitmp2<<"	n-n^2/2:   "<<density-density*density/2.0<<"	difference:  "<<density-density*density/2.0 - trChitmp2<<'\n';
 
}


// Calculate chi(q,iwn) on the finite lattice
double chi0::chiqw_discret(int *q, int n)
{
	int nq0=2*(nbq-1);

	double params[4];
	params[1]=q[0]*2*PI/nq0;
	params[2]=q[1]*2*PI/nq0;
	params[3]=2*n*PI*tem;
	
	double kx, chitmp=0;
	int j,l;
	for (j=-nq0/2; j<nq0/2; j++ )
	{
		kx=j*2*PI/nq0;
		for (l=-nq0/2; l<nq0/2; l++)		
		{
			params[0]=l*2*PI/nq0;
			chitmp-=2*chi0IntegRe(kx, params)/(nq0*nq0);
		}
	}

return chitmp;
}


// return integrand of chi0
double chi0::chi0IntegRe(double kx, double params[])
{
	double ky=params[0], qx=params[1], qy=params[2], wn=params[3];
	double integ;
	double ektmp, ekq;

	ektmp=dispk(kx, ky);
	ekq=dispk(kx+qx, ky+qy);

	double delta_ekq=ektmp-ekq;
	if (qx==0.0 && qy==0.0) delta_ekq=0.0;
	double expk=exp((ektmp-mu0)/tem), expkm;

	double eps=1e-12;
	
	if ( wn )
	{
		if (delta_ekq)
			integ=(1.0/( expk + 1.0 )  - 1.0/( exp((ekq-mu0)/tem) + 1.0 ) )*delta_ekq/(wn*wn+delta_ekq*delta_ekq);
		else
			integ=0;
	}
	else if (fabs(delta_ekq)>eps)
		integ=(1.0/( expk + 1.0 )  - 1.0/(exp((ekq-mu0)/tem)+1.0) )/delta_ekq;
	else
	{
		if ((ektmp-mu0)<0)
			integ=-(1.0/tem)*expk/( ( expk + 1.0 )*( expk + 1.0 ) );
		else
		{
			expkm=exp(-(ektmp-mu0)/tem);
			integ=-(1.0/tem)*expkm/( ( expkm + 1.0 )*( expkm + 1.0 ) );
		}
	}
	

	return integ;
}
