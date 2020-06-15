
/*
 *  chi.cpp
 *  TPSC
 *
 *  Created by Dominic Bergeron
 *
 */

#include "includeDef.h"
#include "chi.h"
#include <cstdlib>


//! Constructor
chi::chi(int nk0, int Nw0, double n, double T, double  *params):chi0(nk0, Nw0, n, T, params)
{
	if (n>1.0) nsr=2.0-n;
	else nsr=n;
	
	Usp=0.0;
	Uch=0.0;
	dblOcc=0.0;
	
	vertex_rt_array=NULL;

}

//! Destructor
chi::~chi()
{
	if (vertex_rt_array)
		delete [] vertex_rt_array;
}

//! calculate V(r,t)=3*Usp*chisp(r,tau)+Uch*chich(r,tau)
void chi::calc_chirt()
{
	cout<<"calcul de vertex_rt_array\n";
	
	if (!keep_chirt || Usp==0.0 || Uch==0.0) 
	{
		keep_chirt=true;
		find_Usp_Uch();
	}

	int nw=nwmax+1;
	int nbq2=nbq*nbq;

	int nq0=2*(nbq-1);
	int nw0=2*(nw-1);

	int nk8th=(nbq*(nbq+1))/2;

	long int j, l, m;

	if (vertex_rt_array)	delete [] vertex_rt_array;

	long int array_size=(long int)nw*nk8th;
	
//	vertex_rt_array=new double[nw*nk8th];
	vertex_rt_array=new double[array_size];	

	if (!vertex_rt_array)
	{
		cout<<"echec d'allocation de memoire de vertex_rt_array\n";
		return;
	}

	double *vertex_r_tmp=new double[nbq2];

	fftw_plan fftplan_r, fftplan_t;
	fftw_r2r_kind kind = FFTW_REDFT00;
	unsigned flag = FFTW_MEASURE;

	fftplan_r=fftw_plan_r2r_2d(nbq, nbq, vertex_r_tmp, vertex_r_tmp, kind, kind, flag);

//	fftplan_t=fftw_plan_many_r2r(1,&nw,1,vertex_rt_array,NULL,nk8th,1,vertex_rt_array,NULL,nk8th,1,&kind,FFTW_MEASURE);
	
	delete [] vertex_r_tmp;

	double normFact=tem/(nq0*nq0);

	double chi0tmp, chisptmp, chichtmp;
#pragma omp parallel private(j,l,m,vertex_r_tmp,chi0tmp,chisptmp,chichtmp)
	{
		vertex_r_tmp=new double[nbq2];
#pragma omp for	
		for (j=0; j<nw; j++)
		{		
			for (l=0; l<nbq; l++)
			{
				for (m=0; m<=l; m++)
				{
					chi0tmp=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
	
					chisptmp=chi0tmp/(1.0-Usp*chi0tmp/2.0);
					chichtmp=chi0tmp/(1.0+Uch*chi0tmp/2.0);
	
	//				vertex_r_tmp[m + nbq*l]=normFact*chi0tmp;
	//					vertex_r_tmp[m + nbq*l]=normFact*U*(3.0*Usp*chisptmp + Uch*chichtmp)/8.0;
					vertex_r_tmp[m + nbq*l]=normFact*U*(3.0*Usp*(chisptmp-chi0tmp) + Uch*(chichtmp-chi0tmp))/8.0;
				}
				for (m=0; m<l; m++)
					vertex_r_tmp[l + nbq*m]=vertex_r_tmp[m + nbq*l];
			}
	
//			fftw_execute(fftplan_r);
			fftw_execute_r2r(fftplan_r,vertex_r_tmp,vertex_r_tmp);
	
			for (l=0; l<nbq; l++)
			{
				for (m=0; m<=l; m++)
					vertex_rt_array[m + (l*(l+1))/2 + j*nk8th]=vertex_r_tmp[m + nbq*l];
			}
	
		}
		delete [] vertex_r_tmp;
	}

	fftw_destroy_plan(fftplan_r);
	
	double *vertex_t_tmp=new double[nw];
	
	fftplan_t=fftw_plan_r2r_1d(nw, vertex_t_tmp, vertex_t_tmp, kind, flag);
	
	delete [] vertex_t_tmp;
	
#pragma omp parallel private(vertex_t_tmp, j, l)
	{
		vertex_t_tmp=new double[nw];
#pragma omp for
		for (l=0; l<nk8th; l++)
		{
			for (j=0; j<nw; j++)	
				vertex_t_tmp[j]=vertex_rt_array[l + j*nk8th];
			
			fftw_execute_r2r(fftplan_t,vertex_t_tmp,vertex_t_tmp);
			
			for (j=0; j<nw; j++)	
				vertex_rt_array[l + j*nk8th]=vertex_t_tmp[j];
		}
		delete [] vertex_t_tmp;
	}

//#pragma omp parallel for private(l)
//	for (l=0; l<nk8th; l++)
//	{
//		fftw_execute_r2r(fftplan_t,&vertex_rt_array[l],&vertex_rt_array[l]);
//	}

	fftw_destroy_plan(fftplan_t);

	double chi0rt, Grtp, Grtm, chirt0;
	double TPSC_fact=U*(3.0*Usp + Uch)/8.0;
	
	cout<<setiosflags(ios::left)<<setprecision(12);

//#pragma omp parallel for private(j,l,m,Grtp,Grtm,chi0rt)
#pragma omp parallel for private(j,l,m)
	for (j=0; j<nw; j++)
		for (l=0; l<nbq; l++)
			for (m=0; m<=l; m++)
			{

				vertex_rt_array[m + (l*(l+1))/2 + j*nk8th]+=TPSC_fact*chirt[m + (l*(l+1))/2 + j*nk8th];
//				cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<setw(25)
//						<<vertex_rt_array[m + (l*(l+1))/2 + j*nk8th]<<setw(25)<<chi0rt
//						<<chi0rt-vertex_rt_array[m + (l*(l+1))/2 + j*nk8th]<<'\n';
			}

	free_chirt();
//	free_Grt();

	cout<<"vertex_rt_array calcule\n";

	chirt0=U*(3.0*Usp*(nsr - Usp*nsr*nsr/(2.0*U)) + Uch*(nsr + Usp*nsr*nsr/(2.0*U) - nsr*nsr))/8.0;
	cout<<"vertex_rt(0,0):  "<<vertex_rt_array[0]<<"		valeur exacte:  "<<chirt0<<"		difference:  "
		<<chirt0-vertex_rt_array[0]<<'\n';


/*
	cout<<setiosflags(ios::left)<<setprecision(12);	

	l=0;
	m=0;
	for (j=0; j<100; j++)
//	for (l=0; l<nbq; l++)
//		for (m=0; m<=l; m++)
		{
			cout<<setw(25)<<chirt[m + (l*(l+1))/2+j*nk8th]<<setw(25)<<vertex_rt_array[m + (l*(l+1))/2+j*nk8th]<<chirt[m + (l*(l+1))/2+j*nk8th]-vertex_rt_array[m + (l*(l+1))/2+j*nk8th]<<'\n';
		}
*/
/*	
	double dchi;
	cout<<'\n';
	for (l=0; l<4; l++)
		{
			dchi=dtau_chirt0[l];
			cout<<setw(25)<<TPSC_fact*dchi
                  <<(vertex_rt_array[l + nk8th]-vertex_rt_array[l])*nw0*tem<<'\n';
		}
*/

	if (!keep_chiqw) free_chi_array();
}


//! Compute Usp and Uch using different grids in chi0 to check convergence
void chi::find_Usp_Uch()
{
	cout<<"calcul de Usp et Uch\n";

	int nk0, nw0;
	double difRelsp=0.0, difRelch=0.0;

	int nw=0;
	fstream paramsFile;

					  
	const char *nameForm="./TPSCparams_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";
	char paramName[5];

	sprintf(paramsFileName, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);

/*	
		paramsFile.open(paramsFileName,ios::in);
		
		if (paramsFile)
		{
			paramsFile>>paramName>>Usp>>paramName>>difRelsp>>paramName>>Uch>>paramName>>difRelch>>paramName>>nbq>>paramName>>nw;
			cout<<"Usp:  "<<Usp<<"   difRelsp:  "<<difRelsp<<"  Uch:  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
			cout<<"nbq: "<<nbq<<"   nw: "<<nw<<'\n';

			paramsFile.close();
		}
*/

/*	
	if (nw>0)
	{
		nwmax=nw-1;
		if (nw>5 && nbq>3)
		{
			nw0=2*(nw-1);
			nk0=2*(nbq-1);
//			calc_chiqwFFT(nk0,nw0);
			calc_chiqw_FFT_spline(nk0,nw0);

			double chi0max=chi0max;
			double chispmax=chi0max/(1.0-Usp*chi0max/2.0);

			if (density>1.0) dblOcc=Usp/U*nsr*nsr/4.0+1.0-nsr;
			else dblOcc=Usp/U*nsr*nsr/4.0;

			
			cout<<"chi0max: "<<chi0max<<"  chispmax: "<<chispmax<<"  rapport: "<<chispmax/chi0max<<'\n';
			cout<<"double occupation: "<<dblOcc<<'\n';

			fstream infosFile;

			const char *nameFormInfos="./TPSC_infos_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d.dat";
			char name[100];

			sprintf(name, nameFormInfos, (double)tp, (double)tpp, (double)U, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);

			infosFile.open(name,ios::out);
			infosFile<<"Usp:  "<<Usp<<"   difRelsp:  "<<difRelsp<<'\n';
			infosFile<<"Uch:  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
			infosFile<<"chi0max: "<<chi0max<<"   chispmax: "<<chispmax<<"   rapport: "<<chispmax/chi0max<<'\n';
			infosFile<<"double occupation (ansatz): "<<dblOcc<<'\n';
			infosFile<<"qxmax, qymax:   "<<qmax[0]<<",  "<<qmax[1]<<'\n';

				infosFile.close();

			return;
		}
	}
	else
		paramsFile.clear();
*/

//	cout<<"nwmax: "<<nwmax<<endl;

	int nq0tmp=nbq-1;
	int nw0tmp=nwmax;

	int nbIterMax=10, nbIter=1;
	double Ustmp=0.0, Uctmp=0.0;
	double Usc[]={0.0, 0.0};

	cout<<"nq: "<<nq0tmp<<"    nw: "<<nw0tmp<<'\n';
	save_chi0_q=false;
	find_Usp_Uch(nq0tmp, nw0tmp, Usc);
	Ustmp=Usc[0];
	Uctmp=Usc[1];
	
//	cout<<"nwmax: "<<nwmax<<endl;

	
	nq0tmp=2*nq0tmp;
	nw0tmp=2*nw0tmp;
	cout<<"nq: "<<nq0tmp<<"    nw: "<<nw0tmp<<'\n';
	save_chi0_q=save_chi0_to_file;
	find_Usp_Uch(nq0tmp, nw0tmp, Usc);
	save_chi0_q=false;

	difRelsp=fabs((Usc[0]-Ustmp)/Usc[0]);
	difRelch=fabs((Usc[1]-Uctmp)/Usc[1]);
	cout<<"iteration q: "<<nbIter<<"  difRelsp: "<<difRelsp<<"  difRelch: "<<difRelch<<'\n';

/*
	int Nmax=550000000;
	int nqw;

	nq0tmp=2*nq0tmp;
	nqw=(nq0tmp/2+1)*(nq0tmp/2+2)*(nw0tmp/2+1)/2;
	while ( nbIter<nbIterMax && ( difRelsp>tolq || difRelch>tolq) && nqw<Nmax)
	{
		Ustmp=Usc[0];
		Uctmp=Usc[1];

		find_Usp_Uch(nq0tmp, nw0tmp, Usc);
		nbIter++;
		difRelsp=fabs((Usc[0]-Ustmp)/Usc[0]);
		difRelch=fabs((Usc[1]-Uctmp)/Usc[1]);

		cout<<"iteration q: "<<nbIter<<"  difRelsp: "<<difRelsp<<"  difRelch: "<<difRelch<<'\n';
		nq0tmp=2*nq0tmp;
		nqw=(nq0tmp/2+1)*(nq0tmp/2+2)*(nw0tmp/2+1)/2;
	}

	nbIter=1;

	nq0tmp=nq0tmp/4;
	nw0tmp=2*nw0tmp;
	nqw=(nq0tmp/2+1)*(nq0tmp/2+2)*(nw0tmp/2+1)/2;
	find_Usp_Uch(nq0tmp, nw0tmp, Usc);
	difRelsp=fabs((Usc[0]-Ustmp)/Usc[0]);
	difRelch=fabs((Usc[1]-Uctmp)/Usc[1]);
	cout<<"iteration w: "<<nbIter<<"  difRelsp: "<<difRelsp<<"  difRelch: "<<difRelch<<'\n';

	nw0tmp=2*nw0tmp;
	nqw=(nq0tmp/2+1)*(nq0tmp/2+2)*(nw0tmp/2+1)/2;

	while ( nbIter<nbIterMax && ( difRelsp>tolw || difRelch>tolw) && nqw<Nmax)
	{
		Ustmp=Usc[0];
		Uctmp=Usc[1];

		find_Usp_Uch(nq0tmp, nw0tmp, Usc);
		nbIter++;
		difRelsp=fabs((Usc[0]-Ustmp)/Usc[0]);
		difRelch=fabs((Usc[1]-Uctmp)/Usc[1]);

		cout<<"iteration w: "<<nbIter<<"  difRelsp: "<<difRelsp<<"  difRelch: "<<difRelch<<'\n';
		nw0tmp=2*nw0tmp;
		nqw=(nq0tmp/2+1)*(nq0tmp/2+2)*(nw0tmp/2+1)/2;
	}

	errUs=difRelsp;
	errUc=difRelch;
*/
	
	Usp=Usc[0];
	Uch=Usc[1];

	double chispmax=chi0max/(1.0-Usp*chi0max/2.0);

	if (density>1.0) dblOcc=Usp/U*nsr*nsr/4.0+1.0-nsr;
	else dblOcc=Usp/U*nsr*nsr/4.0;

//	save_chi();

	fstream infosFile;

	const char *nameFormInfos="./TPSC_infos_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";
	char name[100];

	sprintf(name, nameFormInfos, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);

	infosFile<<setprecision(16);
	infosFile.open(name,ios::out);
	infosFile<<"Usp:  "<<Usp<<"   difRelsp:  "<<difRelsp<<'\n';
	infosFile<<"Uch:  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
	infosFile<<"chi0max: "<<chi0max<<"   chispmax: "<<chispmax<<"   rapport: "<<chispmax/chi0max<<'\n';
	infosFile<<"double occupation (ansatz): "<<dblOcc<<'\n';
	infosFile<<"qxmax, qymax:   "<<qmax[0]<<",  "<<qmax[1]<<'\n';

	cout<<"chi0max: "<<chi0max<<"  chispmax: "<<chispmax<<"  rapport: "<<chispmax/chi0max<<'\n';
	cout<<"double occupation: "<<dblOcc<<'\n';

	infosFile.close();

	paramsFile.open(paramsFileName,ios::out);
	if (!paramsFile)
	{
		cout<<"echec de creation du fichier de parametres\n";
		return;
	}

	paramsFile<<setprecision(16);
	paramsFile<<"Usp  "<<Usp<<"   difRelsp:  "<<difRelsp<<'\n';
	paramsFile<<"Uch  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
	paramsFile<<"nbq  "<<nbq<<'\n';
	paramsFile<<"nbw  "<<nwmax+1<<'\n';
	paramsFile<<"dblOcc_ansatz  "<<dblOcc<<'\n';
	paramsFile<<"qxmax  "<<qmax[0]<<'\n';
	paramsFile<<"qymax  "<<qmax[1]<<'\n';
	paramsFile<<"chispmax_over_chi0max  "<<chispmax/chi0max<<'\n';
	paramsFile<<"mu0  "<<mu0<<'\n';

	cout<<"nbq: "<<nbq<<"   nw: "<<nwmax+1<<'\n';

	paramsFile.close();

}

//! Calculate the Fermi velocity at the hot spots (where the Fermi surface cross the magnetic Brillouin zone boundary)
void chi::vF_hot_spots()
{
	double kx1, ky1, kx2, ky2;
	
	kx1=FS[0];
	ky1=FS[1];
	kx2=FS[2];
	ky2=FS[3];
	
	int j=1;
	while ((kx1+ky1-PI)*(kx2+ky2-PI)>0 && j<nkFS-1)
	{
		j++;
		kx1=kx2;
		ky1=ky2;
		kx2=FS[2*j];
		ky2=FS[2*j+1];
	}
	
//	cout<<setw(20)<<kx1+ky1-PI<<kx2+ky2-PI<<endl;
	
	double ktmp[]={(kx1+kx2)/2,(ky1+ky2)/2};
	double gr_ek[2];
	
	grad_ek(ktmp, gr_ek);
	
	vF=sqrt(gr_ek[0]*gr_ek[0]+gr_ek[1]*gr_ek[1]);
}


//! obtain the spin correlation length
void chi::sp_corr_length()
{
	vF_hot_spots();

	if ( qmax[0]!=(nbq-1) || qmax[1]!=(nbq-1) )
		cout<<"chisp incommensurable (formule asymptotique non valable)\n";

	double chisp_max=chi0max/(1.0-Usp*chi0max/2.0);

	int nk8th=(nbq*(nbq+1))/2;

	int l=qmax[0], m=0;
	double *chi0tmp;
	chi0tmp=&( chiqw_array[(l*(l+1))/2] );

	double chisptmp=chi0tmp[0]/(1.0-Usp*chi0tmp[0]/2.0);
	double chisptmp2;

	while (chisptmp<(chisp_max/2.0) && m<=qmax[1])
	{
		m++;
		chisptmp2=chisptmp;
		chisptmp=chi0tmp[m]/(1.0-Usp*chi0tmp[m]/2.0);
	}
	double q0=(m-1)*PI/(nbq-1);
	double qm=(qmax[1])*PI/(nbq-1);

	double qHM=PI/(nbq-1)*(chisp_max/2.0 - chisptmp2)/(chisptmp-chisptmp2)+q0;

	spCorrLgth_y=1.0/(qm-qHM);

	fstream infosFile;

	const char *nameFormInfos="./TPSC_infos_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat";
	char name[100];

	sprintf(name, nameFormInfos, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);
	infosFile.open(name,ios::out| ios::app);

	infosFile<<setprecision(16);
	cout<<"vitesse de Fermi aux points chauds:  "<<vF<<endl;
	infosFile<<"vitesse de Fermi aux points chauds:  "<<vF<<endl;
	
	cout<<"longueur de correlation de spin (1/(larg.mi-haut.)) (axe y):  "<<spCorrLgth_y<<'\n';
	cout<<"xi_sp/xi_th:  "<<spCorrLgth_y*PI*tem/vF<<'\n';
	cout<<"xi_sp*T:  "<<spCorrLgth_y*tem<<'\n';
	infosFile<<"longueur de correlation de spin (1/(larg.mi-haut.)) (axe y):  "<<spCorrLgth_y<<'\n';
	infosFile<<"xi_sp/xi_th:  "<<spCorrLgth_y*PI*tem/vF<<'\n';
	infosFile<<"xi_sp*T:  "<<spCorrLgth_y*tem<<'\n';
	
	
	m=qmax[1];
	l=m;
	chisptmp=chiqw_array[m+(l*(l+1))/2]/(1.0-Usp*chiqw_array[m+(l*(l+1))/2]/2.0);
	if (chisptmp<(chisp_max/2.0))
	{
		while (chisptmp<(chisp_max/2.0) && l<=qmax[0])
		{
			l++;
			chisptmp2=chisptmp;
			chisptmp=chiqw_array[m+(l*(l+1))/2]/(1.0-Usp*chiqw_array[m+(l*(l+1))/2]/2.0);
		}
		q0=(l-1)*PI/(nbq-1);
		qm=(qmax[0])*PI/(nbq-1);
		qHM=PI/(nbq-1)*(chisp_max/2.0 - chisptmp2)/(chisptmp-chisptmp2)+q0;
		
		spCorrLgth_x=1.0/(qm-qHM);
	}
	else
		spCorrLgth_x=spCorrLgth_y;
	
	cout<<"longueur de correlation de spin (1/(larg.mi-haut.)) (axe x):  "<<spCorrLgth_x<<'\n';
	cout<<"xi_sp/xi_th:  "<<spCorrLgth_x*PI*tem/vF<<'\n';
	cout<<"xi_sp*T:  "<<spCorrLgth_x*tem<<'\n';
	infosFile<<"longueur de correlation de spin (1/(larg.mi-haut.)) (axe x):  "<<spCorrLgth_x<<'\n';
	infosFile<<"xi_sp/xi_th:  "<<spCorrLgth_x*PI*tem/vF<<'\n';
	infosFile<<"xi_sp*T:  "<<spCorrLgth_x*tem<<'\n';
		
	double dchi2;
	if ( qmax[0]==(nbq-1) && qmax[1]==(nbq-1) )
	{
		dchi2=2.0*(chi0tmp[nbq-2]-chi0max)*(nbq-1)*(nbq-1)/(PI*PI);
		xi0_x=sqrt(-0.5*dchi2/chi0max);
		spCorrLgthAs_x=xi0_x*sqrt(Usp/(2.0/chi0max-Usp));
//		double SigAs=-3.0*U*tem*spCorrLgthAs/(16*xi0_x*xi0_x);

		xi0_y=xi0_x;
		spCorrLgthAs_y=spCorrLgthAs_x;
	}
	else
	{
		
//		infosFile<<"chisp incommensurable \n";
//		infosFile.close();
		
//		l=qmax[0];
//		m=qmax[1];
//		if (m>0)
//		{
//			if (m<l)
//				dchi2=(chi0tmp[m+1]-2*chi0tmp[m]+chi0tmp[m-1])*(nbq-1)*(nbq-1)/(PI*PI);
//			else
//				dchi2=(chi0tmp[m+1]-2*chi0tmp[m]+chi0tmp[m-1])*(nbq-1)*(nbq-1)/(PI*PI);
//		}
//		else
//			dchi2=2*(chi0tmp[m+1]-chi0tmp[m])*(nbq-1)*(nbq-1)/(PI*PI);
		
//		dchi2=(chi0tmp[m+1]-2*chi0tmp[m]+chi0tmp[m-1])*(nbq-1)*(nbq-1)/(PI*PI);
//		xi0_y=sqrt(-0.5*dchi2/chi0max);
//		double spCorrLgthAstmp=xi0_y*sqrt(Usp/(2.0/chi0max-Usp));		
//		cout<<"spCorrLgthAstmp:  "<<spCorrLgthAstmp<<endl;
		
		int x, y, indtmp;
		double chitmpl, chitmpm, chitmpr;
		
		x=qmax[0];
		y=qmax[1]-1;
		if (y<0)	y=qmax[1]+1;
		chitmpl=chiqw_array[y+(x*(x+1))/2];
		y=qmax[1];
		chitmpm=chiqw_array[y+(x*(x+1))/2];
		y=qmax[1]+1;
		if (y>x)
		{
			x=qmax[1]+1;
			y=qmax[0];
		}
		chitmpr=chiqw_array[y+(x*(x+1))/2];

		dchi2=(chitmpr-2*chitmpm+chitmpl)*(nbq-1)*(nbq-1)/(PI*PI);
		xi0_y=sqrt(-0.5*dchi2/chi0max);
		spCorrLgthAs_y=xi0_y*sqrt(Usp/(2.0/chi0max-Usp));
		
		
		x=qmax[0]-1;
		y=qmax[1];
		if (y>x)
		{
			x=qmax[1];
			y=qmax[0]-1;
			if (y<0) y=qmax[0]+1;
		}
		chitmpl=chiqw_array[y+(x*(x+1))/2];
		x=qmax[0];
		y=qmax[1];
		chitmpm=chiqw_array[y+(x*(x+1))/2];
		x=qmax[0]+1;
		if (x>(nbq-1)) 	x=qmax[0]-1;
		chitmpr=chiqw_array[y+(x*(x+1))/2];
		
		dchi2=(chitmpr-2*chitmpm+chitmpl)*(nbq-1)*(nbq-1)/(PI*PI);
		xi0_x=sqrt(-0.5*dchi2/chi0max);
		spCorrLgthAs_x=xi0_x*sqrt(Usp/(2.0/chi0max-Usp));
		
	}
	
	cout<<"longueur de correlation de spin (formule asymptotique) (axe y):  "<<spCorrLgthAs_y<<'\n';
	cout<<"xi_sp_As/xi_th:  "<<spCorrLgthAs_y*PI*tem/vF<<'\n';
	cout<<"xi_sp_As*T:  "<<spCorrLgthAs_y*tem<<'\n';
	
	infosFile<<"longueur de correlation de spin (formule asymptotique) (axe y):  "<<spCorrLgthAs_y<<'\n';
	infosFile<<"xi_sp_As/xi_th:  "<<spCorrLgthAs_y*PI*tem/vF<<'\n';
	infosFile<<"xi_sp_As*T:  "<<spCorrLgthAs_y*tem<<'\n';
	
	cout<<"longueur de correlation de spin (formule asymptotique) (axe x):  "<<spCorrLgthAs_x<<'\n';
	cout<<"xi_sp_As/xi_th:  "<<spCorrLgthAs_x*PI*tem/vF<<'\n';
	cout<<"xi_sp_As*T:  "<<spCorrLgthAs_x*tem<<'\n';
	
	infosFile<<"longueur de correlation de spin (formule asymptotique) (axe x):  "<<spCorrLgthAs_x<<'\n';
	infosFile<<"xi_sp_As/xi_th:  "<<spCorrLgthAs_x*PI*tem/vF<<'\n';
	infosFile<<"xi_sp_As*T:  "<<spCorrLgthAs_x*tem<<'\n';
	
	cout<<"(max(chi_sp)/max(chi_0))^0.5:  "<<sqrt(chisp_max/chi0max)<<endl;
	
/*
	int j, w, Nq=4*(nbq-1)*(nbq-1);
	double *chisp_q=new double[nk8th];
	double wn, chi0val;
	for (j=0; j<nk8th; j++) chisp_q[j]=0;
	
	for (l=0; l<nbq; l++)
		for (m=0; m<=l; m++)
		{
			chi0val=chiqw_array[m + (l*(l+1))/2];
			chisp_q[m + (l*(l+1))/2]+=tem*( chi0val/(1-Usp*chi0val/2.0) );
			
			for (j=nwmax; j>0; j--)
			{	
				chi0val=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
				chisp_q[m + (l*(l+1))/2]+=2*tem*( chi0val/(1-Usp*chi0val/2.0) );
			}
			
			for (j=2*nwmax; j>nwmax; j--)
			{	
				wn=2*j*PI*tem;
				chi0val=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
				chisp_q[m + (l*(l+1))/2]+=2*tem*( chi0val/(1-Usp*chi0val/2.0) );
			}
		}
	
	const char *nameFormStruct="./spin_Struct_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d.dat";
	
	fstream struct_file;
	
	sprintf(name, nameFormStruct, (double)tp, (double)tpp, (double)U, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);
	struct_file.open(name,ios::out);
	
	struct_file<<setiosflags(ios::left)<<setprecision(12);
	
	for (l=0; l<nbq; l++)
		for (m=0; m<=l; m++)
		{
			struct_file<<setw(10)<<l<<setw(10)<<m<<chisp_q[m + (l*(l+1))/2]<<endl;
		}
	
	struct_file.close();
	
	l=qmax[0], m=qmax[1];
	double chisp_q_max=chisp_q[m + (l*(l+1))/2];
	
	m=0;
	chisptmp=chisp_q[m + (l*(l+1))/2];
	while (chisptmp<(chisp_q_max/2.0) && m<=qmax[1])
	{
		m++;
		chisptmp2=chisptmp;
		chisptmp=chisp_q[m + (l*(l+1))/2];
	}
	q0=(m-1)*PI/(nbq-1);
	qm=(qmax[1])*PI/(nbq-1);
	qHM=PI/(nbq-1)*(chisp_q_max/2.0 - chisptmp2)/(chisptmp-chisptmp2)+q0;
	spCorrLgthS_y=1.0/(qm-qHM);
	
	cout<<"longueur de correlation de spin d'apres S_sp(q) (axe y):  "<<spCorrLgthS_y<<'\n';
	cout<<"xi_sp/xi_th:  "<<spCorrLgthS_y*PI*tem/vF<<'\n';
	infosFile<<"longueur de correlation de spin d'apres S_sp(q) (axe y):  "<<spCorrLgthS_y<<'\n';
	infosFile<<"xi_sp/xi_th:  "<<spCorrLgthS_y*PI*tem/vF<<'\n';
	
	m=qmax[1];
	l=m;
	chisptmp=chisp_q[m + (l*(l+1))/2];
	if (chisptmp<(chisp_q_max/2.0))
	{
		while (chisptmp<(chisp_q_max/2.0) && l<=qmax[0])
		{
			l++;
			chisptmp2=chisptmp;
			chisptmp=chisp_q[m + (l*(l+1))/2];
		}
		q0=(l-1)*PI/(nbq-1);
		qm=(qmax[0])*PI/(nbq-1);
		qHM=PI/(nbq-1)*(chisp_q_max/2.0 - chisptmp2)/(chisptmp-chisptmp2)+q0;
		spCorrLgthS_x=1.0/(qm-qHM);
	}
	else
		spCorrLgthS_x=spCorrLgthS_y;
	
	cout<<"longueur de correlation de spin d'apres S_sp(q) (axe x):  "<<spCorrLgthS_x<<'\n';
	cout<<"xi_sp/xi_th:  "<<spCorrLgthS_x*PI*tem/vF<<'\n';
	infosFile<<"longueur de correlation de spin d'apres S_sp(q) (axe x):  "<<spCorrLgthS_x<<'\n';
	infosFile<<"xi_sp/xi_th:  "<<spCorrLgthS_x*PI*tem/vF<<'\n';
 
	delete [] chisp_q;
*/	
	
	infosFile.close();
	
}

//return a value of chisp(q,iw), params=(qy,w)
double chi::chisp(double qx, double params[])
{
	double chi0tmp=chiqw(qx,params);

	return chi0tmp/(1-Usp*chi0tmp/2.0);
}

//return a value of chich(q,iw),  params=(qy,w)
double chi::chich(double qx, double params[])
{
	double chi0tmp=chiqw(qx,params);

	return chi0tmp/(1+Uch*chi0tmp/2.0);
}


//! trace of chisp(q,iw)-chi0(q,iw) as a function of Usp=Us
double chi::traceDiffChispChi0(double Us)
{
 double chi0tmp, trChitmp=0, trChitmp1=0;

 long int j, l, m;
 int w;
 int Nq=4*(nbq-1)*(nbq-1);
 int nk8th=(nbq*(nbq+1))/2;
 
 double wn;
 
 
 #pragma omp parallel for private(j,l,m,w,chi0tmp,wn) reduction(+:trChitmp1)
 for (j=2*nwmax; j>nwmax; j--)
	 for (l=0; l<nbq; l++)
		 for (m=0; m<=l; m++)
		 {
			 w=16;
			 if (j==0) w=w/2;
			 if (l==0 || l==nbq-1) w=w/2;
			 if (m==0 || m==nbq-1) w=w/2;
			 if (m==l) w=w/2;

			 wn=2*j*PI*tem;
			 chi0tmp=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);
			 trChitmp1+=w*tem/Nq*( chi0tmp/(1-Us*chi0tmp/2.0) - chi0tmp);
		 }
		 
 
#pragma omp parallel for private(j,l,m,w,chi0tmp) reduction(+:trChitmp)
 for (j=nwmax; j>=0; j--)
	 for (l=0; l<nbq; l++)
		 for (m=0; m<=l; m++)
		 {
			 w=16;
			 if (j==0) w=w/2;
			 if (l==0 || l==nbq-1) w=w/2;
			 if (m==0 || m==nbq-1) w=w/2;
			 if (m==l) w=w/2;

			 chi0tmp=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
			 trChitmp+=w*tem/Nq*( chi0tmp/(1-Us*chi0tmp/2.0) - chi0tmp);
		 }
 
 return trChitmp+trChitmp1;
}


//! trace of chich(q,iw)-chi0(q,iw) as a function of Uch=Uc
double chi::traceDiffChichChi0(double Uc)
{
 double chi0tmp, trChitmp=0, trChitmp1=0;

 long int j, l, m;
 int w;
 int Nq=4*(nbq-1)*(nbq-1);
 int nk8th=(nbq*(nbq+1))/2;

 
 double wn;
 #pragma omp parallel for private(j,l,m,w,chi0tmp,wn) reduction(+:trChitmp1)
 for (j=2*nwmax; j>nwmax; j--)
	 for (l=0; l<nbq; l++)
		 for (m=0; m<=l; m++)
		 {
			 w=16;
			 if (j==0) w=w/2;
			 if (l==0 || l==nbq-1) w=w/2;
			 if (m==0 || m==nbq-1) w=w/2;
			 if (m==l) w=w/2;

			 wn=2*j*PI*tem;
			 chi0tmp=chi0_inf[m + (l*(l+1))/2]/(wn*wn)+chi0_inf2[m + (l*(l+1))/2]/(wn*wn*wn*wn);			 
			 trChitmp1+=w*tem/Nq*( chi0tmp/(1+Uc*chi0tmp/2.0) - chi0tmp);
		 }
 
#pragma omp parallel for private(j,l,m,w,chi0tmp) reduction(+:trChitmp)
 for (j=nwmax; j>=0; j--)
	 for (l=0; l<nbq; l++)
		 for (m=0; m<=l; m++)
		 {
			 w=16;
			 if (j==0) w=w/2;
			 if (l==0 || l==nbq-1) w=w/2;
			 if (m==0 || m==nbq-1) w=w/2;
			 if (m==l) w=w/2;

			 chi0tmp=chiqw_array[m + (l*(l+1))/2 + j*nk8th];
			 trChitmp+=w*tem/Nq*( chi0tmp/(1+Uc*chi0tmp/2.0) - chi0tmp);
		 }
 
 return trChitmp+trChitmp1;
}

//! function of which the root is Usp
double chi::sum_rule_sp(double Us, double par[])
{
	return traceDiffChispChi0(Us) + density-density*density/2.0 - nsr + Us*nsr*nsr/(2.0*U);
}


//! function of which the root is Uch
double chi::sum_rule_ch(double Uc, double par[])
{
	return traceDiffChichChi0(Uc) + density-density*density/2.0 - (nsr + Usp*nsr*nsr/(2.0*U) - nsr*nsr);
}

//! find Usp and Uch for chi0 defined on a given grid
void chi::find_Usp_Uch(int nq0, int nw0, double Usc[])
{	
	calc_chiqw_FFT_spline(nq0,nw0);

	fctPtr Ptr=static_cast<fctPtr> (&chi::sum_rule_sp);

	double root[]={0.0};
	double init[]={U/2.0, 0.0, 2.0/chi0max-2*EPSILON};

	if (find_zero(Ptr, init, NULL, root))
	{
		Usc[0]=root[0];
		Usp=Usc[0];
		cout<<"Usp: "<<Usc[0]<<'\n';
	}
	else
		cout<<"Usp n'a pas ete trouve\n";

	Ptr=static_cast<fctPtr> (&chi::sum_rule_ch);

	root[0]=0.0;
	init[0]=0;
	init[2]=1000.0*U;

	if (find_zero(Ptr, init, NULL, root))
	{
		Usc[1]=root[0];
		cout<<"Uch: "<<Usc[1]<<'\n';
	}
	else
		cout<<"Uch n'a pas ete trouve\n";

/*	
	fstream infosFile;
	
	const char *nameFormInfos="./TPSC_infos_%6.4f_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d.dat";
	char name[100];
	
	sprintf(name, nameFormInfos, (double)tp, (double)tpp, (double)U, (double)density, (double)tem, 2*(nbq-1), 2*nwmax);
	
	infosFile.open(name,ios::out);
	infosFile<<"Usp:  "<<Usp<<"   difRelsp:  "<<difRelsp<<'\n';
	infosFile<<"Uch:  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
	infosFile<<"chi0max: "<<chi0max<<"   chispmax: "<<chispmax<<"   rapport: "<<chispmax/chi0max<<'\n';
	infosFile<<"double occupation (ansatz): "<<dblOcc<<'\n';
	infosFile<<"qxmax, qymax:   "<<qmax[0]<<",  "<<qmax[1]<<'\n';
	
	cout<<"chi0max: "<<chi0max<<"  chispmax: "<<chispmax<<"  rapport: "<<chispmax/chi0max<<'\n';
	cout<<"double occupation: "<<dblOcc<<'\n';
	
	infosFile.close();
	
	paramsFile.open(paramsFileName,ios::out);
	if (!paramsFile)
	{
		cout<<"echec de creation du fichier de parametres\n";
		return;
	}
	
	paramsFile<<setprecision(15);
	paramsFile<<"Usp  "<<Usp<<"   difRelsp:  "<<difRelsp<<'\n';
	paramsFile<<"Uch  "<<Uch<<"   difRelch:  "<<difRelch<<'\n';
	paramsFile<<"nbq  "<<nbq<<'\n';
	paramsFile<<"nbw  "<<nwmax+1<<'\n';
	paramsFile<<"dblOcc_ansatz  "<<dblOcc<<'\n';
	paramsFile<<"qxmax  "<<qmax[0]<<'\n';
	paramsFile<<"qymax  "<<qmax[1]<<'\n';
	paramsFile<<"chispmax_over_chi0max  "<<chispmax/chi0max<<'\n';
	paramsFile<<"mu0  "<<mu0<<'\n';
	
	cout<<"nbq: "<<nbq<<"   nw: "<<nwmax+1<<'\n';
	
	paramsFile.close();	
*/
}

//print the parameters on screen
void chi::get_params() const
{
//	cout<<setiosflags(ios::left)<<"U: "<<setw(10)<<U;
	green0::get_params();
}


//! calculate the correlation length for a set of densities and temperature, nD is the number of densities in the vector dens and nT, the number of temperatures in the vector Temp
void chi::calc_corr_length(double dens[], int nD, double Temp[], int nT)
{
	
	fstream file;
	
	const char *nameForm="./corrLgth_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%7.5f_%6.4f_%6.4f.dat";
	char name[100];
	
	int j, l;
	
	double T, T0=Temp[0], Tf=Temp[1];
	double dT=(Tf-T0)/(nT-1);
	
	double n, n0=dens[0], nf=dens[1];
	double dn=(nf-n0)/(nD-1);
	
	sprintf(name, nameForm, (double)U, (double)tp, (double)tpp, (double)t3, (double)n0, (double)nf, (double)T0, (double)Tf);
	file.open(name,ios::out);
	file<<setiosflags(ios::left)<<setprecision(15);
	
	for (l=0; l<nD; l++)
	{
		n=n0+l*dn;
		set_density(n);
		
		for (j=0; j<nT; j++)
		{
			T=T0+j*dT;
			set_tem(T);
			
			find_Usp_Uch();
			
			//			vF_hot_spots();
			
			sp_corr_length();
			
			
			file<<setw(25)<<n<<setw(25)<<T<<setw(25)<<mu0;
			file<<setw(25)<<vF<<setw(25)<<Usp<<setw(25)<<errUs;
			file<<setw(25)<<Uch<<setw(25)<<errUc<<setw(25)<<dblOcc;
			file<<setw(25)<<spCorrLgth_x<<setw(25)<<spCorrLgth_y;
			file<<setw(25)<<spCorrLgthAs_x<<setw(25)<<spCorrLgthAs_y;
			file<<setw(25)<<xi0_x<<setw(25)<<xi0_y<<setw(25)<<nbq;
			file<<setw(25)<<qmax[0]<<setw(25)<<qmax[1];
			file<<chi0max<<'\n';
			
		}
	}
	
	file.close();
}

void chi::chi_Re_w(int nk, int Nwn, double dens, double T, double *w, int Nwr, int N0, double *eta, int Neta, int m, char *rep_in, char *rep_out, char *nom, int type_fic)
{
	fstream file, file_info;
	int j, n;
	char name[300];
	
	int nk0=2*(nbq-1);
	int nbw=2*nwmax;
	
	if (nk0!=nk || nwmax!=Nwn || density!=dens || tem!=T)
	{
		cout<<"chi_Re_w(): incoherence de parametres\n";
		return;
/*		
		nk0=nk;
		nbq=nk0/2+1;
		nbw=Nwn;
		nwmax=Nwn/2;
		density=dens;
		tem=T;
*/ 
	}
	
	char nameForm_TPSCparams[300];
	char nameForm_chiqn[300];
	char nameForm_chiwRe[300];
	char nom_param[10];
	
	double flparam;
	int qxmax, qymax, intparam;
	
	strcpy(nameForm_TPSCparams,rep_in);
	strcat(nameForm_TPSCparams,"/TPSCparams_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d.dat");
	sprintf(name,nameForm_TPSCparams,(double)U,(double)tp,(double)tpp,(double)t3,(double)density,(double)tem,nk0,nbw);
	file.open(name, ios::in);
	if (!file)
	{
		cout<<"chi_Re_w(): fichier "<<name<<" inexistant:\n";
		return;
	}
	file>>nom_param>>Usp>>nom_param>>flparam;
	file>>nom_param>>Uch>>nom_param>>flparam;
	file>>nom_param>>intparam>>nom_param>>intparam;
	file>>nom_param>>flparam;
	file>>nom_param>>qxmax;
	file>>nom_param>>qymax;
	file.close();
	
	cout<<"Usp:  "<<Usp<<endl;
	cout<<"qxmax:  "<<qxmax<<"    qymax:  "<<qymax<<endl;
	
	
	strcpy(nameForm_chiqn,rep_in);
	strcat(nameForm_chiqn, "/chi0_qmax_wn_%6.4f_%6.4f_%6.4f_%6.4f_%1d_%1d_%1d_%1d.dat");
	sprintf(name, nameForm_chiqn, (double)tp, (double)tpp, (double)density, (double)tem, nk0, nbw, qxmax, qymax);
	
	file.open(name, ios::in );
	if (!file)
	{
		cout<<"chi_Re_w(): fichier "<<name<<" inexistant:\n";
		return;
	}
	
	double *chiqn=new double[nbw/2+1];
	
	for (j=0; j<=nbw/2; j++)
	{
		file>>n>>chiqn[j];
	}
	file.close();
	
	
	dcomplex *coef=new dcomplex[N0];
	dcomplex *z0=new dcomplex[N0];
	dcomplex *func=new dcomplex[N0];
	int *n_pade=new int[N0];
	
	int Np=N0;
	for (j=0; j<Np; j++) n_pade[j]=j;
	
	for (j=0; j<Np; j++) z0[j]=dcomplex(0,2*n_pade[j]*PI*tem);
	
	for (j=0; j<Np; j++) func[j]=chiqn[n_pade[j]]/(1.0-Usp*chiqn[n_pade[j]]/2.0);
	
	pade_cont_frac_coef_rec(func, z0, Np, coef);
	
	strcpy(nameForm_chiwRe,rep_out);
	strcat(nameForm_chiwRe, "/chi0_qmax_wRe_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%d_%d_%d_%d.dat");
	sprintf(name, nameForm_chiwRe, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, qxmax, qymax);
	
	file.open(name, ios::out );
	if (!file)
	{
		cout<<"chi_Re_w(): fichier "<<name<<" inexistant:\n";
		return;
	}
	
	cout<<"chi_Re_w():  ecriture dans le fichier\n";
	dcomplex chitmp;
	file<<setiosflags(ios::left)<<setprecision(12);
	for (j=0; j<Nwr; j++)
	{
		file<<setw(25)<<w[j];
		for (n=0; n<Neta; n++)
		{
			chitmp=pade(dcomplex(w[j],eta[n]), Np, z0, coef);
			file<<setw(25)<<chitmp.real()<<setw(35)<<chitmp.imag();
		}
		file<<endl;
	}
	file.close();
	
	delete [] n_pade;
	delete [] coef;
	delete [] z0;
	delete [] func;
	
}
