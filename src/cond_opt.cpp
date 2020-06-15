
/*! \file cond_opt.cpp
(C) Dominic Bergeron
 */

#include "includeDef.h"
#include "cond_opt.h"


cond_opt::cond_opt(int nk0, int Nw0, double n, double T, double eta, double  params[]): green(nk0, Nw0, n, T, eta, params)
{
	
	eta0=eta;
	
	get_params();
		
	sigma_array=NULL;
	V_array=NULL;
	ImChijj0=NULL;
	chijj1_qn=NULL;
	chijj1_qn_corr=NULL;
	chijj1_qn4=NULL;
	chijj2_qn=NULL;
	chijj2_qn_corr=NULL;
	chijj2_qn6=NULL;
	chijj3_qn=NULL;
	chijj3_qn_corr=NULL;
	chijj3_qn4=NULL;
	
	chijj_bulle_qn=NULL;
	
	chijj_w0=0.0;
}


cond_opt::~cond_opt()
{	
	if (sigma_array)
		delete [] sigma_array;
	
	if (ImChijj0) delete [] ImChijj0;
	
	if (chijj1_qn) delete [] chijj1_qn;
	
	if (chijj1_qn_corr) delete [] chijj1_qn_corr;
	
	if (chijj1_qn4) delete [] chijj1_qn4;
	
	if (chijj2_qn) delete [] chijj2_qn;
	
	if (chijj2_qn_corr) delete [] chijj2_qn_corr;
	
	if (chijj2_qn6) delete [] chijj2_qn6;
	
	if (chijj3_qn) delete [] chijj3_qn;
	
	if (chijj3_qn_corr) delete [] chijj3_qn_corr;
	
	if (chijj3_qn4) delete [] chijj3_qn4;
	
	if (chijj_bulle_qn) delete [] chijj_bulle_qn;
}

void cond_opt::calc_chijj_bulle_1()
{
	cout<<"calcul de chijj1 par calc_chijj_bulle()\n";
	
	int j, l, m;
	
	int nw=nbw/2;
	int nk0=2*(nbk-1);
	int nk4th=nbk*(nbk-2);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (!Self_kw_array)
	{
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}
	
	if (!mu_calcule)
	{
		cout<<"calc_chijj_bulle(): mu non calcule\n";
		return;
	}
	
	find_FS_inter();
	
	double normFact=1.0/(nk0*nk0);
	
	double *dek=new double[nbk*(nbk-2)];
	double *d2ek=new double[nbk*nbk];
	int indHss[]={0,0};
	double vk[2];
	
	double kx, ky;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			d2ek[m+l*nbk]=Hessian_ek(vk, indHss);
		}
	}
	
	
	for (l=0; l<nbk-2; l++)
	{
		kx=(l+1)*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			dek[m+l*nbk]=grad_ek(vk, 0);
		}
	}
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex *Gtmp=new dcomplex[nbw];
	fftw_plan fftplan=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (Gtmp), reinterpret_cast<fftw_complex*> (Gtmp), FFTW_FORWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	//	int k0[]={0,0};
	//	int wn0=0;
	
	//	int x1, y1, tmp1;
	//	double ektmp, kn, sgn;
	//	 dcomplex Gtmp1, HGsum=0;
	//	int l0=k0[0], m0=k0[1], n0=wn0;
	
	
	Gtmp=new dcomplex[nbw];
	fftw_plan fftplan_t_back=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*>(Gtmp), reinterpret_cast<fftw_complex*>(Gtmp), FFTW_BACKWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	if (chijj1_qn) delete [] chijj1_qn;
	
	chijj1_qn=new double[nw+1];
	chijj1_qn4=new double[nw];
	
	for (j=0; j<nw; j++)
	{
		chijj1_qn[j]=0;
		chijj1_qn4[j]=0;
	}
	chijj1_qn[nw]=0;
	
	chijj1_w0=0;
	chijj1_inf2=0;
	
	int NS0=nbw/2+1;
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex k0sum=0;
	
	//	l0=k0[0], m0=k0[1];
	
#pragma omp parallel private(m,j)
	{
		int x, y, tmp, p;
		double wt, t, kn, G0t, Jinf;
		dcomplex k0, ph;
		
		//		dcomplex *J0tmp2=new dcomplex[nbw];
		
		//		dcomplex Jinf_tmp, Jbeta;
		//		dcomplex *J0tmp=new dcomplex[nbw];
		//		dcomplex *Jtmp=new dcomplex[nbw];
		//		dcomplex *HGtmp=new dcomplex[nbw];
		//		double CJ1, CJ2, CJ3, d2J0, d2Jbeta, CJ01, CJ02, CJ03;
		
		dcomplex *Gtmp1=new dcomplex[nbw];
		dcomplex *GGtmp1=new dcomplex[nbw];
		dcomplex Gtmp2, Gtmp0, Gbeta;
		
		double *coeffs=new double[4*(NS0-1)];
		double *tau=new double[NS0];
		double *ft_tmp=new double[NS0];
		double wn, integ, a, b, c, d, x1, x2;
		
		double *tau2=new double[nbw+1];
		double *ft_tmp2=new double[nbw+1];
		double *coeffs2=new double[4*nbw];
		
		dcomplex k0sum_loc=0, chijj2_qn0_loc;
		double *chijj1_qn_loc=new double[nw+1];
		double *chijj1_qn4_loc=new double[nw];
		//		double *chijj2_qn_loc=new double[nw+1];
		//		double *chijj2_qn6_loc=new double[nw];
		for (j=0; j<nw; j++)
		{
			chijj1_qn_loc[j]=0;
			chijj1_qn4_loc[j]=0;
			//			chijj2_qn_loc[j]=0;
			//			chijj2_qn6_loc[j]=0;
		}
		chijj1_qn_loc[nw]=0;
		//		chijj2_qn_loc[nw]=0;
		
		//		dcomplex chijj3_qn0_loc=0;
		//		dcomplex chijj3_qn0_tmp;
		
		double chijj2_inf2_loc=0, chijj2_inf4_loc=0, chijj1_inf2_loc=0, chijj1_w0_loc=0;
		
		//		double dtau_J0, dtau_Jbeta;
		double dtau_Gkp, dtau_Gkm, dtau_chi;
		double C2, C3;
		dcomplex selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex exp1, exp2, exp3, exp4;
		dcomplex f1p, f2p, f3p, f4p, f1m, f2m, f3m, f4m;
		
		dcomplex exptmp;
		dcomplex b1, b2, b3, b4;
		dcomplex denj1, denj2, denj3, denj4;
		dcomplex expj1, expj2, expj3, expj4;
		dcomplex fj1, fj2, fj3, fj4, fj1b, fj2b, fj3b, fj4b;
		
		//		double theta_k;
		//		int ind_theta;
		//		double *chijj1_k_loc=new double[(nw+1)*N_theta_k];
		//		double *chijj2_k_loc=new double[(nw+1)*N_theta_k];
		
		//		for (j=0; j<(nw+1)*N_theta_k; j++)
		//		{
		//			chijj1_k_loc[j]=0;
		//			chijj2_k_loc[j]=0;
		//		}
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
		
#pragma omp for
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				k0=0;
				dtau_Gkp=0;
				//				chijj2_qn0_loc=0;
				//dtau_J0=0;
				//				chijj3_qn0_tmp=0;
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp2=tem/dcomplex(-ek[m+(l+1)*nbk]+mu, kn+Sigma_inf/kn);
//					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);

					Gtmp2=tem/dcomplex(-ek[m+(l+1)*nbk]+mu, kn+Sigma_inf/kn);
//					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				
				// calcul du G(k,t) asymptotique et du J(k,t) asymtotique
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				find_roots_2pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				
				a1=r1/(r1-r2);
				a2=r2/(r2-r1);
				
				if (real(r1)>0)
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
				}
				else
				{
					exptmp=exp(r1/tem);
					den1=exptmp+1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
				}
				else
				{
					exptmp=exp(r2/tem);
					den2=exptmp+1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
				}
				
				k0+=real(f1m+f2m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p);
				dtau_Gkm+=real(r1*f1m+r2*f2m);
				
/*
 				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0)
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
				}
				else
				{
					exptmp=exp(r1/tem);
					den1=exptmp+1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
				}
				else
				{
					exptmp=exp(r2/tem);
					den2=exptmp+1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
				}
				if (real(r3)>0)
				{
					exptmp=exp(-r3/tem);
					den3=exptmp+1.0;
					f3p=a3/den3;
					f3m=a3*exptmp/den3;
				}
				else
				{
					exptmp=exp(r3/tem);
					den3=exptmp+1.0;
					f3p=a3*exptmp/den3;
					f3m=a3/den3;
				}
				if (real(r4)>0)
				{
					exptmp=exp(-r4/tem);
					den4=exptmp+1.0;
					f4p=a4/den4;
					f4m=a4*exptmp/den4;
				}
				else
				{
					exptmp=exp(r4/tem);
					den4=exptmp+1.0;
					f4p=a4*exptmp/den4;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p+r3*f3p+r4*f4p);
				dtau_Gkm+=real(r1*f1m+r2*f2m+r3*f3m+r4*f4m);
				*/
				
				//TF de G(k,ikn) dans le temps
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				Gbeta=Gtmp1[0];
				
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					
					if (real(r1)>0)
					{
						exp1=exp(-t*r1);
					}
					else
					{
						exp1=exp((1.0/tem-t)*r1);
					}
					if (real(r2)>0)
					{
						exp2=exp(-t*r2);
					}
					else
					{
						exp2=exp((1.0/tem-t)*r2);
					}
/*
					if (real(r3)>0)
					{
						exp3=exp(-t*r3);
					}
					else
					{
						exp3=exp((1.0/tem-t)*r3);
					}
					if (real(r4)>0)
					{
						exp4=exp(-t*r4);
					}
					else
					{
						exp4=exp((1.0/tem-t)*r4);
					}
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
*/
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2);
					
					ph=exp(dcomplex(0.0, (j*PI*(nbw-1))/nbw));
					
					Gtmp1[j]=ph*Gtmp1[j]+G0t;
					
				}
				
				t=1.0/tem;
				
				if (real(r1)>0)
				{
					exp1=exp(-t*r1);
				}
				else
				{
					exp1=1.0;
				}
				if (real(r2)>0)
				{
					exp2=exp(-t*r2);
				}
				else
				{
					exp2=1.0;
				}
/*
				if (real(r3)>0)
				{
					exp3=exp(-t*r3);
				}
				else
				{
					exp3=1.0;
				}
				if (real(r4)>0)
				{
					exp4=exp(-t*r4);
				}
				else
				{
					exp4=1.0;
				}
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
*/
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2);
				
				ph=exp(dcomplex(0.0, PI*(nbw-1)));
				
				Gbeta=ph*Gbeta+G0t;
				
/*
				 if ((l+1)%4==0)
				 {
				 cout<<"erreur rel dtau_Gkp:  "<<setw(20)<<(real(Gtmp1[1]-Gtmp1[0])*tem*nbw-dtau_Gkp)/dtau_Gkp<<'\n';
				 cout<<"erreur rel dtau_Gkm:  "<<setw(20)<<(real(-Gtmp1[nbw-1]-Gtmp1[0]-1.0)*tem*nbw-dtau_Gkm)/dtau_Gkm<<'\n';
//				 cout<<"dtau_J0:  "<<setw(20)<<real(Jtmp[1]-Jtmp[0])*tem*nbw<<setw(20)<<dtau_J0<<(real(Jtmp[1]-Jtmp[0])*tem*nbw-dtau_J0)/dtau_J0<<endl;
//				 cout<<"dtau_Jbeta:  "<<setw(20)<<real(Jbeta-Jtmp[nbw-1])*tem*nbw<<setw(20)<<dtau_Jbeta<<(real(Jbeta-Jtmp[nbw-1])*tem*nbw-dtau_Jbeta)/dtau_Jbeta<<endl;
				 }
*/				
				
				// Definition de -G(k,t)G(k,-t)
				GGtmp1[0]=Gtmp1[0]*(Gtmp1[0]+1.0);
				for (j=1; j<nbw; j++)
					GGtmp1[j]=-Gtmp1[j]*Gtmp1[nbw-j];
				
				// calcul du spline de -G(k,t)G(k,-t)
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nbw*tem);
					ft_tmp[j]=real(GGtmp1[j]);
				}
				dtau_chi=dtau_Gkp*(real(Gtmp1[0])+1.0)+real(Gtmp1[0])*dtau_Gkm;
				
				//				if ((l+1)%4==0)
				//					cout<<"dtau_chi:  "<<setw(20)<<(real(GGtmp1[1]-GGtmp1[0])*(tem*nbw)-dtau_chitmp)/dtau_chitmp<<setw(20)<<(dtau_chi-dtau_chitmp)/dtau_chitmp<<'\n';
				
				coeffs[0]=dtau_chi;
				coeffs[1]=0;
				spline_coeffs(tau, ft_tmp, NS0, coeffs);
				
				// TF sur t de -G(k,t)G(k,-t)
				for (j=0; j<nbw/2; j++)
					Gtmp1[j]=6.0*coeffs[4*j];
				for (j=nbw/2; j<nbw; j++)
					Gtmp1[j]=-Gtmp1[nbw-j-1];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				integ=0;
				for (j=0; j<nbw/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j];
					x2=tau[j+1];
					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
				
				Gtmp1[0]=2*integ;
				
				wt=4;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+(l+1)*nbk]*k0;
				chijj1_qn_loc[0]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_w0_loc-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_inf2_loc+=4.0*dtau_chi*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk];
				
				for (j=1; j<=nw; j++)
				{
					wn=2*j*PI*tem;
					
					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real((1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]);
					
					Gtmp1[j]=-2.0*dtau_chi/(wn*wn)+real((1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]/(wn*wn*wn*wn));
					
					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
				}
				
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
			
			for (j=0; j<nw; j++)
			{
				chijj1_qn[j]+=chijj1_qn_loc[j];
				chijj1_qn4[j]+=chijj1_qn4_loc[j];
			}
			chijj1_qn[nw]+=chijj1_qn_loc[nw];
			
			chijj1_w0+=chijj1_w0_loc;
			chijj1_inf2+=chijj1_inf2_loc;
		}
		delete [] chijj1_qn4_loc;
		delete [] chijj1_qn_loc;
		delete [] Gtmp1;
		delete [] GGtmp1;
		delete [] tau;
		delete [] ft_tmp;
		delete [] coeffs;
		delete [] tau2;
		delete [] ft_tmp2;
		delete [] coeffs2;
	}
	
#pragma omp parallel private(m,j)
	{
		int x, y, tmp;
		double wt, kn;
		dcomplex k0;
		dcomplex Gtmp2, Gtmp0;
		
		dcomplex k0sum_loc=0;
		
		double C2, C3;
		dcomplex z, selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex f1m, f2m, f3m, f4m;
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk; l+=nk0/2)
		{
			for (m=0; m<nbk; m++)
			{
				x=l;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				//				C2=G1->Sigma_inf2[y+(x*(x+1))/2];
				//				C3=G1->Sigma_inf3[y+(x*(x+1))/2];
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				k0=0;
				
				/*
				 for (j=0; j<nw; j++)
				 {
				 kn=(2.0*(j-nbw)+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th]));
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;
				 }
				 for (j=nbw-1; j>=nw; j--)
				 {
				 kn=(2.0*j+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+j*nk8th]);
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;
				 }
				 */
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0)
				{
					den1=exp(-r1/tem)+1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else
				{
					den1=exp(r1/tem)+1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					den2=exp(-r2/tem)+1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else
				{
					den2=exp(r2/tem)+1.0;
					f2m=a2/den2;
				}
				if (real(r3)>0)
				{
					den3=exp(-r3/tem)+1.0;
					f3m=a3*exp(-r3/tem)/den3;
				}
				else
				{
					den3=exp(r3/tem)+1.0;
					f3m=a3/den3;
				}
				if (real(r4)>0)
				{
					den4=exp(-r4/tem)+1.0;
					f4m=a4*exp(-r4/tem)/den4;
				}
				else
				{
					den4=exp(r4/tem)+1.0;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				wt=2;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+l*nbk]*k0;
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
		}
	}
	
	//	cout<<"chitmp3:  ";
	//	calc_chijj_vertex_corr2_wn(0, k0, 0, Htmp, NULL, NULL, NULL, NULL);
	
	delete [] d2ek;
	
	//	for (j=0; j<4; j++)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	//	for (j=nbw/8; j<=nbw/2; j+=nbw/8)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	
	fstream file;
	char name[200];
	
	kx0=real(k0sum);
	sprintf(name, "kx0_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat", U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out);
	file<<setprecision(15)<<setiosflags(ios::left)<<kx0<<'\n';
	file.close();
	
	sprintf(name, "kx0_bin_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat", U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&kx0, sizeof(kx0));
	file.close();
	
	for (j=0; j<=nw; j++)
	{
		if (chijj1_qn[j]<0.0)
		{
			cout<<"attention! valeur non-physique de chi_jj\n";
			cout<<setw(30)<<chijj1_qn[j]<<'\n';
			break;
		}
	}
	
	
	//pour enregistrer chijj(iq_n) en binaire
	
	
	const char *nameForm_chijj1_qn_bin="./chijj1_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn_bin, U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	
	
	j=0;
	file.write((char*)&j,sizeof(j));
	file.write((char*)&(chijj1_qn[j]),sizeof(double));
	for (j=1; j<=nbw/2; j++)
	{
		file.write((char*)&j,sizeof(j));
		file.write((char*)&(chijj1_qn[j]),sizeof(double));
	}
	file.close();
	// fin de l'enregistrement en binaire
	
	// pour enregistrer chijj(iq_n) en ascii
	
	const char *nameForm_chijj1_qn="./chijj1_qn_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn, U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left);
	
	j=0;
	file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	for (j=1; j<=nbw/2; j++)
	{
		file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	}
	file.close();
	// fin de l'enregistrement en ascii
	
	cout<<"-k0:  "<<-kx0<<'\n';
	cout<<"chijj1(iqn=0):  "<<chijj1_qn[0]<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+kx0)/kx0<<'\n';
	
	
	fftw_destroy_plan(fftplan);
	fftw_destroy_plan(fftplan_t_back);
	
	delete [] dek;
	
	if (chijj_bulle_qn) delete [] chijj_bulle_qn;
	chijj_bulle_qn=new double[nw+1];
	
	for (j=0; j<=nw; j++) chijj_bulle_qn[j]=chijj1_qn[j];
}

void cond_opt::calc_chijj_bulle()
{
	cout<<"calcul de chijj1 par calc_chijj_bulle()\n";
	
	int j, l, m;
	
	int nw=nbw/2;
	int nk0=2*(nbk-1);
	int nk4th=nbk*(nbk-2);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (!Self_kw_array)
	{
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}
	
	if (!mu_calcule)
	{
		cout<<"calc_chijj_bulle(): mu non calcule\n";
		return;
	}
	
	find_FS_inter();
	
	double normFact=1.0/(nk0*nk0);
	
	double *dek=new double[nbk*(nbk-2)];
	double *d2ek=new double[nbk*nbk];
	int indHss[]={0,0};
	double vk[2];
	
	double kx, ky;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			d2ek[m+l*nbk]=Hessian_ek(vk, indHss);
		}
	}
	
	
	for (l=0; l<nbk-2; l++)
	{
		kx=(l+1)*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			dek[m+l*nbk]=grad_ek(vk, 0);
		}
	}
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex *Gtmp=new dcomplex[nbw];
	fftw_plan fftplan=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (Gtmp), reinterpret_cast<fftw_complex*> (Gtmp), FFTW_FORWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	//	int k0[]={0,0};
	//	int wn0=0;
	
	//	int x1, y1, tmp1;
	//	double ektmp, kn, sgn;
	//	 dcomplex Gtmp1, HGsum=0;
	//	int l0=k0[0], m0=k0[1], n0=wn0;
	
	
	Gtmp=new dcomplex[nbw];
	fftw_plan fftplan_t_back=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*>(Gtmp), reinterpret_cast<fftw_complex*>(Gtmp), FFTW_BACKWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	if (chijj1_qn) delete [] chijj1_qn;
	
	chijj1_qn=new double[nw+1];
	chijj1_qn4=new double[nw];
	
	for (j=0; j<nw; j++)
	{
		chijj1_qn[j]=0;
		chijj1_qn4[j]=0;
	}
	chijj1_qn[nw]=0;
	
	chijj1_w0=0;
	chijj1_inf2=0;
	
	int NS0=nbw/2+1;
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex k0sum=0;
	
	//	l0=k0[0], m0=k0[1];
	
#pragma omp parallel private(m,j)
	{
		int x, y, tmp, p;
		double wt, t, kn, G0t, Jinf;
		dcomplex k0, ph;
		
		//		dcomplex *J0tmp2=new dcomplex[nbw];
		
		//		dcomplex Jinf_tmp, Jbeta;
		//		dcomplex *J0tmp=new dcomplex[nbw];
		//		dcomplex *Jtmp=new dcomplex[nbw];
		//		dcomplex *HGtmp=new dcomplex[nbw];
		//		double CJ1, CJ2, CJ3, d2J0, d2Jbeta, CJ01, CJ02, CJ03;
		
		dcomplex *Gtmp1=new dcomplex[nbw];
		dcomplex *GGtmp1=new dcomplex[nbw];
		dcomplex Gtmp2, Gtmp0, Gbeta;
		
		double *coeffs=new double[4*(NS0-1)];
		double *tau=new double[NS0];
		double *ft_tmp=new double[NS0];
		double wn, integ, a, b, c, d, x1, x2;
		
		double *tau2=new double[nbw+1];
		double *ft_tmp2=new double[nbw+1];
		double *coeffs2=new double[4*nbw];
		
		dcomplex k0sum_loc=0, chijj2_qn0_loc;
		double *chijj1_qn_loc=new double[nw+1];
		double *chijj1_qn4_loc=new double[nw];
		//		double *chijj2_qn_loc=new double[nw+1];
		//		double *chijj2_qn6_loc=new double[nw];
		for (j=0; j<nw; j++)
		{
			chijj1_qn_loc[j]=0;
			chijj1_qn4_loc[j]=0;
			//			chijj2_qn_loc[j]=0;
			//			chijj2_qn6_loc[j]=0;
		}
		chijj1_qn_loc[nw]=0;
		//		chijj2_qn_loc[nw]=0;
		
		//		dcomplex chijj3_qn0_loc=0;
		//		dcomplex chijj3_qn0_tmp;
		
		double chijj2_inf2_loc=0, chijj2_inf4_loc=0, chijj1_inf2_loc=0, chijj1_w0_loc=0;
		
		//		double dtau_J0, dtau_Jbeta;
		double dtau_Gkp, dtau_Gkm, dtau_chi;
		double C2, C3;
		dcomplex selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex exp1, exp2, exp3, exp4;
		dcomplex f1p, f2p, f3p, f4p, f1m, f2m, f3m, f4m;
		
		dcomplex exptmp;
		dcomplex b1, b2, b3, b4;
		dcomplex denj1, denj2, denj3, denj4;
		dcomplex expj1, expj2, expj3, expj4;
		dcomplex fj1, fj2, fj3, fj4, fj1b, fj2b, fj3b, fj4b;
		
		//		double theta_k;
		//		int ind_theta;
		//		double *chijj1_k_loc=new double[(nw+1)*N_theta_k];
		//		double *chijj2_k_loc=new double[(nw+1)*N_theta_k];
		
		//		for (j=0; j<(nw+1)*N_theta_k; j++)
		//		{
		//			chijj1_k_loc[j]=0;
		//			chijj2_k_loc[j]=0;
		//		}
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
		
#pragma omp for
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				k0=0;
				dtau_Gkp=0;
				//				chijj2_qn0_loc=0;
				//				dtau_J0=0;
				//				chijj3_qn0_tmp=0;
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				
				// calcul du G(k,t) asymptotique et du J(k,t) asymtotique
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0)
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
				}
				else
				{
					exptmp=exp(r1/tem);
					den1=exptmp+1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
				}
				else
				{
					exptmp=exp(r2/tem);
					den2=exptmp+1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
				}
				if (real(r3)>0)
				{
					exptmp=exp(-r3/tem);
					den3=exptmp+1.0;
					f3p=a3/den3;
					f3m=a3*exptmp/den3;
				}
				else
				{
					exptmp=exp(r3/tem);
					den3=exptmp+1.0;
					f3p=a3*exptmp/den3;
					f3m=a3/den3;
				}
				if (real(r4)>0)
				{
					exptmp=exp(-r4/tem);
					den4=exptmp+1.0;
					f4p=a4/den4;
					f4m=a4*exptmp/den4;
				}
				else
				{
					exptmp=exp(r4/tem);
					den4=exptmp+1.0;
					f4p=a4*exptmp/den4;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p+r3*f3p+r4*f4p);
				dtau_Gkm+=real(r1*f1m+r2*f2m+r3*f3m+r4*f4m);
				
				//TF de G(k,ikn) dans le temps
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				Gbeta=Gtmp1[0];
				
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					
					if (real(r1)>0)
					{
						exp1=exp(-t*r1);
					}
					else
					{
						exp1=exp((1.0/tem-t)*r1);
					}
					if (real(r2)>0)
					{
						exp2=exp(-t*r2);
					}
					else
					{
						exp2=exp((1.0/tem-t)*r2);
					}
					if (real(r3)>0)
					{
						exp3=exp(-t*r3);
					}
					else
					{
						exp3=exp((1.0/tem-t)*r3);
					}
					if (real(r4)>0)
					{
						exp4=exp(-t*r4);
					}
					else
					{
						exp4=exp((1.0/tem-t)*r4);
					}
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
					
					ph=exp(dcomplex(0.0, (j*PI*(nbw-1))/nbw));
					
					Gtmp1[j]=ph*Gtmp1[j]+G0t;
					
				}
				
				t=1.0/tem;
				
				if (real(r1)>0)
				{
					exp1=exp(-t*r1);
				}
				else
				{
					exp1=1.0;
				}
				if (real(r2)>0)
				{
					exp2=exp(-t*r2);
				}
				else
				{
					exp2=1.0;
				}
				if (real(r3)>0)
				{
					exp3=exp(-t*r3);
				}
				else
				{
					exp3=1.0;
				}
				if (real(r4)>0)
				{
					exp4=exp(-t*r4);
				}
				else
				{
					exp4=1.0;
				}
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
				
				ph=exp(dcomplex(0.0, PI*(nbw-1)));
				
				Gbeta=ph*Gbeta+G0t;
				
				/*
				 if ((l+1)%4==0)
				 {
				 cout<<"dtau_Gkp:  "<<setw(20)<<(real(Gtmp1[1]-Gtmp1[0])*tem*nbw-dtau_Gkp)/dtau_Gkp<<'\n';
				 cout<<"dtau_Gkm:  "<<setw(20)<<(real(-Gtmp1[nbw-1]-Gtmp1[0]-1.0)*tem*nbw-dtau_Gkm)/dtau_Gkm<<'\n';
				 cout<<"dtau_J0:  "<<setw(20)<<real(Jtmp[1]-Jtmp[0])*tem*nbw<<setw(20)<<dtau_J0<<(real(Jtmp[1]-Jtmp[0])*tem*nbw-dtau_J0)/dtau_J0<<endl;
				 cout<<"dtau_Jbeta:  "<<setw(20)<<real(Jbeta-Jtmp[nbw-1])*tem*nbw<<setw(20)<<dtau_Jbeta<<(real(Jbeta-Jtmp[nbw-1])*tem*nbw-dtau_Jbeta)/dtau_Jbeta<<endl;
				 }
				 */
				
				// Definition de -G(k,t)G(k,-t)
				GGtmp1[0]=Gtmp1[0]*(Gtmp1[0]+1.0);
				for (j=1; j<nbw; j++)
					GGtmp1[j]=-Gtmp1[j]*Gtmp1[nbw-j];
				
				// calcul du spline de -G(k,t)G(k,-t)
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nbw*tem);
					ft_tmp[j]=real(GGtmp1[j]);
				}
				dtau_chi=dtau_Gkp*(real(Gtmp1[0])+1.0)+real(Gtmp1[0])*dtau_Gkm;
				
				//				if ((l+1)%4==0)
				//					cout<<"dtau_chi:  "<<setw(20)<<(real(GGtmp1[1]-GGtmp1[0])*(tem*nbw)-dtau_chitmp)/dtau_chitmp<<setw(20)<<(dtau_chi-dtau_chitmp)/dtau_chitmp<<'\n';
				
				coeffs[0]=dtau_chi;
				coeffs[1]=0;
				spline_coeffs(tau, ft_tmp, NS0, coeffs);
				
				// TF sur t de -G(k,t)G(k,-t)
				for (j=0; j<nbw/2; j++)
					Gtmp1[j]=6.0*coeffs[4*j];
				for (j=nbw/2; j<nbw; j++)
					Gtmp1[j]=-Gtmp1[nbw-j-1];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				integ=0;
				for (j=0; j<nbw/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j];
					x2=tau[j+1];
					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
				
				Gtmp1[0]=2*integ;
				
				wt=4;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+(l+1)*nbk]*k0;
				chijj1_qn_loc[0]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_w0_loc-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_inf2_loc+=4.0*dtau_chi*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk];
				
				for (j=1; j<=nw; j++)
				{
					wn=2*j*PI*tem;
					
					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real((1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]);
					
					Gtmp1[j]=-2.0*dtau_chi/(wn*wn)+real((1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]/(wn*wn*wn*wn));
					
					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
				}
				
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
			
			for (j=0; j<nw; j++)
			{
				chijj1_qn[j]+=chijj1_qn_loc[j];
				chijj1_qn4[j]+=chijj1_qn4_loc[j];
			}
			chijj1_qn[nw]+=chijj1_qn_loc[nw];
			
			chijj1_w0+=chijj1_w0_loc;
			chijj1_inf2+=chijj1_inf2_loc;
		}
		delete [] chijj1_qn4_loc;
		delete [] chijj1_qn_loc;
		delete [] Gtmp1;
		delete [] GGtmp1;
		delete [] tau;
		delete [] ft_tmp;
		delete [] coeffs;
		delete [] tau2;
		delete [] ft_tmp2;
		delete [] coeffs2;
	}
	
#pragma omp parallel private(m,j)
	{
		int x, y, tmp;
		double wt, kn;
		dcomplex k0;
		dcomplex Gtmp2, Gtmp0;
		
		dcomplex k0sum_loc=0;
		
		double C2, C3;
		dcomplex z, selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex f1m, f2m, f3m, f4m;
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk; l+=nk0/2)
		{
			for (m=0; m<nbk; m++)
			{
				x=l;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				//				C2=G1->Sigma_inf2[y+(x*(x+1))/2];
				//				C3=G1->Sigma_inf3[y+(x*(x+1))/2];
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				k0=0;
				
				/*
				 for (j=0; j<nw; j++)
				 {
				 kn=(2.0*(j-nbw)+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th]));
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;
				 }
				 for (j=nbw-1; j>=nw; j--)
				 {
				 kn=(2.0*j+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+j*nk8th]);
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;
				 }
				 */
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0)
				{
					den1=exp(-r1/tem)+1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else
				{
					den1=exp(r1/tem)+1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					den2=exp(-r2/tem)+1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else
				{
					den2=exp(r2/tem)+1.0;
					f2m=a2/den2;
				}
				if (real(r3)>0)
				{
					den3=exp(-r3/tem)+1.0;
					f3m=a3*exp(-r3/tem)/den3;
				}
				else
				{
					den3=exp(r3/tem)+1.0;
					f3m=a3/den3;
				}
				if (real(r4)>0)
				{
					den4=exp(-r4/tem)+1.0;
					f4m=a4*exp(-r4/tem)/den4;
				}
				else
				{
					den4=exp(r4/tem)+1.0;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				wt=2;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+l*nbk]*k0;
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
		}
	}
	
	//	cout<<"chitmp3:  ";
	//	calc_chijj_vertex_corr2_wn(0, k0, 0, Htmp, NULL, NULL, NULL, NULL);
	
	delete [] d2ek;
	
	//	for (j=0; j<4; j++)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	//	for (j=nbw/8; j<=nbw/2; j+=nbw/8)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	
	fstream file;
	char name[200];
	
	kx0=real(k0sum);
	sprintf(name, "kx0_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat", U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out);
	file<<setprecision(15)<<setiosflags(ios::left)<<kx0<<'\n';
	file.close();
	
	sprintf(name, "kx0_bin_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat", U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&kx0, sizeof(kx0));
	file.close();
	
	for (j=0; j<=nw; j++)
	{
		if (chijj1_qn[j]<0.0)
		{
			cout<<"attention! valeur non-physique de chi_jj\n";
			cout<<setw(30)<<chijj1_qn[j]<<'\n';
			break;
		}
	}
	
	
	//pour enregistrer chijj(iq_n) en binaire
	
	
	const char *nameForm_chijj1_qn_bin="./chijj1_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn_bin, U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	
	
	j=0;
	file.write((char*)&j,sizeof(j));
	file.write((char*)&(chijj1_qn[j]),sizeof(double));
	for (j=1; j<=nbw/2; j++)
	{
		file.write((char*)&j,sizeof(j));
		file.write((char*)&(chijj1_qn[j]),sizeof(double));
	}
	file.close();
	// fin de l'enregistrement en binaire
	
	// pour enregistrer chijj(iq_n) en ascii
	
	const char *nameForm_chijj1_qn="./chijj1_qn_%5.3f_%5.3f_%5.3f_%5.3f_%5.3f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn, U, tp, tpp, density, tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left);
	
	j=0;
	file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	for (j=1; j<=nbw/2; j++)
	{
		file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	}
	file.close();
	// fin de l'enregistrement en ascii
	
	cout<<"-k0:  "<<-kx0<<'\n';
	cout<<"chijj1(iqn=0):  "<<chijj1_qn[0]<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+kx0)/kx0<<'\n';
	
	
	fftw_destroy_plan(fftplan);
	fftw_destroy_plan(fftplan_t_back);
	
	delete [] dek;
	
	if (chijj_bulle_qn) delete [] chijj_bulle_qn;
	chijj_bulle_qn=new double[nw+1];
	
	for (j=0; j<=nw; j++) chijj_bulle_qn[j]=chijj1_qn[j];
}

void cond_opt::calc_chijj_bulle(int N_theta_k, int *k_chijj, int nk_chijj)
{	
	cout<<"calcul de chijj1 par calc_chijj_bulle()\n";
	
	set_keep_chiqw(false);
	
	long int j, l, m;
	
	int nw=nbw/2;
	int nk0=2*(nbk-1);
	int nk4th=nbk*(nbk-2);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (!Self_kw_array)
	{
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}
	
	if (!mu_calcule)
	{
		cout<<"calc_chijj_bulle(): mu non calcule\n";
		return;
	}
	
//	save_self_ikn(k_chijj, nk_chijj);
	find_FS_inter();
	
	double normFact=1.0/(nk0*nk0);
	
	double *dek=new double[nbk*(nbk-2)];
	double *d2ek=new double[nbk*nbk];
	int indHss[]={0,0};
	double vk[2];
	
	double kx, ky;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			d2ek[m+l*nbk]=Hessian_ek(vk, indHss);
		}
	}
	
	
	for (l=0; l<nbk-2; l++)
	{
		kx=(l+1)*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			dek[m+l*nbk]=grad_ek(vk, 0);
		}
	}
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex *Gtmp=new dcomplex[nbw];
	fftw_plan fftplan=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (Gtmp), reinterpret_cast<fftw_complex*> (Gtmp), FFTW_FORWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	if (chijj1_qn) delete [] chijj1_qn;
	
	chijj1_qn=new double[nw+1];
	chijj1_qn4=new double[nw];
	
	for (j=0; j<nw; j++)	
	{
		chijj1_qn[j]=0;
		chijj1_qn4[j]=0;
	}
	chijj1_qn[nw]=0;
	
	chijj1_w0=0;
	chijj1_inf2=0;
	
	int NS0=nbw/2+1;
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex k0sum=0;
	
	double *chijj1_k=new double[(nw+1)*N_theta_k];
	
	for (j=0; j<(nw+1)*N_theta_k; j++)
	{
		chijj1_k[j]=0;
	}
	
	double Dtheta=PI/(2.0*N_theta_k);
	
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp, p;
		double wt, t, kn, G0t, Jinf;
		dcomplex ph, dcplx;
		dcomplex k0;
				
		dcomplex *Gtmp1=new dcomplex[nbw];
		dcomplex *GGtmp1=new dcomplex[nbw];
		dcomplex Gtmp0;
		dcomplex Gtmp2;
		
		double *coeffs=new double[4*(NS0-1)];
		double *tau=new double[NS0];
		double *ft_tmp=new double[NS0];
		double wn, integ, a, b, c, d, x1, x2;
		
		double *tau2=new double[nbw+1];
		double *ft_tmp2=new double[nbw+1];
		double *coeffs2=new double[4*nbw];
		
		dcomplex k0sum_loc=0;
		dcomplex chijj2_qn0_loc;
		double *chijj1_qn_loc=new double[nw+1];
		double *chijj1_qn4_loc=new double[nw];
		for (j=0; j<nw; j++)	
		{
			chijj1_qn_loc[j]=0;
			chijj1_qn4_loc[j]=0;
		}
		chijj1_qn_loc[nw]=0;
		
		double chijj2_inf2_loc=0, chijj2_inf4_loc=0, chijj1_inf2_loc=0, chijj1_w0_loc=0;
		
		double dtau_Gkp, dtau_Gkm, dtau_chi;
		double C2, C3;
		dcomplex selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex exp1, exp2, exp3, exp4;
		dcomplex f1p, f2p, f3p, f4p, f1m, f2m, f3m, f4m;
		
		dcomplex exptmp;
		dcomplex b1, b2, b3, b4;
		dcomplex denj1, denj2, denj3, denj4;
		dcomplex expj1, expj2, expj3, expj4;
		dcomplex fj1, fj2, fj3, fj4, fj1b, fj2b, fj3b, fj4b;
		
		double theta_k;
		int ind_theta;
		double *chijj1_k_loc=new double[(nw+1)*N_theta_k];
		
		for (j=0; j<(nw+1)*N_theta_k; j++)
		{
			chijj1_k_loc[j]=0;
		}
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
		
#pragma omp for
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
								
				k0=0;
				dtau_Gkp=0;
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					dcplx=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					Gtmp1[j]=dcomplex(dcplx.real(), dcplx.imag());
										
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));

					dcplx=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp2=dcomplex(dcplx.real(), dcplx.imag());
//					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
										
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					dcplx=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					Gtmp1[j]=dcomplex(dcplx.real(), dcplx.imag());
										
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);

					dcplx=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp2=dcomplex(dcplx.real(), dcplx.imag());
//					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
										
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
				}
				
				// calcul du G(k,t) asymptotique	
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];				
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0) 
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
				}
				else 
				{
					exptmp=exp(r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
				}
				if (real(r2)>0) 
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
				}
				else 
				{	
					exptmp=exp(r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
				}
				if (real(r3)>0) 
				{
					exptmp=exp(-r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3/den3;
					f3m=a3*exptmp/den3;
				}
				else 
				{
					exptmp=exp(r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3*exptmp/den3;
					f3m=a3/den3;
				}
				if (real(r4)>0) 
				{
					exptmp=exp(-r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4/den4;
					f4m=a4*exptmp/den4;
				}
				else 
				{
					exptmp=exp(r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4*exptmp/den4;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p+r3*f3p+r4*f4p);
				dtau_Gkm+=real(r1*f1m+r2*f2m+r3*f3m+r4*f4m);		
				
				//TF de G(k,ikn) dans le temps
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
//				Gbeta=Gtmp1[0];
				
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					
					if (real(r1)>0) 
					{
						exp1=exp(-t*r1);
					}
					else 
					{
						exp1=exp((1.0/tem-t)*r1);
					}
					if (real(r2)>0) 
					{
						exp2=exp(-t*r2);
					}
					else 
					{
						exp2=exp((1.0/tem-t)*r2);
					}
					if (real(r3)>0) 
					{
						exp3=exp(-t*r3);
					}
					else 
					{
						exp3=exp((1.0/tem-t)*r3);
					}
					if (real(r4)>0) 
					{
						exp4=exp(-t*r4);
					}
					else 
					{
						exp4=exp((1.0/tem-t)*r4);
					}
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
					
					ph=exp(dcomplex(0.0, (j*PI*(nbw-1))/nbw));
					
					Gtmp1[j]=ph*Gtmp1[j]+G0t;					
					
				}			
				
				t=1.0/tem;
				
				if (real(r1)>0) 
				{
					exp1=exp(-t*r1);
				}
				else 
				{
					exp1=1.0;
				}
				if (real(r2)>0) 
				{
					exp2=exp(-t*r2);
				}
				else 
				{
					exp2=1.0;
				}
				if (real(r3)>0) 
				{
					exp3=exp(-t*r3);
				}
				else 
				{
					exp3=1.0;
				}
				if (real(r4)>0) 
				{
					exp4=exp(-t*r4);
				}
				else 
				{
					exp4=1.0;
				}
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
				
				ph=exp(dcomplex(0.0, PI*(nbw-1)));
				
//				Gbeta=ph*Gbeta+G0t;
				
				/*								
				 if ((l+1)%4==0)
				 {				
				 cout<<"dtau_Gkp:  "<<setw(20)<<(real(Gtmp1[1]-Gtmp1[0])*tem*nbw-dtau_Gkp)/dtau_Gkp<<'\n';
				 cout<<"dtau_Gkm:  "<<setw(20)<<(real(-Gtmp1[nbw-1]-Gtmp1[0]-1.0)*tem*nbw-dtau_Gkm)/dtau_Gkm<<'\n';
				 cout<<"dtau_J0:  "<<setw(20)<<real(Jtmp[1]-Jtmp[0])*tem*nbw<<setw(20)<<dtau_J0<<(real(Jtmp[1]-Jtmp[0])*tem*nbw-dtau_J0)/dtau_J0<<endl;
				 cout<<"dtau_Jbeta:  "<<setw(20)<<real(Jbeta-Jtmp[nbw-1])*tem*nbw<<setw(20)<<dtau_Jbeta<<(real(Jbeta-Jtmp[nbw-1])*tem*nbw-dtau_Jbeta)/dtau_Jbeta<<endl;
				 }
				 */
								
				// Definition de -G(k,t)G(k,-t)				
				GGtmp1[0]=Gtmp1[0]*(Gtmp1[0]+(double)1.0);
				for (j=1; j<nbw; j++)
					GGtmp1[j]=-Gtmp1[j]*Gtmp1[nbw-j];
				
				// calcul du spline de -G(k,t)G(k,-t)			
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nbw*tem);
					ft_tmp[j]=real(GGtmp1[j]);
				}
				dtau_chi=dtau_Gkp*(real(Gtmp1[0])+1.0)+real(Gtmp1[0])*dtau_Gkm;
				
				//				if ((l+1)%4==0)
				//					cout<<"dtau_chi:  "<<setw(20)<<(real(GGtmp1[1]-GGtmp1[0])*(tem*nbw)-dtau_chitmp)/dtau_chitmp<<setw(20)<<(dtau_chi-dtau_chitmp)/dtau_chitmp<<'\n';
				
				coeffs[0]=dtau_chi;
				coeffs[1]=0;				
//				spline_coeffs(tau, ft_tmp, NS0, coeffs);
				spline_coeffs_rel(tau, ft_tmp, NS0, coeffs);				
				
				// TF sur t de -G(k,t)G(k,-t)			
				for (j=0; j<nbw/2; j++)
					Gtmp1[j]=6.0*coeffs[4*j];
				for (j=nbw/2; j<nbw; j++)
					Gtmp1[j]=-Gtmp1[nbw-j-1];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));

				integ=0;
				for (j=0; j<nbw/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j+1]-tau[j];
					integ+=a*x1*x1*x1*x1/4.0+b*x1*x1*x1/3.0+c*x1*x1/2.0+d*x1;
				}				
/*				
				integ=0;
				for (j=0; j<nbw/2; j++)
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
				Gtmp1[0]=2*integ;
				
				wt=4;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+(l+1)*nbk]*k0;
				chijj1_qn_loc[0]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_w0_loc-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_inf2_loc+=4.0*dtau_chi*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk];
								
				theta_k=atan( (1.0*(nk0/2-m))/(nk0/2-l-1) );
				ind_theta=(int)floor(theta_k/Dtheta);
								
				chijj1_k_loc[ind_theta]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				
				for (j=1; j<=nw; j++)
				{
					wn=2*j*PI*tem;
					
					dcplx=dcomplex(real(Gtmp1[j]), imag(Gtmp1[j]));

					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*dcplx);					
//					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]);

					dcplx=-2.0*dtau_chi/(wn*wn)+real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*dcplx)/(wn*wn*wn*wn);					
//					Gtmp1[j]=-2.0*dtau_chi/(wn*wn)+real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j])/(wn*wn*wn*wn);

					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(dcplx);					
//					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);

					chijj1_k_loc[ind_theta+j*N_theta_k]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(dcplx);					
//					chijj1_k_loc[ind_theta+j*N_theta_k]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
				}				
				
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
			
			for (j=0; j<nw; j++)
			{
				chijj1_qn[j]+=chijj1_qn_loc[j];
				chijj1_qn4[j]+=chijj1_qn4_loc[j];
				for (p=0; p<N_theta_k; p++)
				{
					chijj1_k[p+j*N_theta_k]+=chijj1_k_loc[p+j*N_theta_k];
				}
			}
			chijj1_qn[nw]+=chijj1_qn_loc[nw];
			for (p=0; p<N_theta_k; p++)
			{
				chijj1_k[p+nw*N_theta_k]+=chijj1_k_loc[p+nw*N_theta_k];
			}
			
			chijj1_w0+=chijj1_w0_loc;
			chijj1_inf2+=chijj1_inf2_loc;
		}
		delete [] chijj1_qn4_loc;
		delete [] chijj1_qn_loc;
		delete [] Gtmp1;
		delete [] GGtmp1;
		delete [] tau;
		delete [] ft_tmp;
		delete [] coeffs;
		delete [] tau2;
		delete [] ft_tmp2;
		delete [] coeffs2;
		delete [] chijj1_k_loc;
	}
	
	double *chijj1_k_tot=new double[nw+1];
	
	for (j=0; j<=nw; j++)
	{
		chijj1_k_tot[j]=0;
	}
	
	int p;
	
	for (j=0; j<=nw; j++)
		for  (p=0; p<N_theta_k; p++)
		{
			chijj1_k_tot[j]+=chijj1_k[p + j*N_theta_k];
		}
	
	fstream file, file2, file_info;
	char name[200];

	const char *nameForm_chijj1_k="./chijj1_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj1_k_tot="./chijj1_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	
	sprintf(name, nameForm_chijj1_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left);
	
	sprintf(name, nameForm_chijj1_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file2.open(name, ios::out );
	file2<<setprecision(16)<<setiosflags(ios::left);
		
	for (j=0; j<=nw; j++)
	{
		file2<<setw(10)<<j<<setw(30)<<chijj1_k_tot[j]<<setw(30)<<chijj1_qn[j]<<chijj1_k_tot[j]/chijj1_qn[j]<<'\n';
		for  (p=0; p<N_theta_k; p++)
		{
			file<<setw(10)<<j<<setw(10)<<p+1<<chijj1_k[p + j*N_theta_k]<<'\n';
		}
	}
	file.close();
	file2.close();
	
	{
		const char *nameForm_self_k="./self_k_ikn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
		sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
		file.open(name, ios::out );
		file<<setprecision(16)<<setiosflags(ios::left);
		
		dcomplex selftmp;
		long int x, y, tmp;
		double C2, C3, kn;
		for (j=0; j<=nw/16; j++)
		{
			kn=(2.0*j+1)*PI*tem;
			for (p=0; p<nk_chijj; p++)
			{
				for (l=-1; l<=1; l++)
				{
					x=k_chijj[2*p]+l;
					y=k_chijj[2*p+1]+l;
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
	}
	
/*	
	{
		const char *nameForm_self_k="./self_k_ikn_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
		sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
		file.open(name, ios::out );
		file<<setprecision(16)<<setiosflags(ios::left);
		
		dcomplex selftmp;
		long int x, y, tmp;
		double C2, C3, kn;
		for (j=0; j<=nw/16; j++)
		{
			kn=(2.0*j+1)*PI*tem;
			for (p=0; p<nk_chijj; p++)
			{
				x=k_chijj[2*p];
				y=k_chijj[2*p+1];
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
				
				file<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
			}
		}
		file.close();
	}
*/	
	delete [] chijj1_k;
	delete [] chijj1_k_tot;
	
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp;
		double wt, kn;
		dcomplex k0;
		dcomplex Gtmp2, Gtmp0;
		
		dcomplex k0sum_loc=0;
		
		double C2, C3;
		dcomplex z, selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex f1m, f2m, f3m, f4m;
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk; l+=nk0/2)
		{
			for (m=0; m<nbk; m++)
			{
				x=l;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				//				C2=G1->Sigma_inf2[y+(x*(x+1))/2];
				//				C3=G1->Sigma_inf3[y+(x*(x+1))/2];				
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
				
				k0=0;
				
				/*				
				 for (j=0; j<nw; j++)
				 {
				 kn=(2.0*(j-nbw)+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th]));
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;
				 }
				 for (j=nbw-1; j>=nw; j--)
				 {
				 kn=(2.0*j+1)*PI*tem;
				 Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+j*nk8th]);
				 Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
				 Gtmp2=Gtmp2-Gtmp0;
				 
				 k0+=Gtmp2;	
				 }
				 */				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0) 
				{
					den1=exp(-r1/tem)+(double)1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else 
				{
					den1=exp(r1/tem)+(double)1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0) 
				{
					den2=exp(-r2/tem)+(double)1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else 
				{	
					den2=exp(r2/tem)+(double)1.0;
					f2m=a2/den2;
				}
				if (real(r3)>0) 
				{
					den3=exp(-r3/tem)+(double)1.0;
					f3m=a3*exp(-r3/tem)/den3;
				}
				else 
				{
					den3=exp(r3/tem)+(double)1.0;
					f3m=a3/den3;
				}
				if (real(r4)>0) 
				{
					den4=exp(-r4/tem)+(double)1.0;
					f4m=a4*exp(-r4/tem)/den4;
				}
				else 
				{
					den4=exp(r4/tem)+(double)1.0;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);				
				
				wt=2;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+l*nbk]*k0;
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
		}
	}
		
	delete [] d2ek;
	
	kx0=real(k0sum);
	sprintf(name, "kx0_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out);
	file<<setprecision(15)<<setiosflags(ios::left)<<kx0<<'\n';
	file.close();
	
	sprintf(name, "kx0_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&kx0, sizeof(kx0));
	file.close();
	
	for (j=0; j<=nw; j++)
	{
		if (chijj1_qn[j]<0.0)
		{
			cout<<"attention! valeur non-physique de chi_jj\n";
			cout<<setw(30)<<chijj1_qn[j]<<'\n';
			break;
		}
	}
	
	
	//pour enregistrer chijj(iq_n) en binaire
	
	
	const char *nameForm_chijj1_qn_bin="./chijj1_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	
	
	j=0;
	file.write((char*)&j,sizeof(j));
	file.write((char*)&(chijj1_qn[j]),sizeof(double));
	for (j=1; j<=nbw/2; j++)
	{
		file.write((char*)&j,sizeof(j));
		file.write((char*)&(chijj1_qn[j]),sizeof(double));
	}
	file.close();
	// fin de l'enregistrement en binaire
	
	// pour enregistrer chijj(iq_n) en ascii
	
	const char *nameForm_chijj1_qn="./chijj1_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	j=0;
	file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	for (j=1; j<=nbw/2; j++)
	{
		file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	}
	file.close();
	// fin de l'enregistrement en ascii
	
	cout<<"-k0:  "<<-kx0<<'\n';
	cout<<"chijj1(iqn=0):  "<<chijj1_qn[0]<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+kx0)/kx0<<'\n';
	
	
	fftw_destroy_plan(fftplan);
	
	delete [] dek;
	
	if (chijj_bulle_qn) delete [] chijj_bulle_qn;
	chijj_bulle_qn=new double[nw+1];
	
	for (j=0; j<=nw; j++) chijj_bulle_qn[j]=chijj1_qn[j];
}

//! calculate all the terms in the optical conductivity with FFT and cubic spline. That version uses less memory than calc_chijj_vertex_corr_all() and is more stable because it only uses 2nd degree roots for the asymptotic G
void cond_opt::calc_chijj_vertex_corr_all_optim(bool chi3, int mmin, int N_theta_k, int *k_chijj, int nk_chijj)
{	
	bool calc_chijj2_qn0=false;
	
	cout<<"calcul de chijj par calc_chijj_vertex_corr_all_optim()\n";
	
	long int j, l, m;
	
	int nw=nbw/2;
	int nk0=2*(nbk-1);
	int nk4th=nbk*(nbk-2);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (!Self_kw_array)
	{
		//		G1->calc_Self_FFT_spline();
		//		G1->traceSelfG();
		//		G1->find_mu();
		//		mu=G1->mu;
		
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}	
	
	if (!mu_calcule)
	{
		cout<<"calc_chijj_vertex_corr_all_optim(): mu non calcule\n";
		return;
	}
	
	if (k_chijj && nk_chijj)	save_self_ikn(k_chijj, nk_chijj, nw/16);
	find_FS_inter();
		
	V_array=vertex_rt_array;
	
	double normFact=1.0/(nk0*nk0);
	
	double *dek=new double[nbk*(nbk-2)];
	double *d2ek=new double[nbk*nbk];
	int indHss[]={0,0};
	double vk[2];
	
	double kx, ky;
	
	double *ekmu=new double[nbk*nbk];
	
	for (l=0; l<nbk; l++)
		for(m=0; m<nbk; m++)
			ekmu[m+l*nbk]=ek[m+l*nbk]-mu0;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			d2ek[m+l*nbk]=Hessian_ek(vk, indHss);
		}
	}
	
	
	for (l=0; l<nbk-2; l++)
	{
		kx=(l+1)*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			dek[m+l*nbk]=grad_ek(vk, 0);
		}
	}
	
	double  *fr_tmp=new double[nk4th];
	fftw_plan fftplan_RO=fftw_plan_r2r_2d(nbk-2, nbk, fr_tmp, fr_tmp, FFTW_RODFT00, FFTW_REDFT00, FFTW_MEASURE);
	delete [] fr_tmp;
	
	double *fr_tmp1;
	
	long int array_size_h=(long int)(nbk*(nbk-2)*(nbw+1));
	
	if (contains_NaN(V_array, nw*nk8th)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans V_array \n";
	fflush(stdout);
	
	double *h=new double[array_size_h];	
//	double *h=new double[nbk*(nbk-2)*(nbw+1)];
	
//	if (contains_NaN(h, array_size_h)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans h (ligne 1061)\n";
//	fflush(stdout);

//	cout<<"0\n";
//	fflush(stdout);
	
	for (j=0; j<=nbw; j++)
		for (l=0; l<(nbk-2); l++)
			for (m=0; m<nbk; m++)
				h[m+l*nbk+j*nk4th]=0;
	
//	cout<<"1\n";
//	fflush(stdout);

	if (!h) 
	{
		cout<<"h=NULL\n";
		fflush(stdout);
		return;
	}
	
// calcul de h(r,t) (voir l'annexe sur les techniques de calcul)	
/*	
	dcomplex *fr_c=new dcomplex[nk0*nk0];
	fftw_plan fftplan2D=fftw_plan_dft_2d(nk0, nk0, reinterpret_cast<fftw_complex*> (fr_c) , reinterpret_cast<fftw_complex*> (fr_c), FFTW_BACKWARD, FFTW_MEASURE);
	delete [] fr_c;
	
//#pragma omp parallel private(fr_c,j,l,m)
	{
		fr_c=new dcomplex[nk0*nk0];
		double t;
		
//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			for (l=0; l<nk0; l++)
				for (m=0; m<nk0; m++)
					fr_c[m+l*nk0]=dcomplex(0,0);
			t=j/(tem*nbw);
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
						fr_c[m+(l+1)*nk0]=dcomplex(normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + (double)1.0 ),0);
					else
						fr_c[m+(l+1)*nk0]=dcomplex(normFact*dek[m+l*nbk]*exp( (ek[m+(l+1)*nbk]-mu0)*((double)1.0/tem-t))/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  + (double)1.0 ),0);
					
					if (m)	
					{
						fr_c[nk0-m+(l+1)*nk0]=fr_c[m+(l+1)*nk0];
						fr_c[nk0-m+(nk0-l-1)*nk0]=-fr_c[nk0-m+(l+1)*nk0];
					}
					fr_c[m+(nk0-l-1)*nk0]=-fr_c[m+(l+1)*nk0];
				}
			
			fftw_execute_dft(fftplan2D,reinterpret_cast<fftw_complex*> (fr_c),reinterpret_cast<fftw_complex*> (fr_c));
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					h[m+l*nbk+j*nk4th]=fr_c[m+(l+1)*nk0].imag();
		}
		delete [] fr_c;
	}
*/
	
//	fstream file_test;	
//	file_test.open("file_h_test.dat",ios::out);
//	file_test<<setiosflags(ios::left);
	
//#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[nk4th];
		double t;
		
//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			t=(double) j/(tem*nbw);
//			if (j) cout<<j<<"  1"<<endl;
			for (l=0; l<(nbk-2); l++)
				for (m=0; m<nbk; m++)
				{
//					if (ek[m+(l+1)*nbk]>mu0)
//						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + (double)1.0 );
//					else
//						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp( (ek[m+(l+1)*nbk]-mu0)*((double)(1.0/tem)-t))/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  + (double)1.0 );
					
					if (ekmu[m+(l+1)*nbk]>0)
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp(-ekmu[m+(l+1)*nbk]*t)/(exp(-ekmu[m+(l+1)*nbk]/tem)  + (double)1.0 );
					else
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp(ekmu[m+(l+1)*nbk]*((double)(1.0/tem)-t))/(exp(ekmu[m+(l+1)*nbk]/tem)  + (double)1.0 );
//					if (!(j%1024)) file_test<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<fr_tmp1[m+l*nbk]<<endl;
//					file_test.flush();
//					if (j==0) cout<<setw(10)<<l<<setw(10)<<m<<ek[m+(l+1)*nbk]<<endl;
//					fflush(stdout);
				}
//			if (j) cout<<j<<"  2"<<endl;
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
//			if (j) cout<<j<<"  3"<<endl;
			for (l=0; l<(nbk-2); l++)
				for (m=0; m<nbk; m++)
					h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
			
//			if (j) cout<<j<<"  4"<<endl;
		}
		delete [] fr_tmp1;
	}
//	file_test.close();	
	
//	cout<<"2\n";
//	fflush(stdout);
	
	if (contains_NaN(h, array_size_h)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans h\n";
	fflush(stdout);
	
	// calcul de la derivee de h(r,t) a t=0 et t=beta
	double *dh=new double[(nbk-2)*nbk*2];
	
	fr_tmp1=new double[(nbk-2)*nbk];
	
	
	
	// t=0
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
			else
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (ek[m+(l+1)*nbk]-mu0)/tem)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dh[m+l*nbk]=fr_tmp1[m+l*nbk];
	
	// t=beta	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (mu0-ek[m+(l+1)*nbk])/tem )/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
			else
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dh[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];	
	
	delete [] fr_tmp1;
	
	if (contains_NaN(dh, (long int)(2*(nbk-2)*nbk))) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans dh \n";
	fflush(stdout);
	
	double *f0=new double[array_size_h];	
//	double *f0=new double[(nbk-2)*nbk*(nbw+1)];
	
//	cout<<"3\n";
//	fflush(stdout);
	
// calcul de f_0(r,t) (voir l'annexe sur les techniques de calcul)	
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[(nbk-2)*nbk];
		double t, expk;
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			t=j/(tem*nbw);
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
					{
						expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)*( t/(expk + 1.0) -expk/(tem*(expk+1.0)*(expk+1.0)));
					}
					else
					{
						expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((ek[m+(l+1)*nbk]-mu0)*(1.0/tem-t))*(t/(expk+1.0) - 1.0/(tem*(expk+1.0)*(expk+1.0)));
					}
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					f0[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		
		delete [] fr_tmp1;
	}
		
	if (contains_NaN(f0, array_size_h)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans f0 \n";
	fflush(stdout);
	
//	cout<<"4\n";
//	fflush(stdout);
	
	// calcul de la derivee de  f0(r,t) a t=0 et t=beta
	double *df0=new double[(nbk-2)*nbk*2];
	
	fr_tmp1=new double[(nbk-2)*nbk];
	double expk;
	
	//	t=0
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
			{
				expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
				fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*(  (ek[m+(l+1)*nbk]-mu0)*expk/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk + 1.0));
			}
			else
			{
				expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
				fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*( (ek[m+(l+1)*nbk]-mu0)/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk+1.0));
			}
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			df0[m+l*nbk]=fr_tmp1[m+l*nbk];
		}
	
	// t=beta	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
			{
				expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
				fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*(  (mu0-ek[m+(l+1)*nbk])*( 1.0/(tem*(expk + 1.0)) -expk/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk + 1.0));
			}
			else
			{
				expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
				fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*( (mu0-ek[m+(l+1)*nbk])*(1.0/(tem*(expk+1.0)) - 1.0/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk+1.0));
			}
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			df0[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
		}
	
	if (contains_NaN(df0, (long int)(2*(nbk-2)*nbk))) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans df0 \n";
	fflush(stdout);
	
	double TPSC_fact=U*(3.0*Usp + Uch)/4;
	
	double *dhV=new double[nbk*(nbk-2)*2];
	
//	cout<<"5\n";
//	fflush(stdout);
	
	// TF sur l'espace de d/dt(h(r,t)V(r,t)) a t=0 et t=beta
		
	long int x, y, tmp;
	
	//t=0
	for (l=0; l<nbk-2; l++)
	{
		for (m=0; m<nbk; m++)
		{
			x=l+1;
			if (x>nk0/2) x=nk0-x;
			y=m;
			if (y>nk0/2) y=nk0-y;
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk]*V_array[y+(x*(x+1))/2]+ h[m+l*nbk]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
		}
	}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dhV[m+l*nbk]=fr_tmp1[m+l*nbk];
	
	//t=beta
	for (l=0; l<nbk-2; l++)
	{
		for (m=0; m<nbk; m++)
		{
			x=l+1;
			if (x>nk0/2) x=nk0-x;
			y=m;
			if (y>nk0/2) y=nk0-y;
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- h[m+l*nbk+nbw*nk4th]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
		}
	}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dhV[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
	
	if (contains_NaN(dhV, (long int)(2*(nbk-2)*nbk))) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans dhV \n";
	fflush(stdout);
		
	// TF sur l'espace de d/dt f_0(r,t)V(r,t) a t=0 et t=beta
	double *df0V=new double[nbk*(nbk-2)*2];
	
//	cout<<"6\n";
//	fflush(stdout);
	
	//t=0
	for (l=0; l<nbk-2; l++)
	{
		for (m=0; m<nbk; m++)
		{
			x=l+1;
			if (x>nk0/2) x=nk0-x;
			y=m;
			if (y>nk0/2) y=nk0-y;
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk]*V_array[y+(x*(x+1))/2]+ f0[m+l*nbk]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
		}
	}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			df0V[m+l*nbk]=fr_tmp1[m+l*nbk];
	
	//t=beta
	for (l=0; l<nbk-2; l++)
	{
		for (m=0; m<nbk; m++)
		{
			x=l+1;
			if (x>nk0/2) x=nk0-x;
			y=m;
			if (y>nk0/2) y=nk0-y;
			if (y>x)
			{
				tmp=x;
				x=y;
				y=tmp;
			}
			
			fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- f0[m+l*nbk+nbw*nk4th]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
		}
	}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			df0V[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
	
	delete [] fr_tmp1;	

	if (contains_NaN(df0V, (long int)(2*(nbk-2)*nbk))) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans df0V \n";
	fflush(stdout);
	
//	cout<<"7\n";
//	fflush(stdout);

//	double *hV=new double[array_size_h];
//	double *hV=new double[nbk*(nbk-2)*(nbw+1)];
	
	// TF sur l'espace de h(r,t)V(r,t), ( TF_r(h(r,t)V(r,t))(k) ) le résultat est enregistré dans h qui est écrasé
	// le facteur 2 vient du fait que vertex_array(r,tau) est proportionnel a U/8 alors que V_array(r,tau) est proportionnel a U/4
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		long int x, y, tmp, wn;
		fr_tmp1=new double[(nbk-2)*nbk];
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			wn=j;
			if (wn>nbw/2) wn=nbw-wn;
			
			for (l=0; l<nbk-2; l++)
			{
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					
					fr_tmp1[m+l*nbk]=2.0*h[m+l*nbk+j*nk4th]*V_array[y+(x*(x+1))/2+wn*nk8th];
				}
			}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];					
//					hV[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		delete [] fr_tmp1;
	}
	
	if (contains_NaN(h, array_size_h)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans h \n";
	fflush(stdout);
	
//	cout<<"8\n";
//	fflush(stdout);

	// TF sur l'espace de f_0(r,t)V(r,t)
	// Attention! le resultat est enregistre dans f0 et le contenu de f0 est ecrase!
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		long int x, y, tmp, wn;
		fr_tmp1=new double[(nbk-2)*nbk];
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			wn=j;
			if (wn>nbw/2) wn=nbw-wn;
			
			for (l=0; l<nbk-2; l++)
			{
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					
					fr_tmp1[m+l*nbk]=2.0*f0[m+l*nbk+j*nk4th]*V_array[y+(x*(x+1))/2+wn*nk8th];
				}
			}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					f0[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		delete [] fr_tmp1;
	}
	
	if (contains_NaN(f0, array_size_h)) cout<<"calc_chijj_vertex_corr_all_optim(): NaN dans f0 \n";
	fflush(stdout);
	
	free_vertex_array();
	
//	cout<<"9\n";
//	fflush(stdout);
	
	dcomplex *Gtmp=new dcomplex[nbw];
	fftw_plan fftplan=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (Gtmp), reinterpret_cast<fftw_complex*> (Gtmp), FFTW_FORWARD, FFTW_MEASURE);
	delete [] Gtmp;

	Gtmp=new dcomplex[nbw];
	fftw_plan fftplan_t_back=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*>(Gtmp), reinterpret_cast<fftw_complex*>(Gtmp), FFTW_BACKWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	if (chijj1_qn) delete [] chijj1_qn;
	if (chijj2_qn) delete [] chijj2_qn;
	
	chijj1_qn=new double[nw+1];
	chijj2_qn=new double[nw+1];
	chijj1_qn4=new double[nw];
	chijj2_qn6=new double[nw];
		
	double *chijj1_k=new double[(nw+1)*N_theta_k];
	double *chijj2_k=new double[(nw+1)*N_theta_k];
	
	for (j=0; j<(nw+1)*N_theta_k; j++)
	{
		chijj1_k[j]=0;
		chijj2_k[j]=0;
	}
	
	for (j=0; j<nw; j++)	
	{
		chijj1_qn[j]=0;
		chijj1_qn4[j]=0;
		chijj2_qn[j]=0;
		chijj2_qn6[j]=0;
	}
	chijj1_qn[nw]=0;
	chijj2_qn[nw]=0;
	
	chijj1_w0=0;
	chijj1_inf2=0;
	chijj2_inf2=0;
	chijj2_inf4=0;
	
	//	double Sigma_inf=G1->Sigma_inf;
	
	int NS0=nbw/2+1;
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex k0sum=0;
	
	//	l0=k0[0], m0=k0[1];
	
	double Dtheta=PI/(2.0*N_theta_k);
	
//	cout<<"10\n";
//	fflush(stdout);
		
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp, p;
		double wt, t, kn, G0t, Jinf;
		dcomplex k0, ph;
		
		//		dcomplex *J0tmp2=new dcomplex[nbw];
		
		dcomplex Jinf_tmp, Jbeta;
		dcomplex *J0tmp=new dcomplex[nbw];
		dcomplex *Jtmp=new dcomplex[nbw];
		double CJ1, CJ2, CJ3, d2J0, d2Jbeta, CJ01, CJ02, CJ03;
		
		dcomplex *Gtmp1=new dcomplex[nbw];
		dcomplex *GGtmp1=new dcomplex[nbw];
		dcomplex Gtmp2, Gtmp0, Gbeta;
		
		double *coeffs=new double[4*(NS0-1)];
		double *tau=new double[NS0];
		double *ft_tmp=new double[NS0];
		double wn, integ, a, b, c, d, x1, x2;
		
		double *tau2=new double[nbw+1];
		double *ft_tmp2=new double[nbw+1];
		double *coeffs2=new double[4*nbw];
		
		dcomplex k0sum_loc=0, chijj2_qn0_loc;
		double *chijj1_qn_loc=new double[nw+1];
		double *chijj1_qn4_loc=new double[nw];
		double *chijj2_qn_loc=new double[nw+1];
		double *chijj2_qn6_loc=new double[nw];
		for (j=0; j<nw; j++)	
		{
			chijj1_qn_loc[j]=0;
			chijj1_qn4_loc[j]=0;
			chijj2_qn_loc[j]=0;
			chijj2_qn6_loc[j]=0;
		}
		chijj1_qn_loc[nw]=0;
		chijj2_qn_loc[nw]=0;
				
		double chijj2_inf2_loc=0, chijj2_inf4_loc=0, chijj1_inf2_loc=0, chijj1_w0_loc=0;
		
		double dtau_J0, dtau_Jbeta;
		double dtau_Gkp, dtau_Gkm, dtau_chi;
		double C2, C3;
		dcomplex selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex exp1, exp2, exp3, exp4;
		dcomplex f1p, f2p, f3p, f4p, f1m, f2m, f3m, f4m;
		
		dcomplex exptmp;
		dcomplex b1, b2, b3, b4;
		dcomplex denj1, denj2, denj3, denj4;
		dcomplex expj1, expj2, expj3, expj4;
		dcomplex fj1, fj2, fj3, fj4, fj1b, fj2b, fj3b, fj4b;
		
		double theta_k;
		int ind_theta;
		double *chijj1_k_loc=new double[(nw+1)*N_theta_k];
		double *chijj2_k_loc=new double[(nw+1)*N_theta_k];
		
		for (j=0; j<(nw+1)*N_theta_k; j++)
		{
			chijj1_k_loc[j]=0;
			chijj2_k_loc[j]=0;
		}
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
				
				// calcul du spline de hV(k,t)
				
				for (j=0; j<=nbw; j++)
				{
					tau2[j]=j/(nbw*tem);
					ft_tmp2[j]=h[m+l*nbk+j*nk4th];
//					ft_tmp2[j]=hV[m+l*nbk+j*nk4th];
				}
				
				coeffs2[0]=dhV[m+l*nbk];
				coeffs2[1]=dhV[m+l*nbk+nk4th];				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);				
				
				
				// TF sur t de hV(k,t)		
				for (j=0; j<nbw; j++)
					Gtmp1[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs2[4*j];
				
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));

				CJ1=-h[m+l*nbk]-h[m+l*nbk+nbw*nk4th];				
//				CJ1=-hV[m+l*nbk]-hV[m+l*nbk+nbw*nk4th];
				CJ2=dhV[m+l*nbk]+dhV[m+l*nbk+nk4th];
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
				CJ3=-d2J0-d2Jbeta;
				
				for (j=0; j<nbw/2; j++)
				{
					kn=(2*j+1)*PI*tem;
					Jtmp[j+nbw/2]=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				for (j=nbw/2; j<nbw; j++)
				{
					kn=(2*(j-nbw)+1)*PI*tem;
					Jtmp[j-nbw/2]=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				
				// calcul du spline de f0V(k,t)
				
				for (j=0; j<=nbw; j++)
				{
					tau2[j]=j/(nbw*tem);
					ft_tmp2[j]=f0[m+l*nbk+j*nk4th];
				}
				
				coeffs2[0]=df0V[m+l*nbk];
				coeffs2[1]=df0V[m+l*nbk+nk4th];				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);
				
				// TF sur t de f0V(k,t)		
				for (j=0; j<nbw; j++)
					Gtmp1[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs2[4*j];
				
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				CJ01=-f0[m+l*nbk]-f0[m+l*nbk+nbw*nk4th];
				CJ02=df0V[m+l*nbk]+df0V[m+l*nbk+nk4th];
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
				CJ03=-d2J0-d2Jbeta;
				
				for (j=0; j<nbw/2; j++)
				{
					kn=(2*j+1)*PI*tem;
					J0tmp[j+nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				for (j=nbw/2; j<nbw; j++)
				{
					kn=(2*(j-nbw)+1)*PI*tem;
					J0tmp[j-nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				
				k0=0;
				dtau_Gkp=0;
				chijj2_qn0_loc=0;
				dtau_J0=0;
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					Jtmp[j]=Gtmp1[j]*Jtmp[j];
					
					chijj2_qn0_loc+=Gtmp1[j]*Gtmp1[j]*J0tmp[j]/tem;
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp2=tem/dcomplex(-ek[m+(l+1)*nbk]+mu, kn+Sigma_inf/kn);
				//	Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
				
					Jinf_tmp=-I*CJ1/kn;
				//	Jinf_tmp=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn);
					//					J0tmp2[j]=Jtmp[j]/tem;
					Jtmp[j]=Jtmp[j]-Gtmp2*Jinf_tmp;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
					dtau_J0-=kn*imag(Jtmp[j]);
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					Jtmp[j]=Gtmp1[j]*Jtmp[j];
					
					chijj2_qn0_loc+=Gtmp1[j]*Gtmp1[j]*J0tmp[j]/tem;
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
				
					Gtmp2=tem/dcomplex(-ek[m+(l+1)*nbk]+mu, kn+Sigma_inf/kn);
				//	Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
				
					Jinf_tmp=-I*CJ1/kn;
				//	Jinf_tmp=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn);
					//					J0tmp2[j]=Jtmp[j]/tem;
					Jtmp[j]=Jtmp[j]-Gtmp2*Jinf_tmp;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);	
					dtau_J0-=kn*imag(Jtmp[j]);
				}
				
				// calcul du G(k,t) asymptotique et du J(k,t) asymtotique
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				find_roots_2pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				
				a1=r1/(r1-r2);
				a2=r2/(r2-r1);
				
				
				b1=CJ1/(r1-r2);
				b2=-b1;
				
				if (real(r1)>0)
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
					fj1=r1*b1*exptmp/den1;
					fj1b=r1*b1/den1;
				}
				else
				{
					exptmp=exp(r1/tem);
					den1=exptmp+1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
					fj1=r1*b1/den1;
					fj1b=r1*b1*exptmp/den1;
				}
				if (real(r2)>0)
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
					fj2=r2*b2*exptmp/den2;
					fj2b=r2*b2/den2;
				}
				else
				{
					exptmp=exp(r2/tem);
					den2=exptmp+1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
					fj2=r2*b2/den2;
					fj2b=r2*b2*exptmp/den2;
				}
				
				k0+=real(f1m+f2m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p);
				dtau_Gkm+=real(r1*f1m+r2*f2m);
				
				dtau_Jbeta=-dtau_J0;
				dtau_J0+=real(fj1+fj2);
				dtau_Jbeta+=real(fj1b+fj2b);
				
			/*
				// calcul du G(k,t) asymptotique et du J(k,t) asymtotique	
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				b1=(CJ1*r1*r1+CJ2*r1+CJ3)/((r1-r2)*(r1-r3)*(r1-r4));
				b2=(CJ1*r2*r2+CJ2*r2+CJ3)/((r2-r1)*(r2-r3)*(r2-r4));
				b3=(CJ1*r3*r3+CJ2*r3+CJ3)/((r3-r1)*(r3-r2)*(r3-r4));
				b4=(CJ1*r4*r4+CJ2*r4+CJ3)/((r4-r1)*(r4-r2)*(r4-r3));				
				
				if (real(r1)>0) 
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
					fj1=r1*b1*exptmp/den1;
					fj1b=r1*b1/den1;
				}
				else 
				{
					exptmp=exp(r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
					fj1=r1*b1/den1;
					fj1b=r1*b1*exptmp/den1;
				}
				if (real(r2)>0) 
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
					fj2=r2*b2*exptmp/den2;
					fj2b=r2*b2/den2;
				}
				else 
				{	
					exptmp=exp(r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
					fj2=r2*b2/den2;
					fj2b=r2*b2*exptmp/den2;
				}
				if (real(r3)>0) 
				{
					exptmp=exp(-r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3/den3;
					f3m=a3*exptmp/den3;
					fj3=r3*b3*exptmp/den3;
					fj3b=r3*b3/den3;
				}
				else 
				{
					exptmp=exp(r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3*exptmp/den3;
					f3m=a3/den3;
					fj3=r3*b3/den3;
					fj3b=r3*b3*exptmp/den3;
				}
				if (real(r4)>0) 
				{
					exptmp=exp(-r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4/den4;
					f4m=a4*exptmp/den4;
					fj4=r4*b4*exptmp/den4;
					fj4b=r4*b4/den4;
				}
				else 
				{
					exptmp=exp(r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4*exptmp/den4;
					f4m=a4/den4;
					fj4=r4*b4/den4;
					fj4b=r4*b4*exptmp/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p+r3*f3p+r4*f4p);
				dtau_Gkm+=real(r1*f1m+r2*f2m+r3*f3m+r4*f4m);		
				
				dtau_Jbeta=-dtau_J0;
				dtau_J0+=real(fj1+fj2+fj3+fj4);
				dtau_Jbeta+=real(fj1b+fj2b+fj3b+fj4b);
				*/
				
				//TF de G(k,ikn) dans le temps
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				Gbeta=Gtmp1[0];
				
				// TF de J(k,ikn) dans le temps
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Jtmp),reinterpret_cast<fftw_complex*>(Jtmp));
				Jbeta=Jtmp[0];
				
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					
					if (real(r1)>0) 
					{
						exp1=exp(-t*r1);
						expj1=exp(r1*(t-1.0/tem));
					}
					else 
					{
						exp1=exp((1.0/tem-t)*r1);
						expj1=exp(r1*t);
					}
					if (real(r2)>0) 
					{
						exp2=exp(-t*r2);
						expj2=exp(r2*(t-1.0/tem));
					}
					else 
					{
						exp2=exp((1.0/tem-t)*r2);
						expj2=exp(r2*t);
					}
				/*
					if (real(r3)>0) 
					{
						exp3=exp(-t*r3);
						expj3=exp(r3*(t-1.0/tem));
					}
					else 
					{
						exp3=exp((1.0/tem-t)*r3);
						expj3=exp(r3*t);
					}
					if (real(r4)>0) 
					{
						exp4=exp(-t*r4);
						expj4=exp(r4*(t-1.0/tem));
					}
					else 
					{
						exp4=exp((1.0/tem-t)*r4);
						expj4=exp(r4*t);
					}
				*/
					
				//	G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
				//	Jinf=real(b1*expj1/den1+b2*expj2/den2+b3*expj3/den3+b4*expj4/den4);
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2);
					Jinf=real(b1*expj1/den1+b2*expj2/den2);
					
					ph=exp(dcomplex(0.0, (j*PI*(nbw-1))/nbw));
					
					Gtmp1[j]=ph*Gtmp1[j]+G0t;					
					
					Jtmp[j]=Jtmp[j]/ph+Jinf;
				}			
				
				t=1.0/tem;
				
				if (real(r1)>0) 
				{
					exp1=exp(-t*r1);
					expj1=1.0;
				}
				else 
				{
					exp1=1.0;
					expj1=exp(r1*t);
				}
				if (real(r2)>0) 
				{
					exp2=exp(-t*r2);
					expj2=1.0;
				}
				else 
				{
					exp2=1.0;
					expj2=exp(r2*t);
				}
				/*
				if (real(r3)>0) 
				{
					exp3=exp(-t*r3);
					expj3=1.0;
				}
				else 
				{
					exp3=1.0;
					expj3=exp(r3*t);
				}
				if (real(r4)>0) 
				{
					exp4=exp(-t*r4);
					expj4=1.0;
				}
				else 
				{
					exp4=1.0;
					expj4=exp(r4*t);
				}
				*/
				
			//	G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
			//	Jinf=real(b1*expj1/den1+b2*expj2/den2+b3*expj3/den3+b4*expj4/den4);
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2);
				Jinf=real(b1*expj1/den1+b2*expj2/den2);
				
				ph=exp(dcomplex(0.0, PI*(nbw-1)));
				
				Gbeta=ph*Gbeta+G0t;					
				
				Jbeta=Jbeta/ph+Jinf;
				
				// calcul du spline et de la TF sur t de G(k,t)J(k,-t) 				
				for (j=0; j<nbw; j++)
				{
					tau2[j]=j/(nbw*tem);					
					ft_tmp2[j]=real(Gtmp1[j]*Jtmp[j]);
				}
				tau2[nbw]=1.0/tem;
				ft_tmp2[nbw]=real((-Gtmp1[0]-(double)1.0)*Jbeta);
				
				coeffs2[0]=real(dtau_Gkp*Jtmp[0]+Gtmp1[0]*dtau_J0);
				coeffs2[1]=real(dtau_Gkm*Jbeta - (Gtmp1[0]+(double)1.0)*dtau_Jbeta);
				
				CJ1=ft_tmp2[0]-ft_tmp2[nbw];
				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);
				//				spline_coeffs(tau2, ft_tmp2, nbw+1, coeffs2);
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
				//				d2Jbeta=6.0*coeffs2[4*nbw-4]/tem+2.0*coeffs2[4*nbw-3];
				CJ3=d2J0-d2Jbeta;
				
				for (j=0; j<nbw; j++)
				{
					J0tmp[j]=6.0*coeffs2[4*j];
				}
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(J0tmp),reinterpret_cast<fftw_complex*>(J0tmp));
				
				for (j=1; j<=nbw/2; j++)
				{
					wn=2*j*PI*tem;
					Jtmp[j]=-((double)2.0)*I*CJ1/wn+((double)2.0)*I*CJ3/(wn*wn*wn)+((double)2.0)*I*imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*J0tmp[j]/(wn*wn*wn*wn));
				}
				
				// Definition de -G(k,t)G(k,-t)				
				GGtmp1[0]=Gtmp1[0]*(Gtmp1[0]+(double)1.0);
				for (j=1; j<nbw; j++)
					GGtmp1[j]=-Gtmp1[j]*Gtmp1[nbw-j];
				
				// calcul du spline de -G(k,t)G(k,-t)			
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nbw*tem);
					ft_tmp[j]=real(GGtmp1[j]);
				}
				dtau_chi=dtau_Gkp*(real(Gtmp1[0])+1.0)+real(Gtmp1[0])*dtau_Gkm;
				
				//				if ((l+1)%4==0)
				//					cout<<"dtau_chi:  "<<setw(20)<<(real(GGtmp1[1]-GGtmp1[0])*(tem*nbw)-dtau_chitmp)/dtau_chitmp<<setw(20)<<(dtau_chi-dtau_chitmp)/dtau_chitmp<<'\n';
				
				coeffs[0]=dtau_chi;
				coeffs[1]=0;				
				spline_coeffs_rel(tau, ft_tmp, NS0, coeffs);				
				//				spline_coeffs(tau, ft_tmp, NS0, coeffs);
				
				// TF sur t de -G(k,t)G(k,-t)			
				for (j=0; j<nbw/2; j++)
					Gtmp1[j]=6.0*coeffs[4*j];
				for (j=nbw/2; j<nbw; j++)
					Gtmp1[j]=-Gtmp1[nbw-j-1];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				integ=0;
				for (j=0; j<nbw/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j+1]-tau[j];
					integ+=a*x1*x1*x1*x1/4.0+b*x1*x1*x1/3.0+c*x1*x1/2.0+d*x1;
					
					//					x1=tau[j];
					//					x2=tau[j+1];
					//					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
				
				Gtmp1[0]=2*integ;
				
				wt=4;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+(l+1)*nbk]*k0;
				chijj1_qn_loc[0]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_w0_loc-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_inf2_loc+=4.0*dtau_chi*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk];
				chijj2_inf2_loc+=normFact*wt*dek[m+l*nbk]*2.0*CJ1;
				chijj2_inf4_loc-=normFact*wt*dek[m+l*nbk]*2.0*CJ3;
				
				chijj2_qn_loc[0]-=normFact*wt*dek[m+l*nbk]*real(chijj2_qn0_loc);
				
				theta_k=atan( (1.0*(nk0/2-m))/(nk0/2-l-1) );
				ind_theta=(int)floor(theta_k/Dtheta);
				
				//				cout<<setiosflags(ios::left)<<setprecision(10);
				//				cout<<setw(10)<<l+1<<setw(10)<<m<<"theta_k:  "<<setw(20)<<theta_k<<"ind_theta:  "<<ind_theta<<endl;
				
				chijj1_k_loc[ind_theta]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj2_k_loc[ind_theta]-=normFact*wt*dek[m+l*nbk]*real(chijj2_qn0_loc);
				
				for (j=1; j<=nw; j++)
				{
					wn=2*j*PI*tem;
					
					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]);
					
					chijj2_qn_loc[j]-=normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
					
					chijj2_k_loc[ind_theta+j*N_theta_k]-=normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
					
					//					if (p<nk_chijj)	chijj2_k[p + j*nk_chijj]=-normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
					
					chijj2_qn6_loc[j-1]-=normFact*wt*dek[m+l*nbk]*2.0*imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*J0tmp[j])/wn;
					
					Gtmp1[j]=-2.0*dtau_chi/(wn*wn)+real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]/(wn*wn*wn*wn));
					
					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
					
					chijj1_k_loc[ind_theta+j*N_theta_k]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
					
					//					if (p<nk_chijj)	chijj1_k[p + j*nk_chijj]=-2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
				}				
				
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
			
			for (j=0; j<nw; j++)
			{
				chijj1_qn[j]+=chijj1_qn_loc[j];
				chijj1_qn4[j]+=chijj1_qn4_loc[j];
				chijj2_qn[j]+=chijj2_qn_loc[j];
				chijj2_qn6[j]+=chijj2_qn6_loc[j];
				for (p=0; p<N_theta_k; p++)
				{
					chijj1_k[p+j*N_theta_k]+=chijj1_k_loc[p+j*N_theta_k];
					chijj2_k[p+j*N_theta_k]+=chijj2_k_loc[p+j*N_theta_k];
				}
			}
			chijj1_qn[nw]+=chijj1_qn_loc[nw];
			chijj2_qn[nw]+=chijj2_qn_loc[nw];
			for (p=0; p<N_theta_k; p++)
			{
				chijj1_k[p+nw*N_theta_k]+=chijj1_k_loc[p+nw*N_theta_k];
				chijj2_k[p+nw*N_theta_k]+=chijj2_k_loc[p+nw*N_theta_k];
			}
			
			chijj1_w0+=chijj1_w0_loc;
			chijj1_inf2+=chijj1_inf2_loc;
			chijj2_inf2+=chijj2_inf2_loc;
			chijj2_inf4+=chijj2_inf4_loc;
		}
		delete [] chijj1_qn4_loc;
		delete [] chijj1_qn_loc;
		delete [] chijj2_qn6_loc;
		delete [] chijj2_qn_loc;
		delete [] Gtmp1;
		delete [] GGtmp1;
		delete [] tau;
		delete [] ft_tmp;
		delete [] coeffs;
		delete [] tau2;
		delete [] ft_tmp2;
		delete [] coeffs2;
		delete [] Jtmp;
		delete [] J0tmp;
		delete [] chijj1_k_loc;
		delete [] chijj2_k_loc;		
	}
	
//	cout<<"11\n";
//	fflush(stdout);
		
//	delete [] h;
	delete [] f0;
	delete [] dhV;
	delete [] df0V;
	
//	delete [] hV;	
//	delete [] df0;	
	
//	delete [] H;
//	delete [] dtauGH;
//	delete [] dtauf0G;
//	delete [] dtauH;
	
	double *chijj1_k_tot=new double[nw+1];
	double *chijj2_k_tot=new double[nw+1];
	
	for (j=0; j<=nw; j++)
	{
		chijj1_k_tot[j]=0;
		chijj2_k_tot[j]=0;
	}
	
	int p;
	
	for (j=0; j<=nw; j++)
		for  (p=0; p<N_theta_k; p++)
		{
			chijj1_k_tot[j]+=chijj1_k[p + j*N_theta_k];
			chijj2_k_tot[j]+=chijj2_k[p + j*N_theta_k];
		}
	
	fstream file, file2, file3, file4, file_info;
	char name[200];
	
	const char *nameForm_chijj1_k="./chijj1_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj1_k_tot="./chijj1_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj2_k="./chijj2_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj2_k_tot="./chijj2_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	
	sprintf(name, nameForm_chijj1_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj1_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file2.open(name, ios::out );
	file2<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj2_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file3.open(name, ios::out );
	file3<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj2_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out );
	file4<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	for (j=0; j<=nw; j++)
	{
		file2<<setw(10)<<j<<setw(30)<<chijj1_k_tot[j]<<setw(30)<<chijj1_qn[j]<<chijj1_k_tot[j]/chijj1_qn[j]<<'\n';
		file4<<setw(10)<<j<<setw(30)<<chijj2_k_tot[j]<<setw(30)<<chijj2_qn[j]<<chijj2_k_tot[j]/chijj2_qn[j]<<'\n';
		for  (p=0; p<N_theta_k; p++)
		{
			file<<setw(10)<<j<<setw(10)<<p+1<<chijj1_k[p + j*N_theta_k]<<'\n';
			file3<<setw(10)<<j<<setw(10)<<p+1<<chijj2_k[p + j*N_theta_k]<<'\n';			
			//			file<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj1_k[p + j*nk_chijj]<<'\n';
			//			file3<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj2_k[p + j*nk_chijj]<<'\n';
		}
	}
	file.close();
	file2.close();
	file3.close();
	file4.close();
	
	delete [] chijj1_k_tot;
	delete [] chijj2_k_tot;
	
	delete [] chijj1_k;
	delete [] chijj2_k;
	
//	cout<<"12\n";
//	fflush(stdout);

/*	
	{
		const char *nameForm_self_k="./self_k_ikn_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
		sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
		file.open(name, ios::out );
		file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
		
		dcomplex selftmp;
		long int x, y, tmp;
		double C2, C3, kn;
		for (j=0; j<=nw/4; j++)
		{
			kn=(2.0*j+1)*PI*tem;
			for (p=0; p<nk_chijj; p++)
			{
				x=k_chijj[2*p];
				y=k_chijj[2*p+1];
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
				
				file<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
			}
		}
		file.close();
	}
*/
	
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp;
		double wt, kn;
		dcomplex k0;
		dcomplex Gtmp2, Gtmp0;
		
		dcomplex k0sum_loc=0;
		
		double C2, C3;
		dcomplex z, selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex f1m, f2m, f3m, f4m;
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk; l+=nk0/2)
		{
			for (m=0; m<nbk; m++)
			{
				x=l;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				//				C2=G1->Sigma_inf2[y+(x*(x+1))/2];
				//				C3=G1->Sigma_inf3[y+(x*(x+1))/2];				
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
				
				k0=0;
				
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
				
					Gtmp0=tem/dcomplex(-ek[m+l*nbk]+mu, kn+Sigma_inf/kn);
				//	Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
				
					Gtmp0=tem/dcomplex(-ek[m+l*nbk]+mu, kn+Sigma_inf/kn);
				//	Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				find_roots_2pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				
				a1=r1/(r1-r2);
				a2=r2/(r2-r1);
				
				if (real(r1)>0)
				{
					den1=exp(-r1/tem)+(double)1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else
				{
					den1=exp(r1/tem)+(double)1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0)
				{
					den2=exp(-r2/tem)+(double)1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else
				{
					den2=exp(r2/tem)+(double)1.0;
					f2m=a2/den2;
				}
				
				k0+=real(f1m+f2m);
				
				/*
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0) 
				{
					den1=exp(-r1/tem)+(double)1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else 
				{
					den1=exp(r1/tem)+(double)1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0) 
				{
					den2=exp(-r2/tem)+(double)1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else 
				{	
					den2=exp(r2/tem)+(double)1.0;
					f2m=a2/den2;
				}
				if (real(r3)>0) 
				{
					den3=exp(-r3/tem)+(double)1.0;
					f3m=a3*exp(-r3/tem)/den3;
				}
				else 
				{
					den3=exp(r3/tem)+(double)1.0;
					f3m=a3/den3;
				}
				if (real(r4)>0) 
				{
					den4=exp(-r4/tem)+(double)1.0;
					f4m=a4*exp(-r4/tem)/den4;
				}
				else 
				{
					den4=exp(r4/tem)+(double)1.0;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);				
				*/
				
				wt=2;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+l*nbk]*k0;
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
		}
	}
	
//	cout<<"13\n";
//	fflush(stdout);
	
	//	cout<<"chitmp3:  ";
	//	calc_chijj_vertex_corr2_wn(0, k0, 0, Htmp, NULL, NULL, NULL, NULL);
	
	delete [] d2ek;
	
	//	for (j=0; j<4; j++)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	//	for (j=nbw/8; j<=nbw/2; j+=nbw/8)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	
	
	kx0=real(k0sum);
	sprintf(name, "kx0_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out);
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific)<<kx0<<'\n';
	file.close();
	
	sprintf(name, "kx0_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&kx0, sizeof(kx0));
	file.close();
	
	for (j=0; j<=nw; j++)
	{
		if (chijj1_qn[j]<0.0)
		{
			cout<<"attention! valeur non-physique de chi_jj\n";
			cout<<setw(30)<<chijj1_qn[j]<<'\n';
			break;
		}
	}
	
	double *chijj_qn=new double[nbw/2+1];
	
	for (j=0; j<=nbw/2; j++)
		chijj_qn[j]=chijj1_qn[j]+chijj2_qn[j];
	
	/*
	//pour enregistrer chijj(iq_n) en binaire
	
	const char *nameForm_chijj1_qn_bin="./chijj1_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	
	const char *nameForm_chijj2_qn_bin="./chijj2_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj2_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out | ios::binary);
	
	
	j=0;
	file.write((char*)&j,sizeof(j));
	file.write((char*)&(chijj1_qn[j]),sizeof(double));
	file4.write((char*)&j,sizeof(j));
	file4.write((char*)&(chijj2_qn[j]),sizeof(double));
	for (j=1; j<=nbw/2; j++)
	{
		file.write((char*)&j,sizeof(j));
		file.write((char*)&(chijj1_qn[j]),sizeof(double));
		file4.write((char*)&j,sizeof(j));
		file4.write((char*)&(chijj2_qn[j]),sizeof(double));
	}
	file.close();
	file4.close();
	// fin de l'enregistrement en binaire
	*/
	
	// pour enregistrer chijj(iq_n) en ascii
	
	const char *nameForm_chijj1_qn="./chijj1_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	const char *nameForm_chijj_vc1_qn="./chijj_with_vc1_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj_vc1_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file2.open(name, ios::out );
	file2<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	const char *nameForm_chijj2_qn="./chijj2_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj2_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out );
	file4<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);	
	
	j=0;
	file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	file4<<setw(10)<<j<<chijj2_qn[j]<<'\n';
	file2<<setw(10)<<j<<chijj_qn[j]<<'\n';
	for (j=1; j<=nbw/2; j++)
	{
		file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
		file4<<setw(10)<<j<<chijj2_qn[j]<<'\n';
		file2<<setw(10)<<j<<chijj_qn[j]<<'\n';
	}
	file.close();
	file2.close();
	file4.close();
	// fin de l'enregistrement en ascii
	
	cout<<"-k0:  "<<-kx0<<'\n';
	cout<<"chijj1(iqn=0):  "<<chijj1_qn[0]<<'\n';
	cout<<"chijj2(iqn=0):  "<<chijj2_qn[0]<<'\n';
	cout<<"chijj1(iqn=0)+chijj2(iqn=0):  "<<chijj1_qn[0]+chijj2_qn[0]<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+chijj2_qn[0]+kx0)/kx0<<'\n';
	
	
	long int nn;
	int N0;
	
	if (chi3)
	{
		cout<<"calcul de la seconde correction de vertex\n";
		
		if (chijj3_qn) delete [] chijj3_qn;
		chijj3_qn=new double[nw+1];
		for (j=0; j<=nw; j++)	chijj3_qn[j]=0;		
		double *chijj3_k=new double[(nw+1)*N_theta_k];
		for (j=0; j<(nw+1)*N_theta_k; j++)	chijj3_k[j]=0;	
		dcomplex chijj3_qn0=0;
		
//		calcul de h(r,t) (voir l'annexe sur les techniques de calcul)
		
#pragma omp parallel private(fr_tmp1,j,l,m)
		{
			fr_tmp1=new double[(nbk-2)*nbk];
			double t;
			
#pragma omp for
			for (j=0; j<=nbw; j++)
			{
				t=j/(tem*nbw);
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						if (ek[m+(l+1)*nbk]>mu0)
							fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
						else
							fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp( (ek[m+(l+1)*nbk]-mu0)*(1.0/tem-t))/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
			}
			delete [] fr_tmp1;
		}
		
/*		
		// calcul de la derivee de h(r,t) a t=0 et t=beta
		
		double *dh=new double[(nbk-2)*nbk*2];
		
		// t=0	
		fr_tmp1=new double[(nbk-2)*nbk];
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				if (ek[m+(l+1)*nbk]>mu0)
					fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
				else
					fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (ek[m+(l+1)*nbk]-mu0)/tem)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				dh[m+l*nbk]=fr_tmp1[m+l*nbk];
		
		// t=beta	
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				if (ek[m+(l+1)*nbk]>mu0)
					fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (mu0-ek[m+(l+1)*nbk])/tem )/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
				else
					fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				dh[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];	
		
		delete [] fr_tmp1;
*/		
		fr_tmp=new double[nbk*nbk];
		fftw_plan fftplan_RE=fftw_plan_r2r_2d(nbk, nbk, fr_tmp, fr_tmp, FFTW_REDFT00, FFTW_REDFT00, FFTW_MEASURE);
		delete [] fr_tmp;
		
		long int array_size_G=(long int)nk8th*(nbw+1);
		
		double *Grt=new double[array_size_G];
		//	double *Grt=new double[nk8th*(nbw+1)];
		double *dtauGrt=new double[nk8th*2];
		
		// calcul de G(r,-t)
#pragma omp parallel private(fr_tmp1,j,l,m)
		{
			fr_tmp1=new double[nbk*nbk];
			double t;
			
#pragma omp for
			for (j=0; j<=nbw; j++)
			{
				t=j/(tem*nbw);
				for (l=0; l<nbk; l++)
					for (m=0; m<=l; m++)
					{
						if (ek[m+l*nbk]<mu0)
							fr_tmp1[m+l*nbk]=normFact*exp((ek[m+l*nbk]-mu0)*t)/( exp((ek[m+l*nbk]-mu0)/tem) + 1.0 );
						else
							fr_tmp1[m+l*nbk]=normFact*exp((ek[m+l*nbk]-mu0)*(t-1.0/tem))/(exp( -(ek[m+l*nbk]-mu0)/tem )  +1.0 );
						fr_tmp1[l+m*nbk]=fr_tmp1[m+l*nbk];
						
						//					Gkt_tmp[m+(l*(l+1))/2+j*nk8th]=fr_tmp1[m+l*nbk]/normFact;
					}
				
				fftw_execute_r2r(fftplan_RE,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk; l++)
					for (m=0; m<=l; m++)
						Grt[m+(l*(l+1))/2+j*nk8th]=fr_tmp1[m+l*nbk];
			}
			delete [] fr_tmp1;
		}
		
		fr_tmp1=new double[nbk*nbk];
		
		// calcul de d/dt G(r,-t) a t=0
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				if (ek[m+l*nbk]<mu0)
					fr_tmp1[m+l*nbk]=normFact*(ek[m+l*nbk]-mu0)/( exp((ek[m+l*nbk]-mu0)/tem) + 1.0 );
				else
					fr_tmp1[m+l*nbk]=normFact*(ek[m+l*nbk]-mu0)*exp(-(ek[m+l*nbk]-mu0)/tem)/(exp( -(ek[m+l*nbk]-mu0)/tem )  +1.0 );
				fr_tmp1[l+m*nbk]=fr_tmp1[m+l*nbk];
			}
		
		fftw_execute_r2r(fftplan_RE,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				dtauGrt[m+(l*(l+1))/2]=fr_tmp1[m+l*nbk];
				//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauGrt[m+(l*(l+1))/2]<<tem*nbw*(Grt[m+(l*(l+1))/2+nk8th]-Grt[m+(l*(l+1))/2])<<endl;
				dtauGrt[m+(l*(l+1))/2+nk8th]=0;
			}
		// 	d/dt G(r,-t) a t=beta
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
				dtauGrt[m+(l*(l+1))/2+nk8th]=er[m+(l*(l+1))/2];
		dtauGrt[nk8th]=-mu0;
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				dtauGrt[m+(l*(l+1))/2+nk8th]-=dtauGrt[m+(l*(l+1))/2];
				//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauGrt[m+(l*(l+1))/2+nk8th]<<tem*nbw*(Grt[m+(l*(l+1))/2+nbw*nk8th]-Grt[m+(l*(l+1))/2+(nbw-1)*nk8th])<<endl;
			}
		
		delete [] fr_tmp1;
		
		fftw_destroy_plan(fftplan_RE);
		

		if (calc_chijj2_qn0)
		{	
			f0=new double[array_size_h];	
			//	double *f0=new double[(nbk-2)*nbk*(nbw+1)];
			
			// calcul de f_0(r,t) (voir l'annexe sur les techniques de calcul)
		#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				double t, expk;
				
			#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					t=j/(tem*nbw);
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							if (ek[m+(l+1)*nbk]>mu0)
							{
								expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
								fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)*( t/(expk + 1.0) -expk/(tem*(expk+1.0)*(expk+1.0)));
							}
							else
							{
								expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
								fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((ek[m+(l+1)*nbk]-mu0)*(1.0/tem-t))*(t/(expk+1.0) - 1.0/(tem*(expk+1.0)*(expk+1.0)));
							}
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							f0[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
				}
				delete [] fr_tmp1;
			}
			
/*				
			// calcul de la derivee de  f0(r,t) a t=0 et t=beta
			double *df0=new double[(nbk-2)*nbk*2];

			fr_tmp1=new double[(nbk-2)*nbk];

			//	t=0
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
					{
						expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*(  (ek[m+(l+1)*nbk]-mu0)*expk/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk + 1.0));
					}
					else
					{
						expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*( (ek[m+(l+1)*nbk]-mu0)/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk+1.0));
					}
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					df0[m+l*nbk]=fr_tmp1[m+l*nbk];
				}
			
			// t=beta	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
					{
						expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*(  (mu0-ek[m+(l+1)*nbk])*( 1.0/(tem*(expk + 1.0)) -expk/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk + 1.0));
					}
					else
					{
						expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*( (mu0-ek[m+(l+1)*nbk])*(1.0/(tem*(expk+1.0)) - 1.0/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk+1.0));
					}
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					df0[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
				}

			delete [] fr_tmp1;
*/
			
		//calcul de H_0(q,iqn)
		
			double *H=new double[array_size_h];	
			//	double *H=new double[nbk*(nbk-2)*(nbw+1)];
			double *dtauf0G=new double[nbk*(nbk-2)*2];
			
			
			//TF dans l'espace de f_0(r,t)G(r,-t)
			#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
					
			#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=f0[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
							
							//					f0G_tmp[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							H[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
				}
				delete [] fr_tmp1;
			}
				
			fr_tmp1=new double[(nbk-2)*nbk];
			
			//TF dans l'espace de d/dt f_0(r,t)G(r,-t) a t=0 et t=beta
			//t=0
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=df0[m+l*nbk]*Grt[y+(x*(x+1))/2]+f0[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauf0G[m+l*nbk]=fr_tmp1[m+l*nbk];
					//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauf0G[m+l*nbk]<<tem*nbw*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<endl;
				}
			
			//  d/dt f_0(r,t)G(r,-t)  a t=beta	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=df0[m+l*nbk+nk4th]*Grt[y+(x*(x+1))/2+nbw*nk8th]+f0[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauf0G[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
					//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauf0G[m+l*nbk+nk4th]<<tem*nbw*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
				}
			
			delete [] fr_tmp1;
			
			delete [] f0;
			
			double *fg1=new double[nbk*(nbk-2)];
			double *fg3=new double[nbk*(nbk-2)];
			
			N0=nbw+1;
			
			// TF dans le temps de (f0*G)(k,tau) et multiplication par le prefacteur fonction de 1.0/(1.0-Usp*chi0/2) et 1.0/(1.0+Uch*chi0/2)
			// (on prend seulement la partie imaginaire, voir la section "Techniques de calcul")
//			#pragma omp parallel private(l,m,j)
			{
				double chitmp, chisp_tmp, chich_tmp;
				long int x, y, tmp, wn;
				dcomplex *ft_tmp1=new dcomplex[nbw];
				double *tau=new double[nbw+1];
				double *ftau=new double[nbw+1];
				double *coeffs=new double[4*nbw];
				double qn, d2S0, d2Sbeta, Htau0, Htau_beta;
				
//			#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=H[m+l*nbk+j*nk4th];
							//					if (l==2 && m==1 ) cout<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
						}
						
						coeffs[0]=dtauf0G[m+l*nbk];
						coeffs[1]=dtauf0G[m+l*nbk+nk4th];				
						
						//				spline_coeffs(tau, ftau, N0, coeffs);
						spline_coeffs_rel(tau, ftau, N0, coeffs);
						
						for (j=0; j<nbw; j++)
							ft_tmp1[j]=6.0*coeffs[4*j];
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						d2S0=2*coeffs[1];
						d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];				
						//				d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
						
						fg1[m+l*nbk]=H[m+l*nbk]-H[m+l*nbk+nbw*nk4th];
						fg3[m+l*nbk]=d2S0-d2Sbeta;
						
						//				cout<<setw(10)<<l<<setw(10)<<m<<setw(20)<<H[m+l*nbk]<<setw(20)<<H[m+l*nbk+nbw*nk4th]<<setw(20)<<fg1[m+l*nbk]<<setw(20)<<d2S0<<setw(20)<<d2Sbeta<<setw(20)<<fg3[m+l*nbk]<<endl;
						//				cout<<setw(10)<<l<<setw(10)<<m<<fg1[m+l*nbk]<<endl;
						
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						Htau0=H[m+l*nbk];
						Htau_beta=H[m+l*nbk+nbw*nk4th];
						for (j=1; j<nbw/2; j++)
						{
							
							wn=j;
							if (wn>nbw/2) wn=nbw-wn;
							chitmp=chiqw_array[y+(x*(x+1))/2+wn*nk8th];
							
							chisp_tmp=1.0/(1.0-Usp*chitmp/2);
							chich_tmp=1.0/(1.0+Uch*chitmp/2);
							chitmp=3*Usp*chisp_tmp*chisp_tmp+Uch*chich_tmp*chich_tmp;
							
							qn=2*j*PI*tem;
							H[m+l*nbk+j*nk4th]=-2*chitmp*(-(Htau0-Htau_beta)/qn + (d2S0-d2Sbeta)/(qn*qn*qn) + imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j])/(qn*qn*qn*qn));					
							//					Htmp[m+l*nbk+j*nk4th]=-I*(H[m+l*nbk]-H[m+l*nbk+nbw*nk4th])/qn - (dtauf0G[m+l*nbk]-dtauf0G[m+l*nbk+nk4th])/(qn*qn)+ I*(d2S0-d2Sbeta)/(qn*qn*qn)+ ((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]/(qn*qn*qn*qn);					
							//					if (l==nbk-3 && m==nbk-1 && j<20) cout<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<endl;
							//					if (j<10) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<-2*imag(Htmp[m+l*nbk+j*nk4th])<<endl;
							//					cout<<setw(10)<<l+1<<setw(10)<<m<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
						}
						j=nbw/2;
						H[m+l*nbk+j*nk4th]=0;
						for (j=nbw/2+1; j<nbw; j++)
						{
							H[m+l*nbk+j*nk4th]=-H[m+l*nbk+(nbw-j)*nk4th];
							//					qn=2*(j-nbw)*PI*tem;
							//					Htmp[m+l*nbk+j*nk4th]=-I*(H[m+l*nbk]-H[m+l*nbk+nbw*nk4th])/qn - (dtauf0G[m+l*nbk]-dtauf0G[m+l*nbk+nk4th])/(qn*qn)+ I*(d2S0-d2Sbeta)/(qn*qn*qn)+ ((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]/(qn*qn*qn*qn);
						}
						H[m+l*nbk]=0;
					}
				delete [] ft_tmp1;
				delete [] tau;
				delete [] ftau;
				delete [] coeffs;
			}
			
			delete [] dtauf0G;
			
			double *dtauH=new double[nbk*(nbk-2)];
			
			// TF dans le temps de H0(q,iq_m) et calcul de la derivee de H0(q,tau) a tau=0 et tau=beta
			//#pragma omp parallel private(l,m,j)
			{
				dcomplex *ft_tmp1=new dcomplex[nbw];
				long int wn, x, y, tmp;
				double t, qn, tmpsp, tmpch;
				dcomplex tmpnum;
				double zsp, zch;
				double cosbt, cosb, cost, sint, sinb, sinbt;
				double et, eb, ebtp, ebtm;
				double H0inf, dtauHtmp, dtauHinf;
				dcomplex Hinfqn;
				double chitmp, chisp_tmp, chich_tmp;
				
				//#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}	
						
						dtauHtmp=0;
						ft_tmp1[0]=0;
						for (j=1; j<nbw/2; j++)
						{					
							qn=2*PI*j*tem;
							tmpsp=(-qn*qn+Usp*chi0_inf[y+(x*(x+1))/2]/2);
							tmpch=(-qn*qn-Uch*chi0_inf[y+(x*(x+1))/2]/2);
							tmpnum=-I*qn*qn*qn*fg1[m+l*nbk]+I*qn*fg3[m+l*nbk];
							Hinfqn=-6.0*Usp*tmpnum/(tmpsp*tmpsp)-2.0*Uch*tmpnum/(tmpch*tmpch);
							ft_tmp1[j]=tem*( I*H[m+l*nbk+j*nk4th]-Hinfqn);
							dtauHtmp+=2*qn*imag(ft_tmp1[j]);
							//					dtauHtmp+=2*qn*tem*H[m+l*nbk+j*nk4th];
							
							//					cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;					
							//					if (j>(nbw/2-20)) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;
							//					if (l==nbk-3 && m==nbk-1) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;
							//					Htmp[m+l*nbk+j*nk4th]=imag(Hinfqn);
							//					Htmp[m+l*nbk+j*nk4th]=H[m+l*nbk+j*nk4th];
						}
						ft_tmp1[nbw/2]=0;
						for (j=nbw/2+1; j<nbw; j++)
							ft_tmp1[j]=-ft_tmp1[nbw-j];
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						zsp=sqrt(Usp*chi0_inf[y+(x*(x+1))/2]/2);
						zch=sqrt(Uch*chi0_inf[y+(x*(x+1))/2]/2);
						
						cosb=cos(zsp/(2.0*tem));
						sinb=sin(zsp/(2.0*tem));
						eb=exp(-zch/tem);
						for (j=0; j<nbw; j++)
						{
							t=j/(nbw*tem);
							sint=sin(zsp*t);
							sinbt=sin((1.0/(2.0*tem)-t)*zsp);
							cosbt=cos((1.0/(2.0*tem)-t)*zsp);
							et=exp(-t*zch);
							ebtp=exp(-(1.0/tem+t)*zch);
							ebtm=exp(-(1.0/tem-t)*zch);
							tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -t*cosbt/sinb +sint/(2.0*tem*sinb*sinb) )/(4.0*zsp)
							-(fg1[m+l*nbk]/2)*sinbt/sinb;
							tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( t*(ebtm+et)/(1.0-eb) + (ebtp-ebtm)/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
							+(fg1[m+l*nbk]/2)*(ebtm-et)/(1.0-eb);
							H0inf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
							H[m+l*nbk+j*nk4th]=real(ft_tmp1[j])+H0inf;	
							//					if (l==nbk-3 && m==nbk-1 && j%(nbw/16)==0) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<endl;
							//					H[m+l*nbk+j*nk4th]=H0inf;
						}
						j=nbw;
						t=1.0/tem;
						sint=sin(zsp*t);
						sinbt=sin((1.0/(2.0*tem)-t)*zsp);
						cosbt=cos((1.0/(2.0*tem)-t)*zsp);
						et=exp(-t*zch);
						ebtp=exp(-(1.0/tem+t)*zch);
						ebtm=exp(-(1.0/tem-t)*zch);
						tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -t*cosbt/sinb + sint/(2.0*tem*sinb*sinb) )/(4.0*zsp)
						-(fg1[m+l*nbk]/2)*sinbt/sinb;
						tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( t*(ebtm+et)/(1.0-eb) + (ebtp-ebtm)/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
						+(fg1[m+l*nbk]/2)*(ebtm-et)/(1.0-eb);
						H0inf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
						H[m+l*nbk+j*nk4th]=H0inf;
						
						// d/dtau H(q,tau) a tau=0 (les derivees sont egales a tau=0 et tau=beta);
						tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -cosb/sinb+zsp/(2.0*tem*sinb*sinb) )/(4.0*zsp) 
						+ zsp*(fg1[m+l*nbk]/2)*cosb/sinb;
						tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( (eb+1.0)/(1.0-eb)-2.0*zch*eb/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
						+ zch*(fg1[m+l*nbk]/2)*(eb+1.0)/(1.0-eb);
						dtauHinf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
						dtauH[m+l*nbk]=dtauHtmp+dtauHinf;
						
						//				dtauH[m+l*nbk]=dtauHinf;
						//				cout<<setw(25)<<dtauHtmp<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
						//				cout<<setw(25)<<dtauH[m+l*nbk]<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
						//				cout<<setw(25)<<dtauHinf<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;				
					}
			}
			
			
			delete [] fg1;
			delete [] fg3;
			
			//TF dans l'espace de H0	
			#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				
				//TF dans l'espace de H_0(q,tau)		
			#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							fr_tmp1[m+l*nbk]=normFact*H[m+l*nbk+j*nk4th];
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							H[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
				}
				
				delete [] fr_tmp1;
			}
			
			fr_tmp1=new double[(nbk-2)*nbk];
			
			//TF dans l'espace de d/dtau H_0(q,tau) a tau=0		
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					fr_tmp1[m+l*nbk]=normFact*dtauH[m+l*nbk];
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					dtauH[m+l*nbk]=fr_tmp1[m+l*nbk];
				
			delete [] fr_tmp1;
				
			//TF dans l'espace de H_0(r,tau)G(r,-tau)
			#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
			#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=H[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							H[m+l*nbk+j*nk4th]=-fr_tmp1[m+l*nbk];
						}
					//			l=2;
					//			m=1;
					//			cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<H[m+l*nbk+j*nk4th]<<endl;
				}
				delete [] fr_tmp1;
			}
				
			double *dtauGH=new double[(nbk-2)*nbk*2];
				
			fr_tmp1=new double[(nbk-2)*nbk];
			//TF dans l'espace de d/dtau (H_0(r,tau)G(r,-tau)) a tau=0
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dtauH[m+l*nbk]*Grt[y+(x*(x+1))/2]+H[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauGH[m+l*nbk]=-fr_tmp1[m+l*nbk];
					//				cout<<setw(25)<<dtauGH[m+l*nbk]<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<endl;
				}
			
			//TF dans l'espace de d/dtau (H_0(r,tau)G(r,-tau)) a tau=beta
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dtauH[m+l*nbk]*Grt[y+(x*(x+1))/2+nbw*nk8th]+H[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauGH[m+l*nbk+nk4th]=-fr_tmp1[m+l*nbk];
					//				cout<<setw(25)<<dtauGH[m+l*nbk+nk4th]<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
				}
			
			delete [] fr_tmp1;
			
			delete [] dtauH;
			
			#pragma omp parallel private(j,m)
			{
				long int x, y, tmp; //, p;
				double wt, kn;		//, t, G0t, Jinf
				
				dcomplex *HGtmp=new dcomplex[nbw];
				double CJ1, CJ2, CJ3, d2J0, d2Jbeta, CJ01, CJ02, CJ03;
				
				double *tau2=new double[nbw+1];
				double *ft_tmp2=new double[nbw+1];
				double *coeffs2=new double[4*nbw];
				
				dcomplex chijj3_qn0_loc=0;
				dcomplex chijj3_qn0_tmp;
				
				double C2, C3;
				dcomplex selftmp;
				dcomplex *Gtmp1=new dcomplex[nbw];
									
			#pragma omp for
				for (l=0; l<nbk-2; l++)
				{
					for (m=0; m<nbk; m++)
					{
						
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						C2=Sigma_inf2[y+(x*(x+1))/2];
						C3=Sigma_inf3[y+(x*(x+1))/2];
						
						// calcul du spline de HG(k,t)
						
						for (j=0; j<=nbw; j++)
						{
							tau2[j]=j/(nbw*tem);
							ft_tmp2[j]=H[m+l*nbk+j*nk4th];
						}
						
						coeffs2[0]=dtauGH[m+l*nbk];
						coeffs2[1]=dtauGH[m+l*nbk+nk4th];				
						spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);
						
						// TF sur t de HG(k,t)		
						for (j=0; j<nbw; j++)
							Gtmp1[j]=exp(dcomplex(0,-j*PI/nbw))*((double)6.0)*coeffs2[4*j];
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
						
						CJ01=H[m+l*nbk]+H[m+l*nbk+nbw*nk4th];
						CJ02=dtauGH[m+l*nbk]+dtauGH[m+l*nbk+nk4th];
						
						d2J0=2*coeffs2[1];
						d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
						CJ03=d2J0+d2Jbeta;
						
						for (j=0; j<nbw/2; j++)
						{
							kn=(2*j+1)*PI*tem;
							HGtmp[j+nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
						}
						for (j=nbw/2; j<nbw; j++)
						{
							kn=(2*(j-nbw)+1)*PI*tem;
							HGtmp[j-nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,-(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
						}
						
						chijj3_qn0_tmp=0;
						
						for (j=0; j<nw; j++)
						{
							kn=(2.0*(j-nw)+1)*PI*tem;
							
							selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
							
							Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							chijj3_qn0_tmp+=Gtmp1[j]*Gtmp1[j]*HGtmp[j]/tem;
						}
						for (j=nbw-1; j>=nw; j--)
						{
							kn=(2.0*(j-nw)+1)*PI*tem;
							
							selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
							
							Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							chijj3_qn0_tmp+=Gtmp1[j]*Gtmp1[j]*HGtmp[j]/tem;
						}
						
						wt=4;
						if (m==0 || m==nbk-1) wt=wt/2;
												
						chijj3_qn0_loc+=U*dek[m+l*nbk]*normFact*wt*chijj3_qn0_tmp/((double)2.0);				
						
					}
				}
			#pragma omp critical
				{					
					chijj3_qn0+=chijj3_qn0_loc;
				}
				
				delete [] Gtmp1;
				delete [] tau2;
				delete [] ft_tmp2;
				delete [] coeffs2;
				delete [] HGtmp;
			}
			
			delete [] H;
			delete [] dtauGH;
			
			cout<<"chijj3(iqn=0) calcule:  "<<chijj3_qn0<<endl;
				
		}
		
		N0=nbw/2+1;
		
		m=mmin;
		int N=m*((N0-1)/((int)pow(2.0,(int)m+1)))+(N0-1)/((int)pow(2.0,(int)m))+1;
		
		cout<<"2^M:  "<<pow(2.0,(int)m)<<"  (N0-1)/((int)pow(2.0,M)):  "<<(N0-1)/((int)pow(2.0,(int)m))<<"   N:  "<<N<<'\n';
		
		long int *ind_freq=new long int[N];
		
		int j0, jf, d;
		
		jf=(N0-1)/((int)pow(2.0,(int)m));
		for (j=0; j<=jf; j++)
		{
			ind_freq[j]=j;
			//				cout<<setw(10)<<j<<ind_freq[j]<<'\n';
		}
		
		int Nj=(N0-1)/((int)pow(2.0,(int)m+1));
		
		long int n=jf;
		d=1;
		for (l=0; l<m; l++)
		{
			j0=jf+1;
			jf=j0+Nj-1;
			d=2*d;
			for (j=j0; j<=jf; j++)
			{
				n=n+d;
				ind_freq[j]=n;
				//				cout<<setw(10)<<j<<ind_freq[j]<<'\n';				
			}
		}
		
		
		double *dtauhG=new double[(nbk-2)*nbk*2];
		
		//TF dans l'espace de d/dtau (h(r,t)G(r,-t))
		{
			fr_tmp1=new double[(nbk-2)*nbk];
			long int x,y,tmp;
			
			//  d/dt h(r,t)G(r,-t)  a t=0	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dh[m+l*nbk]*Grt[y+(x*(x+1))/2]+h[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauhG[m+l*nbk]=fr_tmp1[m+l*nbk];
				}
			
			//  d/dt h(r,t)G(r,-t)  a t=beta	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dh[m+l*nbk+nk4th]*Grt[y+(x*(x+1))/2+nbw*nk8th]+h[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauhG[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
				}
			delete [] fr_tmp1;
		}
		
		//TF dans l'espace de h(r,t)G(r,-t) et de d/dtau (h(r,t)G(r,-t)), le résultat est enregistré dans h dont le contenu est écrasé
		#pragma omp parallel private(fr_tmp1,j,l,m)
		{
			fr_tmp1=new double[(nbk-2)*nbk];
			long int x,y,tmp;	
			
		#pragma omp for
			for (j=0; j<=nbw; j++)
			{
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=h[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
			}
			delete [] fr_tmp1;
		}
		
		//		cout<<" tem: "<<tem<<"\n nbw: "<<nbw<<"\n nbk: "<<nbk<<"\n nk4th: "<<nk4th<<endl;
		//verifier si les derivees sont bonnes		
		//		for (l=0; l<nbk-2; l++)
		//			for (m=0; m<nbk; m++)
		//			{
		//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauhG[m+l*nbk]<<tem*nbw*(h[m+l*nbk+nk4th]-h[m+l*nbk])<<endl;
		//			}
		//		for (l=0; l<nbk-2; l++)
		//			for (m=0; m<nbk; m++)
		//			{
		//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauhG[m+l*nbk+nk4th]<<tem*nbw*(h[m+l*nbk+nbw*nk4th]-h[m+l*nbk+(nbw-1)*nk4th])<<endl;
		//			}
		
		
		double *hg2=new double[nbk*(nbk-2)];
		
		//TF dans le temps de (h(r,t)G(r,-t))(q) (voir les notes "techniques de calcul")
#pragma omp parallel private(l,m,j)
		{
			double chitmp, chisp_tmp1, chisp_tmp2, chich_tmp1, chich_tmp2;
			long int x, y, tmp, wn;
			dcomplex *ft_tmp1=new dcomplex[nbw];
			double *tau=new double[nbw+1];
			double *ftau=new double[nbw+1];
			double *coeffs=new double[4*nbw];
			double integ, a, b, c, d, x1, x2;
			
#pragma omp for
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					for (j=0; j<=nbw; j++)
					{
						tau[j]=j/(nbw*tem);					
						ftau[j]=h[m+l*nbk+j*nk4th];
						//					if (l==2 && m==1 ) cout<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
					}
					
					coeffs[0]=dtauhG[m+l*nbk];
					coeffs[1]=dtauhG[m+l*nbk+nk4th];				
					spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
					//					spline_coeffs(tau, ftau, nbw+1, coeffs);
					
					for (j=0; j<nbw; j++)
						ft_tmp1[j]=6.0*coeffs[4*j];
					
					fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
					
					integ=0;
					for (j=0; j<nbw; j++)
					{
						a=coeffs[4*j];
						b=coeffs[4*j+1];
						c=coeffs[4*j+2];
						d=coeffs[4*j+3];
						x1=tau[j+1]-tau[j];
						integ+=a*x1*x1*x1*x1/4.0+b*x1*x1*x1/3.0+c*x1*x1/2.0+d*x1;
						
						//						x1=tau[j];
						//						x2=tau[j+1];
						//						integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
					}
					h[m+l*nbk]=2*integ;
					
					//					cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<0<<h[m+l*nbk]<<endl;
					
					hg2[m+l*nbk]=2*(dtauhG[m+l*nbk]-dtauhG[m+l*nbk+nk4th]);
					for (j=1; j<nbw/2; j++)
					{
						h[m+l*nbk+j*nk4th]=2*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]);
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<h[m+l*nbk+j*nk4th]<<endl;
					}
					
				}
			delete [] ft_tmp1;
			delete [] tau;
			delete [] ftau;
			delete [] coeffs;
		}
		
		delete [] dtauhG;

		
		dcomplex *HG=new dcomplex[array_size_h];
		dcomplex *HG2=new dcomplex[array_size_h];		
//		dcomplex *HG=new dcomplex[nbk*(nbk-2)*(nbw+1)];
//		dcomplex *HG2=new dcomplex[nbk*(nbk-2)*(nbw+1)];
		dcomplex *dtauHn=new dcomplex[nbk*(nbk-2)*2];
		dcomplex *dtauHn2=new dcomplex[nbk*(nbk-2)*2];
				
		double qn;
		int nc;
		
		dcomplex *chijj3tmp=new dcomplex[N];
		for (j=0; j<N; j++) chijj3tmp[j]=0;
		
		//		n=1;
		//		for (n=1; n<N; n+=N/8)
		for (n=1; n<N; n++)
		{
			nn=ind_freq[n];
			qn=2.0*nn*PI*tem;
			
			nc=nn/2;
			//			nc=nn/4;
			//			cout<<"nc:  "<<nc<<endl;
			
			//definition de I_n(q,iqm)
#pragma omp parallel private(l,m,j)
			{
				double HGtmp1, HGtmp2, wn;
				double chi0tmp, chi0n, chisp_tmp,chich_tmp,chisp_n,chich_n;
				double D1, D2;
				long int x, y, tmp, p, q;
				
#pragma omp for				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						//						cout<<-hg2[m+l*nbk]/(chi1->chi_0->chi0_inf[y+(x*(x+1))/2])<<endl;
						
						D1=0;
						D2=0;
						for (j=0; j<nbw/2-nc; j++)
						{	
							p=j;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
								chi0tmp=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							else
							{
								//								chi0tmp=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0tmp=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							chisp_tmp=chi0tmp/(1.0-Usp*chi0tmp/2);
							chich_tmp=chi0tmp/(1.0+Uch*chi0tmp/2);
							
							if (p<nbw/2 && p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp2=h[m+l*nbk];							
							
							p=j+nn;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
							{
								chi0n=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
								//								chi0n=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0n=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							
							if (p<nbw/2 && p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp1=h[m+l*nbk];
							
							HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							//								HG[m+l*nbk+j*nk4th]=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<HG[m+l*nbk+j*nk4th]<<endl;
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{	
							p=j-nbw;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
							{
								chi0tmp=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
								//								chi0tmp=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0tmp=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							chisp_tmp=chi0tmp/(1.0-Usp*chi0tmp/2);
							chich_tmp=chi0tmp/(1.0+Uch*chi0tmp/2);
							
							if (p<nbw/2 && p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp2=h[m+l*nbk];
							
							p=j-nbw+nn;
							wn=2*PI*p*tem;
							if (p<0) p=-p;
							if (p<=nbw/2)
							{
								chi0n=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
								//								chi0n=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0n=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							
							if (p<nbw/2 && p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp1=h[m+l*nbk];
							
							if (nn%2==0 && (j-nbw)!=-nn/2)
							{
								D2=D1;
								D1=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							}
							
							if ( (j-nbw)!=-nn/2 || nn%2 )
								HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							//								HG[m+l*nbk+j*nk4th]=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);								
							else
								HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(4*D1-D2)/3.0;
							//								HG[m+l*nbk+j*nk4th]=(4*D1-D2)/3.0;
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j-nbw<<HG[m+l*nbk+j*nk4th]<<endl;
						}
						
					}
			}
			
			//			dcomplex *Htmp=new dcomplex[nbk*(nbk-2)*(nbw+1)];
			
			//TF dans le temps de I_n(q,iqm) et calcul de d/dtau I_n(q,tau) 
#pragma omp parallel private(l,m,j)
			{
				long int x, y, tmp;
				double rttmp, wn, wn1, t;
				dcomplex zsp0, zsp1;
				double zch0, zch1;
				double chisp, chich, Hinf, cq;
				dcomplex *ft_tmp1=new dcomplex[nbw];
				dcomplex dtauH0, dtauH_inf, Htau_inf, tmp_sp, tmp_ch, exp0, exp1;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						dtauH0=0;
						//						cq=-chi1->chi_0->chi0_inf[y+(x*(x+1))/2];
						cq=-chi0_inf[y+(x*(x+1))/2];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=2*j*PI*tem;
							chisp=1.0/(-wn*wn-Usp*cq/2.0);
							chich=1.0/(-wn*wn+Uch*cq/2.0);
							Hinf=hg2[m+l*nbk]*(3.0*Usp*chisp+Uch*chich);
							ft_tmp1[j]=tem*(HG[m+l*nbk+j*nk4th]-Hinf);
							dtauH0-=wn*ft_tmp1[j];
							
							//							Htmp[m+l*nbk+j*nk4th]=Hinf;
							//							Htmp[m+l*nbk+j*nk4th]=HG[m+l*nbk+j*nk4th];
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<setw(30)<<Hinf<<real(HG[m+l*nbk+j*nk4th]-Hinf)/real(HG[m+l*nbk+j*nk4th])<<endl;
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=2*(j-nbw)*PI*tem;
							chisp=1.0/(-wn*wn-Usp*cq/2.0);
							chich=1.0/(-wn*wn+Uch*cq/2.0);
							Hinf=hg2[m+l*nbk]*(3.0*Usp*chisp+Uch*chich);
							ft_tmp1[j]=tem*(HG[m+l*nbk+j*nk4th]-Hinf);
							dtauH0-=wn*ft_tmp1[j];
							
							//							Htmp[m+l*nbk+j*nk4th]=Hinf;
							//							Htmp[m+l*nbk+j*nk4th]=HG[m+l*nbk+j*nk4th];
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG[m+l*nbk+j*nk4th]<<setw(30)<<Hinf<<real(HG[m+l*nbk+j*nk4th]-Hinf)/real(HG[m+l*nbk+j*nk4th])<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						//						rttmp=sqrt(Usp*chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/2.0);
						rttmp=sqrt(Usp*chi0_inf[y+(x*(x+1))/2]/2.0);
						zsp0=dcomplex(0,rttmp);
						zsp1=-zsp0;
						//						zch0=sqrt(Uch*chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/2.0);
						zch0=sqrt(Uch*chi0_inf[y+(x*(x+1))/2]/2.0);
						zch1=-zch0;
						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<(zsp0*zsp0-Usp*cq/2.0)<<setw(40)<<(zsp1*zsp1-Usp*cq/2.0)<<endl;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<(zch0*zch0+Uch*cq/2.0)<<setw(40)<<(zch1*zch1+Uch*cq/2.0)<<endl;
						
						for (j=0; j<nbw; j++)
						{
							t=j/(nbw*tem);
							tmp_sp=exp(-zsp0*t)/((exp(-zsp0/tem)-(double)1.0)*(zsp0-zsp1))							
							+exp(-zsp1*t)/((exp(-zsp1/tem)-(double)1.0)*(zsp1-zsp0));
							tmp_ch=exp(-zch0*t)/((exp(-zch0/tem)-1.0)*(zch0-zch1))
							+exp(zch1*(1.0/tem-t))/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
							Htau_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
							HG[m+l*nbk+j*nk4th]=ft_tmp1[j]+Htau_inf;
							//							HG[m+l*nbk+j*nk4th]=Htau_inf;
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<ft_tmp1[j]<<setw(38)<<Htau_inf<<HG[m+l*nbk+j*nk4th]<<endl;
							//							Htmp[m+l*nbk+j*nk4th]=Htau_inf;
						}
						j=nbw;
						exp0=exp(zsp0/tem);
						exp1=exp(zsp1/tem);
						tmp_sp=((double)1.0)/(((double)1.0-exp0)*(zsp0-zsp1))
						+((double)1.0)/(((double)1.0-exp1)*(zsp1-zsp0));
						exp0=exp(-zch0/tem);
						tmp_ch=exp0/((exp0-(double)1.0)*(zch0-zch1))
						+1.0/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
						Htau_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						HG[m+l*nbk+nbw*nk4th]=ft_tmp1[0]+Htau_inf;
						//						HG[m+l*nbk+nbw*nk4th]=Htau_inf;						
						
						//						Htmp[m+l*nbk+nbw*nk4th]=Htau_inf;
						
						//						if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<ft_tmp1[0]<<setw(38)<<Htau_inf<<HG[m+l*nbk+j*nk4th]<<endl;
						
						// d/dtau Hn(q,tau) a tau=0	
						tmp_sp=-zsp0/((exp(-zsp0/tem)-(double)1.0)*(zsp0-zsp1))							
						-zsp1/((exp(-zsp1/tem)-(double)1.0)*(zsp1-zsp0));
						exp1=exp(zch1/tem);
						tmp_ch=-zch0/((exp(-zch0/tem)-1.0)*(zch0-zch1))
						-zch1*exp1/(((double)1.0-exp1)*(zch1-zch0));
						dtauH_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						dtauHn[m+l*nbk]=I*dtauH0+dtauH_inf;
						//						dtauHn[m+l*nbk]=dtauH_inf;						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<setw(38)<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauH_inf<<nbw*tem*(Htmp[m+l*nbk+nk4th]-Htmp[m+l*nbk])<<endl;
						
						// d/dtau Hn(q,tau) a tau=beta
						exp0=exp(-zsp0/tem);
						exp1=exp(-zsp1/tem);
						tmp_sp=-zsp0*exp0/((exp0-(double)1.0)*(zsp0-zsp1))
						-zsp1*exp1/((exp1-(double)1.0)*(zsp1-zsp0));
						exp0=exp(-zch0/tem);
						tmp_ch=-zch0*exp0/((exp0-(double)1.0)*(zch0-zch1))
						-zch1/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
						dtauH_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						dtauHn[m+l*nbk+nk4th]=I*dtauH0+dtauH_inf;						
						//						dtauHn[m+l*nbk+nk4th]=dtauH_inf;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauH_inf<<nbw*tem*(Htmp[m+l*nbk+nbw*nk4th]-Htmp[m+l*nbk+(nbw-1)*nk4th])<<endl;
					}
				delete [] ft_tmp1;
			}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
                
				//TF dans l'espace de I_n(q,tau)				
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							fr_tmp1[m+l*nbk]=normFact*real(HG[m+l*nbk+j*nk4th]);
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                        {
                            HG[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG[m+l*nbk+j*nk4th].imag());
//							HG[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
                        }
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							fr_tmp1[m+l*nbk]=normFact*imag(HG[m+l*nbk+j*nk4th]);
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                        {
                            HG[m+l*nbk+j*nk4th]=dcomplex(HG[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
                        }
				}
				
				// TF dans l'espace de d/dtau I_n(q,tau) a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*real(dtauHn[m+l*nbk]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk].imag());
//						dtauHn[m+l*nbk].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*imag(dtauHn[m+l*nbk]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
                        dtauHn[m+l*nbk]=dcomplex(dtauHn[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
					}
				
				// TF dans l'espace de d/dtau I_n(q,tau) a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*real(dtauHn[m+l*nbk+nk4th]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk+nk4th].imag());
//						dtauHn[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*imag(dtauHn[m+l*nbk+nk4th]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
                        dtauHn[m+l*nbk+nk4th]=dcomplex(dtauHn[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;
					}
				delete [] fr_tmp1;
			}
			
			for (j=0; j<=nbw; j++)
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						HG2[m+l*nbk+j*nk4th]=HG[m+l*nbk+(nbw-j)*nk4th];
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauHn2[m+l*nbk]=-dtauHn[m+l*nbk+nk4th];
					dtauHn2[m+l*nbk+nk4th]=-dtauHn[m+l*nbk];
				}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
				//TF dans l'espace de d/dtau I_n(r,tau)G_n(r,-tau)
				// a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+real(dtauHn[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk].imag());
//						dtauHn[m+l*nbk].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+imag(dtauHn[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn[m+l*nbk]=dcomplex(dtauHn[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
                    }
				
				// a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+real(dtauHn[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk+nk4th].imag());
//						dtauHn[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+imag(dtauHn[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk+nk4th]=dcomplex(dtauHn[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
				
				
				//TF dans l'espace de I_n(r,tau)G_n(r,-tau)
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=real(HG[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                        {
                            HG[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG[m+l*nbk+j*nk4th].imag());
//							HG[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
                        }
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                        {
                            HG[m+l*nbk+j*nk4th]=dcomplex(HG[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
                        }
				}
				
				//				for (l=0; l<nbk-2; l++)
				//					for (m=0; m<nbk; m++)
				//					{
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//					}
				
				delete [] fr_tmp1;
			}
			
			//			l=2, m=4;
			//			for (j=0; j<=nbw; j++)
			//				cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
			
			
			
			//TF dans le temps de (I_n(r,tau)G_n(r,-tau))(k)
#pragma omp parallel private(l,m,j)
			{
				dcomplex *ft_tmp1=new dcomplex[nbw];
				double *tau=new double[nbw+1];
				double *ftau=new double[nbw+1];
				double *coeffs=new double[4*nbw];
				double wn;
				dcomplex coefwn, d2S0, d2Sbeta;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=real(HG[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=real(dtauHn[m+l*nbk]);
						coeffs[1]=real(dtauHn[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
						//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0=2*coeffs[1];
						d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];
						//						d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
						
						for (j=0; j<nbw; j++)
							ft_tmp1[j]=exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
						
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=imag(HG[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=imag(dtauHn[m+l*nbk]);
						coeffs[1]=imag(dtauHn[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
						//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0+=I*((double)2.0)*coeffs[1];
						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3]);
						//						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3]);
						
						for (j=0; j<nbw; j++)
						{
							ft_tmp1[j]+=I*exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
							//							if (l==0 && m==0)  cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<ft_tmp1[j]<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						coefwn=HG[m+l*nbk]+HG[m+l*nbk+nbw*nk4th];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=(2*(j+nn)+1)*PI*tem;
							if ((j+nn)<nbw)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn]/(wn*wn*wn*wn);
							else
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
						}						
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=(2*(j+nn-nbw)+1)*PI*tem;
							if ((j+nn-nbw)<0)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn]/(wn*wn*wn*wn);
							else if ((j+nn-nbw)<nbw)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn)
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn-nbw]/(wn*wn*wn*wn);
							else
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
						}						
						
					}
				delete [] ft_tmp1;
				delete [] tau;
				delete [] ftau;
				delete [] coeffs;
			}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
				//TF dans l'espace de d/dtau I_n(r,-tau)G_n(r,-tau)
				// a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+real(dtauHn2[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn2[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn2[m+l*nbk].imag());
//						dtauHn2[m+l*nbk].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+imag(dtauHn2[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn2[m+l*nbk]=dcomplex(dtauHn2[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn2[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
                    }
				
				// a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+real(dtauHn2[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn2[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn2[m+l*nbk+nk4th].imag());
//						dtauHn2[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
                    }
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+imag(dtauHn2[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                    {
                        dtauHn2[m+l*nbk+nk4th]=dcomplex(dtauHn2[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn2[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
                    }
				
				
				//TF dans l'espace de I_n(r,-tau)G_n(r,-tau)
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG2[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG2[m+l*nbk+j*nk4th].imag());
//							HG2[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG2[m+l*nbk+j*nk4th]=dcomplex(HG2[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG2[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
				}
				
				//				for (l=0; l<nbk-2; l++)
				//					for (m=0; m<nbk; m++)
				//					{
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn2[m+l*nbk]<<nbw*tem*(HG2[m+l*nbk+nk4th]-HG2[m+l*nbk])<<endl;
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn2[m+l*nbk+nk4th]<<nbw*tem*(HG2[m+l*nbk+nbw*nk4th]-HG2[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//					}
				
				delete [] fr_tmp1;
			}
			
			
			//TF dans le temps de (I_n(r,-tau)G_n(r,-tau))(k)
#pragma omp parallel private(l,m,j)
			{
				dcomplex *ft_tmp1=new dcomplex[nbw];
				double *tau=new double[nbw+1];
				double *ftau=new double[nbw+1];
				double *coeffs=new double[4*nbw];
				double wn;
				dcomplex coefwn, d2S0, d2Sbeta;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=real(HG2[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=real(dtauHn2[m+l*nbk]);
						coeffs[1]=real(dtauHn2[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
						//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0=2*coeffs[1];
						d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];
						//						d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
						
						for (j=0; j<nbw; j++)
							ft_tmp1[j]=exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
						
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=imag(HG2[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=imag(dtauHn2[m+l*nbk]);
						coeffs[1]=imag(dtauHn2[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
						//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0+=I*((double)2.0)*coeffs[1];
						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3]);
						//						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3]);
						
						for (j=0; j<nbw; j++)
						{
							ft_tmp1[j]+=I*exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
							//							if (l==0 && m==0)  cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<ft_tmp1[j]<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						coefwn=HG2[m+l*nbk]+HG2[m+l*nbk+nbw*nk4th];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=(2*j+1)*PI*tem;
							HG2[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn2[m+l*nbk]+dtauHn2[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
							+ ((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*ft_tmp1[j]/(wn*wn*wn*wn);							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG2[m+l*nbk+j*nk4th]<<endl;
						}						
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=(2*(j-nbw)+1)*PI*tem;
							HG2[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn2[m+l*nbk]+dtauHn2[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
							+ ((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*ft_tmp1[j]/(wn*wn*wn*wn);							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG2[m+l*nbk+j*nk4th]<<endl;
						}						
						
					}
				delete [] ft_tmp1;
				delete [] tau;
				delete [] ftau;
				delete [] coeffs;
			}
			
			//			dcomplex *Gr1=new dcomplex[nbk*(nbk-2)*nbw];
			//			dcomplex *Gr2=new dcomplex[nbk*(nbk-2)*nbw];
			
			// somme finale sur k et ik_n			
#pragma omp parallel private(l,m,j)
			{
				long int x, y, tmp, p;
				double C2, C3, kn, wt;
				dcomplex Gtmp, Gtmpn, selftmp;
				dcomplex chijj3_qn_loc=0;
				
				double theta_k;
				int ind_theta;
				double *chijj3_k_loc=new double[(nw+1)*N_theta_k];
				for (j=0; j<(nw+1)*N_theta_k; j++)	chijj3_k_loc[j]=0;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
				{
					for (m=0; m<nbk; m++)
					{
						
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						//						C2=G1->Sigma_inf2[y+(x*(x+1))/2];
						//						C3=G1->Sigma_inf3[y+(x*(x+1))/2];
						
						C2=Sigma_inf2[y+(x*(x+1))/2];
						C3=Sigma_inf3[y+(x*(x+1))/2];						
						
						wt=4;
						if (m==0 || m==nbk-1) wt=wt/2;
						
						//						for (p=0; p<nk_chijj; p++)
						//						{
						//							if (k_chijj[2*p]==(l+1) && k_chijj[2*p+1]==m)	break;
						//						}
						
						//						if (p<nk_chijj) chijj3_k[p + n*nk_chijj]=0;
						
						theta_k=atan( (1.0*(nk0/2-m))/(nk0/2-l-1) );
						ind_theta=(int)floor(theta_k/Dtheta);
						
						for (j=0; j<nbw/2-nc; j++)
						{
							kn=(2*j+1)*PI*tem;
							
							selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), 
											  -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
							
							Gtmp=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr1[m+l*nbk+j*nk4th]=Gtmp;
							//							Gr1[m+l*nbk+j*nk4th]=selftmp;
							//							Gr1[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							kn=(2*(j+nn)+1)*PI*tem;
							if ((j+nn)<nbw/2)
								selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j+nn)*nk8th])/(kn*kn*kn*kn), 
												  -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j+nn)*nk8th])/(kn*kn*kn*kn));
							else
								selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							//								selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							
							Gtmpn=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr2[m+l*nbk+j*nk4th]=Gtmpn;
							//							Gr2[m+l*nbk+j*nk4th]=selftmp;
							//							Gr2[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							chijj3_qn_loc-=I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
							
							//							if (p<nk_chijj)	
							//								chijj3_k[p + n*nk_chijj]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							chijj3_k_loc[ind_theta + n*N_theta_k]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							//							cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<chijj3_qn_loc<<endl;
							//							HGsum[m+l*nbk+j*nk4th]=-I*U*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{
							kn=(2*(j-nbw)+1)*PI*tem;
							
							if ((nbw-j-1)<nbw/2)
								selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th])/(kn*kn*kn*kn), 
												  -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th])/(kn*kn*kn*kn));
							else
								selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							//								selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							
							Gtmp=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr1[m+l*nbk+j*nk4th]=Gtmp;
							//							Gr1[m+l*nbk+j*nk4th]=selftmp;
							//							Gr1[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							kn=(2*(j+nn-nbw)+1)*PI*tem;
							if ((j+nn-nbw)<0)
							{
								if ((nbw-j-nn-1)<nbw/2)
									selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nbw-j-nn-1)*nk8th])/(kn*kn*kn*kn), 
													  -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nbw-j-nn-1)*nk8th])/(kn*kn*kn*kn));
								else
									selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
								//									selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							}
							else
							{
								if ((j+nn-nbw)<nbw/2)
									selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j+nn-nbw)*nk8th])/(kn*kn*kn*kn), 
													  -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j+nn-nbw)*nk8th])/(kn*kn*kn*kn));
								else
									selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
								//									selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							}		
							
							Gtmpn=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr2[m+l*nbk+j*nk4th]=Gtmpn;
							//							Gr2[m+l*nbk+j*nk4th]=selftmp;
							//							Gr2[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							chijj3_qn_loc-=I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
							
							//							if (p<nk_chijj)	
							//								chijj3_k[p + n*nk_chijj]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							chijj3_k_loc[ind_theta + n*N_theta_k]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							//							cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<chijj3_qn_loc<<endl;
							//							HGsum[m+l*nbk+j*nk4th]=-I*U*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
						}					
					}
				}
#pragma omp critical
				{
					chijj3tmp[n]+=chijj3_qn_loc;
					chijj3_qn[n]+=real(chijj3_qn_loc);
					for (p=0; p<N_theta_k; p++) chijj3_k[p + n*N_theta_k]+=chijj3_k_loc[p + n*N_theta_k];
				}
				delete [] chijj3_k_loc;
			}
			
			cout<<setw(10)<<nn<<setw(50)<<chijj3tmp[n]<<endl;
			
			//			cout<<setw(10)<<nn<<setw(50)<<chijj3tmp[n];	
			//			if (n%8==0 || n==1)
			//				calc_chijj_vertex_corr2_wn(nn, k0, wn0, NULL, HG0tmp, NULL, NULL, NULL);
			//			else
			//				cout<<endl;
			
			//			int k0[]={0,0};
			//			k0[0]=1;
			//			k0[1]=2;
			//			int wn0=-20;
			//			calc_chijj_vertex_corr2_wn(nn, k0, wn0, NULL, HG0tmp, HGsum, Gr1, Gr2);
			//			if (wn0<0)
			//				cout<<HGsum[k0[1]+(k0[0]-1)*nbk+(wn0+nbw)*nk4th]<<endl;
			//			else
			//				cout<<HGsum[k0[1]+(k0[0]-1)*nbk+wn0*nk4th]<<endl;
			
			
		}
		
		delete [] chijj3tmp;
		delete [] dtauHn;
		delete [] dtauHn2;
		delete [] HG;
		delete [] HG2;		
		
		
		/*
		// Pour enregistrer chijj3(iq_n) en binaire
		
		const char *nameForm_chijj3_qn_bin="./chijj3_qn_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameForm_chijj3_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
		file3.open(name, ios::out | ios::binary);
		
		for (j=0; j<N; j++)
		{
			file3.write((char*)&(ind_freq[j]),sizeof(ind_freq[j]));
			file3.write((char*)&(chijj3_qn[j]),sizeof(chijj3_qn[j]));
		}
		file3.close();
		// fin de l'enregistrement en binaire
		*/
		
		// Pour enregistrer chijj3(iq_n) en ascii
		const char *nameForm_chijj3_qn="./chijj3_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameForm_chijj3_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
		file3.open(name, ios::out );
		file3<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
		
		for (j=0; j<N; j++)
		{
			file3<<setw(10)<<ind_freq[j]<<chijj3_qn[j]<<'\n';
		}
		file3.close();
		// fin de l'enregistrement en ascii
		
		for (j=0; j<=nbw/2; j++) chijj_qn[j]=0;
		
		chijj_qn[0]=chijj1_qn[0]+chijj2_qn[0];
		
		for (j=1; j<N; j++)
			chijj_qn[j]=chijj1_qn[ind_freq[j]]+chijj2_qn[ind_freq[j]]+chijj3_qn[j];
		
		// Pour enregistrer chijj(iq_n) en ascii
		const char *nameForm_chijj_qn="./chijj_with_all_vc_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameForm_chijj_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
		file.open(name, ios::out );
		file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
		
		for (j=0; j<N; j++)
		{
			file<<setw(10)<<ind_freq[j]<<chijj_qn[j]<<'\n';
		}
		file.close();
		// fin de l'enregistrement en ascii
		
		if (N_theta_k>1)
		{
			// sommer chijj3_k pour tous les k
			
			double *chijj3_k_tot=new double[nw+1];
			
			for (j=0; j<=nw; j++)
				chijj3_k_tot[j]=0;
			
			for (j=0; j<=nw; j++)
				for  (p=0; p<N_theta_k; p++)
					chijj3_k_tot[j]+=chijj3_k[p + j*N_theta_k];
			
			//enregistrer chijj3_k et chijj3_k_tot et la comparaison de chijj3_k_tot avec chijj3_qn dans un fichier
			
			const char *nameForm_chijj3_k="./chijj3_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
			const char *nameForm_chijj3_k_tot="./chijj3_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
			
			sprintf(name, nameForm_chijj3_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
			file.open(name, ios::out );
			file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
			
			sprintf(name, nameForm_chijj3_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
			file2.open(name, ios::out );
			file2<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
			
			for (j=1; j<N; j++)
			{
				file2<<setw(10)<<ind_freq[j]<<setw(30)<<chijj3_k_tot[j]<<setw(30)<<chijj3_qn[j]<<chijj3_k_tot[j]/chijj3_qn[j]<<'\n';
				for  (p=0; p<N_theta_k; p++)
					file<<setw(10)<<ind_freq[j]<<setw(10)<<(p+1)<<chijj3_k[p + j*N_theta_k]<<'\n';
				//				file<<setw(10)<<ind_freq[j]<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj3_k[p + j*nk_chijj]<<'\n';
			}
			file.close();
			file2.close();
			
			delete [] chijj3_k_tot;
		}
		
		delete [] ind_freq;
		delete [] hg2;
		
		delete [] Grt;
		delete [] dtauGrt;
		
		delete [] chijj3_k;
	}
	
	delete [] h;
	delete [] df0;
	delete [] dh;
	
	delete [] chijj_qn;	
	
	fftw_destroy_plan(fftplan);
	fftw_destroy_plan(fftplan_t_back);
	fftw_destroy_plan(fftplan_RO);
	
	delete [] dek;
	
}


//! calculate all the terms in the optical conductivity with FFT and cubic spline (unstable at high temperature because it uses 4th degree roots for asymptotic G)
void cond_opt::calc_chijj_vertex_corr_all(bool chi3, int mmin, int N_theta_k, int *k_chijj, int nk_chijj)
{	
	cout<<"calcul de chijj par calc_chijj_vertex_corr_all()\n";
	
	long int j, l, m;
	
	int nw=nbw/2;
	int nk0=2*(nbk-1);
	int nk4th=nbk*(nbk-2);
	int nk8th=(nbk*(nbk+1))/2;
	
	if (!Self_kw_array)
	{		
		calc_Self_FFT_spline();
		traceSelfG();
		green::find_mu();
	}
	
	if (!mu_calcule)
	{
		cout<<"calc_chijj_vertex_corr_all(): mu non calcule\n";
		return;
	}
	
	if (k_chijj && nk_chijj) save_self_ikn(k_chijj, nk_chijj, nw/16);
	find_FS_inter();
		
	V_array=vertex_rt_array;
		
	double normFact=1.0/(nk0*nk0);
	
	double *dek=new double[nbk*(nbk-2)];
	double *d2ek=new double[nbk*nbk];
	int indHss[]={0,0};
	double vk[2];
	
	double kx, ky;
	
	for (l=0; l<nbk; l++)
	{
		kx=l*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			d2ek[m+l*nbk]=Hessian_ek(vk, indHss);
		}
	}
	
	
	for (l=0; l<nbk-2; l++)
	{
		kx=(l+1)*PI/(nbk-1);
		vk[0]=kx;
		for (m=0; m<nbk; m++)
		{
			ky=m*PI/(nbk-1);
			vk[1]=ky;
			dek[m+l*nbk]=grad_ek(vk, 0);
		}
	}
	
	double  *fr_tmp=new double[(nbk-2)*nbk];
	fftw_plan fftplan_RO=fftw_plan_r2r_2d(nbk-2, nbk, fr_tmp, fr_tmp, FFTW_RODFT00, FFTW_REDFT00, FFTW_MEASURE);
	delete [] fr_tmp;
	
	double *fr_tmp1;
	
	long int array_size_h=(long int) nbk*(nbk-2)*(nbw+1);

	double *h=new double[array_size_h];	
//	double *h=new double[nbk*(nbk-2)*(nbw+1)];
	
	// calcul de h(r,t) (voir l'annexe sur les techniques de calcul)
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[(nbk-2)*nbk];
		double t;
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			t=j/(tem*nbw);
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
					else
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp( (ek[m+(l+1)*nbk]-mu0)*(1.0/tem-t))/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		delete [] fr_tmp1;
	}
	
	// calcul de la derivee de h(r,t) a t=0 et t=beta
	
	double *dh=new double[(nbk-2)*nbk*2];
	
	// t=0	
	fr_tmp1=new double[(nbk-2)*nbk];
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
			else
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (ek[m+(l+1)*nbk]-mu0)/tem)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dh[m+l*nbk]=fr_tmp1[m+l*nbk];
	
	// t=beta	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
		{
			if (ek[m+(l+1)*nbk]>mu0)
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)*exp( (mu0-ek[m+(l+1)*nbk])/tem )/(exp( (mu0-ek[m+(l+1)*nbk])/tem )  + 1.0 );
			else
				fr_tmp1[m+l*nbk]=-normFact*dek[m+l*nbk]*(ek[m+(l+1)*nbk]-mu0)/(exp( (ek[m+(l+1)*nbk]-mu0)/tem )  +1.0 );
		}
	
	fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
	
	for (l=0; l<nbk-2; l++)
		for (m=0; m<nbk; m++)
			dh[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];	
	
	delete [] fr_tmp1;
	
	double *f0=new double[array_size_h];	
//	double *f0=new double[(nbk-2)*nbk*(nbw+1)];
	double *df0=new double[(nbk-2)*nbk*2];
	
	// calcul de f_0(r,t) (voir l'annexe sur les techniques de calcul)
	//#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[(nbk-2)*nbk];
		double t, expk;
		
		//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			t=j/(tem*nbw);
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					if (ek[m+(l+1)*nbk]>mu0)
					{
						expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((mu0-ek[m+(l+1)*nbk])*t)*( t/(expk + 1.0) -expk/(tem*(expk+1.0)*(expk+1.0)));
					}
					else
					{
						expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
						fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*exp((ek[m+(l+1)*nbk]-mu0)*(1.0/tem-t))*(t/(expk+1.0) - 1.0/(tem*(expk+1.0)*(expk+1.0)));
					}
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					f0[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		
		// calcul de la derivee de  f0(r,t) a t=0 et t=beta
		
		//	t=0
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				if (ek[m+(l+1)*nbk]>mu0)
				{
					expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
					fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*(  (ek[m+(l+1)*nbk]-mu0)*expk/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk + 1.0));
				}
				else
				{
					expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
					fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*( (ek[m+(l+1)*nbk]-mu0)/(tem*(expk+1.0)*(expk+1.0))+1.0/(expk+1.0));
				}
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				df0[m+l*nbk]=fr_tmp1[m+l*nbk];
			}
		
		// t=beta	
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				if (ek[m+(l+1)*nbk]>mu0)
				{
					expk=exp((mu0-ek[m+(l+1)*nbk])/tem);
					fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*expk*(  (mu0-ek[m+(l+1)*nbk])*( 1.0/(tem*(expk + 1.0)) -expk/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk + 1.0));
				}
				else
				{
					expk=exp((ek[m+(l+1)*nbk]-mu0)/tem);
					fr_tmp1[m+l*nbk]=normFact*dek[m+l*nbk]*( (mu0-ek[m+(l+1)*nbk])*(1.0/(tem*(expk+1.0)) - 1.0/(tem*(expk+1.0)*(expk+1.0)))+1.0/(expk+1.0));
				}
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				df0[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
			}
		
		delete [] fr_tmp1;
	}
	
	double *hV=new double[array_size_h];
//	double *hV=new double[nbk*(nbk-2)*(nbw+1)];
	
	// TF sur l'espace de h(r,t)V(r,t), ( TF_r(h(r,t)V(r,t))(k) )
	// le facteur 2 vient du fait que vertex_array(r,tau) est proportionnel a U/8 alors que V_array(r,tau) est proportionnel a U/4
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		long int x, y, tmp, wn;
		fr_tmp1=new double[(nbk-2)*nbk];
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			wn=j;
			if (wn>nbw/2) wn=nbw-wn;
			
			for (l=0; l<nbk-2; l++)
			{
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					
					fr_tmp1[m+l*nbk]=2.0*h[m+l*nbk+j*nk4th]*V_array[y+(x*(x+1))/2+wn*nk8th];
				}
			}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					hV[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		delete [] fr_tmp1;
	}
	
	double TPSC_fact=U*(3.0*Usp + Uch)/4;
	
	double *dhV=new double[nbk*(nbk-2)*2];
	double *df0V=new double[nbk*(nbk-2)*2];
	
	// TF sur l'espace de d/dt(h(r,t)V(r,t)) a t=0 et t=beta
	{	
		long int x, y, tmp, wn;
		fr_tmp1=new double[(nbk-2)*nbk];
		
		//t=0
		
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
//				fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk]*V_array[y+(x*(x+1))/2]+ h[m+l*nbk]*TPSC_fact*G1->chi_0->dtau_chirt0[y+(x*(x+1))/2];
				fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk]*V_array[y+(x*(x+1))/2]+ h[m+l*nbk]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
			}
		}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				dhV[m+l*nbk]=fr_tmp1[m+l*nbk];
		
		//t=beta
		
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
//				fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- h[m+l*nbk+nbw*nk4th]*TPSC_fact*G1->chi_0->dtau_chirt0[y+(x*(x+1))/2];
				fr_tmp1[m+l*nbk]=2.0*dh[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- h[m+l*nbk+nbw*nk4th]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
			}
		}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				dhV[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
		
		
		// TF sur l'espace de d/dt f_0(r,t)V(r,t) a t=0 et t=beta
		
		//t=0
		
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
//				fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk]*V_array[y+(x*(x+1))/2]+ f0[m+l*nbk]*TPSC_fact*G1->chi_0->dtau_chirt0[y+(x*(x+1))/2];
				fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk]*V_array[y+(x*(x+1))/2]+ f0[m+l*nbk]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
			}
		}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				df0V[m+l*nbk]=fr_tmp1[m+l*nbk];
		
		//t=beta
		
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
//				fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- f0[m+l*nbk+nbw*nk4th]*TPSC_fact*G1->chi_0->dtau_chirt0[y+(x*(x+1))/2];
				fr_tmp1[m+l*nbk]=2.0*df0[m+l*nbk+nk4th]*V_array[y+(x*(x+1))/2]- f0[m+l*nbk+nbw*nk4th]*TPSC_fact*dtau_chirt0[y+(x*(x+1))/2];				
			}
		}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				df0V[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
		
		delete [] fr_tmp1;	
	}	
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	// Partie du calcul de chi_jj3 (voir les notes "Techniques de calcul")
	
	fr_tmp=new double[nbk*nbk];
	
	fftw_plan fftplan_RE=fftw_plan_r2r_2d(nbk, nbk, fr_tmp, fr_tmp, FFTW_REDFT00, FFTW_REDFT00, FFTW_MEASURE);
	
	delete [] fr_tmp;

	long int array_size_G=(long int)nk8th*(nbw+1);
	
	double *Grt=new double[array_size_G];
//	double *Grt=new double[nk8th*(nbw+1)];
	double *dtauGrt=new double[nk8th*2];
	
	//	double *Gkt_tmp=new double[nk8th*(nbw+1)];
	
	// calcul de G(r,-t) et d/dt G(r,-t)
	//#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[nbk*nbk];
		double t;
		
		// 	G(r,-t)
		//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			t=j/(tem*nbw);
			for (l=0; l<nbk; l++)
				for (m=0; m<=l; m++)
				{
					if (ek[m+l*nbk]<mu0)
						fr_tmp1[m+l*nbk]=normFact*exp((ek[m+l*nbk]-mu0)*t)/( exp((ek[m+l*nbk]-mu0)/tem) + 1.0 );
					else
						fr_tmp1[m+l*nbk]=normFact*exp((ek[m+l*nbk]-mu0)*(t-1.0/tem))/(exp( -(ek[m+l*nbk]-mu0)/tem )  +1.0 );
					fr_tmp1[l+m*nbk]=fr_tmp1[m+l*nbk];
					
					//					Gkt_tmp[m+(l*(l+1))/2+j*nk8th]=fr_tmp1[m+l*nbk]/normFact;
				}
			
			fftw_execute_r2r(fftplan_RE,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk; l++)
				for (m=0; m<=l; m++)
					Grt[m+(l*(l+1))/2+j*nk8th]=fr_tmp1[m+l*nbk];
		}
		
		// d/dt G(r,-t) a t=0
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				if (ek[m+l*nbk]<mu0)
					fr_tmp1[m+l*nbk]=normFact*(ek[m+l*nbk]-mu0)/( exp((ek[m+l*nbk]-mu0)/tem) + 1.0 );
				else
					fr_tmp1[m+l*nbk]=normFact*(ek[m+l*nbk]-mu0)*exp(-(ek[m+l*nbk]-mu0)/tem)/(exp( -(ek[m+l*nbk]-mu0)/tem )  +1.0 );
				fr_tmp1[l+m*nbk]=fr_tmp1[m+l*nbk];
			}
		
		fftw_execute_r2r(fftplan_RE,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				dtauGrt[m+(l*(l+1))/2]=fr_tmp1[m+l*nbk];
//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauGrt[m+(l*(l+1))/2]<<tem*nbw*(Grt[m+(l*(l+1))/2+nk8th]-Grt[m+(l*(l+1))/2])<<endl;
				dtauGrt[m+(l*(l+1))/2+nk8th]=0;
			}
// 	d/dt G(r,-t) a t=beta
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
				dtauGrt[m+(l*(l+1))/2+nk8th]=er[m+(l*(l+1))/2];
		dtauGrt[nk8th]=-mu0;
		for (l=0; l<nbk; l++)
			for (m=0; m<=l; m++)
			{
				dtauGrt[m+(l*(l+1))/2+nk8th]-=dtauGrt[m+(l*(l+1))/2];
//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauGrt[m+(l*(l+1))/2+nk8th]<<tem*nbw*(Grt[m+(l*(l+1))/2+nbw*nk8th]-Grt[m+(l*(l+1))/2+(nbw-1)*nk8th])<<endl;
			}
		
		delete [] fr_tmp1;
	}
	
	fftw_destroy_plan(fftplan_RE);
	
	
	//calcul de H_0(q,iqn)

	double *H=new double[array_size_h];	
//	double *H=new double[nbk*(nbk-2)*(nbw+1)];
	double *dtauf0G=new double[nbk*(nbk-2)*2];
	
	//	dcomplex *f0G_tmp=new dcomplex[nbk*(nbk-2)*(nbw+1)];
	
	
	//TF dans l'espace de f_0(r,t)G(r,-t) et d/dt f_0(r,t)G(r,-t) a t=0 et t=beta
	//#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[(nbk-2)*nbk];
		long int x, y, tmp;
		
		//	f_0(r,t)G(r,-t)	
		//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=f0[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
					
					//					f0G_tmp[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					H[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		
		//  d/dt f_0(r,t)G(r,-t)  a t=0	
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				fr_tmp1[m+l*nbk]=df0[m+l*nbk]*Grt[y+(x*(x+1))/2]+f0[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				dtauf0G[m+l*nbk]=fr_tmp1[m+l*nbk];
				//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauf0G[m+l*nbk]<<tem*nbw*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<endl;
			}
		
		//  d/dt f_0(r,t)G(r,-t)  a t=beta	
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				fr_tmp1[m+l*nbk]=df0[m+l*nbk+nk4th]*Grt[y+(x*(x+1))/2+nbw*nk8th]+f0[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				dtauf0G[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
				//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauf0G[m+l*nbk+nk4th]<<tem*nbw*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
			}
		
		delete [] fr_tmp1;
	}
	
	
	// TF sur l'espace de f_0(r,t)V(r,t)
	// Attention! le resultat est enregistre dans f0 et le contenu de f0 est ecrase!
#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		long int x, y, tmp, wn;
		fr_tmp1=new double[(nbk-2)*nbk];
		
#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			wn=j;
			if (wn>nbw/2) wn=nbw-wn;
			
			for (l=0; l<nbk-2; l++)
			{
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					
					fr_tmp1[m+l*nbk]=2.0*f0[m+l*nbk+j*nk4th]*V_array[y+(x*(x+1))/2+wn*nk8th];
				}
			}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					f0[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		delete [] fr_tmp1;
	}
	
	free_vertex_array();
	
	dcomplex *Gtmp=new dcomplex[nbw];
	fftw_plan fftplan=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*> (Gtmp), reinterpret_cast<fftw_complex*> (Gtmp), FFTW_FORWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	double *fg1=new double[nbk*(nbk-2)];
	double *fg3=new double[nbk*(nbk-2)];
	
	int N0=nbw+1;
	
	
	// TF dans le temps de (f0*G)(k,tau) et multiplication par le prefacteur fonction de 1.0/(1.0-Usp*chi0/2) et 1.0/(1.0+Uch*chi0/2)
	// (on prend seulement la partie imaginaire, voir la section "Techniques de calcul")
	//#pragma omp parallel private(l,m,j)
	{
		double chitmp, chisp_tmp, chich_tmp;
		long int x, y, tmp, wn;
		dcomplex *ft_tmp1=new dcomplex[nbw];
		double *tau=new double[nbw+1];
		double *ftau=new double[nbw+1];
		double *coeffs=new double[4*nbw];
		double qn, d2S0, d2Sbeta, Htau0, Htau_beta;
		
		//#pragma omp for
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				for (j=0; j<=nbw; j++)
				{
					tau[j]=j/(nbw*tem);					
					ftau[j]=H[m+l*nbk+j*nk4th];
					//					if (l==2 && m==1 ) cout<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
				}
			    
				coeffs[0]=dtauf0G[m+l*nbk];
				coeffs[1]=dtauf0G[m+l*nbk+nk4th];				
				
//				spline_coeffs(tau, ftau, N0, coeffs);
				spline_coeffs_rel(tau, ftau, N0, coeffs);
				
				for (j=0; j<nbw; j++)
					ft_tmp1[j]=6.0*coeffs[4*j];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
				
				d2S0=2*coeffs[1];
				d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];				
//				d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
				
				fg1[m+l*nbk]=H[m+l*nbk]-H[m+l*nbk+nbw*nk4th];
				fg3[m+l*nbk]=d2S0-d2Sbeta;
				
				//				cout<<setw(10)<<l<<setw(10)<<m<<setw(20)<<H[m+l*nbk]<<setw(20)<<H[m+l*nbk+nbw*nk4th]<<setw(20)<<fg1[m+l*nbk]<<setw(20)<<d2S0<<setw(20)<<d2Sbeta<<setw(20)<<fg3[m+l*nbk]<<endl;
				//				cout<<setw(10)<<l<<setw(10)<<m<<fg1[m+l*nbk]<<endl;
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}

				Htau0=H[m+l*nbk];
				Htau_beta=H[m+l*nbk+nbw*nk4th];
				for (j=1; j<nbw/2; j++)
				{
					
					wn=j;
					if (wn>nbw/2) wn=nbw-wn;
					chitmp=chiqw_array[y+(x*(x+1))/2+wn*nk8th];
					
					chisp_tmp=1.0/(1.0-Usp*chitmp/2);
					chich_tmp=1.0/(1.0+Uch*chitmp/2);
					chitmp=3*Usp*chisp_tmp*chisp_tmp+Uch*chich_tmp*chich_tmp;
					
					qn=2*j*PI*tem;
					H[m+l*nbk+j*nk4th]=-2*chitmp*(-(Htau0-Htau_beta)/qn + (d2S0-d2Sbeta)/(qn*qn*qn) + imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j])/(qn*qn*qn*qn));					
					//					Htmp[m+l*nbk+j*nk4th]=-I*(H[m+l*nbk]-H[m+l*nbk+nbw*nk4th])/qn - (dtauf0G[m+l*nbk]-dtauf0G[m+l*nbk+nk4th])/(qn*qn)+ I*(d2S0-d2Sbeta)/(qn*qn*qn)+ ((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]/(qn*qn*qn*qn);					
					//					if (l==nbk-3 && m==nbk-1 && j<20) cout<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<endl;
					//					if (j<10) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<-2*imag(Htmp[m+l*nbk+j*nk4th])<<endl;
					//					cout<<setw(10)<<l+1<<setw(10)<<m<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
				}
				j=nbw/2;
				H[m+l*nbk+j*nk4th]=0;
				for (j=nbw/2+1; j<nbw; j++)
				{
					H[m+l*nbk+j*nk4th]=-H[m+l*nbk+(nbw-j)*nk4th];
					//					qn=2*(j-nbw)*PI*tem;
					//					Htmp[m+l*nbk+j*nk4th]=-I*(H[m+l*nbk]-H[m+l*nbk+nbw*nk4th])/qn - (dtauf0G[m+l*nbk]-dtauf0G[m+l*nbk+nk4th])/(qn*qn)+ I*(d2S0-d2Sbeta)/(qn*qn*qn)+ ((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]/(qn*qn*qn*qn);
				}
				H[m+l*nbk]=0;
			}
		delete [] ft_tmp1;
		delete [] tau;
		delete [] ftau;
		delete [] coeffs;
	}
	
	//	double *Htmp=new double[nbk*(nbk-2)*nbw];	
	//	for (l=0; l<nbk-2; l++)
	//		for (m=0; m<nbk; m++)
	//			for (j=0; j<nbw; j++)
	//				Htmp[m+l*nbk+j*nk4th]=H[m+l*nbk+j*nk4th];
	
	
	double *dtauH=new double[nbk*(nbk-2)];
	
	// TF dans le temps de H0(q,iq_m) et calcul de la derivee de H0(q,tau) a tau=0 et tau=beta
	//#pragma omp parallel private(l,m,j)
	{
		dcomplex *ft_tmp1=new dcomplex[nbw];
		long int wn, x, y, tmp;
		double t, qn, tmpsp, tmpch;
		dcomplex tmpnum;
		double zsp, zch;
		double cosbt, cosb, cost, sint, sinb, sinbt;
		double et, eb, ebtp, ebtm;
		double H0inf, dtauHtmp, dtauHinf;
		dcomplex Hinfqn;
		double chitmp, chisp_tmp, chich_tmp;
		//#pragma omp for
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}	
				
				dtauHtmp=0;
				ft_tmp1[0]=0;
				for (j=1; j<nbw/2; j++)
				{					
					qn=2*PI*j*tem;
					tmpsp=(-qn*qn+Usp*chi0_inf[y+(x*(x+1))/2]/2);
					tmpch=(-qn*qn-Uch*chi0_inf[y+(x*(x+1))/2]/2);
					tmpnum=-I*qn*qn*qn*fg1[m+l*nbk]+I*qn*fg3[m+l*nbk];
					Hinfqn=-6.0*Usp*tmpnum/(tmpsp*tmpsp)-2.0*Uch*tmpnum/(tmpch*tmpch);
					ft_tmp1[j]=tem*( I*H[m+l*nbk+j*nk4th]-Hinfqn);
					dtauHtmp+=2*qn*imag(ft_tmp1[j]);
					//					dtauHtmp+=2*qn*tem*H[m+l*nbk+j*nk4th];
					
					//					cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;					
					//					if (j>(nbw/2-20)) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;
					//					if (l==nbk-3 && m==nbk-1) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<setw(25)<<imag(Hinfqn)<<ft_tmp1[j]/tem<<endl;
					//					Htmp[m+l*nbk+j*nk4th]=imag(Hinfqn);
					//					Htmp[m+l*nbk+j*nk4th]=H[m+l*nbk+j*nk4th];
				}
				ft_tmp1[nbw/2]=0;
				for (j=nbw/2+1; j<nbw; j++)
					ft_tmp1[j]=-ft_tmp1[nbw-j];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
				
				zsp=sqrt(Usp*chi0_inf[y+(x*(x+1))/2]/2);
				zch=sqrt(Uch*chi0_inf[y+(x*(x+1))/2]/2);
				
				cosb=cos(zsp/(2.0*tem));
				sinb=sin(zsp/(2.0*tem));
				eb=exp(-zch/tem);
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					sint=sin(zsp*t);
					sinbt=sin((1.0/(2.0*tem)-t)*zsp);
					cosbt=cos((1.0/(2.0*tem)-t)*zsp);
					et=exp(-t*zch);
					ebtp=exp(-(1.0/tem+t)*zch);
					ebtm=exp(-(1.0/tem-t)*zch);
					tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -t*cosbt/sinb +sint/(2.0*tem*sinb*sinb) )/(4.0*zsp)
					-(fg1[m+l*nbk]/2)*sinbt/sinb;
					tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( t*(ebtm+et)/(1.0-eb) + (ebtp-ebtm)/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
					+(fg1[m+l*nbk]/2)*(ebtm-et)/(1.0-eb);
					H0inf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
					H[m+l*nbk+j*nk4th]=real(ft_tmp1[j])+H0inf;	
					//					if (l==nbk-3 && m==nbk-1 && j%(nbw/16)==0) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(25)<<H[m+l*nbk+j*nk4th]<<endl;
					//					H[m+l*nbk+j*nk4th]=H0inf;
				}
				j=nbw;
				t=1.0/tem;
				sint=sin(zsp*t);
				sinbt=sin((1.0/(2.0*tem)-t)*zsp);
				cosbt=cos((1.0/(2.0*tem)-t)*zsp);
				et=exp(-t*zch);
				ebtp=exp(-(1.0/tem+t)*zch);
				ebtm=exp(-(1.0/tem-t)*zch);
				tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -t*cosbt/sinb + sint/(2.0*tem*sinb*sinb) )/(4.0*zsp)
				-(fg1[m+l*nbk]/2)*sinbt/sinb;
				tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( t*(ebtm+et)/(1.0-eb) + (ebtp-ebtm)/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
				+(fg1[m+l*nbk]/2)*(ebtm-et)/(1.0-eb);
				H0inf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
				H[m+l*nbk+j*nk4th]=H0inf;
				
				// d/dtau H(q,tau) a tau=0 (les derivees sont egales a tau=0 et tau=beta);
				tmpsp=(-zsp*zsp*fg1[m+l*nbk]+fg3[m+l*nbk])*( -cosb/sinb+zsp/(2.0*tem*sinb*sinb) )/(4.0*zsp) 
				+ zsp*(fg1[m+l*nbk]/2)*cosb/sinb;
				tmpch=(zch*zch*fg1[m+l*nbk]+fg3[m+l*nbk])*( (eb+1.0)/(1.0-eb)-2.0*zch*eb/(tem*(1.0-eb)*(1.0-eb)) )/(4.0*zch)
				+ zch*(fg1[m+l*nbk]/2)*(eb+1.0)/(1.0-eb);
				dtauHinf=-6.0*Usp*tmpsp-2.0*Uch*tmpch;
				dtauH[m+l*nbk]=dtauHtmp+dtauHinf;
				
				//				dtauH[m+l*nbk]=dtauHinf;
				//				cout<<setw(25)<<dtauHtmp<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//				cout<<setw(25)<<dtauH[m+l*nbk]<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//				cout<<setw(25)<<dtauHinf<<setw(20)<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;				
			}
	}
	
	
	delete [] fg1;
	delete [] fg3;
	
	double *dtauGH=new double[(nbk-2)*nbk*2];
	
	//TF dans l'espace de H0, H0*G et de leur derivee par rapport a tau	
	//#pragma omp parallel private(fr_tmp1,j,l,m)
	{
		fr_tmp1=new double[(nbk-2)*nbk];
		long int x, y, tmp;
		
		//TF dans l'espace de H_0(q,tau)		
		//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					fr_tmp1[m+l*nbk]=normFact*H[m+l*nbk+j*nk4th];
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
					H[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
		}
		
		//TF dans l'espace de d/dtau H_0(q,tau) a tau=0		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				fr_tmp1[m+l*nbk]=normFact*dtauH[m+l*nbk];
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
				dtauH[m+l*nbk]=fr_tmp1[m+l*nbk];
		
		
		//TF dans l'espace de H_0(r,tau)G(r,-tau)		
		//#pragma omp for
		for (j=0; j<=nbw; j++)
		{
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=H[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					H[m+l*nbk+j*nk4th]=-fr_tmp1[m+l*nbk];
				}
			//			l=2;
			//			m=1;
			//			cout<<setw(10)<<j<<setw(10)<<l<<setw(10)<<m<<H[m+l*nbk+j*nk4th]<<endl;
		}
		
		//TF dans l'espace de d/dtau (H_0(r,tau)G(r,-tau)) a tau=0
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				fr_tmp1[m+l*nbk]=dtauH[m+l*nbk]*Grt[y+(x*(x+1))/2]+H[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				dtauGH[m+l*nbk]=-fr_tmp1[m+l*nbk];
				//				cout<<setw(25)<<dtauGH[m+l*nbk]<<nbw*tem*(H[m+l*nbk+nk4th]-H[m+l*nbk])<<endl;
			}
		
		//TF dans l'espace de d/dtau (H_0(r,tau)G(r,-tau)) a tau=beta
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				fr_tmp1[m+l*nbk]=dtauH[m+l*nbk]*Grt[y+(x*(x+1))/2+nbw*nk8th]+H[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
			}
		
		fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
		
		for (l=0; l<nbk-2; l++)
			for (m=0; m<nbk; m++)
			{
				dtauGH[m+l*nbk+nk4th]=-fr_tmp1[m+l*nbk];
				//				cout<<setw(25)<<dtauGH[m+l*nbk+nk4th]<<nbw*tem*(H[m+l*nbk+nbw*nk4th]-H[m+l*nbk+(nbw-1)*nk4th])<<endl;
			}
		
		
		delete [] fr_tmp1;
	}	
	
	Gtmp=new dcomplex[nbw];
	fftw_plan fftplan_t_back=fftw_plan_dft_1d(nbw, reinterpret_cast<fftw_complex*>(Gtmp), reinterpret_cast<fftw_complex*>(Gtmp), FFTW_BACKWARD, FFTW_MEASURE);
	delete [] Gtmp;
	
	
	if (chijj1_qn) delete [] chijj1_qn;
	if (chijj2_qn) delete [] chijj2_qn;
	if (chijj3_qn) delete [] chijj3_qn;
	
	chijj1_qn=new double[nw+1];
	chijj2_qn=new double[nw+1];
	chijj1_qn4=new double[nw];
	chijj2_qn6=new double[nw];
	chijj3_qn=new double[nw+1];
	
	double *chijj1_k=new double[(nw+1)*N_theta_k];
	double *chijj2_k=new double[(nw+1)*N_theta_k];
	double *chijj3_k=new double[(nw+1)*N_theta_k];
	
	for (j=0; j<(nw+1)*N_theta_k; j++)
	{
		chijj1_k[j]=0;
		chijj2_k[j]=0;
		chijj3_k[j]=0;
	}
	
	for (j=0; j<nw; j++)	
	{
		chijj1_qn[j]=0;
		chijj1_qn4[j]=0;
		chijj2_qn[j]=0;
		chijj2_qn6[j]=0;
		chijj3_qn[j]=0;
	}
	chijj1_qn[nw]=0;
	chijj2_qn[nw]=0;
	chijj3_qn[nw]=0;
	
	chijj1_w0=0;
	chijj1_inf2=0;
	chijj2_inf2=0;
	chijj2_inf4=0;
	
	dcomplex chijj3_qn0=0;
	
//	double Sigma_inf=G1->Sigma_inf;
	
	int NS0=nbw/2+1;
	
	cout<<setiosflags(ios::left)<<setprecision(10);
	
	dcomplex k0sum=0;
	
//	l0=k0[0], m0=k0[1];
	
	double Dtheta=PI/(2.0*N_theta_k);
	
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp, p;
		double wt, t, kn, G0t, Jinf;
		dcomplex k0, ph;
		
		//		dcomplex *J0tmp2=new dcomplex[nbw];
		
		dcomplex Jinf_tmp, Jbeta;
		dcomplex *J0tmp=new dcomplex[nbw];
		dcomplex *Jtmp=new dcomplex[nbw];
		dcomplex *HGtmp=new dcomplex[nbw];
		double CJ1, CJ2, CJ3, d2J0, d2Jbeta, CJ01, CJ02, CJ03;
		
		dcomplex *Gtmp1=new dcomplex[nbw];
		dcomplex *GGtmp1=new dcomplex[nbw];
		dcomplex Gtmp2, Gtmp0, Gbeta;
		
		double *coeffs=new double[4*(NS0-1)];
		double *tau=new double[NS0];
		double *ft_tmp=new double[NS0];
		double wn, integ, a, b, c, d, x1, x2;
		
		double *tau2=new double[nbw+1];
		double *ft_tmp2=new double[nbw+1];
		double *coeffs2=new double[4*nbw];
		
		dcomplex k0sum_loc=0, chijj2_qn0_loc;
		double *chijj1_qn_loc=new double[nw+1];
		double *chijj1_qn4_loc=new double[nw];
		double *chijj2_qn_loc=new double[nw+1];
		double *chijj2_qn6_loc=new double[nw];
		for (j=0; j<nw; j++)	
		{
			chijj1_qn_loc[j]=0;
			chijj1_qn4_loc[j]=0;
			chijj2_qn_loc[j]=0;
			chijj2_qn6_loc[j]=0;
		}
		chijj1_qn_loc[nw]=0;
		chijj2_qn_loc[nw]=0;
		
		dcomplex chijj3_qn0_loc=0;
		dcomplex chijj3_qn0_tmp;
		
		double chijj2_inf2_loc=0, chijj2_inf4_loc=0, chijj1_inf2_loc=0, chijj1_w0_loc=0;
		
		double dtau_J0, dtau_Jbeta;
		double dtau_Gkp, dtau_Gkm, dtau_chi;
		double C2, C3;
		dcomplex selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex exp1, exp2, exp3, exp4;
		dcomplex f1p, f2p, f3p, f4p, f1m, f2m, f3m, f4m;
		
		dcomplex exptmp;
		dcomplex b1, b2, b3, b4;
		dcomplex denj1, denj2, denj3, denj4;
		dcomplex expj1, expj2, expj3, expj4;
		dcomplex fj1, fj2, fj3, fj4, fj1b, fj2b, fj3b, fj4b;
		
		double theta_k;
		int ind_theta;
		double *chijj1_k_loc=new double[(nw+1)*N_theta_k];
		double *chijj2_k_loc=new double[(nw+1)*N_theta_k];
		
		for (j=0; j<(nw+1)*N_theta_k; j++)
		{
			chijj1_k_loc[j]=0;
			chijj2_k_loc[j]=0;
		}
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
		
#pragma omp for
		for (l=0; l<nbk-2; l++)
		{
			for (m=0; m<nbk; m++)
			{
				
				x=l+1;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
				
				// calcul du spline de hV(k,t)
				
				for (j=0; j<=nbw; j++)
				{
					tau2[j]=j/(nbw*tem);
					ft_tmp2[j]=hV[m+l*nbk+j*nk4th];
				}
				
				coeffs2[0]=dhV[m+l*nbk];
				coeffs2[1]=dhV[m+l*nbk+nk4th];				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);				
//				spline_coeffs(tau2, ft_tmp2, nbw+1, coeffs2);
				
				
				// TF sur t de hV(k,t)		
				for (j=0; j<nbw; j++)
					Gtmp1[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs2[4*j];
				
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				CJ1=-hV[m+l*nbk]-hV[m+l*nbk+nbw*nk4th];
				CJ2=dhV[m+l*nbk]+dhV[m+l*nbk+nk4th];
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
//				d2Jbeta=6.0*coeffs2[4*nbw-4]/tem+2.0*coeffs2[4*nbw-3];
				CJ3=-d2J0-d2Jbeta;
				
				for (j=0; j<nbw/2; j++)
				{
					kn=(2*j+1)*PI*tem;
					Jtmp[j+nbw/2]=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				for (j=nbw/2; j<nbw; j++)
				{
					kn=(2*(j-nbw)+1)*PI*tem;
					Jtmp[j-nbw/2]=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				
				// calcul du spline de f0V(k,t)
				
				for (j=0; j<=nbw; j++)
				{
					tau2[j]=j/(nbw*tem);
					ft_tmp2[j]=f0[m+l*nbk+j*nk4th];
				}
				
				coeffs2[0]=df0V[m+l*nbk];
				coeffs2[1]=df0V[m+l*nbk+nk4th];				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);				
//				spline_coeffs(tau2, ft_tmp2, nbw+1, coeffs2);
				
				// TF sur t de f0V(k,t)		
				for (j=0; j<nbw; j++)
					Gtmp1[j]=exp(dcomplex(0,j*PI/nbw))*((double)6.0)*coeffs2[4*j];
				
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				CJ01=-f0[m+l*nbk]-f0[m+l*nbk+nbw*nk4th];
				CJ02=df0V[m+l*nbk]+df0V[m+l*nbk+nk4th];
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];				
//				d2Jbeta=6.0*coeffs2[4*nbw-4]/tem+2.0*coeffs2[4*nbw-3];
				CJ03=-d2J0-d2Jbeta;
				
				for (j=0; j<nbw/2; j++)
				{
					kn=(2*j+1)*PI*tem;
					J0tmp[j+nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				for (j=nbw/2; j<nbw; j++)
				{
					kn=(2*(j-nbw)+1)*PI*tem;
					J0tmp[j-nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				
				// calcul du spline de HG(k,t)
				
				for (j=0; j<=nbw; j++)
				{
					tau2[j]=j/(nbw*tem);
					ft_tmp2[j]=H[m+l*nbk+j*nk4th];
				}
				
				coeffs2[0]=dtauGH[m+l*nbk];
				coeffs2[1]=dtauGH[m+l*nbk+nk4th];				
				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);
//				spline_coeffs(tau2, ft_tmp2, nbw+1, coeffs2);
				
				// TF sur t de HG(k,t)		
				for (j=0; j<nbw; j++)
					Gtmp1[j]=exp(dcomplex(0,-j*PI/nbw))*((double)6.0)*coeffs2[4*j];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				CJ01=H[m+l*nbk]+H[m+l*nbk+nbw*nk4th];
				CJ02=dtauGH[m+l*nbk]+dtauGH[m+l*nbk+nk4th];
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
//				d2Jbeta=6.0*coeffs2[4*nbw-4]/tem+2.0*coeffs2[4*nbw-3];
				CJ03=d2J0+d2Jbeta;
				
				for (j=0; j<nbw/2; j++)
				{
					kn=(2*j+1)*PI*tem;
					HGtmp[j+nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				for (j=nbw/2; j<nbw; j++)
				{
					kn=(2*(j-nbw)+1)*PI*tem;
					HGtmp[j-nbw/2]=-I*CJ01/kn-CJ02/(kn*kn)+I*CJ03/(kn*kn*kn)+((double)1.0-exp(dcomplex(0,-(2*(j-nbw)+1)*PI/nbw)))*Gtmp1[j]/(kn*kn*kn*kn);
				}
				
				//				if (l==l0-1 && m==m0)
				//				{
				//					for (j=0; j<nbw/2; j++)
				//						cout<<setw(40)<<U*HGtmp[wn0+nbw/2]/2.0<<endl;					
				//				}
				
				k0=0;
				dtau_Gkp=0;
				chijj2_qn0_loc=0;
				dtau_J0=0;
				chijj3_qn0_tmp=0;
				
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					chijj3_qn0_tmp+=Gtmp1[j]*Gtmp1[j]*HGtmp[j]/tem;
					
					Jtmp[j]=Gtmp1[j]*Jtmp[j];
					
					chijj2_qn0_loc+=Gtmp1[j]*Gtmp1[j]*J0tmp[j]/tem;
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					Jinf_tmp=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn);
					//					J0tmp2[j]=Jtmp[j]/tem;
					Jtmp[j]=Jtmp[j]-Gtmp2*Jinf_tmp;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);
					dtau_J0-=kn*imag(Jtmp[j]);
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
					
					chijj3_qn0_tmp+=Gtmp1[j]*Gtmp1[j]*HGtmp[j]/tem;
					
					Jtmp[j]=Gtmp1[j]*Jtmp[j];
					
					chijj2_qn0_loc+=Gtmp1[j]*Gtmp1[j]*J0tmp[j]/tem;
					
					//					Gtmp1[j]=tem/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp2=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+(l+1)*nbk]+mu+C2/(kn*kn));
					Gtmp1[j]=Gtmp1[j]-Gtmp2;
					
					Jinf_tmp=-I*CJ1/kn-CJ2/(kn*kn)+I*CJ3/(kn*kn*kn);
					//					J0tmp2[j]=Jtmp[j]/tem;
					Jtmp[j]=Jtmp[j]-Gtmp2*Jinf_tmp;
					
					k0+=Gtmp1[j];
					dtau_Gkp+=kn*imag(Gtmp1[j]);	
					dtau_J0-=kn*imag(Jtmp[j]);
				}
				
				// calcul du G(k,t) asymptotique et du J(k,t) asymtotique	
				coeff_pol[1]=-ek[m+(l+1)*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];				
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				b1=(CJ1*r1*r1+CJ2*r1+CJ3)/((r1-r2)*(r1-r3)*(r1-r4));
				b2=(CJ1*r2*r2+CJ2*r2+CJ3)/((r2-r1)*(r2-r3)*(r2-r4));
				b3=(CJ1*r3*r3+CJ2*r3+CJ3)/((r3-r1)*(r3-r2)*(r3-r4));
				b4=(CJ1*r4*r4+CJ2*r4+CJ3)/((r4-r1)*(r4-r2)*(r4-r3));				
				
				if (real(r1)>0) 
				{
					exptmp=exp(-r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1/den1;
					f1m=a1*exptmp/den1;
					fj1=r1*b1*exptmp/den1;
					fj1b=r1*b1/den1;
				}
				else 
				{
					exptmp=exp(r1/tem);
					den1=exptmp+(double)1.0;
					f1p=a1*exptmp/den1;
					f1m=a1/den1;
					fj1=r1*b1/den1;
					fj1b=r1*b1*exptmp/den1;
				}
				if (real(r2)>0) 
				{
					exptmp=exp(-r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2/den2;
					f2m=a2*exptmp/den2;
					fj2=r2*b2*exptmp/den2;
					fj2b=r2*b2/den2;
				}
				else 
				{	
					exptmp=exp(r2/tem);
					den2=exptmp+(double)1.0;
					f2p=a2*exptmp/den2;
					f2m=a2/den2;
					fj2=	r2*b2/den2;
					fj2b=r2*b2*exptmp/den2;
				}
				if (real(r3)>0) 
				{
					exptmp=exp(-r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3/den3;
					f3m=a3*exptmp/den3;
					fj3=r3*b3*exptmp/den3;
					fj3b=r3*b3/den3;
				}
				else 
				{
					exptmp=exp(r3/tem);
					den3=exptmp+(double)1.0;
					f3p=a3*exptmp/den3;
					f3m=a3/den3;
					fj3=r3*b3/den3;
					fj3b=r3*b3*exptmp/den3;
				}
				if (real(r4)>0) 
				{
					exptmp=exp(-r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4/den4;
					f4m=a4*exptmp/den4;
					fj4=r4*b4*exptmp/den4;
					fj4b=r4*b4/den4;
				}
				else 
				{
					exptmp=exp(r4/tem);
					den4=exptmp+(double)1.0;
					f4p=a4*exptmp/den4;
					f4m=a4/den4;
					fj4=r4*b4/den4;
					fj4b=r4*b4*exptmp/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);
				
				dtau_Gkm=-dtau_Gkp;
				dtau_Gkp+=real(r1*f1p+r2*f2p+r3*f3p+r4*f4p);
				dtau_Gkm+=real(r1*f1m+r2*f2m+r3*f3m+r4*f4m);		
				
				dtau_Jbeta=-dtau_J0;
				dtau_J0+=real(fj1+fj2+fj3+fj4);
				dtau_Jbeta+=real(fj1b+fj2b+fj3b+fj4b);
				
				
				//TF de G(k,ikn) dans le temps
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				Gbeta=Gtmp1[0];
				
				// TF de J(k,ikn) dans le temps
				fftw_execute_dft(fftplan_t_back,reinterpret_cast<fftw_complex*>(Jtmp),reinterpret_cast<fftw_complex*>(Jtmp));
				Jbeta=Jtmp[0];
				
				for (j=0; j<nbw; j++)
				{
					t=j/(nbw*tem);
					
					if (real(r1)>0) 
					{
						exp1=exp(-t*r1);
						expj1=exp(r1*(t-1.0/tem));
					}
					else 
					{
						exp1=exp((1.0/tem-t)*r1);
						expj1=exp(r1*t);
					}
					if (real(r2)>0) 
					{
						exp2=exp(-t*r2);
						expj2=exp(r2*(t-1.0/tem));
					}
					else 
					{
						exp2=exp((1.0/tem-t)*r2);
						expj2=exp(r2*t);
					}
					if (real(r3)>0) 
					{
						exp3=exp(-t*r3);
						expj3=exp(r3*(t-1.0/tem));
					}
					else 
					{
						exp3=exp((1.0/tem-t)*r3);
						expj3=exp(r3*t);
					}
					if (real(r4)>0) 
					{
						exp4=exp(-t*r4);
						expj4=exp(r4*(t-1.0/tem));
					}
					else 
					{
						exp4=exp((1.0/tem-t)*r4);
						expj4=exp(r4*t);
					}
					
					G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
					Jinf=real(b1*expj1/den1+b2*expj2/den2+b3*expj3/den3+b4*expj4/den4);
					
					ph=exp(dcomplex(0.0, (j*PI*(nbw-1))/nbw));
					
					Gtmp1[j]=ph*Gtmp1[j]+G0t;					
					
					Jtmp[j]=Jtmp[j]/ph+Jinf;
				}			
				
				t=1.0/tem;
				
				if (real(r1)>0) 
				{
					exp1=exp(-t*r1);
					expj1=1.0;
				}
				else 
				{
					exp1=1.0;
					expj1=exp(r1*t);
				}
				if (real(r2)>0) 
				{
					exp2=exp(-t*r2);
					expj2=1.0;
				}
				else 
				{
					exp2=1.0;
					expj2=exp(r2*t);
				}
				if (real(r3)>0) 
				{
					exp3=exp(-t*r3);
					expj3=1.0;
				}
				else 
				{
					exp3=1.0;
					expj3=exp(r3*t);
				}
				if (real(r4)>0) 
				{
					exp4=exp(-t*r4);
					expj4=1.0;
				}
				else 
				{
					exp4=1.0;
					expj4=exp(r4*t);
				}
				
				G0t=real(-a1*exp1/den1-a2*exp2/den2-a3*exp3/den3-a4*exp4/den4);
				Jinf=real(b1*expj1/den1+b2*expj2/den2+b3*expj3/den3+b4*expj4/den4);
				
				ph=exp(dcomplex(0.0, PI*(nbw-1)));
				
				Gbeta=ph*Gbeta+G0t;					
				
				Jbeta=Jbeta/ph+Jinf;
				
				// calcul du spline et de la TF sur t de G(k,t)J(k,-t) 				
				for (j=0; j<nbw; j++)
				{
					tau2[j]=j/(nbw*tem);					
					ft_tmp2[j]=real(Gtmp1[j]*Jtmp[j]);
				}
				tau2[nbw]=1.0/tem;
				ft_tmp2[nbw]=real((-Gtmp1[0]-(double)1.0)*Jbeta);
				
				coeffs2[0]=real(dtau_Gkp*Jtmp[0]+Gtmp1[0]*dtau_J0);
				coeffs2[1]=real(dtau_Gkm*Jbeta - (Gtmp1[0]+(double)1.0)*dtau_Jbeta);
				
				CJ1=ft_tmp2[0]-ft_tmp2[nbw];

				spline_coeffs_rel(tau2, ft_tmp2, nbw+1, coeffs2);
//				spline_coeffs(tau2, ft_tmp2, nbw+1, coeffs2);
				
				d2J0=2*coeffs2[1];
				d2Jbeta=6.0*coeffs2[4*nbw-4]/(nbw*tem)+2.0*coeffs2[4*nbw-3];
//				d2Jbeta=6.0*coeffs2[4*nbw-4]/tem+2.0*coeffs2[4*nbw-3];
				CJ3=d2J0-d2Jbeta;
				
				for (j=0; j<nbw; j++)
				{
					J0tmp[j]=6.0*coeffs2[4*j];
				}
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(J0tmp),reinterpret_cast<fftw_complex*>(J0tmp));
				
				for (j=1; j<=nbw/2; j++)
				{
					wn=2*j*PI*tem;
					Jtmp[j]=-((double)2.0)*I*CJ1/wn+((double)2.0)*I*CJ3/(wn*wn*wn)+((double)2.0)*I*imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*J0tmp[j]/(wn*wn*wn*wn));
				}
				
				// Definition de -G(k,t)G(k,-t)				
				GGtmp1[0]=Gtmp1[0]*(Gtmp1[0]+(double)1.0);
				for (j=1; j<nbw; j++)
					GGtmp1[j]=-Gtmp1[j]*Gtmp1[nbw-j];
				
				// calcul du spline de -G(k,t)G(k,-t)			
				for (j=0; j<NS0; j++)
				{
					tau[j]=j/(nbw*tem);
					ft_tmp[j]=real(GGtmp1[j]);
				}
				dtau_chi=dtau_Gkp*(real(Gtmp1[0])+1.0)+real(Gtmp1[0])*dtau_Gkm;
				
				//				if ((l+1)%4==0)
				//					cout<<"dtau_chi:  "<<setw(20)<<(real(GGtmp1[1]-GGtmp1[0])*(tem*nbw)-dtau_chitmp)/dtau_chitmp<<setw(20)<<(dtau_chi-dtau_chitmp)/dtau_chitmp<<'\n';
				
				coeffs[0]=dtau_chi;
				coeffs[1]=0;				
				spline_coeffs_rel(tau, ft_tmp, NS0, coeffs);				
//				spline_coeffs(tau, ft_tmp, NS0, coeffs);
				
				// TF sur t de -G(k,t)G(k,-t)			
				for (j=0; j<nbw/2; j++)
					Gtmp1[j]=6.0*coeffs[4*j];
				for (j=nbw/2; j<nbw; j++)
					Gtmp1[j]=-Gtmp1[nbw-j-1];
				
				fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(Gtmp1),reinterpret_cast<fftw_complex*>(Gtmp1));
				
				integ=0;
				for (j=0; j<nbw/2; j++)
				{
					a=coeffs[4*j];
					b=coeffs[4*j+1];
					c=coeffs[4*j+2];
					d=coeffs[4*j+3];
					x1=tau[j+1]-tau[j];
					integ+=a*x1*x1*x1*x1/4.0+b*x1*x1*x1/3.0+c*x1*x1/2.0+d*x1;
					
//					x1=tau[j];
//					x2=tau[j+1];
//					integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
				}
				
				Gtmp1[0]=2*integ;
				
				wt=4;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+(l+1)*nbk]*k0;
				chijj1_qn_loc[0]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_w0_loc-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj1_inf2_loc+=4.0*dtau_chi*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk];
				chijj2_inf2_loc+=normFact*wt*dek[m+l*nbk]*2.0*CJ1;
				chijj2_inf4_loc-=normFact*wt*dek[m+l*nbk]*2.0*CJ3;
				
				chijj3_qn0_loc+=U*dek[m+l*nbk]*normFact*wt*chijj3_qn0_tmp/((double)2.0);
				
				chijj2_qn_loc[0]-=normFact*wt*dek[m+l*nbk]*real(chijj2_qn0_loc);
				
				theta_k=atan( (1.0*(nk0/2-m))/(nk0/2-l-1) );
				ind_theta=(int)floor(theta_k/Dtheta);
				
//				cout<<setiosflags(ios::left)<<setprecision(10);
//				cout<<setw(10)<<l+1<<setw(10)<<m<<"theta_k:  "<<setw(20)<<theta_k<<"ind_theta:  "<<ind_theta<<endl;
				
				chijj1_k_loc[ind_theta]-=2.0*wt*normFact*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[0]);
				chijj2_k_loc[ind_theta]-=normFact*wt*dek[m+l*nbk]*real(chijj2_qn0_loc);
				
				for (j=1; j<=nw; j++)
				{
					wn=2*j*PI*tem;
					
					chijj1_qn4_loc[j-1]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]);
					
					chijj2_qn_loc[j]-=normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
					
					chijj2_k_loc[ind_theta+j*N_theta_k]-=normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
							 
//					if (p<nk_chijj)	chijj2_k[p + j*nk_chijj]=-normFact*wt*dek[m+l*nbk]*imag(Jtmp[j])/wn;
					
					chijj2_qn6_loc[j-1]-=normFact*wt*dek[m+l*nbk]*2.0*imag(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*J0tmp[j])/wn;
					
					Gtmp1[j]=-2.0*dtau_chi/(wn*wn)+real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*Gtmp1[j]/(wn*wn*wn*wn));
					
					chijj1_qn_loc[j]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
					
					chijj1_k_loc[ind_theta+j*N_theta_k]-=2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
					
//					if (p<nk_chijj)	chijj1_k[p + j*nk_chijj]=-2.0*normFact*wt*dek[m+l*nbk]*dek[m+l*nbk]*real(Gtmp1[j]);
				}				

			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
			
			for (j=0; j<nw; j++)
			{
				chijj1_qn[j]+=chijj1_qn_loc[j];
				chijj1_qn4[j]+=chijj1_qn4_loc[j];
				chijj2_qn[j]+=chijj2_qn_loc[j];
				chijj2_qn6[j]+=chijj2_qn6_loc[j];
				for (p=0; p<N_theta_k; p++)
				{
					chijj1_k[p+j*N_theta_k]+=chijj1_k_loc[p+j*N_theta_k];
					chijj2_k[p+j*N_theta_k]+=chijj2_k_loc[p+j*N_theta_k];
				}
			}
			chijj1_qn[nw]+=chijj1_qn_loc[nw];
			chijj2_qn[nw]+=chijj2_qn_loc[nw];
			for (p=0; p<N_theta_k; p++)
			{
				chijj1_k[p+nw*N_theta_k]+=chijj1_k_loc[p+nw*N_theta_k];
				chijj2_k[p+nw*N_theta_k]+=chijj2_k_loc[p+nw*N_theta_k];
			}
			
			chijj1_w0+=chijj1_w0_loc;
			chijj1_inf2+=chijj1_inf2_loc;
			chijj2_inf2+=chijj2_inf2_loc;
			chijj2_inf4+=chijj2_inf4_loc;
			chijj3_qn0+=chijj3_qn0_loc;
		}
		delete [] chijj1_qn4_loc;
		delete [] chijj1_qn_loc;
		delete [] chijj2_qn6_loc;
		delete [] chijj2_qn_loc;
		delete [] Gtmp1;
		delete [] GGtmp1;
		delete [] tau;
		delete [] ft_tmp;
		delete [] coeffs;
		delete [] tau2;
		delete [] ft_tmp2;
		delete [] coeffs2;
		delete [] Jtmp;
		delete [] J0tmp;
		delete [] HGtmp;
		delete [] chijj1_k_loc;
		delete [] chijj2_k_loc;		
	}
	
	delete [] hV;
	delete [] f0;
	delete [] H;
	delete [] dhV;
	delete [] df0V;
	delete [] dtauGH;
	delete [] dtauf0G;
	delete [] df0;
	delete [] dtauH;
	
	double *chijj1_k_tot=new double[nw+1];
	double *chijj2_k_tot=new double[nw+1];
	
	for (j=0; j<=nw; j++)
	{
		chijj1_k_tot[j]=0;
		chijj2_k_tot[j]=0;
	}
	
	int p;
	
	for (j=0; j<=nw; j++)
		for  (p=0; p<N_theta_k; p++)
		{
			chijj1_k_tot[j]+=chijj1_k[p + j*N_theta_k];
			chijj2_k_tot[j]+=chijj2_k[p + j*N_theta_k];
		}
	
	fstream file, file2, file3, file4, file_info;
	char name[200];
	
	const char *nameForm_chijj1_k="./chijj1_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj1_k_tot="./chijj1_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj2_k="./chijj2_k_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	const char *nameForm_chijj2_k_tot="./chijj2_k_tot_qn_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	
	sprintf(name, nameForm_chijj1_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj1_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file2.open(name, ios::out );
	file2<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj2_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file3.open(name, ios::out );
	file3<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	sprintf(name, nameForm_chijj2_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out );
	file4<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	for (j=0; j<=nw; j++)
	{
		file2<<setw(10)<<j<<setw(30)<<chijj1_k_tot[j]<<setw(30)<<chijj1_qn[j]<<chijj1_k_tot[j]/chijj1_qn[j]<<'\n';
		file4<<setw(10)<<j<<setw(30)<<chijj2_k_tot[j]<<setw(30)<<chijj2_qn[j]<<chijj2_k_tot[j]/chijj2_qn[j]<<'\n';
		for  (p=0; p<N_theta_k; p++)
		{
			file<<setw(10)<<j<<setw(10)<<p+1<<chijj1_k[p + j*N_theta_k]<<'\n';
			file3<<setw(10)<<j<<setw(10)<<p+1<<chijj2_k[p + j*N_theta_k]<<'\n';			
//			file<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj1_k[p + j*nk_chijj]<<'\n';
//			file3<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj2_k[p + j*nk_chijj]<<'\n';
		}
	}
	file.close();
	file2.close();
	file3.close();
	file4.close();
	
	delete [] chijj1_k_tot;
	delete [] chijj2_k_tot;
	
/*	
	{
		const char *nameForm_self_k="./self_k_ikn_%4.2f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
		sprintf(name, nameForm_self_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
		file.open(name, ios::out );
		file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
		dcomplex selftmp;
		long int x, y, tmp;
		double C2, C3, kn;
		for (j=0; j<=nw/4; j++)
		{
			kn=(2.0*j+1)*PI*tem;
			for (p=0; p<nk_chijj; p++)
			{
				x=k_chijj[2*p];
				y=k_chijj[2*p+1];
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];
				
				selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
				
				file<<setw(10)<<j<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<setw(25)<<real(selftmp)<<imag(selftmp)<<endl;
			}
		}
		file.close();
	}
*/
	
#pragma omp parallel private(m,j)
	{
		long int x, y, tmp;
		double wt, kn;
		dcomplex k0;
		dcomplex Gtmp2, Gtmp0;
		
		dcomplex k0sum_loc=0;
		
		double C2, C3;
		dcomplex z, selftmp;
		dcomplex coeff_pol[5];
		dcomplex roots[4];	
		dcomplex r1,r2,r3,r4;
		dcomplex a1, a2, a3, a4;
		dcomplex den1, den2, den3, den4;
		dcomplex f1m, f2m, f3m, f4m;
		
		coeff_pol[0]=1.0;
		coeff_pol[2]=-Sigma_inf;
		
#pragma omp for
		for (l=0; l<nbk; l+=nk0/2)
		{
			for (m=0; m<nbk; m++)
			{
				x=l;
				if (x>nk0/2) x=nk0-x;
				y=m;
				if (y>nk0/2) y=nk0-y;
				if (y>x)
				{
					tmp=x;
					x=y;
					y=tmp;
				}
				
//				C2=G1->Sigma_inf2[y+(x*(x+1))/2];
//				C3=G1->Sigma_inf3[y+(x*(x+1))/2];				
				
				C2=Sigma_inf2[y+(x*(x+1))/2];
				C3=Sigma_inf3[y+(x*(x+1))/2];				
				
				k0=0;
				
							
				for (j=0; j<nw; j++)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-conj(Self_kw_array[y+(x*(x+1))/2+(nw-j-1)*nk8th]));
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				for (j=nbw-1; j>=nw; j--)
				{
					kn=(2.0*(j-nw)+1)*PI*tem;
					
					selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th])/(kn*kn*kn*kn));
					
					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-selftmp);
					
					//					Gtmp2=tem/(dcomplex(0.0, kn)-ek[m+l*nbk]+mu-Self_kw_array[y+(x*(x+1))/2+(j-nw)*nk8th]);
					
					Gtmp0=tem/(dcomplex(0.0, kn+Sigma_inf/kn-C3/(kn*kn*kn))-ek[m+l*nbk]+mu+C2/(kn*kn));
					Gtmp2=Gtmp2-Gtmp0;
					
					k0+=Gtmp2;
				}
				
				coeff_pol[1]=-ek[m+l*nbk]+mu;
				coeff_pol[3]=-C2;
				coeff_pol[4]=-C3;
				find_roots_4pol(coeff_pol, roots);
				r1=roots[0];
				r2=roots[1];
				r3=roots[2];
				r4=roots[3];
				
				a1=r1*r1*r1/((r1-r2)*(r1-r3)*(r1-r4));
				a2=r2*r2*r2/((r2-r1)*(r2-r3)*(r2-r4));
				a3=r3*r3*r3/((r3-r1)*(r3-r2)*(r3-r4));
				a4=r4*r4*r4/((r4-r1)*(r4-r2)*(r4-r3));
				
				if (real(r1)>0) 
				{
					den1=exp(-r1/tem)+(double)1.0;
					f1m=a1*exp(-r1/tem)/den1;
				}
				else 
				{
					den1=exp(r1/tem)+(double)1.0;
					f1m=a1/den1;
				}
				if (real(r2)>0) 
				{
					den2=exp(-r2/tem)+(double)1.0;
					f2m=a2*exp(-r2/tem)/den2;
				}
				else 
				{	
					den2=exp(r2/tem)+(double)1.0;
					f2m=a2/den2;
				}
				if (real(r3)>0) 
				{
					den3=exp(-r3/tem)+(double)1.0;
					f3m=a3*exp(-r3/tem)/den3;
				}
				else 
				{
					den3=exp(r3/tem)+(double)1.0;
					f3m=a3/den3;
				}
				if (real(r4)>0) 
				{
					den4=exp(-r4/tem)+(double)1.0;
					f4m=a4*exp(-r4/tem)/den4;
				}
				else 
				{
					den4=exp(r4/tem)+(double)1.0;
					f4m=a4/den4;
				}
				
				k0+=real(f1m+f2m+f3m+f4m);				
				
				wt=2;
				if (m==0 || m==nbk-1) wt=wt/2;
				
				k0sum_loc-=2.0*wt*normFact*d2ek[m+l*nbk]*k0;
			}
		}
#pragma omp critical
		{
			k0sum+=k0sum_loc;
		}
	}
	
	//	cout<<"chitmp3:  ";
	//	calc_chijj_vertex_corr2_wn(0, k0, 0, Htmp, NULL, NULL, NULL, NULL);
	
	delete [] d2ek;
	
	//	for (j=0; j<4; j++)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	//	for (j=nbw/8; j<=nbw/2; j+=nbw/8)
	//	{
	//		cout<<setw(10)<<j<<setw(40)<<chijj2_qn[j];
	//		calc_chijj_vertex_corr1_wn(j);
	//	}
	
	
	kx0=real(k0sum);
	sprintf(name, "kx0_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out);
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific)<<kx0<<'\n';
	file.close();
	
	sprintf(name, "kx0_bin_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat", (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	file.write((char*)&kx0, sizeof(kx0));
	file.close();
	
	for (j=0; j<=nw; j++)
	{
		if (chijj1_qn[j]<0.0)
		{
			cout<<"attention! valeur non-physique de chi_jj\n";
			cout<<setw(30)<<chijj1_qn[j]<<'\n';
			break;
		}
	}
	

	//pour enregistrer chijj(iq_n) en binaire
	
	
	const char *nameForm_chijj1_qn_bin="./chijj1_qn_bin_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out | ios::binary);
	
	const char *nameForm_chijj2_qn_bin="./chijj2_qn_bin_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj2_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out | ios::binary);
	
	
	j=0;
	file.write((char*)&j,sizeof(j));
	file.write((char*)&(chijj1_qn[j]),sizeof(double));
	file4.write((char*)&j,sizeof(j));
	file4.write((char*)&(chijj2_qn[j]),sizeof(double));
	for (j=1; j<=nbw/2; j++)
	{
		file.write((char*)&j,sizeof(j));
		file.write((char*)&(chijj1_qn[j]),sizeof(double));
		file4.write((char*)&j,sizeof(j));
		file4.write((char*)&(chijj2_qn[j]),sizeof(double));
	}
	file.close();
	file4.close();
	// fin de l'enregistrement en binaire
	
	// pour enregistrer chijj(iq_n) en ascii
	
	const char *nameForm_chijj1_qn="./chijj1_qn_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj1_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file.open(name, ios::out );
	file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
	
	const char *nameForm_chijj2_qn="./chijj2_qn_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d.dat";
	sprintf(name, nameForm_chijj2_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw);
	file4.open(name, ios::out );
	file4<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);	
	
	j=0;
	file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
	file4<<setw(10)<<j<<chijj2_qn[j]<<'\n';
	for (j=1; j<=nbw/2; j++)
	{
		file<<setw(10)<<j<<chijj1_qn[j]<<'\n';
		file4<<setw(10)<<j<<chijj2_qn[j]<<'\n';
	}
	file.close();
	file2.close();
	file3.close();
	file4.close();
	// fin de l'enregistrement en ascii
	
	chijj3_qn[0]=real(chijj3_qn0);
	
	cout<<"-k0:  "<<-kx0<<'\n';
	cout<<"chijj1(iqn=0):  "<<chijj1_qn[0]<<'\n';
	cout<<"chijj2(iqn=0):  "<<chijj2_qn[0]<<'\n';
	cout<<"chijj3(iqn=0):  "<<chijj3_qn0<<endl;
	cout<<"chijj1(iqn=0)+chijj2(iqn=0):  "<<chijj1_qn[0]+chijj2_qn[0]<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+chijj2_qn[0]+kx0)/kx0<<'\n';
	cout<<"chijj1(iqn=0)+chijj2(iqn=0)+chijj3(iqn=0):  "<<chijj1_qn[0]+chijj2_qn[0]+real(chijj3_qn0)<<'\n';
	cout<<"difference relative avec -k0:  "<<(chijj1_qn[0]+chijj2_qn[0]+real(chijj3_qn0)+kx0)/kx0<<'\n';
	
	
	N0=nbw/2+1;
	
	double *chijj_qn=new double[nbw/2+1];
	
	for (j=0; j<=nbw/2; j++)
		chijj_qn[j]=chijj1_qn[j]+chijj2_qn[j];
	
	
	long int nn;
	
	if (chi3)
	{
		m=mmin;
		int N=m*((N0-1)/((int)pow(2.0,(int)m+1)))+(N0-1)/((int)pow(2.0,(int)m))+1;
		
		cout<<"2^mmin:  "<<pow(2.0,(int)m)<<"  (N0-1)/((int)pow(2.0,m)):  "<<(N0-1)/((int)pow(2.0,(int)m))<<"   N:  "<<N<<'\n';
		
		long int *ind_freq=new long int[N];
		
		int j0, jf, d;
		
		jf=(N0-1)/((int)pow(2.0,(int)m));
		for (j=0; j<=jf; j++)
		{
			ind_freq[j]=j;
			//				cout<<setw(10)<<j<<ind_freq[j]<<'\n';
		}
		
		int Nj=(N0-1)/((int)pow(2.0,(int)m+1));
		
		long int n=jf;
		d=1;
		for (l=0; l<m; l++)
		{
			j0=jf+1;
			jf=j0+Nj-1;
			d=2*d;
			for (j=j0; j<=jf; j++)
			{
				n=n+d;
				ind_freq[j]=n;
				//				cout<<setw(10)<<j<<ind_freq[j]<<'\n';				
			}
		}
		
		
		double *dtauhG=new double[(nbk-2)*nbk*2];
		
		//TF dans l'espace de h(r,t)G(r,-t) et de d/dtau (h(r,t)G(r,-t))
		{
			fr_tmp1=new double[(nbk-2)*nbk];
			long int x,y,tmp;
			
			//  d/dt h(r,t)G(r,-t)  a t=0	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dh[m+l*nbk]*Grt[y+(x*(x+1))/2]+h[m+l*nbk]*dtauGrt[y+(x*(x+1))/2];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauhG[m+l*nbk]=fr_tmp1[m+l*nbk];
				}
			
			//  d/dt h(r,t)G(r,-t)  a t=beta	
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					x=l+1;
					if (x>nk0/2) x=nk0-x;
					y=m;
					if (y>nk0/2) y=nk0-y;
					if (y>x)
					{
						tmp=x;
						x=y;
						y=tmp;
					}
					fr_tmp1[m+l*nbk]=dh[m+l*nbk+nk4th]*Grt[y+(x*(x+1))/2+nbw*nk8th]+h[m+l*nbk+nbw*nk4th]*dtauGrt[y+(x*(x+1))/2+nk8th];
				}
			
			fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
			
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauhG[m+l*nbk+nk4th]=fr_tmp1[m+l*nbk];
				}
			delete [] fr_tmp1;
		}
		
		//#pragma omp parallel private(fr_tmp1,j,l,m)
		{
			fr_tmp1=new double[(nbk-2)*nbk];
			long int x,y,tmp;			
			//#pragma omp for
			for (j=0; j<=nbw; j++)
			{
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=h[m+l*nbk+j*nk4th]*Grt[y+(x*(x+1))/2+j*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						h[m+l*nbk+j*nk4th]=fr_tmp1[m+l*nbk];
			}
			delete [] fr_tmp1;
		}
		
//		cout<<" tem: "<<tem<<"\n nbw: "<<nbw<<"\n nbk: "<<nbk<<"\n nk4th: "<<nk4th<<endl;
//verifier si les derivees sont bonnes		
//		for (l=0; l<nbk-2; l++)
//			for (m=0; m<nbk; m++)
//			{
//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauhG[m+l*nbk]<<tem*nbw*(h[m+l*nbk+nk4th]-h[m+l*nbk])<<endl;
//			}
//		for (l=0; l<nbk-2; l++)
//			for (m=0; m<nbk; m++)
//			{
//				cout<<setw(10)<<l<<setw(10)<<m<<setw(25)<<dtauhG[m+l*nbk+nk4th]<<tem*nbw*(h[m+l*nbk+nbw*nk4th]-h[m+l*nbk+(nbw-1)*nk4th])<<endl;
//			}
		 	
		
		double *hg2=new double[nbk*(nbk-2)];
		
		
		//TF dans le temps de (h(r,t)G(r,-t))(q) (voir les notes "techniques de calcul")
#pragma omp parallel private(l,m,j)
		{
			double chitmp, chisp_tmp1, chisp_tmp2, chich_tmp1, chich_tmp2;
			long int x, y, tmp, wn;
			dcomplex *ft_tmp1=new dcomplex[nbw];
			double *tau=new double[nbw+1];
			double *ftau=new double[nbw+1];
			double *coeffs=new double[4*nbw];
			double integ, a, b, c, d, x1, x2;
			
#pragma omp for
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					for (j=0; j<=nbw; j++)
					{
						tau[j]=j/(nbw*tem);					
						ftau[j]=h[m+l*nbk+j*nk4th];
						//					if (l==2 && m==1 ) cout<<setw(10)<<j<<H[m+l*nbk+j*nk4th]<<endl;
					}
					
					coeffs[0]=dtauhG[m+l*nbk];
					coeffs[1]=dtauhG[m+l*nbk+nk4th];				
					spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
//					spline_coeffs(tau, ftau, nbw+1, coeffs);
					
					for (j=0; j<nbw; j++)
						ft_tmp1[j]=6.0*coeffs[4*j];
					
					fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
					
					integ=0;
					for (j=0; j<nbw; j++)
					{
						a=coeffs[4*j];
						b=coeffs[4*j+1];
						c=coeffs[4*j+2];
						d=coeffs[4*j+3];
						x1=tau[j+1]-tau[j];
						integ+=a*x1*x1*x1*x1/4.0+b*x1*x1*x1/3.0+c*x1*x1/2.0+d*x1;
						
//						x1=tau[j];
//						x2=tau[j+1];
//						integ+=a*(x2*x2*x2*x2-x1*x1*x1*x1)/4.0+b*(x2*x2*x2-x1*x1*x1)/3.0+c*(x2*x2-x1*x1)/2.0+d*(x2-x1);
					}
					h[m+l*nbk]=2*integ;
					
					//					cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<0<<h[m+l*nbk]<<endl;
					
					hg2[m+l*nbk]=2*(dtauhG[m+l*nbk]-dtauhG[m+l*nbk+nk4th]);
					for (j=1; j<nbw/2; j++)
					{
						h[m+l*nbk+j*nk4th]=2*real(((double)1.0-exp(dcomplex(0,-2*j*PI/nbw)))*ft_tmp1[j]);
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<h[m+l*nbk+j*nk4th]<<endl;
					}
					
				}
			delete [] ft_tmp1;
			delete [] tau;
			delete [] ftau;
			delete [] coeffs;
		}
		
		delete [] dtauhG;
		
		//		double *Sigma_inf4=G1->Sigma_inf4;
		
		dcomplex *HG=new dcomplex[array_size_h];
		dcomplex *HG2=new dcomplex[array_size_h];		
//		dcomplex *HG=new dcomplex[nbk*(nbk-2)*(nbw+1)];
//		dcomplex *HG2=new dcomplex[nbk*(nbk-2)*(nbw+1)];
		dcomplex *dtauHn=new dcomplex[nbk*(nbk-2)*2];
		dcomplex *dtauHn2=new dcomplex[nbk*(nbk-2)*2];
		
		//		dcomplex *HG0tmp=new dcomplex[nbk*(nbk-2)*(nbw+1)];
		//		dcomplex *HGsum=new dcomplex[nbk*(nbk-2)*(nbw+1)];
		
		double qn;
		int nc;
		
		dcomplex *chijj3tmp=new dcomplex[N];
		for (j=0; j<N; j++) chijj3tmp[j]=0;
		
		//		n=1;
		//		for (n=1; n<N; n+=N/8)
		for (n=1; n<N; n++)
		{
			nn=ind_freq[n];
			qn=2.0*nn*PI*tem;
			
			nc=nn/2;
			//			nc=nn/4;
			//			cout<<"nc:  "<<nc<<endl;
			
			//definition de I_n(q,iqm)
#pragma omp parallel private(l,m,j)
			{
				double HGtmp1, HGtmp2, wn;
				double chi0tmp, chi0n, chisp_tmp,chich_tmp,chisp_n,chich_n;
				double D1, D2;
				long int x, y, tmp, p, q;
				
#pragma omp for				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
//						cout<<-hg2[m+l*nbk]/(chi1->chi_0->chi0_inf[y+(x*(x+1))/2])<<endl;
						
						D1=0;
						D2=0;
						for (j=0; j<nbw/2-nc; j++)
						{	
							p=j;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
								chi0tmp=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							else
							{
//								chi0tmp=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0tmp=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							chisp_tmp=chi0tmp/(1.0-Usp*chi0tmp/2);
							chich_tmp=chi0tmp/(1.0+Uch*chi0tmp/2);
							
							if (p<nbw/2 && p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp2=h[m+l*nbk];							
							
							p=j+nn;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
							{
								chi0n=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
//								chi0n=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0n=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							
							if (p<nbw/2 && p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp1=h[m+l*nbk];
							
							HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							//								HG[m+l*nbk+j*nk4th]=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<HG[m+l*nbk+j*nk4th]<<endl;
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{	
							p=j-nbw;
							if (p<0) p=-p;
							wn=2*PI*p*tem;
							if (p<=nbw/2)
							{
								chi0tmp=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
//								chi0tmp=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0tmp=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							chisp_tmp=chi0tmp/(1.0-Usp*chi0tmp/2);
							chich_tmp=chi0tmp/(1.0+Uch*chi0tmp/2);
							
							if (p<nbw/2 && p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp2=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp2=h[m+l*nbk];
							
							p=j-nbw+nn;
							wn=2*PI*p*tem;
							if (p<0) p=-p;
							if (p<=nbw/2)
							{
								chi0n=chiqw_array[y+(x*(x+1))/2+p*nk8th];
							}
							else
							{
//								chi0n=chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi1->chi_0->chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
								chi0n=chi0_inf[y+(x*(x+1))/2]/(wn*wn)+chi0_inf2[y+(x*(x+1))/2]/(wn*wn*wn*wn);
							}
							
							if (p<nbw/2 && p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+p*nk4th]/(wn*wn*wn*wn);
							else if (p)
								HGtmp1=-hg2[m+l*nbk]/(wn*wn)+h[m+l*nbk+(nbw/2-1)*nk4th]/(wn*wn*wn*wn);
							else
								HGtmp1=h[m+l*nbk];
							
							if (nn%2==0 && (j-nbw)!=-nn/2)
							{
								D2=D1;
								D1=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							}
							
							if ( (j-nbw)!=-nn/2 || nn%2 )
								HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(HGtmp2-HGtmp1)/(chi0tmp-chi0n);
							//								HG[m+l*nbk+j*nk4th]=(HGtmp2-HGtmp1)/(chi0tmp-chi0n);								
							else
								HG[m+l*nbk+j*nk4th]=(3*Usp*chisp_tmp+Uch*chich_tmp)*(4*D1-D2)/3.0;
							//								HG[m+l*nbk+j*nk4th]=(4*D1-D2)/3.0;
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j-nbw<<HG[m+l*nbk+j*nk4th]<<endl;
						}
						
					}
			}
			
			//			dcomplex *Htmp=new dcomplex[nbk*(nbk-2)*(nbw+1)];
			
			//TF dans le temps de I_n(q,iqm) et calcul de d/dtau I_n(q,tau) 
#pragma omp parallel private(l,m,j)
			{
				long int x, y, tmp;
				double rttmp, wn, wn1, t;
				dcomplex zsp0, zsp1;
				double zch0, zch1;
				double chisp, chich, Hinf, cq;
				dcomplex *ft_tmp1=new dcomplex[nbw];
				dcomplex dtauH0, dtauH_inf, Htau_inf, tmp_sp, tmp_ch, exp0, exp1;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
						dtauH0=0;
//						cq=-chi1->chi_0->chi0_inf[y+(x*(x+1))/2];
						cq=-chi0_inf[y+(x*(x+1))/2];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=2*j*PI*tem;
							chisp=1.0/(-wn*wn-Usp*cq/2.0);
							chich=1.0/(-wn*wn+Uch*cq/2.0);
							Hinf=hg2[m+l*nbk]*(3.0*Usp*chisp+Uch*chich);
							ft_tmp1[j]=tem*(HG[m+l*nbk+j*nk4th]-Hinf);
							dtauH0-=wn*ft_tmp1[j];
							
							//							Htmp[m+l*nbk+j*nk4th]=Hinf;
							//							Htmp[m+l*nbk+j*nk4th]=HG[m+l*nbk+j*nk4th];
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<setw(30)<<Hinf<<real(HG[m+l*nbk+j*nk4th]-Hinf)/real(HG[m+l*nbk+j*nk4th])<<endl;
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=2*(j-nbw)*PI*tem;
							chisp=1.0/(-wn*wn-Usp*cq/2.0);
							chich=1.0/(-wn*wn+Uch*cq/2.0);
							Hinf=hg2[m+l*nbk]*(3.0*Usp*chisp+Uch*chich);
							ft_tmp1[j]=tem*(HG[m+l*nbk+j*nk4th]-Hinf);
							dtauH0-=wn*ft_tmp1[j];
							
							//							Htmp[m+l*nbk+j*nk4th]=Hinf;
							//							Htmp[m+l*nbk+j*nk4th]=HG[m+l*nbk+j*nk4th];
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG[m+l*nbk+j*nk4th]<<setw(30)<<Hinf<<real(HG[m+l*nbk+j*nk4th]-Hinf)/real(HG[m+l*nbk+j*nk4th])<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
//						rttmp=sqrt(Usp*chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/2.0);
						rttmp=sqrt(Usp*chi0_inf[y+(x*(x+1))/2]/2.0);
						zsp0=dcomplex(0,rttmp);
						zsp1=-zsp0;
//						zch0=sqrt(Uch*chi1->chi_0->chi0_inf[y+(x*(x+1))/2]/2.0);
						zch0=sqrt(Uch*chi0_inf[y+(x*(x+1))/2]/2.0);
						zch1=-zch0;
						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<(zsp0*zsp0-Usp*cq/2.0)<<setw(40)<<(zsp1*zsp1-Usp*cq/2.0)<<endl;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<(zch0*zch0+Uch*cq/2.0)<<setw(40)<<(zch1*zch1+Uch*cq/2.0)<<endl;
						
						for (j=0; j<nbw; j++)
						{
							t=j/(nbw*tem);
							tmp_sp=exp(-zsp0*t)/((exp(-zsp0/tem)-(double)1.0)*(zsp0-zsp1))							
							+exp(-zsp1*t)/((exp(-zsp1/tem)-(double)1.0)*(zsp1-zsp0));
							tmp_ch=exp(-zch0*t)/((exp(-zch0/tem)-1.0)*(zch0-zch1))
							+exp(zch1*(1.0/tem-t))/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
							Htau_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
							HG[m+l*nbk+j*nk4th]=ft_tmp1[j]+Htau_inf;
							//							HG[m+l*nbk+j*nk4th]=Htau_inf;
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<ft_tmp1[j]<<setw(38)<<Htau_inf<<HG[m+l*nbk+j*nk4th]<<endl;
							//							Htmp[m+l*nbk+j*nk4th]=Htau_inf;
						}
						j=nbw;
						exp0=exp(zsp0/tem);
						exp1=exp(zsp1/tem);
						tmp_sp=((double)1.0)/(((double)1.0-exp0)*(zsp0-zsp1))
						+((double)1.0)/(((double)1.0-exp1)*(zsp1-zsp0));
						exp0=exp(-zch0/tem);
						tmp_ch=exp0/((exp0-(double)1.0)*(zch0-zch1))
						+1.0/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
						Htau_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						HG[m+l*nbk+nbw*nk4th]=ft_tmp1[0]+Htau_inf;
						//						HG[m+l*nbk+nbw*nk4th]=Htau_inf;						
						
						//						Htmp[m+l*nbk+nbw*nk4th]=Htau_inf;
						
						//						if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<ft_tmp1[0]<<setw(38)<<Htau_inf<<HG[m+l*nbk+j*nk4th]<<endl;
						
						// d/dtau Hn(q,tau) a tau=0	
						tmp_sp=-zsp0/((exp(-zsp0/tem)-(double)1.0)*(zsp0-zsp1))							
						-zsp1/((exp(-zsp1/tem)-(double)1.0)*(zsp1-zsp0));
						exp1=exp(zch1/tem);
						tmp_ch=-zch0/((exp(-zch0/tem)-1.0)*(zch0-zch1))
						-zch1*exp1/(((double)1.0-exp1)*(zch1-zch0));
						dtauH_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						dtauHn[m+l*nbk]=I*dtauH0+dtauH_inf;
						//						dtauHn[m+l*nbk]=dtauH_inf;						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<setw(38)<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauH_inf<<nbw*tem*(Htmp[m+l*nbk+nk4th]-Htmp[m+l*nbk])<<endl;
						
						// d/dtau Hn(q,tau) a tau=beta
						exp0=exp(-zsp0/tem);
						exp1=exp(-zsp1/tem);
						tmp_sp=-zsp0*exp0/((exp0-(double)1.0)*(zsp0-zsp1))
						-zsp1*exp1/((exp1-(double)1.0)*(zsp1-zsp0));
						exp0=exp(-zch0/tem);
						tmp_ch=-zch0*exp0/((exp0-(double)1.0)*(zch0-zch1))
						-zch1/(((double)1.0-exp(zch1/tem))*(zch1-zch0));
						dtauH_inf=hg2[m+l*nbk]*(3.0*Usp*tmp_sp+Uch*tmp_ch);
						dtauHn[m+l*nbk+nk4th]=I*dtauH0+dtauH_inf;						
						//						dtauHn[m+l*nbk+nk4th]=dtauH_inf;
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;						
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauH_inf<<nbw*tem*(Htmp[m+l*nbk+nbw*nk4th]-Htmp[m+l*nbk+(nbw-1)*nk4th])<<endl;
					}
				delete [] ft_tmp1;
			}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
				//TF dans l'espace de I_n(q,tau)				
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							fr_tmp1[m+l*nbk]=normFact*real(HG[m+l*nbk+j*nk4th]);
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG[m+l*nbk+j*nk4th].imag());
//							HG[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
							fr_tmp1[m+l*nbk]=normFact*imag(HG[m+l*nbk+j*nk4th]);
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG[m+l*nbk+j*nk4th]=dcomplex(HG[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
				}
				
				// TF dans l'espace de d/dtau I_n(q,tau) a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*real(dtauHn[m+l*nbk]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk].imag());
//						dtauHn[m+l*nbk].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*imag(dtauHn[m+l*nbk]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
                        dtauHn[m+l*nbk]=dcomplex(dtauHn[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
					}
				
				// TF dans l'espace de d/dtau I_n(q,tau) a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*real(dtauHn[m+l*nbk+nk4th]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk+nk4th].imag());
//						dtauHn[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						fr_tmp1[m+l*nbk]=normFact*imag(dtauHn[m+l*nbk+nk4th]);
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
                        dtauHn[m+l*nbk+nk4th]=dcomplex(dtauHn[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
						//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;
					}
				delete [] fr_tmp1;
			}
			
			for (j=0; j<=nbw; j++)
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
						HG2[m+l*nbk+j*nk4th]=HG[m+l*nbk+(nbw-j)*nk4th];
			for (l=0; l<nbk-2; l++)
				for (m=0; m<nbk; m++)
				{
					dtauHn2[m+l*nbk]=-dtauHn[m+l*nbk+nk4th];
					dtauHn2[m+l*nbk+nk4th]=-dtauHn[m+l*nbk];
				}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
				//TF dans l'espace de d/dtau I_n(r,tau)G_n(r,-tau)
				// a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+real(dtauHn[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk].imag());
//						dtauHn[m+l*nbk].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+imag(dtauHn[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk]=dcomplex(dtauHn[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
				
				// a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+real(dtauHn[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn[m+l*nbk+nk4th].imag());
//						dtauHn[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+imag(dtauHn[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn[m+l*nbk+nk4th]=dcomplex(dtauHn[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
				
				
				//TF dans l'espace de I_n(r,tau)G_n(r,-tau)
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=real(HG[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG[m+l*nbk+j*nk4th].imag());
//							HG[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=imag(HG[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG[m+l*nbk+j*nk4th]=dcomplex(HG[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
				}
				
				//				for (l=0; l<nbk-2; l++)
				//					for (m=0; m<nbk; m++)
				//					{
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk]<<nbw*tem*(HG[m+l*nbk+nk4th]-HG[m+l*nbk])<<endl;
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn[m+l*nbk+nk4th]<<nbw*tem*(HG[m+l*nbk+nbw*nk4th]-HG[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//					}
				
				delete [] fr_tmp1;
			}
			
			//			l=2, m=4;
			//			for (j=0; j<=nbw; j++)
			//				cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
			
			
			
			//TF dans le temps de (I_n(r,tau)G_n(r,-tau))(k)
#pragma omp parallel private(l,m,j)
			{
				dcomplex *ft_tmp1=new dcomplex[nbw];
				double *tau=new double[nbw+1];
				double *ftau=new double[nbw+1];
				double *coeffs=new double[4*nbw];
				double wn;
				dcomplex coefwn, d2S0, d2Sbeta;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=real(HG[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=real(dtauHn[m+l*nbk]);
						coeffs[1]=real(dtauHn[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0=2*coeffs[1];
						d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];
//						d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
						
						for (j=0; j<nbw; j++)
							ft_tmp1[j]=exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
						
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=imag(HG[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=imag(dtauHn[m+l*nbk]);
						coeffs[1]=imag(dtauHn[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0+=I*((double)2.0)*coeffs[1];
						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3]);
//						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3]);
						
						for (j=0; j<nbw; j++)
						{
							ft_tmp1[j]+=I*exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
							//							if (l==0 && m==0)  cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<ft_tmp1[j]<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						coefwn=HG[m+l*nbk]+HG[m+l*nbk+nbw*nk4th];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=(2*(j+nn)+1)*PI*tem;
							if ((j+nn)<nbw)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn]/(wn*wn*wn*wn);
							else
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
						}						
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=(2*(j+nn-nbw)+1)*PI*tem;
							if ((j+nn-nbw)<0)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn]/(wn*wn*wn*wn);
							else if ((j+nn-nbw)<nbw)
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn)
								+ ((double)1.0-exp(dcomplex(0,-(2*(j+nn)+1)*PI/nbw)))*ft_tmp1[j+nn-nbw]/(wn*wn*wn*wn);
							else
								HG[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn[m+l*nbk]+dtauHn[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn);
							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG[m+l*nbk+j*nk4th]<<endl;
						}						
						
					}
				delete [] ft_tmp1;
				delete [] tau;
				delete [] ftau;
				delete [] coeffs;
			}
			
			
			//#pragma omp parallel private(fr_tmp1,j,l,m)
			{
				fr_tmp1=new double[(nbk-2)*nbk];
				long int x, y, tmp;
				
				//TF dans l'espace de d/dtau I_n(r,-tau)G_n(r,-tau)
				// a tau=0
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+real(dtauHn2[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn2[m+l*nbk]=dcomplex(fr_tmp1[m+l*nbk],dtauHn2[m+l*nbk].imag());
//						dtauHn2[m+l*nbk].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk])*dtauGrt[y+(x*(x+1))/2]+imag(dtauHn2[m+l*nbk])*Grt[y+(x*(x+1))/2];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn2[m+l*nbk]=dcomplex(dtauHn2[m+l*nbk].real(),fr_tmp1[m+l*nbk]);
//						dtauHn2[m+l*nbk].imag()=fr_tmp1[m+l*nbk];
				
				// a tau=beta
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+real(dtauHn2[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn2[m+l*nbk+nk4th]=dcomplex(fr_tmp1[m+l*nbk],dtauHn2[m+l*nbk+nk4th].imag());
//						dtauHn2[m+l*nbk+nk4th].real()=fr_tmp1[m+l*nbk];
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk+nbw*nk4th])*dtauGrt[y+(x*(x+1))/2+nk8th]+imag(dtauHn2[m+l*nbk+nk4th])*Grt[y+(x*(x+1))/2+nbw*nk8th];
					}
				
				fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
				
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
                        dtauHn2[m+l*nbk+nk4th]=dcomplex(dtauHn2[m+l*nbk+nk4th].real(),fr_tmp1[m+l*nbk]);
//						dtauHn2[m+l*nbk+nk4th].imag()=fr_tmp1[m+l*nbk];
				
				
				//TF dans l'espace de I_n(r,-tau)G_n(r,-tau)
				//#pragma omp for
				for (j=0; j<=nbw; j++)
				{
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=real(HG2[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG2[m+l*nbk+j*nk4th]=dcomplex(fr_tmp1[m+l*nbk],HG2[m+l*nbk+j*nk4th].imag());
//							HG2[m+l*nbk+j*nk4th].real()=fr_tmp1[m+l*nbk];
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
						{
							x=l+1;
							if (x>nk0/2) x=nk0-x;
							y=m;
							if (y>nk0/2) y=nk0-y;
							if (y>x)
							{
								tmp=x;
								x=y;
								y=tmp;
							}
							fr_tmp1[m+l*nbk]=imag(HG2[m+l*nbk+j*nk4th])*Grt[y+(x*(x+1))/2+j*nk8th];
						}
					
					fftw_execute_r2r(fftplan_RO,fr_tmp1,fr_tmp1);
					
					for (l=0; l<nbk-2; l++)
						for (m=0; m<nbk; m++)
                            HG2[m+l*nbk+j*nk4th]=dcomplex(HG2[m+l*nbk+j*nk4th].real(),fr_tmp1[m+l*nbk]);
//							HG2[m+l*nbk+j*nk4th].imag()=fr_tmp1[m+l*nbk];
				}
				
				//				for (l=0; l<nbk-2; l++)
				//					for (m=0; m<nbk; m++)
				//					{
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn2[m+l*nbk]<<nbw*tem*(HG2[m+l*nbk+nk4th]-HG2[m+l*nbk])<<endl;
				//						cout<<setw(10)<<l<<setw(10)<<m<<setw(40)<<dtauHn2[m+l*nbk+nk4th]<<nbw*tem*(HG2[m+l*nbk+nbw*nk4th]-HG2[m+l*nbk+(nbw-1)*nk4th])<<endl;
				//					}
				
				delete [] fr_tmp1;
			}
			
			
			//TF dans le temps de (I_n(r,-tau)G_n(r,-tau))(k)
#pragma omp parallel private(l,m,j)
			{
				dcomplex *ft_tmp1=new dcomplex[nbw];
				double *tau=new double[nbw+1];
				double *ftau=new double[nbw+1];
				double *coeffs=new double[4*nbw];
				double wn;
				dcomplex coefwn, d2S0, d2Sbeta;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
					for (m=0; m<nbk; m++)
					{
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=real(HG2[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=real(dtauHn2[m+l*nbk]);
						coeffs[1]=real(dtauHn2[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0=2*coeffs[1];
						d2Sbeta=6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3];
//						d2Sbeta=6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3];
						
						for (j=0; j<nbw; j++)
							ft_tmp1[j]=exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
						
						for (j=0; j<=nbw; j++)
						{
							tau[j]=j/(nbw*tem);					
							ftau[j]=imag(HG2[m+l*nbk+j*nk4th]);
						}
						
						coeffs[0]=imag(dtauHn2[m+l*nbk]);
						coeffs[1]=imag(dtauHn2[m+l*nbk+nk4th]);
						spline_coeffs_rel(tau, ftau, nbw+1, coeffs);
//						spline_coeffs(tau, ftau, nbw+1, coeffs);
						
						d2S0+=I*((double)2.0)*coeffs[1];
						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/(nbw*tem)+2.0*coeffs[4*nbw-3]);
//						d2Sbeta+=I*(6.0*coeffs[4*nbw-4]/tem+2.0*coeffs[4*nbw-3]);
						
						for (j=0; j<nbw; j++)
						{
							ft_tmp1[j]+=I*exp(dcomplex(0,-PI*j/nbw))*((double)6.0)*coeffs[4*j];
							//							if (l==0 && m==0)  cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<ft_tmp1[j]<<endl;
						}
						
						fftw_execute_dft(fftplan,reinterpret_cast<fftw_complex*>(ft_tmp1),reinterpret_cast<fftw_complex*>(ft_tmp1));
						
						coefwn=HG2[m+l*nbk]+HG2[m+l*nbk+nbw*nk4th];
						for (j=0; j<nbw/2-nc; j++)
						{
							wn=(2*j+1)*PI*tem;
							HG2[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn2[m+l*nbk]+dtauHn2[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
							+ ((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*ft_tmp1[j]/(wn*wn*wn*wn);							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<setw(38)<<HG2[m+l*nbk+j*nk4th]<<endl;
						}						
						for (j=nbw/2-nc; j<nbw; j++)
						{
							wn=(2*(j-nbw)+1)*PI*tem;
							HG2[m+l*nbk+j*nk4th]=-I*coefwn/wn - (dtauHn2[m+l*nbk]+dtauHn2[m+l*nbk+nk4th])/(wn*wn)+I*(d2S0+d2Sbeta)/(wn*wn*wn) 
							+ ((double)1.0-exp(dcomplex(0,-(2*j+1)*PI/nbw)))*ft_tmp1[j]/(wn*wn*wn*wn);							
							//							if (l==2 && m==4) cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<setw(38)<<HG2[m+l*nbk+j*nk4th]<<endl;
						}						
						
					}
				delete [] ft_tmp1;
				delete [] tau;
				delete [] ftau;
				delete [] coeffs;
			}
			
			//			dcomplex *Gr1=new dcomplex[nbk*(nbk-2)*nbw];
			//			dcomplex *Gr2=new dcomplex[nbk*(nbk-2)*nbw];
			
			// somme finale sur k et ik_n			
#pragma omp parallel private(l,m,j)
			{
				long int x, y, tmp, p;
				double C2, C3, kn, wt;
				dcomplex Gtmp, Gtmpn, selftmp;
				dcomplex chijj3_qn_loc=0;
				
				double theta_k;
				int ind_theta;
				double *chijj3_k_loc=new double[(nw+1)*N_theta_k];
				for (j=0; j<(nw+1)*N_theta_k; j++)	chijj3_k_loc[j]=0;
				
#pragma omp for
				for (l=0; l<nbk-2; l++)
				{
					for (m=0; m<nbk; m++)
					{
						
						x=l+1;
						if (x>nk0/2) x=nk0-x;
						y=m;
						if (y>nk0/2) y=nk0-y;
						if (y>x)
						{
							tmp=x;
							x=y;
							y=tmp;
						}
						
//						C2=G1->Sigma_inf2[y+(x*(x+1))/2];
//						C3=G1->Sigma_inf3[y+(x*(x+1))/2];
						
						C2=Sigma_inf2[y+(x*(x+1))/2];
						C3=Sigma_inf3[y+(x*(x+1))/2];						
						
						wt=4;
						if (m==0 || m==nbk-1) wt=wt/2;
						
//						for (p=0; p<nk_chijj; p++)
//						{
//							if (k_chijj[2*p]==(l+1) && k_chijj[2*p+1]==m)	break;
//						}
						
//						if (p<nk_chijj) chijj3_k[p + n*nk_chijj]=0;
						
						theta_k=atan( (1.0*(nk0/2-m))/(nk0/2-l-1) );
						ind_theta=(int)floor(theta_k/Dtheta);
						
						for (j=0; j<nbw/2-nc; j++)
						{
							kn=(2*j+1)*PI*tem;
							
							selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn), 
											 -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+j*nk8th])/(kn*kn*kn*kn));
							
							Gtmp=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr1[m+l*nbk+j*nk4th]=Gtmp;
							//							Gr1[m+l*nbk+j*nk4th]=selftmp;
							//							Gr1[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							kn=(2*(j+nn)+1)*PI*tem;
							if ((j+nn)<nbw/2)
								selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j+nn)*nk8th])/(kn*kn*kn*kn), 
												 -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j+nn)*nk8th])/(kn*kn*kn*kn));
							else
								selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							//								selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							
							Gtmpn=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr2[m+l*nbk+j*nk4th]=Gtmpn;
							//							Gr2[m+l*nbk+j*nk4th]=selftmp;
							//							Gr2[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							chijj3_qn_loc-=I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
							
//							if (p<nk_chijj)	
//								chijj3_k[p + n*nk_chijj]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							chijj3_k_loc[ind_theta + n*N_theta_k]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							//							cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<j<<chijj3_qn_loc<<endl;
							//							HGsum[m+l*nbk+j*nk4th]=-I*U*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
						}
						for (j=nbw/2-nc; j<nbw; j++)
						{
							kn=(2*(j-nbw)+1)*PI*tem;
							
							if ((nbw-j-1)<nbw/2)
								selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th])/(kn*kn*kn*kn), 
												 -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nbw-j-1)*nk8th])/(kn*kn*kn*kn));
							else
								selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							//								selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							
							Gtmp=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr1[m+l*nbk+j*nk4th]=Gtmp;
							//							Gr1[m+l*nbk+j*nk4th]=selftmp;
							//							Gr1[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							kn=(2*(j+nn-nbw)+1)*PI*tem;
							if ((j+nn-nbw)<0)
							{
								if ((nbw-j-nn-1)<nbw/2)
									selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(nbw-j-nn-1)*nk8th])/(kn*kn*kn*kn), 
													 -Sigma_inf/kn+C3/(kn*kn*kn)-imag(Self_kw_array[y+(x*(x+1))/2+(nbw-j-nn-1)*nk8th])/(kn*kn*kn*kn));
								else
									selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
								//									selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							}
							else
							{
								if ((j+nn-nbw)<nbw/2)
									selftmp=dcomplex(-C2/(kn*kn)+real(Self_kw_array[y+(x*(x+1))/2+(j+nn-nbw)*nk8th])/(kn*kn*kn*kn), 
													 -Sigma_inf/kn+C3/(kn*kn*kn)+imag(Self_kw_array[y+(x*(x+1))/2+(j+nn-nbw)*nk8th])/(kn*kn*kn*kn));
								else
									selftmp=dcomplex(-C2/(kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
								//									selftmp=dcomplex(-C2/(kn*kn)+Sigma_inf4[y+(x*(x+1))/2]/(kn*kn*kn*kn), -Sigma_inf/kn+C3/(kn*kn*kn));
							}		
							
							Gtmpn=((double)1.0)/(dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu-selftmp);
							
							//							Gr2[m+l*nbk+j*nk4th]=Gtmpn;
							//							Gr2[m+l*nbk+j*nk4th]=selftmp;
							//							Gr2[m+l*nbk+j*nk4th]=dcomplex(0.0, kn)-ek[m+(l+1)*nbk]+mu;
							
							chijj3_qn_loc-=I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
							
//							if (p<nk_chijj)	
//								chijj3_k[p + n*nk_chijj]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							chijj3_k_loc[ind_theta + n*N_theta_k]-=real(I*U*tem*normFact*wt*dek[m+l*nbk]*Gtmp*Gtmpn*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th]))/(2*qn);
							
							//							cout<<setw(10)<<l<<setw(10)<<m<<setw(10)<<(j-nbw)<<chijj3_qn_loc<<endl;
							//							HGsum[m+l*nbk+j*nk4th]=-I*U*(HG[m+l*nbk+j*nk4th]-HG2[m+l*nbk+j*nk4th])/(2*qn);
						}					
					}
				}
#pragma omp critical
				{
					chijj3tmp[n]+=chijj3_qn_loc;
					chijj3_qn[n]+=real(chijj3_qn_loc);
					for (p=0; p<N_theta_k; p++) chijj3_k[p + n*N_theta_k]+=chijj3_k_loc[p + n*N_theta_k];
				}
				delete [] chijj3_k_loc;
			}
			
			cout<<setw(10)<<nn<<setw(50)<<chijj3tmp[n]<<endl;
			//			cout<<setw(10)<<nn<<setw(50)<<chijj3tmp[n];	
			//			if (n%8==0 || n==1)
			//				calc_chijj_vertex_corr2_wn(nn, k0, wn0, NULL, HG0tmp, NULL, NULL, NULL);
			//			else
			//				cout<<endl;
			
			//			k0[0]=1;
			//			k0[1]=2;
			//			int wn0=-20;
			//			calc_chijj_vertex_corr2_wn(nn, k0, wn0, NULL, HG0tmp, HGsum, Gr1, Gr2);
			//			if (wn0<0)
			//				cout<<HGsum[k0[1]+(k0[0]-1)*nbk+(wn0+nbw)*nk4th]<<endl;
			//			else
			//				cout<<HGsum[k0[1]+(k0[0]-1)*nbk+wn0*nk4th]<<endl;
			
			
		}
		
		delete [] chijj3tmp;
		delete [] dtauHn;
		delete [] dtauHn2;
		delete [] HG;
		delete [] HG2;		
		
		
		// Pour enregistrer chijj3(iq_n) en binaire
		
		const char *nameForm_chijj3_qn_bin="./chijj3_qn_bin_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameForm_chijj3_qn_bin, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
		file3.open(name, ios::out | ios::binary);
		
		for (j=0; j<N; j++)
		{
			file3.write((char*)&(ind_freq[j]),sizeof(ind_freq[j]));
			file3.write((char*)&(chijj3_qn[j]),sizeof(chijj3_qn[j]));
		}
		file3.close();
		// fin de l'enregistrement en binaire
		
		
		// Pour enregistrer chijj3(iq_n) en ascii
		const char *nameForm_chijj3_qn="./chijj3_qn_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
		
		sprintf(name, nameForm_chijj3_qn, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
		file3.open(name, ios::out );
		file3<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
		
		for (j=0; j<N; j++)
		{
			file3<<setw(10)<<ind_freq[j]<<chijj3_qn[j]<<'\n';
		}
		file3.close();
		// fin de l'enregistrement en ascii
		
		for (j=0; j<=nbw/2; j++) chijj_qn[j]=0;
		
		chijj_qn[0]=chijj1_qn[0]+chijj2_qn[0];
		
		for (j=1; j<N; j++)
			chijj_qn[ind_freq[j]]=chijj1_qn[ind_freq[j]]+chijj2_qn[ind_freq[j]]+chijj3_qn[j];
		
		// sommer chijj3_k pour tous les k
		
		if (N_theta_k>1)
		{
			double *chijj3_k_tot=new double[nw+1];
			
			for (j=0; j<=nw; j++)
				chijj3_k_tot[j]=0;
			
			for (j=0; j<=nw; j++)
				for  (p=0; p<N_theta_k; p++)
					chijj3_k_tot[j]+=chijj3_k[p + j*N_theta_k];
			
			//enregistrer chijj3_k et chijj3_k_tot et la comparaison de chijj3_k_tot avec chijj3_qn dans un fichier
			
			const char *nameForm_chijj3_k="./chijj3_k_qn_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
			const char *nameForm_chijj3_k_tot="./chijj3_k_tot_qn_old_%5.3f_%5.3f_%5.3f_%5.3f_%7.5f_%6.4f_%1d_%1d_%1d.dat";
			
			sprintf(name, nameForm_chijj3_k, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
			file.open(name, ios::out );
			file<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
			
			sprintf(name, nameForm_chijj3_k_tot, (double)U, (double)tp, (double)tpp, (double)t3, (double)density, (double)tem, nk0, nbw, mmin);
			file2.open(name, ios::out );
			file2<<setprecision(16)<<setiosflags(ios::left)<<setiosflags(ios::scientific);
			
			for (j=1; j<N; j++)
			{
				file2<<setw(10)<<ind_freq[j]<<setw(30)<<chijj3_k_tot[j]<<setw(30)<<chijj3_qn[j]<<chijj3_k_tot[j]/chijj3_qn[j]<<'\n';
				for  (p=0; p<N_theta_k; p++)
					file<<setw(10)<<ind_freq[j]<<setw(10)<<(p+1)<<chijj3_k[p + j*N_theta_k]<<'\n';
				//				file<<setw(10)<<ind_freq[j]<<setw(10)<<k_chijj[2*p]<<setw(10)<<k_chijj[2*p+1]<<chijj3_k[p + j*nk_chijj]<<'\n';
			}
			file.close();
			file2.close();
			
			delete [] chijj3_k;
		}
		
		
		delete [] ind_freq;
		delete [] hg2;
	}
	
	delete [] chijj_qn;
	delete [] chijj1_k;
	delete [] chijj2_k;
	
	
	fftw_destroy_plan(fftplan);
	fftw_destroy_plan(fftplan_t_back);
	fftw_destroy_plan(fftplan_RO);
	
	delete [] h;
	delete [] dek;
	delete [] Grt;
	delete [] dtauGrt;
	
}

