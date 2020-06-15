
//! tpsc_main.cpp: main function for the two-particle self-consistent computation of one- and two-particle properties for the Hubbard model
//! created by Dominic Bergeron

#include "cond_opt.h"
#include <time.h>

//! note: I used objects of different classes to compute different types of quantities, namely the classes in which the functions that compute those quantities are defined. However, one can simply declare an object of class "cond_opt", which inherits from all the other classes, and call any function in any of the classes defined in the code. The structure is as follow: "cond_opt" inherits from "green", which inherits from "chi", which inherits from "chi0" which inherits from "green0",  which inherits from "hamiltonien", which inherits from "generique". In principle, the structure of the code would allow to consider different types of lattices, but here the sysmetries of the square lattice are used to save  memory. This is necessary because the fast Fourier transform approach used extensively in the code consumes a lot of memory.

int main(int argc, char **argv)
{
	long int j, l, m;
	
 	double tp, tpp, t3, U, T, n, eta=0.1;
	int nk0, Nw0, NPadeMax, mmin, Nfreq_save;
	bool chijj3;	
	
//	Load the parameters of the system
//  T: temperature	
//	n:density
//	nk0: number of sites along one direction, must be a power of 2
//	Nw0: total number of frequencies (+ et -), must be a power of 2
//	chijj3: boolean deciding if the Aslamazov-Larkin (AL) contributions to the current-current correlation function (chijj) are to be calculated
//  mmin: integer controling the number of Matsubara frequencies for which the AL contributions are computed, (parameter m of appendix D of the article, the larger m is, the more sparse the Matsubara grid will be)
	fstream inputFile;
	inputFile.open("../../input.dat",ios::in);
	if (!inputFile)
	{
		cout<<"echec d'ouverture de input.dat"<<'\n';
		return 0;
	}
//	inputFile>>U>>tp>>tpp>>t3>>n>>T>>nk0>>Nw0>>Nfreq_save;
	inputFile>>U>>tp>>tpp>>t3>>n>>T>>nk0>>Nw0>>chijj3>>mmin;
	inputFile.close();
//  system parameters loaded
	
	cout<<"U: "<<U<<"  tp:  "<<tp<<"  tpp:  "<<tpp<<"  t3:  "<<t3<<"  n:  "<<n<<"  T:  "<<T<<"  nk0:  "<<nk0<<"  Nw0:  "<<Nw0<<"  chijj3:  "<<chijj3<<"  mmin:  "<<mmin<<endl;

	double paramsHub[]={tp, tpp, t3, U};

	struct stat file_stat;
	
	char sw_name[200];
	
//	const char name_form_sw[]="EDC_U%g_tp%g_tpp%g/spectral_weight";
//	sprintf(sw_name,name_form_sw,U,tp,tpp);
	
	const char name_form_rep[]="U%g_tp%g_tpp%g/EDC/tem%g";
	sprintf(sw_name,name_form_rep,U,tp,tpp,T);
	
	if (stat(sw_name,&file_stat))
	{
		mkdir(sw_name, S_IRWXU | S_IRWXG | S_IRWXO);
	}
	
	strcat(sw_name,"/spectral_weight");
	cout<<"sw_name: "<<sw_name<<endl;
	
//  Compute only the Green function and one-particle quantities 
//  look in the files green.h and green.cpp to see the other functions that compute one-particle quantities
	green *G1=new green(nk0, Nw0, n, T, eta, paramsHub);
 

//	To put chi0(q, iwn=0) and chi0(qmax,iwn) in files, comment this line:
//	G1->set_save_chi0_to_file(false);
//  print system parameters on screen
	G1->get_params();
//  the matrix containing chi0 does not need to be kept in memory after the self-energy is computed
	G1->set_keep_chiqw(false);
//  compute the self-energy
	G1->calc_Self_FFT_spline();
//  compute the trace of Sigma2*G1	
	G1->traceSelfG();
//  compute the chemical potentiel
	G1->find_mu();
//  find the "Fermi surface" (maxima of A(k,omega=0)) with interactions	
	G1->find_FS_inter();

//	G1->save_Self(Nfreq_save);
	
//	int k[]={nk0/4, nk0/4};
//	G1->save_Green(k);


	
//  get the spectral density versus k at a given energy
//  analytical continuation is done with the values in etap as the imaginary part of the energy in the Green function, 
//  if the last parameter is true, Pade approximants are used, otherwise a polynomial of degree 2 or 3 is used, which yields only qualitative results but is more stable than Pade.
	int NP=10;
	int N_eta=4;
	double eta_p[]={ 0.01, 0.001, 0.0001, 0.00001};
	int *n_qn=new int[Nw0/2];
	for (j=0; j<Nw0/2; j++) n_qn[j]=j;
//   for the parameters, see the function definition in green.h or green.cpp
	
	double w_lim=0.1;
	double Dw=0.02;
	int Nw_Ak=(int)(2*w_lim/Dw)+1;
	double w=0;
	
	G1->MDC(n_qn, NP, w, eta_p, N_eta, true);


//	for (j=0; j<Nw_Ak; j++)
//	{
//		w=-w_lim+j*Dw;
//		G1->MDC(n_qn, NP, w, eta_p, N_eta, true);
//	}

	
//  get the non-interacting density of states 
//	 G1->G0wn_local();
//  get the density of states with interactions		
	G1->Gwn_local();
	

/*
 // compute energy distribution curves
	NP=Nw0;
	
	double wm=10.0;
	Dw=0.01;
	int Nw_EDC=(int)(2*wm/Dw)+1;
	
	double *wv=new double[Nw_EDC];
	for (j=0; j<Nw_EDC; j++)
	{
		wv[j]=-wm+j*Dw;
	}
	
	int dkx=nk0/16;
	int kx0=0;
	int kxf;
	
	int ky0i=nk0/8;
	int ky0f=7*nk0/8;
	
	int *kx=new int[nk0/2+1];
	int *ky=new int[nk0/2+1];
	
	while (kx0<nk0/2)
	{
		kxf=(ky0f+kx0)/2;
		kx[0]=(kx0+ky0i)/2;
		ky[0]=kx[0]-kx0;
		j=0;
		while ((kx[j]+1)<=kxf)
		{
			j++;
			kx[j]=kx[j-1]+1;
			ky[j]=ky[j-1]+1;
		}
		G1->EDC(kx, ky, j+1, wv, Nw_EDC, mmin, NP, eta_p, N_eta, sw_name);
		kx0+=dkx;
	}
*/
	
/*
	while (kx0<nk0/2)
	{
		kx[0]=kx0;
		ky[0]=0;
		j=0;
		while ((kx[j]+1)<=nk0/2)
		{
			j++;
			kx[j]=kx[j-1]+1;
			ky[j]=ky[j-1]+1;
		}
		G1->EDC(kx, ky, j+1, wv, Nw_EDC, mmin, NP, eta_p, N_eta, sw_name);
		kx0+=dkx;
	}
*/
/*
	
	int Nkx0=3;
	int kx0_v[]={3*nk0/8, 7*nk0/16, 15*nk0/32};
	int ky0[]={0, 0, -nk0/32};
	int kxf_v[]={5*nk0/8, 11*nk0/16, 3*nk0/4};
	
	for (l=0; l<Nkx0; l++)
	{
		kx[0]=kx0_v[l];
		ky[0]=ky0[l];
		j=0;
		while ((kx[j]+1)<=kxf_v[l])
		{
			j++;
			kx[j]=kx[j-1]+1;
			ky[j]=ky[j-1]+1;
		}
		G1->EDC(kx, ky, j+1, wv, Nw_EDC, mmin, NP, eta_p, N_eta, sw_name);
	}
 
	delete [] wv;
	delete [] kx;
	delete [] ky;
*/

	
/*
//  get the Matsubara self-energy for a set of wave vectors k given in the file "k_selfk.dat". 
//  The file has the number of k values on the first line and the values kx and ky in two columns on the following lines
	int nk_selfk;
	inputFile.open("k_selfk.dat",ios::in);
	inputFile>>nk_selfk;
//	cout<<"nk_selfk: "<<nk_selfk<<endl;
	int *k_selfk=new int[2*nk_selfk];
	cout<<setiosflags(ios::left);
	for (j=0; j<nk_selfk; j++)
	{
		inputFile>>k_selfk[2*j]>>k_selfk[2*j+1];
//		cout<<setw(10)<<k_selfk[2*j]<<k_selfk[2*j+1]<<endl;
	}
//	G1->save_self_ikn(k_selfk, nk_selfk, Nw0/32);
*/

	delete G1;

/*
//  Compute the current-current correlation function chijj(iwn)
//  the bubble contribution is in the file starting with chijj1_qn
//  the different contributions by angles to the bubble are in the file starting with chijj1_k_qn	
//  the Maki-Thompson vertex correction (MK) is in the file starting with chijj2_qn
//  the different contributions by angles to the MK correction are in the file starting with chijj2_k_qn
//  the Aslamazov-Larkin vertex correction (AL) is in the file starting with chijj3_qn
//  the different contributions by angles to the AL correction are in the file starting with chijj3_k_qn
//  for the bubble and MK contributions, all the Matsubara frequencies up to the cutoff are in the files (Matsubara indices are in the first columns of the files)
//  for the AL contribution, the Matsubara frequencies are those of the grid determined by Nw0 and mmin (see appendix D of the article). 
//  Be careful when adding all three contributions together! To do so, one must first keep only the frequencies of the sparse Matsubara grid for chijj1 and chijj2 
//  and then add chijj3, which is known only on that grid.
//  see the files cond_opt.h and cond_opt.cpp for the description and définition of the functions	
	cond_opt *sigma1;
	sigma1=new cond_opt(nk0, Nw0, n, T, eta, paramsHub);
//	the contributions to chijj(iwn) from N_theta_k angle in the Brillouin zone are kept
//	int N_theta_k=1;
	
//  to compute only the bubble contribution, use this line
//	sigma1->calc_chijj_bulle_1();
//	sigma1->calc_chijj_bulle(N_theta_k, k_selfk, nk_selfk);
//  to compute the bubble and the vertex corrections, use this line (if chijj3=false, only the MK correction will be computed)	
//	sigma1->calc_chijj_vertex_corr_all_optim(chijj3, mmin, N_theta_k, k_selfk, nk_selfk);
	sigma1->calc_chijj_vertex_corr_all_optim(chijj3, mmin, 1, NULL, 0);

//  old version of the function that computes the bubble and the vertex corrections, less memory efficient, but stil works	
//	sigma1->calc_chijj_vertex_corr_all(chijj3, mmin, N_theta_k, k_selfk, nk_selfk);
//	sigma1->calc_chijj_vertex_corr_all(chijj3, mmin, 1, NULL, 0);
	
	delete sigma1;
*/

/*
// Compute the spin and charge suscptibilities and the spin correlation length.
// here the correlation length is computed for a set of values of U
// look in the files chi.h and chi.cpp to see the other available functions of the class chi
	chi *chi1;
	
//	double Uvals[]={6.0};
	int NU=1;
	
	for (j=0; j<NU; j++)
	{
//		paramsHub[3]=Uvals[j];
		
		chi1=new chi(nk0, Nw0, n, T, paramsHub);
		chi1->get_params();
		chi1->set_keep_chirt(false);		
		chi1->find_Usp_Uch();
		chi1->sp_corr_length();
		delete chi1;
	}
*/
	
/*
//  Compute some values of chi0
//  look at the files chi0.h and chi0.cpp to see the available functions of the class chi0
	chi0 *chi0_1;

//  here chi0 is computed in a region delimited by qlims={qx1,qy1,qx2,qy2} in the first eighth of the Brillouin zone. calc_chi0_qlims uses the function chiqw that
//  computes a single value of chi0(q,omega) at a time on the infinite lattice.
	int nk=nk0/2;	
	double qlims[]={PI,0.99*PI,PI,PI};
	int Nq[]={1,10};
	double n1=1.1, n2=1.3, dn=0.05;
	for (n=n1; n<=n2; n+=dn)
	{
		chi0_1=new chi0(nk0, Nw0, n, T, paramsHub);
		chi0_1->calc_chi0_qlims(0, qlims, Nq);
		delete chi0_1;
	}
*/
/*
//  Compute chi0 with the fast Fourier transforms and cubic splines approach.
	chi0 *chi0_1=new chi0(nk0, Nw0, n, T, paramsHub);
	chi0_1->get_params();	
//  no need to keep the real space version of chi0 if the self-energy is not to be calculated
	chi0_1->set_keep_chirt(false);
//  uncomment this line if you do not want chi0(q,iwn=0) and chi0(qmax,iwn) to be saved in files
//	chi0_1->set_save_chi0_q(false);
	chi0_1->calc_chiqw_FFT_spline(nk0, Nw0);
	
	delete chi0_1; 
*/

	cout<<"fin\n";
	
	return 0;

}
