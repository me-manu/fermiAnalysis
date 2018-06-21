// written by Benoit Lott & Stefan Funk based on the note LAT-SAP-TN-OOO1 by Jean Ballet 
// Compile with: .L SensitivityMap.C++// version 3/9/2006
// 
#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLine.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TMinuit.h"
#include <fstream>
int const evtlim=1000000;
// Global variables
double GTime[evtlim],GE[evtlim],Gr[evtlim],Glon[evtlim],Glat[evtlim],Gid[evtlim],Grs[evtlim][10];
int GConvtype[evtlim],Gevt, n_source;
double  glon_s[10],glat_s[10],Index_source[10],Flux_source[10];

//Effective area
TH2F* diffuseMap[34];
TH2F* diffuseMap2[30];
TH2D* exposure_maps[2][18];
TH2D *dta[2];
TH2F* limb_smooth;
double par_bgd[10]; 
double enisot[100],fluxisot[100],enisrat;
int ind_enisot;
double enlimb[100],fluxlimb[100],enlimbrat;
double engal[100],engalrat;
double eexp[100],eexprat;
int ind_enlimb;
Double_t pi2=2*TMath::Pi();

int get_evtlim(){
return evtlim;
}
Double_t thet_68(Double_t E){   
   double pf[5];
if (E>2500) {pf[0]=1.011503;pf[1]=0.023886;pf[2]=2.85e-4;pf[3]=-0.70298;} else {pf[0]=1.262954;pf[1]=0.030774;pf[2]=0;pf[3]=-0.854023;}
  Double_t gam=2.35-.095*log(E/100);
  double t= pow(E/100.,pf[3]);
  // double xx=2.*log10(E/100.)+2;
  double sigma=(pf[0])*sqrt(pow(pf[1]*t,2)+pow(pf[2],2.));
 Double_t v=sigma*sqrt(2*(gam*(pow(0.32,1./(1-gam))-1)));
 return v*180/3.1415927;
 } 

double psf(double E, double r,int i){
 Float_t conv2=28.5/thet_68(E);
 double ind=r*conv2*180./TMath::Pi(); //index from map creation with IRFloader
 if (ind >100) return 0;
 double elog=log10(E);
 double z1,z2;
 int yb=dta[i]->GetYaxis()->FindBin(elog);
 int xb=dta[i]->GetXaxis()->FindBin(ind);
 double x0=dta[i]->GetXaxis()->GetBinCenter(xb);
 if (ind>x0) {z1=dta[i]->GetBinContent(xb,yb);
   z2=dta[i]->GetBinContent(xb+1,yb);
   return ((x0+1-ind)*z1+ (ind-x0)*z2)*conv2*conv2;}
 else {z1=dta[i]->GetBinContent(xb-1,yb);
       z2=dta[i]->GetBinContent(xb,yb);
   return ((x0-ind)*z1+ (ind-x0+1)*z2)*conv2*conv2; } 
}

Double_t Exposure(Double_t E, Double_t l, Double_t b,int i){
   if (E<30) return 0;
  int xBin= exposure_maps[i][0]->GetXaxis()->FindBin(l);
  int yBin= exposure_maps[i][0]->GetYaxis()->FindBin(b);
  double preIndex = log(E/eexp[0])/eexprat;
  int index= static_cast<int>(preIndex); 
  double exposure = 0;
  if (index>=17) exposure = exposure_maps[i][17]->GetBinContent(xBin,yBin);
  else 
  {double Emin=eexp[index];
   double Emax=eexp[index+1]; 
  exposure = ((Emax-E)*exposure_maps[i][index]->GetBinContent(xBin, yBin) 
    + (E-Emin)*exposure_maps[i][index + 1]->GetBinContent(xBin, yBin))/(Emax-Emin);
  }
  return exposure;
}
 
Double_t Np(Double_t *x, Double_t *par){
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100; 
  Double_t l      = par[2];
  Double_t b      = par[3]; 
 Double_t sum=0;
for (int i=0;i<2;i++){  
      sum+=ss*pow(E/100,-alpha)*Exposure(E,l,b,i)*E*log(10);
 }
  return sum;
 }

Double_t Np_front(Double_t *x, Double_t *par){
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100; 
  Double_t l      = par[2];
  Double_t b      = par[3]; 
  return ss*pow(E/100,-alpha)*Exposure(E,l,b,0)*E*log(10);
 }
Double_t Np_back(Double_t *x, Double_t *par){
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100; 
  Double_t l      = par[2];
  Double_t b      = par[3]; 
  return ss*pow(E/100,-alpha)*Exposure(E,l,b,1)*E*log(10);
 }


Double_t Np2(Double_t *x, Double_t *par){
  Double_t E=pow(10,x[0]); 
  Double_t r=x[1]; 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100; 
  Double_t l      = par[2];
  Double_t b      = par[3]; 
 Double_t sum=0;
for (int i=0;i<2;i++){  
      sum+=ss*pow(E/100,-alpha)*psf(E,r,i)*2*TMath::Pi()*r*Exposure(E,l,b,i)*E*log(10);
 }
  return sum;
 }
double ResidualBackground(Double_t E){
  double preIndex=log(E/enisot[0])/enisrat;
  if (preIndex <0 || preIndex>ind_enisot) return 0;
  Int_t index =int(preIndex);
  return (1-preIndex+index)*fluxisot[index]+(preIndex-index)*fluxisot[index+1];
}

double Limb(Double_t E,Double_t l,Double_t b){
  return 0;
  int xBin= limb_smooth->GetXaxis()->FindBin(l);
  int yBin= limb_smooth->GetYaxis()->FindBin(b);
  double preIndex=log(E/enlimb[0])/enlimbrat;
  if (preIndex <0 || preIndex>ind_enlimb) return 0;
  Int_t index =int(preIndex); 
  return ((1-preIndex+index)*fluxlimb[index]+(preIndex-index)*fluxlimb[index+1])*limb_smooth->GetBinContent(xBin,yBin);

}


Double_t BackgroundFlux(Double_t E, Double_t l, Double_t b, Double_t normeg=1, Double_t normgal=1,Double_t index_gal=0){

//   std::cout << "Background Flux " << std::endl;
  double galactic = 1.e-5;
  int xBin= diffuseMap[0]->GetXaxis()->FindBin(l);
  int yBin= diffuseMap[0]->GetYaxis()->FindBin(b);
  Double_t preIndex    = log(E/engal[0])/engalrat;
  Int_t index = static_cast<int> (preIndex);
  if (index < 0) galactic = 1.e-5;
  if (index >= 29) galactic = diffuseMap[29]->GetBinContent(xBin,yBin);
  else galactic= (1-preIndex+index)*diffuseMap[index]->GetBinContent(xBin,yBin)+ 
    (preIndex-index)*diffuseMap[index+1]->GetBinContent(xBin,yBin);
  return normeg*ResidualBackground(E) + pow(E/100.,-index_gal)*galactic*normgal; 
}

void LoadMaps_Gal(std::string filename){
  char ttrc[50];  
  TFile* f = new TFile(filename.c_str(),"r");
  TH1F* h1=(TH1F *) f->Get("en");
  for (int i=0;i<30;i++){
    sprintf(ttrc,"h_%i",i);
    diffuseMap[i] =(TH2F *) f->Get(ttrc);
    engal[i]=h1->GetBinContent(i+1);
    //printf("%lf \n",engal[i]);
  } 
 engalrat=log(engal[1]/engal[0]);  
}

void LoadMaps_Gal_p6(std::string filename){
  char ttrc[50];  
  TFile* f = new TFile(filename.c_str(),"r");
  for (int i=0;i<34;i++){
    sprintf(ttrc,"h_%i",i);
    diffuseMap[i] =(TH2F *) f->Get(ttrc);
  }  
}

void LoadMaps_Aeff(std::string filename){
  char ttrc[50];  
  TFile* f = new TFile(filename.c_str(),"r");
 TH1F* h1=(TH1F *) f->Get("en");
 if (h1!=NULL) {
   for (int i=0;i<18;i++) eexp[i]=h1->GetBinContent(i+1);}
 else { for (int i=0;i<18;i++) eexp[i]=pow(10,1+0.25*i);}
  eexprat=log(eexp[1]/eexp[0]);  
 for (int k=0;k<2;k++){
   for (int i=0;i<18;i++){
   if (k==0) sprintf(ttrc,"exposure_fa_%i",i); 
   if (k==1) sprintf(ttrc,"exposure_ba_%i",i);
    exposure_maps[k][i] =(TH2D*)f->Get(ttrc);
  }  
 }
 //for (int i=0;i<18;i++) printf("%lf \n",eexp[i]);
}

void LoadMaps_PSF(std::string filename){
  char ttrc[50];  
  TFile* f = new TFile(filename.c_str(),"r");
  for (int k=0;k<2;k++){
  sprintf(ttrc,"ta_%d",k);
  dta[k] =(TH2D *) f->Get(ttrc);
  }     
} 

void LoadCoef(std::string filename){
  // loads coefficients generated by bgd_fit.C
 FILE* fp=fopen(filename.c_str(),"r");
   fscanf(fp,"%lf %lf %lf %lf",&par_bgd[0],&par_bgd[1],&par_bgd[2],&par_bgd[3]);
}
void Load_Isot(std::string filename){
 FILE* fp=fopen(filename.c_str(),"r");
 int i=0;
 double errfluxisot;
 while (fscanf(fp,"%lf %lf %lf",&enisot[i],&fluxisot[i],&errfluxisot)>0)
   {fluxisot[i]*=pow(enisot[i]/100.,2.1); i++;}
 enisrat=log(enisot[1]/enisot[0]);
 ind_enisot=i; 
 return;
}
void Load_Limb(std::string filename1,std::string filename2){
 FILE* fp=fopen(filename1.c_str(),"r");
 int i=0;
 double errfluxlimb;
 while (fscanf(fp,"%lf %lf %lf",&enlimb[i],&fluxlimb[i],&errfluxlimb)>0)
   {fluxlimb[i]*=pow(enlimb[i]/100.,2.1); printf("%e \n",fluxlimb[i]);i++; }
 enlimbrat=log(enlimb[1]/enlimb[0]);
 ind_enlimb=i;
 TFile* f2 = new TFile(filename2.c_str());
 limb_smooth =(TH2F *) f2->Get("h_0");
 return;
}


void initialize_TS(std::string path = "/nfs/farm/g/glast/u/mmeyer/adaptive/"){
  LoadMaps_Aeff(path + "expmap_tot_P8R2_SOURCE_72m.root");
  LoadMaps_Gal(path + "ring_P8_SOURCE.root");  
  LoadMaps_PSF(path + "IRF_P8R2_SOURCE_V6.root");
  Load_Isot(path + "iso_P8R2_SOURCE_V6_v06.txt");
  return;
}

void initialize_TS_clean(){
  LoadMaps_Aeff("/data/lightcurves/expmap_tot_P7_12m_clean.root");
  LoadMaps_Gal("../2LAC/ring_P7v6.root");  
  LoadMaps_PSF("IRF_P7V6_clean.root");
  Load_Isot("isotrop_2year_P76_clean_v0.txt");
  Load_Limb("limb_2year_P76_clean_v0_smooth.txt","limb_smooth.root");
  return;
}

Double_t TSKernel(Double_t *x, Double_t *par){  

  // !!! Internally calculate everything with respect to the
  // normalisation at 100 MeV

  Double_t E      = pow(10,x[0]);             // from log Energy to Energy
  Double_t r      = x[1];
  Double_t alpha  = par[0];               //photon index 
  Double_t sFlux  = par[1];            
  Double_t l      = par[2];
  Double_t b      = par[3];
  bool addGalactic = par[4];

  // g = (dN/dE)_100MeV of the source, divided by the dN/dE_100MeV of the bg
  // ss = (alpha-1)*par[1]/100*pow(E, -alpha)
  // sb = par[2] * (E/E0)^{-2.1} <--- this -2.1 is put into the ss for computational reasons
  Double_t kk      = (alpha-1)*sFlux/100;  // differential flux
  Double_t ss     = kk * pow(E/100, -alpha+2.1);   // to go from norm
						    // to differential
  Double_t sb=BackgroundFlux(E,l,b, addGalactic);
  Double_t sum=0;
  for (int i=0;i<2;i++){  
  Double_t g=psf(E,r,i)*ss/sb; //PSF
  Double_t I=(1+g)*log(1+g)-g; 
  Double_t integrandOmega = 2*TMath::Pi()*I*r; //(dOmega == 2 Pi dr)
  Double_t integrandE = Exposure(E, l, b,i)*sb*E*log(10.)*pow(E/100,-2.1);
  sum+=integrandE * integrandOmega;
  }
  return sum; 
}
Double_t Ts_integrand(Double_t *x, Double_t *par){  

  // !!! Internally calculate everything with respect to the
  // normalisation at 100 MeV

  Double_t E      = pow(10,x[0]);             // from log Energy to Energy
  Double_t alpha  = par[0];               //photon index 
  Double_t sFlux  = par[1];            
  Double_t l      = par[2];
  Double_t b      = par[3];
  Double_t normeg = par[4];
  Double_t normgal =par[5];
  Double_t indexgal =par[6];
  Double_t time =par[7];
  
  // g = (dN/dE)_100MeV of the source, divided by the dN/dE_100MeV of the bg
  // ss = (alpha-1)*par[1]/100*pow(E, -alpha)
  // sb = par[2] * (E/E0)^{-2.1} <--- this -2.1 is put into the ss for computational reasons
  Double_t kk      = (alpha-1)*sFlux/100;  // differential flux
  Double_t ss     = kk * pow(E/100, -alpha+2.1);   // to go from norm
						    // to differential
  Double_t sb=BackgroundFlux(E,l,b,normeg,normgal,indexgal);
  int imax=50;
  double sum=0;  
  double tstep=3.*thet_68(E)*TMath::Pi()/(180.*imax);
  for (int mm=0;mm<imax;mm++) { 
     double r=(mm+.5)*tstep; 			//rad
     for (int i=0;i<2;i++){  
        Double_t g=psf(E,r,i)*ss/sb; //PSF
        Double_t I=(1+g)*log(1+g)-g; 
        Double_t integrandOmega = 2*TMath::Pi()*I*r; //(dOmega == 2 Pi dr)
        Double_t integrandE = Exposure(E, l, b,i)*sb*E*log(10.)*pow(E/100,-2.1);
        sum+=integrandE * integrandOmega*tstep;
     }
  }
  return 2.*time*sum*par[8]/(6.*365.25*86400.*25); 
}

double Integr(TF2* f,double ymin,double ymax) {
// FILE *f60 = fopen("essai_10.dat","w");
//  int nstep=50;
  int nstep=50;
  int imax=50;
  double estep=(ymax-ymin)/nstep;
  double sum=0;  
  for (int k=0;k<nstep-1;k++){ 				//loop over log10(energy), min=100 MeV
   double y=ymin+(k+0.5)*estep;
   double E=pow(10,y); 
   double tstep=3.*thet_68(E)*TMath::Pi()/(180.*imax);
   for (int mm=0;mm<imax;mm++) { 
     double theta=(mm+.5)*tstep; 			//rad
     sum+=f->Eval(y,theta)*tstep;
  //     prov+=f->Eval(y,theta)*tstep;
   }
//fprintf(f60,"%e %e \n",y, prov);
 }
//  fclose(f60);
   sum*=estep;
   return sum;
}
double significance(double flux, double index,double l, double b){ 
  // flux: integral flux [E> 100 MeV] ph cm^-2 s^-1
  // index: spectral index (positive)
  // time: accumulation time (s)
  // l: galactic longitude (deg)
  // b: galactic latitude (deg)
  Double_t pf[6];    
  TF2 *f1=new TF2("f1",TSKernel,1,6,0,1.5,6);
  pf[0]=index; //index
  pf[1]=flux; 
  pf[2]=l;
  pf[3]=b; 
  pf[4] =1;
  pf[5]=1;
  f1->SetParameters(&pf[0]); 
  double elim=log10(100.);
  double elmax=log10(2e5);
  double TS=2.*Integr(f1,elim,elmax); 
  return TS;
}
Double_t S_1(Double_t *x, Double_t *par){ 
  Double_t r=x[1];
  Double_t l      = par[2];
  Double_t b      = par[3]; 
  Double_t E=pow(10,x[0]);
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100;//
  bool addGalactic = true;
  Double_t sb=BackgroundFlux(E,l,b, addGalactic); //diffuse bgd
  Double_t sum=0;
  for (int i=0;i<2;i++){
      Double_t g=psf(E,r,i)*ss*pow(E/100,-alpha+2.1)/sb;
  Double_t I=g*g/(1+g); 
  sum+=E*sb*pow(E/100,-2.1)*log(10.)*I*r*pi2*Exposure(E,l,b,i);
  }
  return sum;
}
Double_t S_AG(Double_t *x, Double_t *par){ 
  Double_t r=x[1];
  Double_t l      = par[2];
  Double_t b      = par[3]; 
  Double_t E=pow(10,x[0]);
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100;//
  bool addGalactic = true;
  Double_t sb=BackgroundFlux(E,l,b, addGalactic); //diffuse bgd
  Double_t sum=0;
  for (int i=0;i<2;i++){
 Double_t g=psf(E,r,i)*ss*pow(E/100,-alpha+2.1)/sb;
 Double_t I=g*g/(1+g); 
  sum+=E*sb*pow(E/100,-2.1)*log(10.)*I*r*pi2*Exposure(E,l,b,i);
  }
  return sum*x[0];
}
Double_t S_GG(Double_t *x, Double_t *par){ 
  Double_t r=x[1];
  Double_t l      = par[2];
  Double_t b      = par[3]; 
  Double_t E=pow(10,x[0]);
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*par[1]/100;//
  bool addGalactic = true;
  Double_t sb=BackgroundFlux(E,l,b, addGalactic); //diffuse bgd
  Double_t sum=0;
  for (int i=0;i<2;i++){
      Double_t g=psf(E,r,i)*ss*pow(E/100,-alpha+2.1)/sb;
      Double_t I=g*g/(1+g);
   sum+=E*sb*pow(E/100,-2.1)*log(10.)*I*r*pi2*Exposure(E,l,b,i);
  }
  return sum*log(E/par[4])*log(E/par[4]);
}
 
float significance_2(float flux, float index,float time,float l, float b){ 
  // flux: integral flux [E> 100 MeV] ph cm^-2 s^-1
  // index: spectral index (positive)
  // time: accumulation time (in s)
  // l: galactic longitude (deg)
  // b: galactic latitude (deg)
  Double_t pf[6]; 
  TF2 *f1=new TF2("f1",TSKernel,1,6,0,1.5,6);
  pf[0]=index; //index
  pf[1]=flux; 
  pf[2]=l;
  pf[3]=b; 
  pf[4] =1.;
  pf[5]=1;
  //printf("time %lf\n",time);
  f1->SetParameters(&pf[0]); 
  double elim=log10(100.);
  double elmax=log10(1e5);
  float TS=time*2.*Integr(f1,elim,elmax)/(6.*365.25*86400); 
  return TS;
}

double uncertainty(double flux, double index,double time, double l, double b, double emin=100){
  Double_t pf[6];
  double emax=1e5;
  pf[0]=index; //index
  pf[1]=flux;
  pf[2]=l;
  pf[3]=b;
  double dft=time/(6.*365.25*86400);
  TF2 *f2=new TF2("f2",S_1,2,6,0,1.5,4);
  f2->SetParameters(&pf[0]);
  TF2 *f_AG=new TF2("f_AG",S_AG,2,6,0,1.5,4);
  f_AG->SetParameters(&pf[0]);
  double elim=log10(emin);
  double elax=log10(emax);
  Double_t s1=Integr(f2,elim,elax)*dft;
  Double_t s=Integr(f_AG,elim,elax)*dft;
  Double_t E0=pow(10,s/s1);
  //Double_t E1=E0*exp(-1./(index-1));
  TF2 *f_GG=new TF2("f_GG",S_GG,1,6,0,1.5,5);
  pf[4]=E0;
  f_GG->SetParameters(&pf[0]);
  double sg=Integr(f_GG,elim,elax)*dft;
  double var_gam=1./sg;
  double df100=100*sqrt(1./s1+pow(log(emin/E0)+1./(index-1),2)*var_gam);
  return df100;
}
double uncertainty_opt(double flux, double index,double time, double l, double b, double epivot=100){
  Double_t pf[6];
  double emax=1e5;
  double emin=100;
  pf[0]=index; //index
  pf[1]=flux;
  pf[2]=l;
  pf[3]=b;
  double dft=time/(6.*365.25*86400);
  TF2 *f2=new TF2("f2",S_1,2,6,0,1.5,4);
  f2->SetParameters(&pf[0]);
  TF2 *f_AG=new TF2("f_AG",S_AG,2,6,0,1.5,4);
  f_AG->SetParameters(&pf[0]);
  double elim=log10(emin);
  double elax=log10(emax);
  Double_t s1=Integr(f2,elim,elax)*dft;
  Double_t s=Integr(f_AG,elim,elax)*dft;
  Double_t E0=pow(10,s/s1);
  //Double_t E1=E0*exp(-1./(index-1));
  TF2 *f_GG=new TF2("f_GG",S_GG,1,6,0,1.5,5);
  pf[4]=E0;
  f_GG->SetParameters(&pf[0]);
  double sg=Integr(f_GG,elim,elax)*dft;
  double var_gam=1./sg;
  double df100=100*sqrt(1./s1+pow(log(epivot/E0)+1./(index-1),2)*var_gam);
  return df100;
}
void Npred_ini(double fl100,double gamma,double GLON,double GLAT){
  initialize_TS();
    double pf[4];
    pf[0]=gamma;
    pf[1]=fl100;
    pf[2]=GLON;
    pf[3]=GLAT;
    TF1 *fn=new TF1("fn",Np,2.,6,4);
    fn->SetParameters(&pf[0]);
    TF2 *fn2=new TF2("fn2",Np2,2.,6,0,1,4);
    fn2->SetParameters(&pf[0]);
    printf("Np: %f %f\n",fn->Integral(2,6),Integr(fn2,2,6));
    printf("Np2: %f ratio: %f\n",fn2->Integral(2.,6.,0.,1./57.6),fn2->Integral(2.,6.,0.,1./57.6)/Integr(fn2,2.,6.));
/*TF2 *f1=new TF2("f1",TSKernel,1.,6.,0.,1.5,4);
    for (int i=0;i<10;i++){
    pf[1]=pow(10,-1+i*0.2)*fl100;    
    f1->SetParameters(&pf[0]);
    printf("flux %e TS: %f %f ratio: %f\n",pf[1],f1->Integral(2.,6.,0.,1.5),f1->Integral(2.,6.,0.,1./57.6),1./(f1->Integral(2.,6.,0.,1./57.6)/f1->Integral(2.,6.,0.,1.5)));
}
*/
}

Double_t Np_bgd(Double_t *x, Double_t *par){
  Double_t E=pow(10,x[0]); 
  Double_t l      = par[0];
  Double_t b      = par[1]; 
  Double_t normeg= par[2];
  Double_t normgal= par[3];
  Double_t index_gal= par[4];
  Double_t sum=0;
  for (int i=0;i<2;i++){  
    sum+=pow(E/100.,-2.1)*BackgroundFlux(E,l,b,normeg,normgal,index_gal)*Exposure(E,l,b,i)*E*log(10);
 }
  return sum/(6.*365.25*86400.);
 }
double Npred(double fl100,double gamma,double time,double GLON,double GLAT, double emin=2, double emax=6){
    double pf[4];
    pf[0]=gamma;
    pf[1]=fl100;
    pf[2]=GLON;
    pf[3]=GLAT;
    TF1 *fn=new TF1("fn",Np,2.,6,4);
    fn->SetParameters(&pf[0]);
    return fn->Integral(emin,emax)*time/(6.*365.25*86400.);   

}

double Npred_front(double fl100,double gamma,double time,double GLON,double GLAT){
    double pf[4];
    pf[0]=gamma;
    pf[1]=fl100;
    pf[2]=GLON;
    pf[3]=GLAT;
    TF1 *fn=new TF1("fn",Np_front,2.,6,4);
    fn->SetParameters(&pf[0]);
    return fn->Integral(2,6)*time/(6.*365.25*86400.);   

}
double Npred_back(double fl100,double gamma,double time,double GLON,double GLAT){
    double pf[4];
    pf[0]=gamma;
    pf[1]=fl100;
    pf[2]=GLON;
    pf[3]=GLAT;
    TF1 *fn=new TF1("fn",Np_back,2.,6,4);
    fn->SetParameters(&pf[0]);
    return fn->Integral(2,6)*time/(6.*365.25*86400.);   

}



// ========== Lise ============ //

void Glob(double *vE, double *vTime, double *vConvtype, double *vr,int evt){
  for (int i=0;i<evt;i++){
        double E=vE[i];
        double Time=vTime[i];
        int Convtype=vConvtype[i];
	double r=vr[i];
	GTime[i]=Time;
	GE[i]=E;
	Gr[i]=r;
	GConvtype[i]=Convtype;}
  Gevt=evt;
  return;
	}

void Glob_2(double *vE, double *vTime, double *vConvtype, double *vr,double *vlon,double *vlat,int evt){ 
  for (int i=0;i<evt;i++){
        double E=vE[i];
        double Time=vTime[i];
        int Convtype=vConvtype[i];
	double r=vr[i];
	GTime[i]=Time;
	GE[i]=E;
	Gr[i]=r;
	GConvtype[i]=Convtype;
        Glon[i]=vlon[i]; 
        Glat[i]=vlat[i];

}
  Gevt=evt;
  return;
	}

void Glob_mc(double *vE, double *vTime, double *vConvtype, double *vr,double *vid,int evt){
 
  for (int i=0;i<evt;i++){
        double E=vE[i];
        double Time=vTime[i];
        int Convtype=vConvtype[i];
	double r=vr[i];
	GTime[i]=Time;
	GE[i]=E;
	Gr[i]=r;
	GConvtype[i]=Convtype;
        Gid[i]=vid[i];
}
  Gevt=evt;
  printf("%d evts in Glob \n",Gevt);
  //int i=Gevt-1;
  //printf("%lf %lf %d %lf %d\n",GE[i],GTime[i],GConvtype[i],Gr[i],Gid[i]);
  return;
	}

void Get_alpha(double *valpha, double *vTstart,
double *vTstop, double *vFlux, double *vIndex, double *vTS, double *vNormGal, double *vNormEG, 
double *vIndexGal, double *vbeta, double *vNph, double L, double B, int bin, int evt){
  double alpha0=1; 
  for (int k=0;k<bin;k++){
     double Tstart=vTstart[k];
     double Tstop=vTstop[k]; 
     double Flux=vFlux[k];
     double Index=vIndex[k];
     double NormGal=vNormGal[k];
     double NormEG=vNormEG[k];
     double IndexGal=vIndexGal[k];
     double beta=vbeta[k];
     double Nph=vNph[k];
     double ss=(Index-1)*Flux/100;
     double ll=0;
     double g=0;
     int nn=0;
     for (int i=0;i<evt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	
	///printf("NormEG=%lf\n",NormEG);
	//printf("NormGal=%lf\n",NormGal);
	//printf("beta=%lf\n",beta);
	//printf("\n");
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    nn++;
	    //printf("g=%lf\n",g);
	    //printf("ll=%lf\n",ll);
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*alpha0*(ll-beta*Nph);
     //printf("TSe=%lf ll=%lf beta:%lf Nph:%lf nn=%d\n",TSe,ll,beta,Nph,nn);
     valpha[k]=vTS[k]/TSe;
     //printf("valpha[k]=%lf\n",valpha[k]);	
	}  
  return;   
 }


double Time_targ(double Index, double Flux, double L, double B, double NormEG, double NormGal,
double IndexGal, double alpha, double beta, double np, double Ti, double TStarg, int evt, int varalpha){
  double ll=0,TSe=0,time_targ=0;
  int nn=0;
  //printf("%lf %lf",TStarg,alpha);
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti){
     nn++;
     double Nph_p=np*(GTime[m]-Ti);
     double E=GE[m];
     double ss=(Index-1)*Flux/100;    
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
     if (varalpha==0) {
        ll+=log(1+g);
        TSe=2*alpha*(ll-beta*Nph_p);
       }
     else{
        ll+=log(1+alpha*g);
	TSe=2*(ll-beta*Nph_p);
        }  
     if (TSe>TStarg){
        time_targ=GTime[m+1];
	break;
     }
     }
     }
  //printf("TSe=%lf nn %d %lf %lf\n",TSe,nn,alpha,beta);        
  return time_targ;	   
 } 
double TS_alpha_g(double *x, double *par){
  double ll=0;
  double Tstart=par[0];
  double Tstop=par[1];
  double Flux=par[2];
  double Index=par[3];
  double NormGal=par[4];
  double NormEG=par[5];
  double IndexGal=par[6];
  double beta=par[7];
  double Nph=par[8];
  int evt=par[9],nn=0;
  double L=par[10];
  double B=par[11];
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<evt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+x[0]*g);
	    nn++;
	  
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-beta*Nph);
     printf("TSe: %lf \n", TSe);
  return TSe;
}

void Get_alpha_g(double *valpha, double *vTstart,
double *vTstop, double *vFlux, double *vIndex, double *vTS, double *vNormGal, double *vNormEG, 
double *vIndexGal, double *vbeta, double *vNph, double L, double B, int bin, int evt){
   double par[12];
   par[9]=evt;
   par[10]=L;
   par[11]=B;
   TF1* f1=new TF1("f1",TS_alpha_g,0.5,1.5,12);
   for (int k=0;k<bin;k++){
     par[0]=vTstart[k];
     par[1]=vTstop[k]; 
     par[2]=vFlux[k];
     par[3]=vIndex[k];
     par[4]=vNormGal[k];
     par[5]=vNormEG[k];
     par[6]=vIndexGal[k];
     par[7]=vbeta[k];
     par[8]=vNph[k];
     f1->SetParameters(par);
     //printf("vTS[k]: %lf \n", vTS[k]);
     valpha[k]=f1->GetX(vTS[k],0.8,1.2);    
	}               
   return;      
}   
double Time_targ2(double Index, double Flux, double L, double B, double NormEG, double NormGal,
double IndexGal, double Ti, double TStarg, int evt){
  double TSe=0,time_targ=0;
  int nn=0;
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
     TSe+=2*(log(1+g)-g/(1+g));
     //printf("TSe=%lf\n",TSe);
     if (TSe>TStarg){
        time_targ=GTime[m+1];
        break;
     }
     }
     }
  return time_targ;
 }

double Get_TS(double Index, double Flux, double L, double B, double Ti, double To,double normeg=1,double normgal=1,double indgal=0){
  double TSe=0;
  int nn=0;
  for (int m=0;m<Gevt-1;m++){
   if (GTime[m]>Ti && GTime[m]<To){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     TSe+=2*(log(1+g)-g/(1+g));
     }
     }
  return TSe;
 }

double Get_time_TS(double TS,double Index, double Flux, double L, double B, double Ti, double normeg=1,double normgal=1,double indgal=0){
  double TSe=0;
  int nn=0;
  for (int m=0;m<Gevt-1;m++){
   if (GTime[m]>Ti){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     TSe+=2*(log(1+g)-g/(1+g));
   if (TSe>TS) return GTime[m];
  }
     }
  return -1;
 }


double Get_unc(double Index, double Flux, double L, double B, double Ti, double To, int evt,double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0;
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti && GTime[m]<To){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     SigA+=g*g/((1+g)*(1+g));
     Sln+=log(E)*g*g/((1+g)*(1+g));
     Sln2+=log(E)*log(E)*g*g/((1+g)*(1+g));
     if (g<1e-12) continue;
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(Elowerl/Epivot)+1/(Index-1))*(log(Elowerl/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     }
   if (GTime[m]>To) break;
  }
  return 100.*errK;
 }
double Get_unc_mod(double Index, double Flux, double L, double B, double Ti, double To, int evt,double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0;
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti && GTime[m]<To){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     L=Glon[m];
     B=Glat[m];
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     SigA+=g*g/((1+g)*(1+g));
     Sln+=log(E)*g*g/((1+g)*(1+g));
     Sln2+=log(E)*log(E)*g*g/((1+g)*(1+g));
     if (g<1e-12) continue;
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(Elowerl/Epivot)+1/(Index-1))*(log(Elowerl/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     }
   if (GTime[m]>To) break;
  }
  return 100.*errK;
 }


double Get_time_unc(double unc,double Index, double Flux, double L, double B, double Ti, double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0;
  for (int m=0;m<Gevt;m++){
   if (GTime[m]>Ti){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     SigA+=g*g/((1+g)*(1+g));
     Sln+=log(E)*g*g/((1+g)*(1+g));
     Sln2+=log(E)*log(E)*g*g/((1+g)*(1+g));
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(Elowerl/Epivot)+1/(Index-1))*(log(Elowerl/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     if (errK <unc/100.) return GTime[m];
   }
  }
  return -1;
 }


double plot_unc(double Index, double Flux, double L, double B, double Ti, double To, int evt){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0;
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti && GTime[m]<To){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B);
     SigA+=g*g/((1+g)*(1+g));
     Sln+=log(E)*g*g/((1+g)*(1+g));
     Sln2+=log(E)*log(E)*g*g/((1+g)*(1+g));
     if (g<1e-12) continue;
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(100/Epivot)+1/(Index-1))*(log(100/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     }
     }
  return 100.*errK;
 }
void cut_flare(){
  gROOT->SetStyle("Plain");
  initialize_TS();
  float index;
  float time=1;
  float l=40;
  float b=40;
  double vx[6][101],vy[6][101];
  double TS0=50;
  int i;
  for (int j=0;j<6;j++){
    index=1.8+j*0.2;
   for (i=0;i<101;i++){
    float flux=1e-8*pow(10,0.03*i);
    double TS=significance_2(flux,index,time,l,b);
    vx[j][i]=flux;
    vy[j][i]=TS0/TS;
    //printf("%d %e %e\n",i,flux,TS0/TS);

   }  
  }
  TCanvas* c1 =new TCanvas( "c1", "", 0, 0, 650, 500 );
  TH2F* frame4=new TH2F("frame4","",100,1e-8,1e-5,100,1e3,1e9);
  gPad->SetLogy(); 
  gPad->SetLogx(); 
  frame4->SetStats(0);
  frame4->SetYTitle("Time (s)");
  frame4->SetXTitle("Flux[E>100 MeV] (ph cm^{-2}s^{-1})");
  frame4->GetXaxis()->CenterTitle(1);
  frame4->GetYaxis()->CenterTitle(1);
  frame4->Draw();
  TGraph *grf[6];
  for (int j=0;j<6;j++){
    grf[j]= new TGraph(101,&vx[j][0],&vy[j][0]);
    grf[j]->SetLineWidth(2);
    grf[j]->SetLineColor(2); 
    grf[j]->Draw("L");
  }
  c1->Update();
}
double TS_flux(double *x, double *par){
  double ll=0;
  double Flux=pow(10,x[0]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
  return TSe;
}

double nnp(double *x, double *par){
  double Flux=pow(10,x[0]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  //double ex=par[6];
  double L=par[7];
  double B=par[8];
  //int nn=0;
  double np=0; 
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){

	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    np+=g/(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
  return np;
}
double Get_Np(double index, double flux, double l, double b, double normeg,double normgal,double indgal,double Tstart, double Tstop, double ex){
  double kdum= normeg+normgal-indgal+Tstart-Tstop;
  kdum=0;
  return ex*Npred(flux,index,6.*365.25*86400,l,b)/Exposure(1000,l,b,0);
  /*
   double par[9];
   TF1* f2=new TF1("f2",nnp,-10,-3,9);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=index;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f2->SetParameters(par);
   return f2->Eval(log10(flux));
  */
}

int Get_Np_mc(double Tstart,double Tstop, int src_id){
          int nn=0;
          for (int i=0;i<Gevt;i++){
            double Time=GTime[i];
	    if ((Tstart<=Time) && (Time<Tstop)){
	       if (Gid[i]==src_id) nn++;
	    }
            if (Time>Tstop) break;
	  }
  return nn;
}
int Get_Np_mc_type(double Tstart,double Tstop, int src_id, int itype){
          int nn=0;
          for (int i=0;i<Gevt;i++){
            double Time=GTime[i];
	    if ((Tstart<=Time) && (Time<Tstop)){
	       if (Gid[i]==src_id and GConvtype[i]==itype) nn++;
	    }
            if (Time>Tstop) break;
	  }
  return nn;
}



int Get_Np_tot(double Tstart,double Tstop){
          int nn=0;
          for (int i=0;i<Gevt;i++){
            double Time=GTime[i];
	    if ((Tstart<=Time) && (Time<Tstop)){
	       nn++;
	    }
            if (Time>Tstop) break;
	  }
  return nn;
}

double Get_TS0(double index, double flux,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex){
   double par[9];
   TF1* f1=new TF1("f1",TS_flux,-10,-3,9);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=index;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f1->SetParameters(par);
   return  f1->Eval(log10(flux));      
} 

double Get_flux(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex){
   double par[9];
   TF1* f1=new TF1("f1",TS_flux,-10,-3,9);
   TF1* f2=new TF1("f2",nnp,-10,-3,9);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=index;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f1->SetParameters(par);
   f2->SetParameters(par);
   //for (int i=0;i<9;i++){
   //  printf("%d %e %e %e %d\n",i,log10(5e-7*pow(1.2,i-4)),f1->Eval(log10(5e-7*pow(1.2,i-4))),f2->Eval(log10(5e-7*pow(1.2,i-4))),get_Np_mc(Tstart,Tstop,1001));   }
   //printf("min:%e \n",f1->GetMaximumX(-9,-4));
   return  pow(10.,f1->GetMaximumX(-9,-4)) ;      
} 


double Get_errflux(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex){
   double par[9];
   TF1* f1=new TF1("f1",TS_flux,-9,-6,9);
   par[0]=Tstart;
   par[1]=Tstop;
   par[2]=index;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f1->SetParameters(par);
   //TCanvas* ce2 = new TCanvas("ce2","", 0, 0, 500, 500);
   //ce2->SetLogy();
   //f1->Draw();
   //ce2->Update();
   double xmin=f1->GetMaximumX(-9,-4);
   //printf("flux=%e\n",pow(10.,xmin)) ;
   double TSmin=f1->Eval(xmin);
   //printf("TSmin=%lf\n",TSmin);
   double TSsup=TSmin-1;
   //printf("TSsup=%lf\n",TSsup);
   double fmin=pow(10,f1->GetX(TSsup,-9,xmin));
   double fmax=pow(10,f1->GetX(TSsup,xmin,-4));
   //printf("fmin=%e\n",fmin);
   //printf("fmax=%e\n",fmax);
   double sigf=(fmax-fmin)/2;
   return 100.*sigf/pow(10.,xmin);
}
double TS_flux_lin(double *x, double *par){
  double ll=0;
  double Flux=x[0];
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
  return TSe;
}



double Get_flux_lin(double flux,double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex){
   double par[9];
  
   TF1* f1=new TF1("f1",TS_flux_lin,flux/2.,flux*2,9);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=index;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f1->SetParameters(par);
   return  f1->GetMaximumX(flux/2.,flux*2.);      
} 

double minus_TS_flux(double *x, double *par){
  double ll=0;
  double Flux=pow(10,x[0]); 
  double Index=x[1];
  double Tstart=par[0];
  double Tstop=par[1];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  int nn=0;
  double ss=(Index-1)*Flux/100;  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];;
      	double r=Gr[i]; 
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    nn++;  
          }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
  return -TSe;
}

double minus_TS_flux_3D(double *x, double *par){
  double ll=0;
  double Flux=pow(10,x[0]); 
  double Index=x[1];
  double Tstart=par[0];
  double Tstop=par[1];
  double NormGal=par[3];
  double NormEG=x[2];
  double IndexGal=par[4];
  double ex=par[5];
  double L=par[6];
  double B=par[7];
  int nn=0;
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    nn++;
	  
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
  return -TSe;
}
double Get_flux_index(double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex, double* ind){
   double par[9];
   if (normgal>0) {
   TF2* f2=new TF2("f2",minus_TS_flux,-10,-3,1,3,9);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=0;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f2->SetParameters(par); 
   f2->SetRange(-9.,1., -3.,3.);
   double xmin;    
   f2->GetMinimumXY(xmin,*ind);
   printf("%e %e\n",xmin,*ind);
   return  pow(10,xmin);
   }
   else
     {
       TF3* f3=new TF3("f3",minus_TS_flux_3D,-10,-3,1,3,0,2,8);
   par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=0;
   par[3]=normgal;
   par[4]=indgal;
   par[5]=ex;
   par[6]=l;
   par[7]=b;
   f3->SetParameters(par); 
   f3->SetRange(-9.,1., -3.,3.,0,2);
   double xmin,zmin;     
   f3->GetMinimumXYZ(xmin,*ind,zmin);
   printf("%e %e %e\n",xmin,*ind,zmin);
   return  pow(10,xmin);
     }
}
double Npred_bgd(double l,double b,double normeg,double normgal,double indgal, double emin=2,double emax=5){
   double par[9];
   par[0]=l;
   par[1]=b;
   par[2]=normeg;
   par[3]=normgal;
   par[4]=indgal;
   TF1* fnb=new TF1("fnb",Np_bgd,2,6,5);
   fnb->SetParameters(par);
   return  fnb->Integral(emin,emax)*6.*365.25*86400.;
      }


double get_E1(double flux, double index, double l,double b, double emin=100, double emax=3e5){
  Double_t pf[6];
  pf[0]=index; 
  pf[1]=flux*pow(emin/100.,index-1);
  pf[2]=l;  
  pf[3]=b;
  TF2 *f2=new TF2("f2",S_1,2,6,0,1.5,4);
  f2->SetParameters(&pf[0]);
  TF2 *f_AG=new TF2("f_AG",S_AG,2,6,0,1.5,4);
  f_AG->SetParameters(&pf[0]);
  double elim=log10(emin);
  double elax=log10(emax);
  Double_t s1=Integr(f2,elim,elax);
  Double_t s=Integr(f_AG,elim,elax);
  Double_t E0=pow(10,s/s1);
  return E0*exp(-1./(index-1));
}
double Time_targ_ErF(double Index, double Flux, double L, double B, double NormEG, double NormGal,
double IndexGal, double Ti, double errKtarg, int evt, double El=100){
  printf("Ti begin cumul = %lf\n",Ti);
  printf("evt=%d\n",evt);
  printf("index=%lf\n",Index);
  printf("flux=%e\n",Flux);
  double SigA=0,time_targ_ErF=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0;
  for (int m=0;m<evt-1;m++){
   if (GTime[m]>Ti){
     nn++;
     double E=GE[m];
     double ss=(Index-1)*Flux/100;    
     double g=psf(E,Gr[m],GConvtype[m])*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
     //printf("g=%e %e %e %e\n",g,E,ss,psf(E,Gr[m],GConvtype[m]));
     SigA+=g*g/((1+g)*(1+g));
     //printf("SigA=%e\n",SigA);
     Sln+=log(E)*g*g/((1+g)*(1+g));
     Sln2+=log(E)*log(E)*g*g/((1+g)*(1+g));
     if (g<1e-12) continue;
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(El/Epivot)+1/(Index-1))*(log(El/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     //printf("m=%d\n",m);
     //printf("errKtarg=%lf\n",errKtarg);
     //printf("errK=%lf\n",errK);
     //printf("m-2=%d\n",m-2);
     //printf("evt=%d\n",evt);     
     if (errK<errKtarg){
        time_targ_ErF=GTime[m+1];
	break;
     }
     if (((m==evt-3) or (m==evt-2) or (m==evt-1)) && (errK>errKtarg)){
     time_targ_ErF=-1;
     printf("errF>errFtarg even with all available photons\n");
     printf("GTime[m+1]=%lf\n",GTime[m+1]);
     }
     }
   //else{
   //time_targ_ErF=Ti;
   //}  
    
   //else{
  //printf("errK=%lf errkTarg=%lf time_targ=%lf\n", errK,errKtarg,time_targ_ErF);
  //}
  }
  //printf("Time_targ_ErF=%lf\n",time_targ_ErF);
  return time_targ_ErF;	   
 } 
 
double der_TS(double *x, double *par){
  double ll=0;
  double Flux=pow(10,x[0]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=g/(Flux*(1+g));
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     return ll-ex*Npred(1,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0);
}

void int_diff(){
  //initialize_TS(); 
  //Load_Isot("../COSPAR/isotropic_iem_v02.txt");
  Load_Isot("../lightcurves/isotrop_2year_P76_clean_v0.txt");
  double sum=0;
  double vy[30],vx[30];
  printf("%d\n",ind_enisot);
  for (int i=0;i<ind_enisot-1;i++){
     double x0=enisot[i];
     double x1=enisot[i+1];
     vx[i]=log10(enisot[i]);
     vx[i+1]=log10(enisot[i+1]);
     vy[i]=fluxisot[i]*pow(enisot[i]/100.,-2.1); 
     vy[i+1]=fluxisot[i+1]*pow(enisot[i+1]/100.,-2.1); 
     double y0=fluxisot[i]*pow(enisot[i]/100.,-2.1); 
     double y1=fluxisot[i+1]*pow(enisot[i+1]/100.,-2.1);
     double minus_alpha=log(y1/y0)/log(x1/x0);
     sum+=(y1*x1-x0*y0)/(minus_alpha+1);
  }
  TGraph* gr=new TGraph(ind_enisot,vx,vy);
  gr->Draw("AL");
  printf("integral about %lf MeV: %lf\n",enisot[0],sum*2*pi2);

}

static void fcn_TS(int &npar, double *gin, double &f, double *par, int
iflag){
  double kdum=npar+gin[0]+f+iflag;
  kdum=0; 
  double ll=0;
  double Flux=pow(10,par[9]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  //double radius=5./57.3;
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
            //if (r< radius) ll2+=log(BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal));
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
     //TSe+=2*(ll2-ex*3.1415927*pow(radius,2.)*Npred_bgd(L,B,NormEG,NormGal,IndexGal)/Exposure(1000,L,B,0)); 
  f=-TSe;
  //printf("fff %e %lf %lf %lf %lf %lf %lf %lf %lf %lf %e\n",f,Tstart,Tstop,Index, NormGal,NormEG,IndexGal,ex,L,B,Flux);
}


double Get_flux_mult_minuit(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex,double* out){
   TMinuit minuit(10);
   int ierflg = 0;
   minuit.SetFCN(fcn_TS);
   minuit.mnparm(0, "TSTART",Tstart,0.01,1e8,1e9,ierflg);
   minuit.FixParameter(0); 
   minuit.mnparm(1, "TSTOP",Tstop,0.01,1e8,1e9,ierflg);
   minuit.FixParameter(1);
   minuit.mnparm(2, "Index",index,0.01,1,5,ierflg);
   minuit.FixParameter(2);
   minuit.mnparm(3, "normgal",normgal,0.01,0,5,ierflg);
   minuit.FixParameter(3);
   minuit.mnparm(4, "normeg",normeg,0.01,0,5,ierflg);
   minuit.FixParameter(4); 
   minuit.mnparm(5, "indgal",indgal,0.01,-2,2,ierflg);
   minuit.FixParameter(5);
   minuit.mnparm(6, "exposure",ex,0.01,0,1e11,ierflg);
   minuit.FixParameter(6); 
   minuit.mnparm(7, "L",l,0.01,0,360,ierflg);
   minuit.FixParameter(7);
   minuit.mnparm(8, "B",b,0.01,-90,90,ierflg);
   minuit.FixParameter(8);
   minuit.mnparm(9, "Flux", -7,0.01,-9,-4,ierflg);
   minuit.SetErrorDef(0.5);
   minuit.SetPrintLevel(-1);
   //minuit.SetMaxIteration(500);
   minuit.Migrad();  
   double flux,errflux;
   //,errindex,errnormeg,errnormgal;
   minuit.GetParameter(9,flux,errflux); 
   //minuit.GetParameter(2,index,errindex);
   //minuit.GetParameter(4,normeg,errnormeg);
   //minuit.GetParameter(3,normgal,errnormgal);
   //printf("index:%lf\n",index); 
   //printf("normeg:%e\n",normeg); 
   //printf("normgal:%e\n",normgal); 
   *out=normeg; 
   return pow(10,flux);    
}
double TS_bgd(double *x, double *par){
  double ll2=0;
  double Tstart=par[0];
  double Tstop=par[1];
  double NormGal=par[3];
  double NormEG=x[0];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  //int nn=0;
  double radius=5./57.3;
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        //int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){ 
	  if (r< radius) {ll2+=log(BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal));
	    // printf("%e %e\n",ll2, BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal));  

}
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll2-ex*3.1415927*pow(radius,2.)*Npred_bgd(L,B,NormEG,NormGal,IndexGal)/Exposure(1000,L,B,0)); 
     //printf("%lf %lf %lf %lf %lf %lf,%lf %lf\n",ll2, TSe, Npred_bgd(L,B,NormEG,NormGal,IndexGal),L,B,NormEG,NormGal,IndexGal);
     printf("%lf \n",ex*3.1415927*pow(radius,2.)*Npred_bgd(L,B,1.47,NormGal,IndexGal)/Exposure(1000,L,B,0));
  return TSe;
}

double TS_flux_bgd(double *x, double *par){
  double ll=0,ll2=0; 
  double Flux=pow(10,x[0]);
  double Tstart=par[0];
  double Tstop=par[1]; 
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=x[1];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  //int nn=0;
  double radius=5./57.3;
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
            if (r< radius) ll2+=log(BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal));
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
     TSe+=2*(ll2-ex*3.1415927*pow(radius,2.)*Npred_bgd(L,B,NormEG,NormGal,IndexGal)/Exposure(1000,L,B,0)); 
  return TSe;
}
double TS_flux_bgd_1D(double *x, double *par){
  double ll=0,ll2=0; 
  double Flux=par[4];
  double Tstart=par[0];
  double Tstop=par[1]; 
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=x[0];
  double IndexGal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  //int nn=0;
  double radius=5./57.3;
  double ss=(Index-1)*Flux/100;  
  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	    double g=psf(E,r,Convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
            if (r< radius) ll2+=log(BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal));
	    }
	if (Time>Tstop){
	    break;
	    }
	    }
     double TSe=2*(ll-ex*Npred(Flux,Index,6.*365.25*86400,L,B)/Exposure(1000,L,B,0));
     TSe+=2*(ll2-ex*3.1415927*pow(radius,2.)*Npred_bgd(L,B,NormEG,NormGal,IndexGal)/Exposure(1000,L,B,0)); 
  return TSe;
}

double Get_bgd(double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex){
  double par[9];
  TF1* f1=new TF1("f1",TS_bgd,0,5,9);
 par[0]=Tstart;
   par[1]=Tstop; 
   par[2]=0;
   par[3]=normgal;
   par[4]=normeg;
   par[5]=indgal;
   par[6]=ex;
   par[7]=l;
   par[8]=b;
   f1->SetParameters(par);
   f1->Draw(); 
return f1->GetMaximumX(0,2);
} 
//*****************************

void Glob_mult(double *vE, double *vTime, double *vConvtype, double *vr,int evt,int N_source=1){
  n_source=N_source;
  printf("n_source %d \n",n_source);
  for (int i=0;i<evt;i++){
        double E=vE[i];
        double Time=vTime[i];
        int Convtype=vConvtype[i];
	GTime[i]=Time;
	GE[i]=E;
        for (int is=0;is<n_source;is++) Grs[i][is]=vr[i*n_source+is];
        if (n_source==1) Gr[i]=vr[i];
	GConvtype[i]=Convtype;}
  Gevt=evt;
  return;
	}

void feed_sources(double *vglon,double *vglat, double *vindex, double *vflux,int N_source)
{
    n_source=N_source;
    for (int is=0;is<n_source;is++) {
      glon_s[is]=vglon[is];
      glat_s[is]=vglat[is];
      Index_source[is]=vindex[is];
      Flux_source[is]=vflux[is];
      printf("%lf %lf %lf %lf\n", glon_s[is],glat_s[is], Index_source[is],Flux_source[is]);
      }
}


double Get_unc_mult(double Index, double Flux, double L, double B, double Ti, double To, int evt,double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0,ss,g,r,sgprim;
  for (int m=0;m<evt-1;m++){
    //if (Gid[m]==0) continue;
   if (GTime[m]>Ti && GTime[m]<To){
     nn++;
     double E=GE[m];
     int Convtype=GConvtype[m]; 
     sgprim=0; 
     for (int is=1;is<n_source;is++){ 
               r=Grs[m][is];
               ss=(Index_source[is]-1)*Flux_source[is]*pow(E/100,-Index_source[is]+2.1)/100;
               
	       sgprim+=psf(E,r,Convtype)*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal);}
     r=Grs[m][0]; 
     ss=(Index-1)*Flux*pow(E/100,-Index+2.1)/100;
     g=psf(E,r,Convtype)*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     //sgprim=0;
     double gprim=g+sgprim;
     SigA+=g*g/((1+gprim)*(1+gprim));
     Sln+=log(E)*g*g/((1+gprim)*(1+gprim));
     Sln2+=log(E)*log(E)*g*g/((1+gprim)*(1+gprim));
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(Elowerl/Epivot)+1/(Index-1))*(log(Elowerl/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
   }
   if (GTime[m]>To) break;
}
  return 100.*errK;
 }


double Get_time_unc_mult(double unc,double Index, double Flux, double L, double B, double Ti, double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
  double SigA=0;
  int nn=0;
  double Sln=0,Sln2=0,Epivot=0,Hgg=0,deltaG=0,correc=0,errK=0,ss,g,r,sgprim;
  for (int m=0;m<Gevt;m++){ 
    //if (Gid[m]==0) continue;
   if (GTime[m]>Ti){
     nn++;
     double E=GE[m];
     int Convtype=GConvtype[m]; 
     sgprim=0;
     for (int is=1;is<n_source;is++){ 
               r=Grs[m][is];
               ss=(Index_source[is]-1)*Flux_source[is]*pow(E/100,-Index_source[is]+2.1)/100;
               
	       sgprim+=psf(E,r,Convtype)*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal);}
    
     r=Grs[m][0]; 
     ss=(Index-1)*Flux*pow(E/100,-Index+2.1)/100;
     g=psf(E,r,Convtype)*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal);
     //sgprim=0;
     double gprim=g+sgprim;
     SigA+=g*g/((1+gprim)*(1+gprim));
     Sln+=log(E)*g*g/((1+gprim)*(1+gprim));
     Sln2+=log(E)*log(E)*g*g/((1+gprim)*(1+gprim));
     if (SigA<1e-19) continue;
     Epivot=exp(Sln/SigA);
     Hgg=Sln2+log(Epivot)*log(Epivot)*SigA-2*log(Epivot)*Sln;
     if (Hgg<1e-19) continue;
     deltaG=1/Hgg;
     correc=(log(Elowerl/Epivot)+1/(Index-1))*(log(Elowerl/Epivot)+1/(Index-1))*deltaG;
     errK=sqrt(1/SigA+correc);
     if (errK <unc/100.) return GTime[m];
   }
  }
  return -1;
}
static void fcn_TS_mult(int &npar, double *gin, double &f, double *par, int
iflag){
  double kdum=npar+gin[0]+f+iflag;
  kdum=0; 
  double ll=0;
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double normgal=par[3];
  double normeg=par[4];
  double indgal=par[5];
  double ex=par[6];
  double L=par[7];
  double B=par[8];
  Index_source[0]=Index;
  for (int k=0;k<n_source;k++) Flux_source[k]=pow(10,par[9+k]);
  double ss,g,r,sgprim;
     for (int i=0;i<Gevt;i++){
        //if (Gid[i]==0) continue;
        double Time=GTime[i];
        if ((Tstart<=Time) && (Time<Tstop)){
           double E=GE[i]; 
           sgprim=0; 
           for (int is=1;is<n_source;is++){ 
               r=Grs[i][is];
               ss=(Index_source[is]-1)*Flux_source[is]*pow(E/100,-Index_source[is]+2.1)/100;
               
	       sgprim+=psf(E,r,GConvtype[i])*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal);}
          r=Grs[i][0]; 
          ss=(Index-1)*Flux_source[0]*pow(E/100,-Index+2.1)/100;
          g=psf(E,r,GConvtype[i])*ss/BackgroundFlux(E,L,B,normeg,normgal,indgal)+sgprim;
	  ll+=log(1+g);
	}
	if (Time>Tstop){
	    break;
	    }
     }
	double Ntot=0;

        for (int k=0;k<n_source;k++) Ntot+=Npred(Flux_source[k],Index_source[k],6.*365.25*86400,L,B);
     double TSe=2*(ll-ex*Ntot/Exposure(1000,L,B,0));
     f=-TSe;
  return;
}

double Get_flux_mult_src_minuit(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop,double ex,double* out){
   char ttrc[50];
   *out=0;
   TMinuit minuit(10);
   int ierflg = 0;
   minuit.SetFCN(fcn_TS_mult);
   //minuit.SetFCN(fcn_TS);
   minuit.mnparm(0, "TSTART",Tstart,0.01,1e8,1e9,ierflg);
   minuit.FixParameter(0); 
   minuit.mnparm(1, "TSTOP",Tstop,0.01,1e8,1e9,ierflg);
   minuit.FixParameter(1);
   minuit.mnparm(2, "Index",index,0.01,1,5,ierflg);
   minuit.FixParameter(2);
   minuit.mnparm(3, "normgal",normgal,0.01,0,5,ierflg);
   minuit.FixParameter(3);
   minuit.mnparm(4, "normeg",normeg,0.01,0,5,ierflg);
   minuit.FixParameter(4); 
   minuit.mnparm(5, "indgal",indgal,0.01,-2,2,ierflg);
   minuit.FixParameter(5);
   minuit.mnparm(6, "exposure",ex,0.01,0,1e11,ierflg);
   minuit.FixParameter(6); 
   minuit.mnparm(7, "L",l,0.01,0,400,ierflg);
   minuit.FixParameter(7);
   minuit.mnparm(8, "B",b,0.01,-100,100,ierflg);
   minuit.FixParameter(8);
   minuit.mnparm(9, "Flux", -7,0.01,-9,-4,ierflg);
   for (int k=1;k<n_source;k++) {
     sprintf(ttrc,"Flux_%i",k);
     minuit.mnparm(9+k,ttrc, -7,0.01,-9,-4,ierflg);}
   minuit.SetErrorDef(0.5);
   minuit.SetPrintLevel(-1);
   //minuit.SetMaxIteration(500);
   minuit.Migrad();  
   double flux,errflux,fluxd,errfluxd;
   minuit.GetParameter(9,flux,errflux); 
   for (int k=1;k<n_source;k++) {
        minuit.GetParameter(9+k,fluxd,errfluxd);
        Flux_source[k]=pow(10.,fluxd);
   }
   //*out=fluxd; 
   return pow(10,flux);     
}

double Get_E0(double flux, double index, double l, double b, double time=1, double emin=100, double emax=1e5){
  Double_t pf[6];
  time=0;
  pf[0]=index; //index
  pf[1]=flux;
  pf[2]=l;
  pf[3]=b;
  TF2 *f2=new TF2("f2",S_1,2,6,0,1.5,4);
  f2->SetParameters(&pf[0]);
  TF2 *f_AG=new TF2("f_AG",S_AG,2,6,0,1.5,4);
  f_AG->SetParameters(&pf[0]);
  double elim=log10(emin);
  double elax=log10(emax);
  Double_t s1=Integr(f2,elim,elax);
  Double_t s=Integr(f_AG,elim,elax);
  Double_t E0=pow(10,s/s1);
  return E0;
}


double Get_g(double E,double r, int convtype,double Index, double Flux, double L, double B, double normeg=1,double normgal=1,double indgal=0,double Elowerl=100){
    double Flux100=Flux*pow(Elowerl/100,Index-1);
    double ss=(Index-1)*Flux100/100;
    double g=psf(E,r,convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,normeg,normgal,indgal);
  return g;
 }

void test_bgd_source(){
  initialize_TS();
  for (int i=0;i<10;i++){
    double E=100*pow(10.,i*0.2);
     printf("%lf %e \n", E,ResidualBackground(E)*pow(E/100.,2.1));
  }
}
void test_bgd_clean(){
  initialize_TS_clean();
  for (int i=0;i<10;i++){
    double E=100*pow(10.,i*0.2);
     printf("%lf %e \n", E,ResidualBackground(E)*pow(E/100.,2.1));
  }
}
double Bgd_integrate(double glon, double glat,double normeg,double normgal,double indexgal,double logE0,double logE1,double rad){
  int nn=10;
  double itstep=(1-cos(rad/57.3))/nn;
  double ipstep=360./(57.3*nn);
  double sume=0;  
  for (int itheta=1; itheta< (nn+1); itheta++){
    double theta=57.3*acos(1-itheta*itstep);
    for (int j=0; j< nn; j++){
      double phi=j*ipstep;
      double glon_i=glon+theta*cos(phi);
      double glat_i=glat+theta*sin(phi);
      sume+=Npred_bgd(glon_i,glat_i,normeg,normgal,indexgal,logE0,logE1);
	      }
	      }
return sume*itstep*ipstep;
}

TGraph *va1,*va2;
void initialize_va(double *Time, double *Aefft1,double *Aefft2, const int len){
  double vx[len],vy1[len],vy2[len];
  for (int i=0;i<len;i++){
    vx[i]=Time[i];
    vy1[i]=Aefft1[i];
    vy2[i]=Aefft2[i];
  }
    va1=new TGraph(len,vx,vy1); 
    va2=new TGraph(len,vx,vy2);
    //printf("%lf \n",va1->Eval(2.59803e+08));
    return;
}

double FT(double t,double t0, double tr, double tf){
  double y=2./(exp((-t+t0)/tr)+exp((t-t0)/tf));
  return y;}

Double_t npt(Double_t *x, Double_t *par){
  double t=par[6];
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*(par[1]*FT(t,par[4],par[5],par[8])+par[7]+par[9]*FT(t,par[10],par[11],par[12]))/100; 
  Double_t sum=0;
  Double_t GLON=par[2]; 
  Double_t GLAT=par[3];
for (int i=0;i<2;i++){
   if (i==0)  {double A1=(va1->Eval(t))*Exposure(E,GLON,GLAT,0)/Exposure(1000,GLON,GLAT,0);sum+=ss*pow(E/100,-alpha)*A1*E*log(10);}
  else {double A2=(va2->Eval(t))*Exposure(E,GLON,GLAT,1)/Exposure(1000,GLON,GLAT,1);sum+=ss*pow(E/100,-alpha)*A2*E*log(10);}
 }
  return sum;
 }

Double_t flare(Double_t *x, Double_t *par){
  double t=x[0]; 
  Double_t ss=par[0]*FT(t,par[1],par[2],par[3])+par[4]+par[5]*FT(t,par[6],par[7],par[8]); 
  return ss;
 }


double NpT(double fl100,double gamma, double GLON,double GLAT, double t0, double tr, double t,double f0, double tf,double Flux1=0,double t0b=1,double trb=1,double tfb=1){
    double pf[13];
    pf[0]=gamma;
    pf[1]=fl100; 
    pf[2]=GLON; 
    pf[3]=GLAT;
    pf[4]=t0; 
    pf[5]=tr; 
    pf[6]=t;
    pf[7]=f0;
    pf[8]=tf;
    pf[9]=Flux1;
    pf[10]=t0b;
    pf[11]=trb;
    pf[12]=tfb;
    TF1 *fnt=new TF1("fnt",npt,2.,6,13);
    fnt->SetParameters(&pf[0]);
    return fnt->Integral(2,6);   
}

 Double_t Npte(Double_t *x, Double_t *par){
  double t=x[1];
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*(par[1]*FT(t,par[4],par[5],par[7])+par[6]+par[8]*FT(t,par[9],par[10],par[11]))/100; 
  Double_t sum=0;
  Double_t GLON=par[2]; 
  Double_t GLAT=par[3];
  for (int i=0;i<2;i++){
    if (i==0)  {double A1=(va1->Eval(t))*Exposure(E,GLON,GLAT,0)/Exposure(1000,GLON,GLAT,0);sum+=ss*pow(E/100,-alpha)*A1*E*log(10);
      //printf("A1:%e t %lf\n",va1->Eval(t),t);
}
  else {double A2=(va2->Eval(t))*Exposure(E,GLON,GLAT,1)/Exposure(1000,GLON,GLAT,1);sum+=ss*pow(E/100,-alpha)*A2*E*log(10);}
 }
  
  return sum;
 }
   
double NpredT(double fl100,double gamma, double GLON,double GLAT, double tmin,double tmax, double t0, double tr, double f0, double tf,double Flux1=0,double t0b=1,double trb=1,double tfb=1){
    double pf[12];
    pf[0]=gamma;
    pf[1]=fl100; 
    pf[2]=GLON; 
    pf[3]=GLAT;
    pf[4]=t0; 
    pf[5]=tr;
    pf[6]=f0;
    pf[7]=tf;
    pf[8]=Flux1;
    pf[9]=t0b;
    pf[10]=trb;
    pf[11]=tfb;
    TF2 *fnt=new TF2("fnt",Npte,2.,6,tmin,tmax,12);
    fnt->SetParameters(&pf[0]);
    int nstepx=50;
    int nstepy=500;
    double xmin=2,xmax=6,ymin=tmin,ymax=tmax;
    double xstep=(xmax-xmin)/nstepx;
    double ystep=(ymax-ymin)/nstepy;
    double sum=0;
    for (int kx=0;kx<nstepx-1;kx++){ 
 	double x=xmin+(kx+0.5)*xstep;
      for (int ky=0;ky<nstepy-1;ky++){ 				//loop over log
        double y=ymin+(ky+0.5)*ystep;
        sum+= fnt->Eval(x,y)*xstep*ystep;
      }
    }
    return sum;  
}

static void fcn_TS_var(int &npar, double *gin, double &f, double *par, int
iflag){
  double kdum=npar+gin[0]+f+iflag;
  kdum=0; 
  double ll=0;
  double Flux=pow(10,par[8]);
  double Flux0=pow(10,par[11]);
  double Flux1=pow(10,par[13]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double L=par[6];
  double B=par[7];
  double t0=par[9];
  double tr=par[10]; 
  double tf=par[12];
  double t0b=par[14];
  double trb=par[15]; 
  double tfb=par[16];
  double ss=(Index-1)*Flux/100;  
  double ss0=(Index-1)*Flux0/100;  
  double ss1=(Index-1)*Flux1/100;  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	  double g=psf(E,r,Convtype)*(ss*FT(Time,t0,tr,tf)+ss0+ss1*FT(Time,t0b,trb,tfb))*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
     }
	double yy=NpredT(Flux,Index,L,B,Tstart,Tstop,t0,tr,Flux0,tf,Flux1,t0b,trb,tfb);
     double TSe=2*(ll-yy); 
  f=-TSe;
  //printf("%lf %lf \n",f,yy);
 }

double Fit_var(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop, double ff0, double* out,double fbgd1=1,double t0bgd=1, double trbgd=1,double tfbgd=1){
   TMinuit minuit(17);
   int ierflg = 0;
   minuit.SetFCN(fcn_TS_var);
   minuit.mnparm(0, "TSTART",Tstart,0.01,0,1e9,ierflg);
   minuit.FixParameter(0); 
   minuit.mnparm(1, "TSTOP",Tstop,0.01,0,1e9,ierflg);
   minuit.FixParameter(1);
   minuit.mnparm(2, "Index",index,0.01,1,5,ierflg);
   minuit.FixParameter(2);
   minuit.mnparm(3, "normgal",normgal,0.01,0,5,ierflg);
   minuit.FixParameter(3);
   minuit.mnparm(4, "normeg",normeg,0.01,0,5,ierflg);
   minuit.FixParameter(4);
   minuit.mnparm(5, "indgal",indgal,0.01,-2,2,ierflg);
   minuit.FixParameter(5);
   minuit.mnparm(6, "L",l,0.01,0,360,ierflg);
   minuit.FixParameter(6);
   minuit.mnparm(7, "B",b,0.01,-90,90,ierflg);
   minuit.FixParameter(7);
   //minuit.mnparm(8, "Flux", log10(5*ff0),0.01,-9,-4,ierflg); 
   minuit.mnparm(8, "Flux", -5.0,0.01,-9,-4,ierflg);
   //minuit.FixParameter(8);
   minuit.mnparm(9, "t0",Tstart+0.55*(-Tstart+Tstop),0.01,Tstart, Tstop,ierflg); 
   //minuit.FixParameter(9);
   minuit.mnparm(10, "tr", 8000,0.01,10,50000,ierflg);
   //minuit.FixParameter(10);
   minuit.mnparm(11, "Flux0", log10(ff0),0.01,-9,-4,ierflg);
   minuit.FixParameter(11);
   minuit.mnparm(12, "tf", 8000,1,100,1e5,ierflg);
   //minuit.FixParameter(12);
   minuit.mnparm(13, "Fbgd1", fbgd1,1e-6,0,10,ierflg);
   minuit.FixParameter(13); 
   minuit.mnparm(14, "t0bgd", t0bgd,0.01,Tstart, Tstop,ierflg);
   minuit.FixParameter(14);
   minuit.mnparm(15, "trbgd", trbgd,0.01,10,50000,ierflg);
   minuit.FixParameter(15); 
   minuit.mnparm(16, "tfbgd", tfbgd,0.01,10,1e6,ierflg);
   minuit.FixParameter(16);
   minuit.SetErrorDef(0.5);
   //minuit.SetPrintLevel(-1);
   //minuit.SetMaxIteration(1);
   minuit.Migrad();  
   double flux,errflux, t0,errt0, tr, errtr,tf,errtf;
   minuit.GetParameter(8,flux,errflux);  
   minuit.GetParameter(9,t0,errt0); 
   minuit.GetParameter(10,tr,errtr);   
   minuit.GetParameter(12,tf,errtf); 
   *out=t0;
   double *out1=out+1; 
   *out1=tr;
   double *out2=out1+1; 
   *out2=tf;
   return pow(10,flux);    
}
double Get_g_F(double E,double r, int convtype,double Index, double Flux, double L, double B, double t0, double tr, double t, double fbgd,double tf){
  double ss=(Index-1)*(Flux*FT(t,t0,tr,tf)+fbgd)/100;
    double g=psf(E,r,convtype)*ss*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,1,1,0);
  return g;
 }
Double_t npt_new(Double_t *x, Double_t *par){
  double t=par[6];
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*(par[1]*FT(t,par[4],par[5],par[8])+par[7])/100; 
  Double_t sum=0;
  Double_t GLON=par[2]; 
  Double_t GLAT=par[3];
for (int i=0;i<2;i++){
   if (i==0)  {double A1=(va1->Eval(t))*Exposure(E,GLON,GLAT,0)/Exposure(1000,GLON,GLAT,0);sum+=ss*pow(E/100,-alpha)*A1*E*log(10);}
  else {double A2=(va2->Eval(t))*Exposure(E,GLON,GLAT,1)/Exposure(1000,GLON,GLAT,1);sum+=ss*pow(E/100,-alpha)*A2*E*log(10);}
 }
  return sum;
 }

double NpT_new(double fl100,double gamma, double GLON,double GLAT, double t0, double tr, double t,double f0, double tf){
    double pf[13];
    pf[0]=gamma;
    pf[1]=fl100; 
    pf[2]=GLON; 
    pf[3]=GLAT;
    pf[4]=t0; 
    pf[5]=tr; 
    pf[6]=t;
    pf[7]=f0;
    pf[8]=tf;
    TF1 *fnt=new TF1("fnt",npt_new,2.,6,9);
    fnt->SetParameters(&pf[0]);
    return fnt->Integral(2,6);   
}

 Double_t Npte_new(Double_t *x, Double_t *par){
  double t=x[1];
  Double_t E=pow(10,x[0]); 
  Double_t alpha=par[0]; //photon index 
  Double_t ss=(alpha-1)*(par[1]*FT(t,par[4],par[5],par[7])+par[6])/100; 
  Double_t sum=0;
  Double_t GLON=par[2]; 
  Double_t GLAT=par[3];
  for (int i=0;i<2;i++){
    if (i==0)  {double A1=(va1->Eval(t))*Exposure(E,GLON,GLAT,0)/Exposure(1000,GLON,GLAT,0);sum+=ss*pow(E/100,-alpha)*A1*E*log(10);
      //printf("A1:%e t %lf\n",va1->Eval(t),t);
}
  else {double A2=(va2->Eval(t))*Exposure(E,GLON,GLAT,1)/Exposure(1000,GLON,GLAT,1);sum+=ss*pow(E/100,-alpha)*A2*E*log(10);}
 }
  
  return sum;
 }
   
double NpredT_new(double fl100,double gamma, double GLON,double GLAT, double tmin,double tmax, double t0, double tr, double f0, double tf){
    double pf[8];
    pf[0]=gamma;
    pf[1]=fl100; 
    pf[2]=GLON; 
    pf[3]=GLAT;
    pf[4]=t0; 
    pf[5]=tr;
    pf[6]=f0;
    pf[7]=tf;
    TF2 *fnt=new TF2("fnt",Npte_new,2.,6,tmin,tmax,8);
    fnt->SetParameters(&pf[0]);
    int nstepx=50;
    int nstepy=500;
    double xmin=2,xmax=6,ymin=tmin,ymax=tmax;
    double xstep=(xmax-xmin)/nstepx;
    double ystep=(ymax-ymin)/nstepy;
    double sum=0;
    for (int kx=0;kx<nstepx-1;kx++){ 
 	double x=xmin+(kx+0.5)*xstep;
      for (int ky=0;ky<nstepy-1;ky++){ 				//loop over log
        double y=ymin+(ky+0.5)*ystep;
        sum+= fnt->Eval(x,y)*xstep*ystep;
      }
    }
    //printf("sum: %lf\n",sum);
    return sum;    
}

static void fcn_TS_var_new(int &npar, double *gin, double &f, double *par, int
iflag){
  double kdum=npar+gin[0]+f+iflag;
  kdum=0;  
  double ll=0;
  double Flux=pow(10,par[8]);
  double Flux0=pow(10,par[11]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double L=par[6];
  double B=par[7];
  double t0=par[9];
  double tr=par[10]; 
  double tf=par[12];
  double ss=(Index-1)*Flux/100;  
  double ss0=(Index-1)*Flux0/100;  
  int nnn=0;
  //printf("%d\n",Gevt);    
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];

	if ((Tstart<=Time) && (Time<Tstop)){
          nnn++; 
	  double g=psf(E,r,Convtype)*(ss*FT(Time,t0,tr,tf)+ss0)*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
     }
	double yy=NpredT_new(Flux,Index,L,B,Tstart,Tstop,t0,tr,Flux0,tf);
     double TSe=2*(ll-yy); 
  f=-TSe;
  //printf("f: %e %d %e %e\n",f,nnn,t0,yy);
  //exit(0);
 }

double Fit_var_new(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop, double ff0i, double t0i, double tri, double tfi, double fbgd, double* out){
   TMinuit minuit(13);
   int ierflg = 0;
   minuit.SetFCN(fcn_TS_var_new);
   minuit.mnparm(0, "TSTART",Tstart,0.01,0,1e9,ierflg);
   minuit.FixParameter(0); 
   minuit.mnparm(1, "TSTOP",Tstop,0.01,0,1e9,ierflg);
   minuit.FixParameter(1);
   minuit.mnparm(2, "Index",index,0.01,1,5,ierflg);
   minuit.FixParameter(2);
   minuit.mnparm(3, "normgal",normgal,0.01,0,5,ierflg);
   minuit.FixParameter(3);
   minuit.mnparm(4, "normeg",normeg,0.01,0,5,ierflg);
   minuit.FixParameter(4);
   minuit.mnparm(5, "indgal",indgal,0.01,-2,2,ierflg);
   minuit.FixParameter(5);
   minuit.mnparm(6, "L",l,0.01,0,360,ierflg);
   minuit.FixParameter(6);
   minuit.mnparm(7, "B",b,0.01,-90,90,ierflg);
   minuit.FixParameter(7);
   minuit.mnparm(8, "Flux", log10(ff0i),0.01,-9,-4,ierflg); 
   //minuit.mnparm(8, "Flux", -5.0,0.01,-9,-4,ierflg);
   //minuit.FixParameter(8);
   minuit.mnparm(9, "t0",t0i,0.01,Tstart, Tstop,ierflg); 
   //minuit.FixParameter(9);
   minuit.mnparm(10, "tr", tri,0.01,10,2e5,ierflg);
   //minuit.FixParameter(10);
   minuit.mnparm(11, "Flux0", log10(fbgd),0.01,-9,-4,ierflg);
   //minuit.FixParameter(11);
   minuit.mnparm(12, "tf", tfi,1,100,2e5,ierflg);
   //minuit.FixParameter(12);
   minuit.SetErrorDef(0.5);
   //minuit.SetPrintLevel(-1);
   //minuit.SetMaxIteration(500);
   minuit.Migrad();
   double flux,errflux, t0,errt0, tr, errtr,tf,errtf;
   minuit.GetParameter(8,flux,errflux);  
   minuit.GetParameter(9,t0,errt0); 
   minuit.GetParameter(10,tr,errtr);   
   minuit.GetParameter(12,tf,errtf); 
   *out=t0;
   double *out1=out+1; 
   *out1=tr;
   double *out2=out1+1; 
   *out2=tf;
   return pow(10,flux);    
}
static void fcn_TS_var_2(int &npar, double *gin, double &f, double *par, int
iflag){
  double kdum=npar+gin[0]+f+iflag;
  kdum=0; 
  double ll=0;
  double Flux=pow(10,par[8]);
  double Flux0=pow(10,par[12]);
  double Flux1=pow(10,par[13]);
  double Tstart=par[0];
  double Tstop=par[1];
  double Index=par[2];
  double NormGal=par[3];
  double NormEG=par[4];
  double IndexGal=par[5];
  double L=par[6];
  double B=par[7];
  double t0=par[9];
  double tr=par[10]; 
  double tf=par[11];
  double t0b=par[14];
  double trb=par[15]; 
  double tfb=par[16];
  double ss=(Index-1)*Flux/100;  
  double ss0=(Index-1)*Flux0/100;  
  double ss1=(Index-1)*Flux1/100;  
     for (int i=0;i<Gevt;i++){
        double E=GE[i];
        double Time=GTime[i];
        int Convtype=GConvtype[i];
	double r=Gr[i];
	    
	if ((Tstart<=Time) && (Time<Tstop)){
	  double g=psf(E,r,Convtype)*(ss*FT(Time,t0,tr,tf)+ss0+ss1*FT(Time,t0b,trb,tfb))*pow(E/100,-Index+2.1)/BackgroundFlux(E,L,B,NormEG,NormGal,IndexGal);
	    ll+=log(1+g);
	    }
	if (Time>Tstop){
	    break;
	    }
     }
	double yy=NpredT(Flux,Index,L,B,Tstart,Tstop,t0,tr,Flux0,tf,Flux1,t0b,trb,tfb);
     double TSe=2*(ll-yy); 
  f=-TSe;
  printf("%lf %lf \n",f,yy);
 }

double Fit_var_2(double index,double l,double b,double normeg,double normgal,double indgal,double Tstart,double Tstop, double ffi, double t00,double tr0,double tf0, double fbgd, double ff1, double t01,double tr1,double tf1, double* out){
   TMinuit minuit(17);
   int ierflg = 0;
   minuit.SetFCN(fcn_TS_var_2);
   minuit.mnparm(0, "TSTART",Tstart,100,2e8,5e8,ierflg);
   minuit.FixParameter(0); 
   minuit.mnparm(1, "TSTOP",Tstop,100,2e8,5e8,ierflg);
   minuit.FixParameter(1);
   minuit.mnparm(2, "Index",index,0.01,1,5,ierflg);
   minuit.FixParameter(2);
   minuit.mnparm(3, "normgal",normgal,0.01,0,5,ierflg);
   minuit.FixParameter(3);
   minuit.mnparm(4, "normeg",normeg,0.01,0,5,ierflg);
   minuit.FixParameter(4);
   minuit.mnparm(5, "indgal",indgal,0.01,-2,2,ierflg);
   minuit.FixParameter(5);
   minuit.mnparm(6, "L",l,0.01,0,360,ierflg);
   minuit.FixParameter(6);
   minuit.mnparm(7, "B",b,0.01,-90,90,ierflg);
   minuit.FixParameter(7);
   minuit.mnparm(8, "Flux", log10(ffi),0.01,-9,-4,ierflg);
   //minuit.FixParameter(8); 
   minuit.mnparm(9, "t0",t00,0.01,Tstart, Tstop,ierflg);
   //minuit.FixParameter(9); 
   minuit.mnparm(10, "tr", tr0,0.01,10,50000,ierflg); 
   //minuit.FixParameter(10); 
   minuit.mnparm(11, "tf", tf0,1,100,1e5,ierflg);
   //minuit.FixParameter(11);
   minuit.mnparm(12, "Flux0", log10(fbgd),0.01,-9,-4,ierflg);
   minuit.FixParameter(12);
   minuit.mnparm(13, "Fbgd1", log10(ff1),0.01,-9,-4,ierflg);
   //minuit.FixParameter(13);
   minuit.mnparm(14, "t0bgd", t01,0.01,Tstart, Tstart+2*(Tstop-Tstart),ierflg);
   //minuit.FixParameter(14);
   minuit.mnparm(15, "trbgd", tr1,0.01,10,50000,ierflg);
   //minuit.FixParameter(15);
   minuit.mnparm(16, "tfbgd", tf1,0.01,10,1e6,ierflg);
   //minuit.FixParameter(16);
   minuit.SetErrorDef(0.5);
   //minuit.SetPrintLevel(-1);
   //minuit.SetMaxIteration(1);
   minuit.Migrad();  
   double flux,errflux, t0,errt0, tr, errtr,tf,errtf, flux1,errflux1, t01f,errt01f, tr1f, errtr1f,tf1f,errtf1f,flux2,errflux2;
   minuit.GetParameter(8,flux,errflux);  
   minuit.GetParameter(9,t0,errt0); 
   minuit.GetParameter(10,tr,errtr);   
   minuit.GetParameter(11,tf,errtf);
   minuit.GetParameter(12,flux2,errflux2);  
   minuit.GetParameter(13,flux1,errflux1);  
   minuit.GetParameter(14,t01f,errt01f); 
   minuit.GetParameter(15,tr1f,errtr1f);   
   minuit.GetParameter(16,tf1f,errtf1f);
   *out=t0;
   double *out1=out+1; 
   *out1=tr;
   double *out2=out+2; 
   *out2=tf;
   double *out3=out+3; 
   *out3=pow(10,flux1);
   double *out4=out+4; 
   *out4=t01f;
   double *out5=out+5; 
   *out5=tr1f;
   double *out6=out+6; 
   *out6=tf1f;
   double *out7=out+7; 
   *out7=pow(10,flux2);
   return pow(10,flux);    
}
