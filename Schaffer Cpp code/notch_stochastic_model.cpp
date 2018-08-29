/***************************************************************************************************
  Stochastic simulation of the Notch chemical reaction system using the Gillespie algorithm

  Version 0.01: Implement Deterministic equations stochastically, print to screen # of each species

***************************************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/****************
  Include files
****************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MersenneTwister.h"
#include <map>

using namespace std;
using std::ofstream;
using std::endl;
using std::cerr;
using std::ios;
using std::cout;

// Species
//Units here are molecules per cell
int Hcm;				//1.6244x10^(-11)
int Hcp;				//5.55879*10^(-10)
int Hnp;				//1.76469*10^(-9)
double H2np;			//
int Rcm;				//2.42085*10^(-11)
int Rcp;				//7.57182*10^(-11)
int Rnp;				//3.27784*10^(-9)
int Nm;				//4.06075*10^(-9)	//Notch mRNA
int Np;				//4.7729*10^(-8)		//Notch protein in the cytoplasm
int Ncp;				//NICD protein in the cytoplasm
int Nnp;				//NICD protein in the nucleus
int delta;			//Delta protein
int NP_free;			//Notch promoter (diploid so start with 2 free promoter regions)
int NP_Rr;
int NP_RNa;
int NP_RNa_Rr;
int NP_Rr2;
int NP_RNa2;
int NP_Hr;
int NP_Rr_Hr;
int NP_RNa_Hr;
int NP_RNa_Rr_Hr;
int NP_Rr2_Hr;
int NP_RNa2_Hr;
int HP_free;			//Hes1 promoter (diploid so start with 2 free promoter regions)
int HP_Hr;
int HP_Hr2;
int HP_Hr3;
int HP_Rr;
int HP_Rr_Hr;
int HP_Rr_Hr2;
int HP_Rr_Hr3;
int HP_RNa;
int HP_RNa_Hr;
int HP_RNa_Hr2;
int HP_RNa_Hr3;
int HP_RNa_Rr;
int HP_RNa_Rr_Hr;
int HP_RNa_Rr_Hr2;
int HP_RNa_Rr_Hr3;
int HP_Rr2;
int HP_Rr2_Hr;
int HP_Rr2_Hr2;
int HP_Rr2_Hr3;
int HP_RNa2;
int HP_RNa2_Hr;
int HP_RNa2_Hr2;
int HP_RNa2_Hr3;
int RP_free;			//RPJK promoter (diploid so start with 2 free promoter regions)
int RP_Hr=0;
int RP_Hr2=0;
int RP_Hr3=0;
int RP_Rr=0;
int RP_Rr_Hr=0;
int RP_Rr_Hr2=0;
int RP_Rr_Hr3=0;
int RP_RNa=0;
int RP_RNa_Hr=0;
int RP_RNa_Hr2=0;
int RP_RNa_Hr3=0;
int RP_Rr2=0;
int RP_Rr2_Hr=0;
int RP_Rr2_Hr2=0;
int RP_Rr2_Hr3=0;
int RP_RNa_Rr=0;
int RP_RNa_Rr_Hr=0;
int RP_RNa_Rr_Hr2=0;
int RP_RNa_Rr_Hr3=0;
int RP_Rr3=0;
int RP_Rr3_Hr=0;
int RP_Rr3_Hr2=0;
int RP_Rr3_Hr3=0;
int RP_RNa_Rr2=0;
int RP_RNa_Rr2_Hr=0;
int RP_RNa_Rr2_Hr2=0;
int RP_RNa_Rr2_Hr3=0;
int RP_RNa2=0;
int RP_RNa2_Hr=0;
int RP_RNa2_Hr2=0;
int RP_RNa2_Hr3=0;
int RP_RNa2_Rr=0;
int RP_RNa2_Rr_Hr=0;
int RP_RNa2_Rr_Hr2=0;
int RP_RNa2_Rr_Hr3=0;
int RP_RNa3=0;
int RP_RNa3_Hr=0;
int RP_RNa3_Hr2=0;
int RP_RNa3_Hr3=0;

//skasdhjfk
map<double, int> delayedRxnList;

bool delayedRxnReady(double timeOfRxn) {

	//merely check if a delayed rxn occurs before or at timeOfRxn
	//if yes return true
	map<double,int>::iterator mi=delayedRxnList.begin();

	if(delayedRxnList.size()==0) {
		return false;
	}

	//for(map<double, int>::iterator iter = delayedRxnList.begin(); iter != delayedRxnList.end(); iter++ ) {
	//	double tim = (*iter).first;
	//	int j = (*iter).second;
	//}

	double firstDelayedRxn = (*mi).first;

	if( firstDelayedRxn <= timeOfRxn) {
		return true;
	}

	return false;

}

double fireRxn() {

	//find RxnNumber
	//update X[] with RxnNumber firing
	//update t to time of firing
	map<double,int>::iterator mi=delayedRxnList.begin();
	double t = (*mi).first;
	int rxn = (*mi).second;

	//cout<<rxn<<endl;
	switch(rxn) {
		//1. 0 --> Nm		0.RfNm
		case 0:
			Nm++;
		break;

		//2. 0 --> Rcm		1.RfRm
		case 1:
			Rcm++;
		break;

		//3. 0 --> Hcm		2.RfHm
		case 2:
			Hcm++;
		break;

		//4. Nm --> Np		3.KtrN
		case 3:
			//Nm--; 	Don't destroy mRNA after transcription
			Np++;
		break;

		//5. Rcm --> Rcp		4.KtrRc
		case 4:
			//Rcm--; 	Don't destroy mRNA after transcription
			Rcp++;
		break;

		//6. Hcm --> Hcp		5.KtrHc
		case 5:
			//Hcm--; 	Don't destroy mRNA after transcription
			Hcp++;
			break;
	}

	delayedRxnList.erase(t);
	return t;

}

void addRxn(int rxn, double timeOfRxn, double delay) {

	//Add reaction rxn to delayed reaction list
	double timeOfRxnFinish = timeOfRxn+delay;
	delayedRxnList[timeOfRxnFinish]=rxn;
}

int main(int argc, char *argv[]) {

	int runTime=3000;
	int numRuns=1;
	//int numSpecies=87;
	int numRxns=299;


	bool delay=true;	//set true here if you want to implement delays

	//Delay times for the 6 delayed reactions
	//units here are minutes delayed
	double TpNc=21;
	double TpRc=4.3;
	double TpHc=2.35;
	double TmNc=70;
	double TmRc=20;
	double TmHc=10;

	//Reaction rates
	//Units here are per minute
	float KdNp=0.017;
	float KdHcp=0.0315;
	float KdHnp=0.0315;
	float KdHcm=0.029;
	float KdH2np=0;		//disable Hes1 dimer degradation
	float KdNcp=0.00385;
	float KdNnp=0.00385;
	float KdNm=0.00035;
	float KdRcp=0.00231;
	float KdRnp=0.00231;
	float KdRcm=0.0075;

	//units here are per minute
	float KtrHc=4.5;
	float KtrN=2;
	float KtrRc=3.2;

	//units here are per minute
	float KniRcp=0.1;
	float KniHcp=0.1;
	float KniNcp=0.1;


	float KfNcp=0.000276;	//7.6*10^8 per molar per minute		//forward rate of NICD binding to Delta
	float KaHp=0.02535;
	float KrHp=0;			// disable disassociation of Hes1 dimer
	float Kr=0.0410;    		//Keqr=8.19e-3
	//float Kr=0.001628;
	float Kr_r=5;
	float Kn=0.0254;    		//Keqn=5.07e-3
	//float Kn=0.00072;
	float Kn_r=5;
	float Ka=0.0127;   		//Keqa=2.54e-3
	//float Ka=0.00036;
	float Ka_r=5;

	//These generation rates are based on molecules generated per cell per minute
	float Vmaxh=197;     //9.98;
	float Vmaxr= 79;     //4.99
	float Vmaxn=21.6;    //1.43
	float Vbh=4.5;		//2.53*10^(-12)M/min
	float Vbr=1.7;		//1.27*10^(-12)M/min
	float Vbn=0.5;		//3.6*10^(-13)M/min
	float rNbox=0.3;
	float rRbox=0.2;
	float tc=0.5;

	int mu=0;
	int i,j,k,tp;


	float r1, r2;

	double t, tau, a0, sum;
	double c[28]={0,0,0,KtrN,KtrRc,KtrHc,KaHp,KrHp,KfNcp,KdNm,KdRcm,KdHcm,KdNnp,KdRnp,KdNp,KdNcp,KdRcp,KdHcp,KdH2np,KniNcp,KniRcp,KniHcp,Kr,Kr_r,Kn,Kn_r,Ka,Ka_r};

	//put species you want to follow here

	for(i=0;i<numRuns;i++)	{

    	//make the A array a size of numRxn+1 so we can start at reaction 1 instead of 0
    	double a[numRxns+1];
    	for(j=0;j<numRxns+1;j++){
    		a[j]=0;
    	}

		long int seed = time(0);

		//random seed
		srand(seed);
		//same seed
		//srand(1);

		 Hcm=1;			//3.64x10^(-12)
         Hcp=35;			//1.25*10^(-10)
         Hnp=109;		//2.77*10^(-9)
         H2np=101;		//2.5662*10^(-9)
         Rcm=1;			//1.2*10^(-12)
         Rcp=11;			//3.74*10^(-11)
         Rnp=447;		//1.13*10^(-8)
         Nm=70;			//2.51*10^(-10)	//Notch mRNA
         Np=8130;		//2.95*10^(-8)	//Notch protein in the cytoplasm
         Ncp=0;			//NICD protein in the cytoplasm
         Nnp=0;			//NICD protein in the nucleus
         delta=0;		//Delta protein
         NP_free=2;		//Notch promoter (diploid so start with 2 free promoter regions)
         NP_Rr=0;
         NP_RNa=0;
         NP_RNa_Rr=0;
         NP_Rr2=0;
         NP_RNa2=0;
         NP_Hr=0;
         NP_Rr_Hr=0;
         NP_RNa_Hr=0;
         NP_RNa_Rr_Hr=0;
         NP_Rr2_Hr=0;
         NP_RNa2_Hr=0;
         HP_free=2;		//Hes1 promoter (diploid so start with 2 free promoter regions)
         HP_Hr=0;
         HP_Hr2=0;
         HP_Hr3=0;
         HP_Rr=0;
         HP_Rr_Hr=0;
         HP_Rr_Hr2=0;
         HP_Rr_Hr3=0;
         HP_RNa=0;
         HP_RNa_Hr=0;
         HP_RNa_Hr2=0;
         HP_RNa_Hr3=0;
         HP_RNa_Rr=0;
         HP_RNa_Rr_Hr=0;
         HP_RNa_Rr_Hr2=0;
         HP_RNa_Rr_Hr3=0;
         HP_Rr2=0;
         HP_Rr2_Hr=0;
         HP_Rr2_Hr2=0;
         HP_Rr2_Hr3=0;
         HP_RNa2=0;
         HP_RNa2_Hr=0;
         HP_RNa2_Hr2=0;
         HP_RNa2_Hr3=0;
         RP_free=2;			//RPJK promoter (diploid so start with 2 free promoter regions)
         RP_Hr=0;
         RP_Hr2=0;
         RP_Hr3=0;
         RP_Rr=0;
         RP_Rr_Hr=0;
         RP_Rr_Hr2=0;
         RP_Rr_Hr3=0;
         RP_RNa=0;
         RP_RNa_Hr=0;
         RP_RNa_Hr2=0;
         RP_RNa_Hr3=0;
         RP_Rr2=0;
         RP_Rr2_Hr=0;
         RP_Rr2_Hr2=0;
         RP_Rr2_Hr3=0;
         RP_RNa_Rr=0;
         RP_RNa_Rr_Hr=0;
         RP_RNa_Rr_Hr2=0;
         RP_RNa_Rr_Hr3=0;
         RP_Rr3=0;
         RP_Rr3_Hr=0;
         RP_Rr3_Hr2=0;
         RP_Rr3_Hr3=0;
         RP_RNa_Rr2=0;
         RP_RNa_Rr2_Hr=0;
         RP_RNa_Rr2_Hr2=0;
         RP_RNa_Rr2_Hr3=0;
         RP_RNa2=0;
         RP_RNa2_Hr=0;
         RP_RNa2_Hr2=0;
         RP_RNa2_Hr3=0;
         RP_RNa2_Rr=0;
         RP_RNa2_Rr_Hr=0;
         RP_RNa2_Rr_Hr2=0;
         RP_RNa2_Rr_Hr3=0;
         RP_RNa3=0;
         RP_RNa3_Hr=0;
         RP_RNa3_Hr2=0;
         RP_RNa3_Hr3=0;
		t=0;
		tp=0;
		cout<<"Run #"<<i<<endl;
		cout<<"t\t"<<"Rcm\t"<<"Hcm\t"<<"Nm\t"<<"Rcp\t"<<"Hcp\t"<<"Np\t"<<"Hnp\t"<<endl;


		//Runs the simulation
		//can stop when we have a certain number of a species or when we reach a simulation time
		while(t<runTime) {

			// finding c values that need to be calculated
			// basically only the formation of the mRNA from the various products
			//Notch promoter "RfNm"
			c[0]=NP_free*Vbn +
			     rRbox*NP_Rr*(Vbn) +
			     NP_RNa*(tc*Vmaxn+Vbn) +
			     rRbox*NP_RNa_Rr*(tc*Vmaxn+Vbn) +
			     rRbox*rRbox*NP_Rr2*(Vbn) +
			     NP_RNa2*(Vmaxn+Vbn) +
			     rNbox*NP_Hr*(Vbn) +
			     rRbox*rNbox*NP_Rr_Hr*(Vbn) +
			     rNbox*NP_RNa_Hr*(tc*Vmaxn+Vbn) +
			     rNbox*rRbox*NP_RNa_Rr_Hr*(tc*Vmaxn+Vbn) +
			     rRbox*rRbox*rNbox*NP_Rr2_Hr*(Vbn) +
			     rNbox*NP_RNa2_Hr*(Vmaxn+Vbn);

			//RBPJ-k promoter "RfRm"
			c[1]=RP_free*(Vbr) +
				rNbox*RP_Hr*(Vbr) +
				rNbox*rNbox*RP_Hr2*(Vbr) +
				rNbox*rNbox*rNbox*RP_Hr3*(Vbr) +
				rRbox*RP_Rr*(Vbr) +
				rRbox*rNbox*RP_Rr_Hr*(Vbr) +
				rRbox*rNbox*rNbox*RP_Rr_Hr2*(Vbr) +
				rRbox*rNbox*rNbox*rNbox*RP_Rr_Hr3*(Vbr) +
				RP_RNa*(tc*tc*Vmaxr+Vbr) +
				rNbox*RP_RNa_Hr*(tc*tc*Vmaxr+Vbr) +
				rNbox*rNbox*RP_RNa_Hr2*(tc*tc*Vmaxr+Vbr) +
				rNbox*rNbox*rNbox*RP_RNa_Hr3*(tc*tc*Vmaxr+Vbr) +
				rRbox*rRbox*RP_Rr2*(Vbr) +
				rRbox*rRbox*rNbox*RP_Rr2_Hr*(Vbr) +
				rRbox*rRbox*rNbox*rNbox*RP_Rr2_Hr2*(Vbr) +
				rRbox*rRbox*rNbox*rNbox*rNbox*RP_Rr2_Hr3*(Vbr) +
				rRbox*RP_RNa_Rr*(tc*tc*Vmaxr+Vbr) +
				rRbox*rNbox*RP_RNa_Rr_Hr*(tc*tc*Vmaxr+Vbr) +
				rRbox*rNbox*rNbox*RP_RNa_Rr_Hr2*(tc*tc*Vmaxr+Vbr) +
				rRbox*rNbox*rNbox*rNbox*RP_RNa_Rr_Hr3*(tc*tc*Vmaxr+Vbr) +
				rRbox*rRbox*rRbox*RP_Rr3*(Vbr) +
				rRbox*rRbox*rRbox*rNbox*RP_Rr3_Hr*(Vbr) +
				rRbox*rRbox*rRbox*rNbox*rNbox*RP_Rr3_Hr2*(Vbr) +
				rRbox*rRbox*rRbox*rNbox*rNbox*rNbox*RP_Rr3_Hr3*(Vbr) +
				rRbox*rRbox*RP_RNa_Rr2*(tc*tc*Vmaxr+Vbr) +
				rRbox*rRbox*rNbox*RP_RNa_Rr2_Hr*(tc*tc*Vmaxr+Vbr) +
				rRbox*rRbox*rNbox*rNbox*RP_RNa_Rr2_Hr2*(tc*tc*Vmaxr+Vbr) +
				rRbox*rRbox*rNbox*rNbox*rNbox*RP_RNa_Rr2_Hr3*(tc*tc*Vmaxr+Vbr) +
				RP_RNa2*(tc*Vmaxr+Vbr) +
				rNbox*RP_RNa2_Hr*(tc*Vmaxr+Vbr) +
				rNbox*rNbox*RP_RNa2_Hr2*(tc*Vmaxr+Vbr) +
				rNbox*rNbox*rNbox*RP_RNa2_Hr3*(tc*Vmaxr+Vbr) +
				rRbox*RP_RNa2_Rr*(tc*Vmaxr+Vbr) +
				rRbox*rNbox*RP_RNa2_Rr_Hr*(tc*Vmaxr+Vbr) +
				rRbox*rNbox*rNbox*RP_RNa2_Rr_Hr2*(tc*Vmaxr+Vbr) +
				rRbox*rNbox*rNbox*rNbox*RP_RNa2_Rr_Hr3*(tc*Vmaxr+Vbr) +
				RP_RNa3*(Vmaxr+Vbr) +
				rNbox*RP_RNa3_Hr*(Vmaxr+Vbr) +
				rNbox*rNbox*RP_RNa3_Hr2*(Vmaxr+Vbr) +
				rNbox*rNbox*rNbox*RP_RNa3_Hr3*(Vmaxr+Vbr);


			//Hes1 promoter "RfHm"
			c[2]=HP_free*Vbh +
			     rNbox*HP_Hr*(Vbh) +
			     rNbox*rNbox*HP_Hr2*(Vbh) +
				 rNbox*rNbox*rNbox*HP_Hr3*(Vbh) +
				 rRbox*HP_Rr*(Vbh) +
				 rRbox*rNbox*HP_Rr_Hr*(Vbh) +
				 rNbox*rNbox*rRbox*HP_Rr_Hr2*(Vbh) +
				 rNbox*rNbox*rNbox*rRbox*HP_Rr_Hr3*(Vbh) +
				 HP_RNa*(tc*Vmaxh+Vbh) +
				 rNbox*HP_RNa_Hr*(tc*Vmaxh+Vbh) +
				 rNbox*rNbox*HP_RNa_Hr2*(tc*Vmaxh+Vbh) +
				 rNbox*rNbox*rNbox*HP_RNa_Hr3*(tc*Vmaxh+Vbh) +
				 rRbox*HP_RNa_Rr*(tc*Vmaxh+Vbh) +
				 rRbox*rNbox*HP_RNa_Rr_Hr*(tc*Vmaxh+Vbh) +
				 rRbox*rNbox*rNbox*HP_RNa_Rr_Hr2*(tc*Vmaxh+Vbh) +
				 rRbox*rNbox*rNbox*rNbox*HP_RNa_Rr_Hr3*(tc*Vmaxh+Vbh) +
				 rRbox*rRbox*HP_Rr2*(Vbh) +
				 rRbox*rRbox*rNbox*HP_Rr2_Hr*(Vbh) +
				 rRbox*rRbox*rNbox*rNbox*HP_Rr2_Hr2*(Vbh) +
				 rRbox*rRbox*rNbox*rNbox*rNbox*HP_Rr2_Hr3*(Vbh) +
				 HP_RNa2*(Vmaxh+Vbh) +
				 rNbox*HP_RNa2_Hr*(Vmaxh+Vbh) +
				 rNbox*rNbox*HP_RNa2_Hr2*(Vmaxh+Vbh) +
				 rNbox*rNbox*rNbox*HP_RNa2_Hr3*(Vmaxh+Vbh);


			H2np = KaHp*Hnp*Hnp;  // THIS MIGHT NEED ADJUSTMENT
			//Propensity functions

			//Generation reactions
			a[0]=c[0];		//1. 0 --> Nm		0.RfNm
			a[1]=c[1];		//2. 0 --> Rcm		1.RfRm
			a[2]=c[2];		//3. 0 --> Hcm		2.RfHm

			//Translation channels (delay these if they occur in delay model)
			a[3]=c[3]*Nm;		//4. Nm --> Np		3.KtrN
			a[4]=c[4]*Rcm;		//5. Rcm --> Rcp		4.KtrRc
			a[5]=c[5]*Hcm;		//6. Hcm --> Hcp		5.KtrHc

			//Other reactions in the cytoplasm
			//a[6]=c[6]*Hnp*Hnp;	//7. Hnp + Hnp --> H2np		6.KaHp
			//a[7]=c[7]*H2np;	//8. H2np --> Hnp + Hnp	7.KrHp
			a[8]=c[8]*delta*Np;	//9. delta + Np --> Ncp		8.KfNcp

			// Effect of GSK3b on NICD half life

			if(Hcp>100)
			{
			  c[12]=c[15]=0.0014;   // KdNnp = KdNcp = 0.0014 - 3hr half life
			}
			else
			{
			  c[12]=c[15]=0.00385;   // KdNnp = KdNcp = 0.00385 - 8hr half life
			}

			//Degeneration channels
			a[9]=c[9]*Nm;		//10. Nm --> 0		9.KdNm
			a[10]=c[10]*Rcm;		//11. Rcm --> 0		10.KdRcm
			a[11]=c[11]*Hcm;		//12. Hcm --> 0		11.KdHcm
			a[12]=c[12]*Nnp;		//13. Nnp --> 0		12.KdNnp
			a[13]=c[13]*Rnp;		//14. Rnp --> 0		13.KdRnp
			a[14]=c[14]*Np;		//15. Np --> 0		14.KdNp
			a[15]=c[15]*Ncp;		//16. Ncp --> 0		15.KdNcp
			a[16]=c[16]*Rcp;		//17. Rcp --> 0		16.KdRcp
			a[17]=c[17]*Hcp;		//18. Hcp --> 0		17.KdHcp
			a[18]=KdHnp*Hnp;		//19. Hnp --> 0		xx.KdHnp
			//a[19]=0;	//20. H2np --> 0	18.KdH2np

			//Protein transport channels cyt --> "nuc only transport into nucleus"
			a[20]=c[19]*Ncp;		//21. Ncp --> Nnp			19.KniNcp
			a[21]=c[20]*Rcp;		//22. Rcp --> Rnp			20.KniRcp
			a[22]=c[21]*Hcp;		//23. Hcp --> Hnp			21.KniHcp

			//Notch Promoter binding (18 reactions)
			a[23]=c[22]*2*NP_free*Rnp;		//24. NP_free      + Rnp  --> NP_Rr
			a[24]=c[24]*NP_free*H2np;		//25. NP_free      + H2np --> NP_Hr
			a[25]=c[24]*NP_Rr*H2np;			//26. NP_Rr        + H2np --> NP_Rr_Hr
			a[26]=c[26]*NP_Rr*Nnp;			//27. NP_Rr        + Nnp  --> NP_RNa
			a[27]=c[22]*NP_Rr*Rnp;			//28. NP_Rr        + Rnp  --> NP_Rr2
			a[28]=c[22]*2*NP_Hr*Rnp;			//29. NP_Hr        + Rnp  --> NP_Rr_Hr
			a[29]=c[22]*NP_Rr_Hr*Rnp;		//30. NP_Rr_Hr     + Rnp  --> NP_Rr2_Hr
			a[30]=c[26]*NP_Rr_Hr*Nnp;		//31. NP_Rr_Hr     + Nnp  --> NP_RNa_Hr
			a[31]=c[22]*NP_RNa*Rnp;			//32. NP_RNa       + Rnp  --> NP_RNa_Rr
			a[32]=c[24]*NP_RNa*H2np;			//33. NP_RNa       + H2np --> NP_RNa_Hr
			a[33]=c[26]*2*NP_Rr2*Nnp;		//34. NP_Rr2       + Nnp  --> NP_RNa_Rr
			a[34]=c[24]*NP_Rr2*H2np;			//35. NP_Rr2       + H2np --> NP_Rr2_Hr
			a[35]=c[26]*2*NP_Rr2_Hr*Nnp;		//36. NP_Rr2_Hr    + Nnp  --> NP_RNA_Rr_Hr
			a[36]=c[22]*NP_RNa_Hr*Rnp;		//37. NP_RNa_Hr    + Rnp  --> NP_RNA_Rr_Hr
			a[37]=c[26]*NP_RNa_Rr*Nnp;		//38. NP_RNa_Rr    + Nnp  --> NP_RNa2
			a[38]=c[24]*NP_RNa_Rr*H2np;	  	//39. NP_RNa_Rr    + H2np --> NP_RNa_RR_Hr
			a[39]=c[26]*NP_RNa_Rr_Hr*Nnp;    //40. NP_RNa_Rr_Hr + Nnp  --> NP_RNa2_Hr
			a[40]=c[24]*NP_RNa2*H2np;		//41. NP_RNa2      + H2np --> NP_RNa2_Hr

			// Notch promoter disassociation (18 reactions)
			a[41]=c[23]*NP_Rr;			//42. NP_Rr        --> NP_free      + Rnp
			a[42]=c[25]*NP_RNa;			//43. NP_RNa       --> NP_Rr        + Nnp
			a[43]=c[25]*NP_RNa_Rr;		//44. NP_RNa_Rr    --> NP_RNa       + Rnp
			a[44]=c[25]*NP_RNa_Rr;		//45. NP_RNa_Rr    --> NP_Rr2       + Nnp
			a[45]=c[25]*2*NP_Rr2;		//46. NP_Rr2       --> NP_Rr        + Rnp
			a[46]=c[25]*NP_RNa2;			//47. NP_RNa2      --> NP_RNa_Rr    + Nnp
			a[47]=c[23]*NP_Hr;			//48. NP_Hr        --> NP_free      + H2np
			a[48]=c[25]*NP_Rr_Hr;		//49. NP_Rr_Hr     --> NP_Rr        + H2np
			a[49]=c[25]*NP_Rr_Hr;		//50. NP_Rr_Hr     --> NP_Hr        + Rnp
			a[50]=c[25]*NP_RNa_Hr;		//51. NP_RNa_Hr    --> NP_Rr_Hr     + Nnp
			a[51]=c[25]*NP_RNa_Hr;		//52. NP_RNa_Hr    --> NP_RNa       + H2np
			a[52]=c[25]*NP_RNa_Rr_Hr;	//53. NP_RNa_Rr_Hr --> NP_RNa_Hr    + Rnp
			a[53]=c[25]*NP_RNa_Rr_Hr;	//54. NP_RNa_Rr_Hr --> NP_Rr2_Hr    + Nnp
			a[54]=c[25]*NP_RNa_Rr_Hr;	//55. NP_RNa_Rr_Hr --> NP_RNa_Rr    + H2np
			a[55]=c[25]*2*NP_Rr2_Hr;		//56. NP_Rr2_Hr    --> NP_Rr_Hr     + Rnp
			a[56]=c[25]*NP_Rr2_Hr;		//57. NP_Rr2_Hr    --> NP_Rr2       + H2np
			a[57]=c[25]*2*NP_RNa2_Hr;	//58. NP_RNa2_Hr   --> NP_RNa_Rr_Hr + Nnp
			a[58]=c[25]*NP_RNa2_Hr;		//59. NP_RNa2_Hr   --> NP_RNa2      + H2np

			// Hes1 promoter binding
			a[59]=c[22]*2*HP_free*Rnp;			//60. HP_free       + Rnp  --> HP_Rr
			a[60]=c[24]*3*HP_free*H2np;			//61. HP_free       + H2np --> HP_Hr
			a[61]=c[22]*2*HP_Hr*Rnp;				//62. HP_Hr         + Rnp  --> HP_Rr_Hr
			a[62]=c[24]*2*HP_Hr*H2np;			//63. HP_Hr         + H2np --> HP_Hr2
			a[63]=c[22]*2*HP_Hr2*Rnp;			//64. HP_Hr2        + Rnp  --> HP_Rr_Hr2
			a[64]=c[24]*HP_Hr2*H2np;				//65. HP_Hr2        + H2np --> HP_Hr3
			a[65]=c[22]*2*HP_Hr3*Rnp;			//66. HP_Hr3        + Rnp  --> HP_Rr_Hr3
			a[66]=c[22]*HP_Rr*Rnp;				//67. HP_Rr         + Rnp  --> HP_Rr2
			a[67]=c[22]*HP_Rr*Nnp;				//68. HP_Rr         + Nnp  --> HP_RNa
			a[68]=c[24]*3*HP_Rr*H2np;			//69. HP_Rr         + H2np --> HP_Rr_Hr
			a[69]=c[22]*HP_Rr_Hr*Rnp;			//70. HP_Rr_Hr      + Rnp  --> HP_Rr2_Hr
			a[70]=c[26]*HP_Rr_Hr*Nnp;			//71. HP_Rr_Hr      + Nnp  --> HP_RNa_Hr
			a[71]=c[24]*2*HP_Rr_Hr*H2np;			//72. HP_Rr_Hr      + H2np --> HP_Rr_Hr2
			a[72]=c[22]*HP_Rr_Hr2*Rnp;			//73. HP_Rr_Hr2     + Rnp  --> HP_Rr2_Hr2
			a[73]=c[26]*HP_Rr_Hr2*Nnp;			//74. HP_Rr_Hr2     + Nnp  --> HP_RNa_Hr2
			a[74]=c[24]*HP_Rr_Hr2*H2np;			//75. HP_Rr_Hr2     + H2np --> HP_Rr_Hr3
			a[75]=c[22]*HP_Rr_Hr3*Rnp;			//76. HP_Rr_Hr3     + Rnp  --> HP_Rr2_Hr3
			a[76]=c[26]*HP_Rr_Hr3*Nnp;			//77. HP_Rr_Hr3     + Nnp  --> HP_RNa_Hr3
			a[77]=c[22]*HP_RNa*Rnp;				//78. HP_RNa        + Rnp  --> HP_RNa_Rr
			a[78]=c[24]*3*HP_RNa*H2np;			//79. HP_RNa        + H2np --> HP_RNa_Hr
			a[79]=c[22]*HP_RNa_Hr*Rnp;			//80. HP_RNa_Hr     + Rnp  --> HP_RNa_Rr_Hr
			a[80]=c[24]*2*HP_RNa_Hr*H2np;		//81. HP_RNa_Hr     + H2np --> HP_RNa_Hr2
			a[81]=c[22]*HP_RNa_Hr2*Rnp;			//82. HP_RNa_Hr2    + Rnp  --> HP_RNa_Rr_Hr2
			a[82]=c[24]*HP_RNa_Hr2*H2np;			//83. HP_RNa_Hr2    + H2np --> HP_RNa_Hr3
			a[83]=c[22]*HP_RNa_Hr3*Rnp;			//84. HP_RNa_Hr3    + Rnp  --> HP_RNa_Rr_Hr3
			a[84]=c[26]*HP_RNa_Rr*Nnp;			//85. HP_RNa_Rr     + Nnp  --> HP_RNa2
			a[85]=c[24]*3*HP_RNa_Rr*H2np;		//86. HP_RNa_Rr     + H2np --> HP_RNa_Rr_Hr
			a[86]=c[26]*HP_RNa_Rr_Hr*Nnp;		//88. HP_RNa_Rr_Hr  + Nnp  --> HP_RNa2_Hr
			a[87]=c[24]*2*HP_RNa_Rr_Hr*H2np;		//89. HP_RNa_Rr_Hr  + H2np --> HP_RNa_Rr_Hr2
			a[88]=c[26]*HP_RNa_Rr_Hr2*Nnp;		//91. HP_RNa_Rr_Hr2 + Nnp  --> HP_RNa2_Hr2
			a[89]=c[24]*HP_RNa_Rr_Hr2*H2np;		//92. HP_RNa_Rr_Hr2 + H2np --> HP_RNa_Rr_Hr3
			a[90]=c[26]*HP_RNa_Rr_Hr3*Nnp;		//94. HP_RNa_Rr_Hr3 + Nnp  --> HP_RNa2_Hr3
			a[91]=c[26]*2*HP_Rr2*Nnp;			//96. HP_Rr2        + Nnp  --> HP_RNa_Rr
			a[92]=c[24]*3*HP_Rr2*H2np;			//97. HP_Rr2        + H2np --> HP_Rr2_Hr
			a[93]=c[26]*2*HP_Rr2_Hr*Nnp;			//99. HP_Rr2_Hr     + Nnp  --> HP_RNa_Rr_Hr
			a[94]=c[24]*2*HP_Rr2_Hr*H2np;		//100. HP_Rr2_Hr     + H2np --> HP_Rr2_Hr2
			a[95]=c[26]*2*HP_Rr2_Hr2*Nnp;		//102. HP_Rr2_Hr2    + Nnp  --> HP_RNa_Rr_Hr2
			a[96]=c[24]*HP_Rr2_Hr2*H2np;			//103. HP_Rr2_Hr2    + H2np --> HP_Rr2_Hr3
			a[97]=c[26]*2*HP_Rr2_Hr3*Nnp;		//105. HP_Rr2_Hr3    + Nnp  --> HP_RNa_Rr_Hr3
			a[98]=c[24]*3*HP_RNa2*H2np;			//106. HP_RNa2       + H2np --> HP_RNa2_Hr
			a[99]=c[24]*2*HP_RNa2_Hr*H2np;		//107. HP_RNa2_Hr    + H2np --> HP_RNa2_Hr2
			a[100]=c[24]*HP_RNa2_Hr2*H2np;		//108. HP_RNa2_Hr2   + H2np --> HP_RNa2_Hr3

			// Hes1 promoter disassociation
			a[101]=c[23]*HP_Hr;				//109. HP_Hr         --> HP_free       + H2np
			a[102]=c[25]*2*HP_Hr2;			//110. HP_Hr2        --> HP_Hr         + H2np
			a[103]=c[25]*3*HP_Hr3;			//111. HP_Hr3        --> HP_Hr2        + H2np
			a[104]=c[23]*HP_Rr;				//112. HP_Rr         --> HP_free       + Rnp
			a[105]=c[25]*HP_Rr_Hr;			//113. HP_Rr_Hr      --> HP_Hr         + Rnp
			a[106]=c[25]*HP_Rr_Hr;			//114. HP_Rr_Hr      --> HP_Rr         + H2np
			a[107]=c[25]*HP_Rr_Hr2;			//115. HP_Rr_Hr2     --> HP_Hr2        + Rnp
			a[108]=c[25]*2*HP_Rr_Hr2;		//116. HP_Rr_Hr2     --> HP_Rr_Hr      + H2np
			a[109]=c[25]*HP_Rr_Hr3;			//117. HP_Rr_Hr3     --> HP_Hr3        + Rnp
			a[110]=c[25]*3*HP_Rr_Hr3;		//118. HP_Rr_Hr3     --> HP_Rr_Hr2     + H2np
			a[111]=c[25]*HP_RNa;				//119. HP_RNa        --> HP_Rr         + Nnp
			a[112]=c[25]*HP_RNa_Hr;			//120. HP_RNa_Hr     --> HP_Rr_Hr      + Nnp
			a[113]=c[25]*HP_RNa_Hr;			//121. HP_RNa_Hr     --> HP_RNa        + H2np
			a[114]=c[25]*HP_RNa_Hr2;			//122. HP_RNa_Hr2    --> HP_Rr_Hr2     + Nnp
			a[115]=c[25]*2*HP_RNa_Hr2;		//123. HP_RNa_Hr2    --> HP_RNa_Hr     + H2np
			a[116]=c[25]*HP_RNa_Hr3;			//124. HP_RNa_Hr3    --> HP_Rr_Hr3     + Nnp
			a[117]=c[25]*3*HP_RNa_Hr3;		//125. HP_RNa_Hr3    --> HP_RNa_Hr2    + H2np
			a[118]=c[25]*HP_RNa_Rr;			//126. HP_RNa_Rr     --> HP_RNa        + Rnp
			a[119]=c[25]*HP_RNa_Rr;			//127. HP_RNa_Rr     --> HP_Rr2        + Nnp
			a[120]=c[25]*HP_RNa_Rr_Hr;		//128. HP_RNa_Rr_Hr  --> HP_RNa_Hr     + Rnp
			a[121]=c[25]*HP_RNa_Rr_Hr;		//129. HP_RNa_Rr_Hr  --> HP_Rr2_Hr     + Nnp
			a[122]=c[25]*HP_RNa_Rr_Hr;		//130. HP_RNa_Rr_Hr  --> HP_RNa_Rr     + H2np
			a[123]=c[25]*HP_RNa_Rr_Hr2;		//131. HP_RNa_Rr_Hr2 --> HP_RNa_Hr2    + Rnp
			a[124]=c[25]*HP_RNa_Rr_Hr2;		//132. HP_RNa_Rr_Hr2 --> HP_Rr2_Hr2    + Nnp
			a[125]=c[25]*2*HP_RNa_Rr_Hr2;	//133. HP_RNa_Rr_Hr2 --> HP_RNa_Rr_Hr  + H2np
			a[126]=c[25]*HP_RNa_Rr_Hr3;		//134. HP_RNa_Rr_Hr3 --> HP_RNa_Hr3    + Rnp
			a[127]=c[25]*HP_RNa_Rr_Hr3;		//135. HP_RNa_Rr_Hr3 --> HP_Rr2_Hr3    + Nnp
			a[128]=c[25]*3*HP_RNa_Rr_Hr3;	//136. HP_RNa_Rr_Hr3 --> HP_RNa_Rr_Hr2 + H2np
			a[129]=c[25]*2*HP_Rr2;			//137. HP_Rr2        --> HP_Rr         + Rnp
			a[130]=c[25]*HP_Rr2_Hr;			//138. HP_Rr2_Hr     --> HP_Rr_Hr      + Rnp
			a[131]=c[25]*HP_Rr2_Hr;			//139. HP_Rr2_Hr     --> HP_Rr2        + H2np
			a[132]=c[25]*HP_Rr2_Hr2;			//140. HP_Rr2_Hr2    --> HP_Rr_Hr2     + Rnp
			a[133]=c[25]*2*HP_Rr2_Hr2;		//141. HP_Rr2_Hr2    --> HP_Rr2_Hr     + H2np
			a[134]=c[25]*HP_Rr2_Hr3;			//142. HP_Rr2_Hr3    --> HP_Rr_Hr3     + Rnp
			a[135]=c[25]*3*HP_Rr2_Hr3;		//143. HP_Rr2_Hr3    --> HP_Rr2_Hr2    + H2np
			a[136]=c[25]*2*HP_RNa2;			//144. HP_RNa2       --> HP_RNa_Rr     + Nnp
			a[137]=c[25]*2*HP_RNa2_Hr;		//145. HP_RNa2_Hr    --> HP_RNa_Rr_Hr  + Nnp
			a[138]=c[25]*HP_RNa2_Hr;			//146. HP_RNa2_Hr    --> HP_RNa2       + H2np
			a[139]=c[25]*2*HP_RNa2_Hr2;		//147. HP_RNa2_Hr2   --> HP_RNa_Rr_Hr2 + Nnp
			a[140]=c[25]*2*HP_RNa2_Hr2;		//148. HP_RNa2_Hr2   --> HP_RNa2_Hr    + H2np
			a[141]=c[25]*2*HP_RNa2_Hr3;		//149. HP_RNa2_Hr3   --> HP_RNa_Rr_Hr3 + Nnp
			a[142]=c[25]*3*HP_RNa2_Hr3;		//150. HP_RNa2_Hr3   --> HP_RNa2_Hr2   + H2np

			//CBF1 promoter binding
			a[143]=c[22]*3*RP_free*Rnp;			//151. RP_free + Rnp         --> RP_Rr
			a[144]=c[24]*3*RP_free*H2np;			//152. RP_free + H2np        --> RP_Hr
			a[145]=c[22]*3*RP_Hr*Rnp;			//153. RP_Hr + Rnp           --> RP_Rr_Hr
			a[146]=c[24]*2*RP_Hr*H2np;			//154. RP_Hr + H2np          --> RP_Hr2
			a[147]=c[22]*3*RP_Hr2*Rnp;			//155. RP_Hr2 + Rnp          --> RP_Rr_Hr2
			a[148]=c[24]*RP_Hr2*H2np;			//156. RP_Hr2 + H2np         --> RP_Hr3
			a[149]=c[22]*3*RP_Hr3*Rnp;			//157. RP_Hr3 + Rnp          --> RP_Rr_Hr3
			a[150]=c[22]*2*RP_Rr*Rnp;			//158. RP_Rr + Rnp           --> RP_Rr2
			a[151]=c[26]*RP_Rr*Nnp;				//159. RP_Rr + Nnp           --> RP_RNa
			a[152]=c[24]*3*RP_Rr*H2np;			//160. RP_Rr + H2np          --> RP_Rr_Hr
			a[153]=c[22]*2*RP_Rr_Hr*Rnp;			//161. RP_Rr_Hr + Rnp        --> RP_Rr2_Hr
			a[154]=c[26]*RP_Rr_Hr*Nnp;			//162. RP_Rr_Hr + Nnp        --> RP_RNa_Hr
			a[155]=c[24]*2*RP_Rr_Hr*H2np;		//163. RP_Rr_Hr + H2np       --> RP_Rr_Hr2
			a[156]=c[22]*2*RP_Rr_Hr2*Rnp;		//164. RP_Rr_Hr2 + Rnp       --> RP_Rr2_Hr2
			a[157]=c[26]*RP_Rr_Hr2*Nnp;			//165. RP_Rr_Hr2 + Nnp       --> RP_RNa_Hr2
			a[158]=c[24]*RP_Rr_Hr2*H2np;			//166. RP_Rr_Hr2 + H2np      --> RP_Rr_Hr3
			a[159]=c[22]*2*RP_Rr_Hr3*Rnp;		//167. RP_Rr_Hr3 + Rnp       --> RP_Rr2_Hr3
			a[160]=c[26]*RP_Rr_Hr3*Nnp;			//168. RP_Rr_Hr3 + Nnp       --> RP_RNa_Hr3
			a[161]=c[22]*2*RP_RNa*Rnp;			//169. RP_RNa + Rnp          --> RP_RNa_Rr
			a[162]=c[24]*3*RP_RNa*H2np;			//170. RP_RNa + H2np         --> RP_RNa_Hr
			a[163]=c[22]*2*RP_RNa_Hr*Rnp;		//171. RP_RNa_Hr + Rnp       --> RP_RNa_Rr_Hr
			a[164]=c[24]*2*RP_RNa_Hr*H2np;		//172. RP_RNa_Hr + H2np      --> RP_RNa_Hr2
			a[165]=c[22]*2*RP_RNa_Hr2*Rnp;		//173. RP_RNa_Hr2 + Rnp      --> RP_RNa_Rr_Hr2
			a[166]=c[24]*RP_RNa_Hr2*H2np;		//174. RP_RNa_Hr2 + H2np     --> RP_RNa_Hr3
			a[167]=c[22]*2*RP_RNa_Hr3*Rnp;		//175. RP_RNa_Hr3 + Rnp      --> RP_RNa_Rr_Hr3
			a[168]=c[22]*RP_Rr2*Rnp;				//176. RP_Rr2 + Rnp          --> RP_Rr3
			a[169]=c[26]*2*RP_Rr2*Nnp;			//177. RP_Rr2 + Nnp          --> RP_RNa_Rr
			a[170]=c[24]*3*RP_Rr2*H2np;			//178. RP_Rr2 + H2np         --> RP_Rr2_Hr
			a[171]=c[24]*RP_Rr2_Hr*Rnp;			//179. RP_Rr2_Hr + Rnp       --> RP_Rr3_Hr
			a[172]=c[26]*2*RP_Rr2_Hr*Nnp;		//180. RP_Rr2_Hr + Nnp       --> RP_RNa_Rr_Hr
			a[173]=c[24]*2*RP_Rr2_Hr*H2np;		//181. RP_Rr2_Hr + H2np      --> RP_Rr2_Hr2
			a[174]=c[22]*RP_Rr2_Hr2*Rnp;			//182. RP_Rr2_Hr2 + Rnp      --> RP_Rr3_Hr2
			a[175]=c[26]*2*RP_Rr2_Hr2*Nnp;		//183. RP_Rr2_Hr2 + Nnp      --> RP_RNa_Rr_Hr2
			a[176]=c[24]*RP_Rr2_Hr2*H2np;		//184. RP_Rr2_Hr2 + H2np     --> RP_Rr2_Hr3
			a[177]=c[22]*RP_Rr2_Hr3*Rnp;			//185. RP_Rr2_Hr3 + Rnp      --> RP_Rr3_Hr3
			a[178]=c[26]*2*RP_Rr2_Hr3*Nnp;		//186. RP_Rr2_Hr3 + Nnp      --> RP_RNa_Rr_Hr3
			a[179]=c[22]*RP_RNa_Rr*Rnp;			//187. RP_RNa_Rr + Rnp       --> RP_RNa_Rr2
			a[180]=c[26]*RP_RNa_Rr*Nnp;			//188. RP_RNa_Rr + Nnp       --> RP_RNa2
			a[181]=c[24]*3*RP_RNa_Rr*H2np;		//189. RP_RNa_Rr + H2np      --> RP_RNa_Rr_Hr
			a[182]=c[22]*RP_RNa_Rr_Hr*Rnp;		//190. RP_RNa_Rr_Hr + Rnp    --> RP_Rr2_Hr
			a[183]=c[26]*RP_RNa_Rr_Hr*Nnp;		//191. RP_RNa_Rr_Hr + Nnp    --> RP_RNa2_Hr
			a[184]=c[24]*2*RP_RNa_Rr_Hr*H2np;	//192. RP_RNa_Rr_Hr + H2np   --> RP_RNa_Rr_Hr2
			a[185]=c[22]*RP_RNa_Rr_Hr2*Rnp;		//193. RP_RNa_Rr_Hr2 + Rnp   --> RP_RNa_Rr2_Hr2
			a[186]=c[26]*RP_RNa_Rr_Hr2*Nnp;		//194. RP_RNa_Rr_Hr2 + Nnp   --> RP_RNa2_Hr2
			a[187]=c[24]*RP_RNa_Rr_Hr2*H2np;		//195. RP_RNa_Rr_Hr2 + H2np  --> RP_RNa_Rr_Hr3
			a[188]=c[22]*RP_RNa_Rr_Hr3*Rnp;		//196. RP_RNa_Rr_Hr3 + Rnp   --> RP_RNa_Rr2_Hr3
			a[189]=c[26]*RP_RNa_Rr_Hr3*Nnp;		//197. RP_RNa_Rr_Hr3 + Nnp   --> RP_RNa2_Hr3
			a[190]=c[26]*3*RP_Rr3*Nnp;			//198. RP_Rr3 + Nnp          --> RP_RNa_Rr2
			a[191]=c[24]*3*RP_Rr3*H2np;			//199. RP_Rr3 + H2np         --> RP_Rr3_Hr
			a[192]=c[26]*3*RP_Rr3_Hr*Nnp;		//200. RP_Rr3_Hr + Nnp       --> RP_RNa_Rr2_Hr
			a[193]=c[24]*2*RP_Rr3_Hr*H2np;		//201. RP_Rr3_Hr + H2np      --> RP_Rr3_Hr2
			a[194]=c[26]*3*RP_Rr3_Hr2*Nnp;		//202. RP_Rr3_Hr2 + Nnp      --> RP_RNa_Rr2_Hr2
			a[195]=c[24]*RP_Rr3_Hr2*H2np;		//203. RP_Rr3_Hr2 + H2np     --> RP_Rr3_Hr3
			a[196]=c[26]*3*RP_Rr3_Hr3*Nnp;		//204. RP_Rr3_Hr3 + Nnp      --> RP_RNa_Rr2_Hr3
			a[197]=c[26]*2*RP_RNa_Rr2*Nnp;		//205. RP_RNa_Rr2 + Nnp      --> RP_RNa2_Rr
			a[198]=c[24]*3*RP_RNa_Rr2*H2np;		//206. RP_RNa_Rr2 + H2np     --> RP_RNa_Rr2_Hr
			a[199]=c[26]*2*RP_RNa_Rr2_Hr*Nnp;	//207. RP_RNa_Rr2_Hr + Nnp   --> RP_RNa2_Rr_Hr
			a[200]=c[24]*2*RP_RNa_Rr2_Hr*H2np;	//208. RP_RNa_Rr2_Hr + H2np  --> RP_RNa_Rr2_Hr2
			a[201]=c[26]*2*RP_RNa_Rr2_Hr2*Nnp;	//209. RP_RNa_Rr2_Hr2 + Nnp  --> RP_RNa2_Rr_Hr2
			a[202]=c[24]*RP_RNa_Rr2_Hr2*H2np;	//210. RP_RNa_Rr2_Hr2 + H2np --> RP_RNa_Rr2_Hr3
			a[203]=c[26]*2*RP_RNa_Rr2_Hr3*Nnp;	//211. RP_RNa_Rr2_Hr3 + Nnp  --> RP_RNa2_Rr_Hr3
			a[204]=c[22]*RP_RNa2*Rnp;			//212. RP_RNa2 + Rnp         --> RP_RNa2_Rr
			a[205]=c[24]*3*RP_RNa2*H2np;			//213. RP_RNa2 + H2np        --> RP_RNa2_Hr
			a[206]=c[22]*RP_RNa2_Hr*Rnp;			//214. RP_RNa2_Hr + Rnp      --> RP_RNa2_Rr_Hr
			a[207]=c[24]*2*RP_RNa2_Hr*H2np;		//215. RP_RNa2_Hr + H2np     --> RP_RNa2_Hr2
			a[208]=c[22]*RP_RNa2_Hr2*Rnp;		//216. RP_RNa2_Hr2 + Rnp     --> RP_RNa2_Rr_Hr2
			a[209]=c[24]*RP_RNa2_Hr2*H2np;		//217. RP_RNa2_Hr2 + H2np    --> RP_RNa2_Hr3
			a[210]=c[22]*RP_RNa2_Hr3*Rnp;		//218. RP_RNa2_Hr3 + Rnp     --> RP_RNa2_Rr_Hr3
			a[211]=c[26]*RP_RNa2_Rr*Nnp;			//219. RP_RNa2_Rr + Nnp      --> RP_RNa3
			a[212]=c[24]*3*RP_RNa2_Rr*H2np;		//220. RP_RNa2_Rr + H2np     --> RP_RNa2_Rr_Hr
			a[213]=c[26]*RP_RNa2_Rr_Hr*Nnp;		//221. RP_RNa2_Rr_Hr + Nnp   --> RP_RNa3_Hr
			a[214]=c[24]*2*RP_RNa2_Rr_Hr*H2np;	//222. RP_RNa2_Rr_Hr + H2np  --> RP_RNa2_Rr_Hr2
			a[215]=c[26]*RP_RNa2_Rr_Hr2*Nnp;		//223. RP_RNa2_Rr_Hr2 + Nnp  --> RP_RNa3_Hr2
			a[216]=c[24]*RP_RNa2_Rr_Hr2*H2np;	//224. RP_RNa2_Rr_Hr2 + H2np --> RP_RNa2_Rr_Hr3
			a[217]=c[26]*RP_RNa2_Rr_Hr3*Nnp;		//225. RP_RNa2_Rr_Hr3 + Nnp  --> RP_RNa3_Hr3
			a[218]=c[24]*3*RP_RNa3*H2np;			//226. RP_RNa3 + H2np        --> RP_RNa3_Hr
			a[219]=c[24]*2*RP_RNa3_Hr*H2np;		//227. RP_RNa3_Hr + H2np     --> RP_RNa3_Hr2
			a[220]=c[24]*RP_RNa3_Hr2*H2np;		//228. RP_RNa3_Hr2 + H2np    --> RP_RNa3_Hr3

			//CBF1 promoter disassociation
			a[221]=c[23]*RP_Hr;				//229. RP_Hr		--> RP_free + H2np
			a[222]=c[25]*2*RP_Hr2;			//230. RP_Hr2		--> RP_Hr + H2np
			a[223]=c[25]*3*RP_Hr3;			//231. RP_Hr3		--> RP_Hr2 + H2np
			a[224]=c[23]*RP_Rr;				//232. RP_Rr		--> RP_free + Rnp
			a[225]=c[25]*RP_Rr_Hr;			//233. RP_Rr_Hr		--> RP_Hr + Rnp
			a[226]=c[25]*RP_Rr_Hr;			//234. RP_Rr_Hr		--> RP_Rr + H2np
			a[227]=c[25]*RP_Rr_Hr2;			//235. RP_Rr_Hr2	--> RP_Hr2 + Rnp
			a[228]=c[25]*2*RP_Rr_Hr2;		//236. RP_Rr_Hr2	--> RP_Rr_Hr + H2np
			a[229]=c[25]*RP_Rr_Hr3;			//237. RP_Rr_Hr3	--> RP_Hr3 + Rnp
			a[230]=c[25]*3*RP_Rr_Hr3;		//238. RP_Rr_Hr3	--> RP_Rr_Hr2 + H2np
			a[231]=c[25]*RP_RNa;				//239. RP_RNa		--> RP_Rr + Nnp
			a[232]=c[25]*RP_RNa_Hr;			//240. RP_RNa_Hr	--> RP_Rr_Hr + Nnp
			a[233]=c[25]*RP_RNa_Hr;			//241. RP_RNa_Hr	--> RP_RNa + H2np
			a[234]=c[25]*RP_RNa_Hr2;			//242. RP_RNa_Hr2	--> RP_Rr_Hr2 + Nnp
			a[235]=c[25]*2*RP_RNa_Hr2;		//243. RP_RNa_Hr2	--> RP_RNa_Hr + H2np
			a[236]=c[25]*RP_RNa_Hr3;			//244. RP_RNa_Hr3	--> RP_Rr_Hr3 + Nnp
			a[237]=c[25]*3*RP_RNa_Hr3;		//245. RP_RNa_Hr3	--> RP_RNa_Hr2 + H2np
			a[238]=c[25]*2*RP_Rr2;			//246. RP_Rr2		--> RP_Rr + Rnp
			a[239]=c[25]*2*RP_Rr2_Hr;		//247. RP_Rr2_Hr	--> RP_Rr_Hr + Rnp
			a[240]=c[25]*RP_Rr2_Hr;			//248. RP_Rr2_Hr	--> RP_Rr2 + H2np
			a[241]=c[25]*2*RP_Rr2_Hr2;		//249. RP_Rr2_Hr2	--> RP_Rr_Hr2 + Rnp
			a[242]=c[25]*2*RP_Rr2_Hr2;		//250. RP_Rr2_Hr2	--> RP_Rr2_Hr + H2np
			a[243]=c[25]*2*RP_Rr2_Hr3;		//251. RP_Rr2_Hr3	--> RP_Rr_Hr3 + Rnp
			a[244]=c[25]*3*RP_Rr2_Hr3;		//252. RP_Rr2_Hr3	--> RP_Rr2_Hr2 + H2np
			a[245]=c[25]*RP_RNa_Rr;			//253. RP_RNa_Rr	--> RP_RNa + Rnp
			a[246]=c[25]*RP_RNa_Rr;			//254. RP_RNa_Rr	--> RP_Rr2 + Nnp
			a[247]=c[25]*RP_RNa_Rr_Hr;		//255. RP_RNa_Rr_Hr	--> RP_RNa_Hr + Rnp
			a[248]=c[25]*RP_RNa_Rr_Hr;		//256. RP_RNa_Rr_Hr	--> RP_Rr2_Hr + Nnp
			a[249]=c[25]*RP_RNa_Rr_Hr;		//257. RP_RNa_Rr_Hr	--> RP_RNa_Rr + H2np
			a[250]=c[25]*RP_RNa_Rr_Hr2;		//258. RP_RNa_Rr_Hr2	--> RP_RNa_Hr2 + Rnp
			a[251]=c[25]*RP_RNa_Rr_Hr2;		//259. RP_RNa_Rr_Hr2	--> RP_Rr2_Hr2 + Nnp
			a[252]=c[25]*2*RP_RNa_Rr_Hr2;	//260. RP_RNa_Rr_Hr2	--> RP_RNa_Rr_Hr+ H2np
			a[253]=c[25]*RP_RNa_Rr_Hr3;		//261. RP_RNa_Rr_Hr3	--> RP_RNa_Hr3 + Rnp
			a[254]=c[25]*RP_RNa_Rr_Hr3;		//262. RP_RNa_Rr_Hr3	--> RP_Rr2_Hr3 + Nnp
			a[255]=c[25]*3*RP_RNa_Rr_Hr3;	//263. RP_RNa_Rr_Hr3	--> RP_RNa_Rr_Hr2 + H2np
			a[256]=c[25]*3*RP_Rr3;			//264. RP_Rr3		--> RP_Rr2 + Rnp
			a[257]=c[25]*3*RP_Rr3_Hr;		//265. RP_Rr3_Hr	--> RP_Rr2_Hr + Rnp
			a[258]=c[25]*RP_Rr3_Hr;			//266. RP_Rr3_Hr	--> RP_Rr3 + H2np
			a[259]=c[25]*3*RP_Rr3_Hr2;		//267. RP_Rr3_Hr2	--> RP_Rr2_Hr2 + Rnp
			a[260]=c[25]*2*RP_Rr3_Hr2;		//268. RP_Rr3_Hr2	--> RP_Rr3_Hr + H2np
			a[261]=c[25]*3*RP_Rr3_Hr3;		//269. RP_Rr3_Hr3	--> RP_Rr2_Hr3 + Rnp
			a[262]=c[25]*3*RP_Rr3_Hr3;		//270. RP_Rr3_Hr3	--> RP_Rr3_Hr2 + H2np
			a[263]=c[25]*2*RP_RNa_Rr2;		//271. RP_RNa_Rr2	--> RP_RNa_Rr + Rnp
			a[264]=c[25]*RP_RNa_Rr2;			//272. RP_RNa_Rr2	--> RP_Rr3 + Nnp
			a[265]=c[25]*2*RP_RNa_Rr2_Hr;	//273. RP_RNa_Rr2_Hr	--> RP_RNa_Rr_Hr + Rnp
			a[266]=c[25]*RP_RNa_Rr2_Hr;		//274. RP_RNa_Rr2_Hr	--> RP_Rr3_Hr + Nnp
			a[267]=c[25]*RP_RNa_Rr2_Hr;		//275. RP_RNa_Rr2_Hr	--> RP_RNa_Rr2 + H2np
			a[268]=c[25]*2*RP_RNa_Rr2_Hr2;	//276. RP_RNa_Rr2_Hr2	--> RP_RNa_Rr_Hr2 + Rnp
			a[269]=c[25]*RP_RNa_Rr2_Hr2;		//277. RP_RNa_Rr2_Hr2	--> RP_Rr3_Hr2 + Nnp
			a[270]=c[25]*2*RP_RNa_Rr2_Hr2;	//278. RP_RNa_Rr2_Hr2	--> RP_RNa_Rr2_Hr + H2np
			a[271]=c[25]*2*RP_RNa_Rr2_Hr3;	//279. RP_RNa_Rr2_Hr3	--> RP_RNa_Rr_Hr3 + Rnp
			a[272]=c[25]*RP_RNa_Rr2_Hr3;		//280. RP_RNa_Rr2_Hr3	--> RP_Rr3_Hr3 + Nnp
			a[273]=c[25]*3*RP_RNa_Rr2_Hr3;	//281. RP_RNa_Rr2_Hr3	--> RP_RNa_Rr2_Hr2 + H2np
			a[274]=c[25]*2*RP_RNa2;			//282. RP_RNa2		--> RP_RNa_Rr + Nnp
			a[275]=c[25]*2*RP_RNa2_Hr;		//283. RP_RNa2_Hr	--> RP_RNa_Rr_Hr + Nnp
			a[276]=c[25]*RP_RNa2_Hr;			//284. RP_RNa2_Hr	--> RP_RNa2 + H2np
			a[277]=c[25]*2*RP_RNa2_Hr2;		//285. RP_RNa2_Hr2	--> RP_RNa_Rr_Hr2 + Nnp
			a[278]=c[25]*2*RP_RNa2_Hr2;		//286. RP_RNa2_Hr2	--> RP_RNa2_Hr + H2np
			a[279]=c[25]*2*RP_RNa2_Hr3;		//287. RP_RNa2_Hr3	--> RP_RNa_Rr_Hr3 + Nnp
			a[280]=c[25]*3*RP_RNa2_Hr3;		//288. RP_RNa2_Hr3	--> RP_RNa2_Hr2 + H2np
			a[281]=c[25]*RP_RNa2_Rr;			//289. RP_RNa2_Rr 	--> RP_RNa2 + Rnp
			a[282]=c[25]*2*RP_RNa2_Rr;		//290. RP_RNa2_Rr 	--> RP_RNa_Rr2 + Nnp
			a[283]=c[25]*RP_RNa2_Rr_Hr;		//291. RP_RNa2_Rr_Hr 	--> RP_RNa2_Hr + Rnp
			a[284]=c[25]*2*RP_RNa2_Rr_Hr;	//292. RP_RNa2_Rr_Hr 	--> RP_RNa_Rr2_Hr + Nnp
			a[285]=c[25]*RP_RNa2_Rr_Hr;		//293. RP_RNa2_Rr_Hr 	--> RP_RNa2_Rr + H2np
			a[286]=c[25]*RP_RNa2_Rr_Hr2;		//294. RP_RNa2_Rr_Hr2	--> RP_RNa2_Hr + Rnp
			a[287]=c[25]*2*RP_RNa2_Rr_Hr2;	//295. RP_RNa2_Rr_Hr2	--> RP_RNa_Rr2_Hr2 + Nnp
			a[288]=c[25]*2*RP_RNa2_Rr_Hr2;	//296. RP_RNa2_Rr_Hr2	--> RP_RNa2_Rr_Hr + H2np
			a[289]=c[25]*RP_RNa2_Rr_Hr3;		//297. RP_RNa2_Rr_Hr3 	--> RP_RNa2_Hr3 + Rnp
			a[290]=c[25]*2*RP_RNa2_Rr_Hr3;	//298. RP_RNa2_Rr_Hr3 	--> RP_RNa_Rr2_Hr3 + Nnp
			a[291]=c[25]*3*RP_RNa2_Rr_Hr3;	//299. RP_RNa2_Rr_Hr3 	--> RP_RNa2_Rr_Hr2 + H2np
			a[292]=c[25]*3*RP_RNa3;			//300. RP_RNa3 		--> RP_RNa2_Rr + Nnp
			a[293]=c[25]*3*RP_RNa3_Hr;		//301. RP_RNa3_Hr 	--> RP_RNa2_Rr_Hr + Nnp
			a[294]=c[25]*RP_RNa3_Hr;			//302. RP_RNa3_Hr 	--> RP_RNa3 + H2np
			a[295]=c[25]*3*RP_RNa3_Hr2;		//303. RP_RNa3_Hr2 	--> RP_RNa2_Rr_Hr2 + Nnp
			a[296]=c[25]*2*RP_RNa3_Hr2;		//304. RP_RNa3_Hr2 	--> RP_RNa3_Hr + H2np
			a[297]=c[25]*3*RP_RNa3_Hr3;		//305. RP_RNa3_Hr3 	--> RP_RNa2_Rr_Hr3 + Nnp
			a[298]=c[25]*3*RP_RNa3_Hr3;		//306. RP_RNa3_Hr3 	--> RP_RNa3_Hr2 + H2np


			// find a0 for weighting the rxn's (sum of the A array)
			a0=0;
			for(k=0;k<numRxns+1;k++) {
				a0+=a[k];
			}

			//Use the Mersenne Twister algorithm to generate the random numbers
			MTRand mtr;
			r1=(float)mtr.rand();
			r2=(float)mtr.rand();

			//Testing random number. need same seed
			//r1=(float)rand()/RAND_MAX;
			//r2=(float)rand()/RAND_MAX;

			// Increment the time step
			tau=1/a0*log(1/r1);

			//Check if we have delayed reactions which would have occurred between t and t+tau
			if(delay && delayedRxnReady(t+tau)) {
				//update system with delayed reaction occuring and to time t when delayed reaction occurs
				t=fireRxn();
			} else {

				// Deciding which reaction occured and updating
				sum=0;
				mu=0;

				//find which reaction fired
				do {
					sum+=a[mu];
					mu++;
				}
				while(sum <= r2*a0);

				t += tau;

				//cout<<mu<<endl;
				switch(mu-1) {

					//delay
					//1. 0 --> Nm		0.RfNm
					case 0:
						if(delay) {addRxn(0,t+tau, TpNc);}
						else {
							Nm++;
						}
						break;
					//delay
					//2. 0 --> Rcm		1.RfRm
					case 1:
						if(delay) {addRxn(1,t+tau, TpRc);}
						else {
							Rcm++;
						}
						break;

					//delay
					//3. 0 --> Hcm		2.RfHm
					case 2:
						if(delay) {addRxn(2,t+tau, TpHc);}
						else {
							Hcm++;
						}
						break;

					//delay
					//4. Nm --> Np		3.KtrN
					case 3:
						if(delay) {addRxn(3,t+tau, TmNc);}
						else {
							//Nm--;
							Np++;
						}
						break;

					//delay
					//5. Rcm --> Rcp		4.KtrRc
					case 4:
						if(delay) {addRxn(4,t+tau, TmRc);}
						else {
							//Rcm--;
							Rcp++;
						}
						break;

					//delay
					//6. Hcm --> Hcp		5.KtrHc
					case 5:
						if(delay) {addRxn(5,t+tau, TmHc);}
						else {
							//Hcm--;
							Hcp++;
						}
						break;

					//7. Hnp + Hnp --> H2np		6.KaHp
					case 6:
						Hnp--;
						Hnp--;
						H2np++;
						break;

					//8. H2np --> Hnp + Hnp		7.KrHp
					case 7:
						H2np--;
						Hnp++;
						Hnp++;
						break;

					//9. delta + Np --> Ncp		8.KfNcp
					case 8:
						delta--;
						Np--;
						Ncp++;
						break;

					//10. Nm --> 0		9.KdNm
					case 9:
						Nm--;
						break;

					//11. Rcm --> 0		10.KdRcm
					case 10:
						Rcm--;
						break;

					//12 .Hcm --> 0		11.KdHcm
					case 11:
						Hcm--;
						break;

					//13. Nnp --> 0		12.KdNnp
					case 12:
						Nnp--;
						break;

					//14. Rnp --> 0		13.KdRnp
					case 13:
						Rnp--;
						break;

					//15. Np --> 0		14.KdNp
					case 14:
						Np--;
						break;

					//16. Ncp --> 0		15.KdNcp
					case 15:
						Ncp--;
						break;

					//17. Rcp --> 0		16.KdRcp
					case 16:
						Rcp--;
						break;

					//18. Hcp --> 0		17.KdHcp
					case 17:
						Hcp--;
						break;

					//19. Hnp --> 0	xx.KdHnp
					case 18:
						Hnp--;
						break;

					//19. H2np --> 0	18.KdH2np
					case 19:
						H2np--;
						break;

					//20. Ncp --> Nnp			19.KniNcp
					case 20:
						Ncp--;
						Nnp++;
						break;

					//21. Rcp --> Rnp			20.KniRcp
					case 21:
						Rcp--;
						Rnp++;
						break;

					//22. Hcp --> Hnp			21.KniHcp
					case 22:
						Hcp--;
						Hnp++;
						break;

					//23. NP_free      + Rnp  --> NP_Rr
					case 23:
						NP_free--;
						Rnp--;
						NP_Rr++;
						break;

					//24. NP_free      + H2np --> NP_Hr
					case 24:
						NP_free--;
						H2np--;
						NP_Hr++;
						break;

					//25. NP_Rr        + H2np --> NP_Rr_Hr
					case 25:
						NP_Rr--;
						H2np--;
						NP_Rr_Hr++;
						break;

					//26. NP_Rr        + Nnp  --> NP_RNa
					case 26:
						NP_Rr--;
						Nnp--;
						NP_RNa++;
						break;

					//27. NP_Rr        + Rnp  --> NP_Rr2
					case 27:
						NP_Rr--;
						Rnp--;
						NP_Rr2++;
						break;

					//28. NP_Hr        + Rnp  --> NP_Rr_Hr
					case 28:
						NP_Hr--;
						Rnp--;
						NP_Rr_Hr++;
						break;

					//29. NP_Rr_Hr     + Rnp  --> NP_Rr2_Hr
					case 29:
						NP_Rr_Hr--;
						Rnp--;
						NP_Rr2_Hr++;
						break;

					//30. NP_Rr_Hr     + Nnp  --> NP_RNa_Hr
					case 30:
						NP_Rr_Hr--;
						Nnp--;
						NP_RNa_Hr++;
						break;

					//31. NP_RNa       + Rnp  --> NP_RNa_Rr
					case 31:
						NP_RNa--;
						Rnp--;
						NP_RNa_Rr++;
						break;

					//32. NP_RNa       + H2np --> NP_RNa_Hr
					case 32:
						NP_RNa--;
						H2np--;
						NP_RNa_Hr++;
						break;

					//33. NP_Rr2       + Nnp  --> NP_RNa_Rr
					case 33:
						NP_Rr2--;
						Nnp--;
						NP_RNa_Rr++;
						break;

					//34. NP_Rr2       + H2np --> NP_Rr2_Hr
					case 34:
						NP_Rr2--;
						H2np--;
						NP_Rr2_Hr++;
						break;

					//35. NP_Rr2_Hr    + Nnp  --> NP_RNa_Rr_Hr
					case 35:
						NP_Rr2_Hr--;
						Nnp--;
						NP_RNa_Rr_Hr++;
						break;

					//36. NP_RNa_Hr    + Rnp  --> NP_RNa_Rr_Hr
					case 36:
						NP_RNa_Hr--;
						Rnp--;
						NP_RNa_Rr_Hr++;
						break;

					//37. NP_RNa_Rr    + Nnp  --> NP_RNa2
					case 37:
						NP_RNa_Rr--;
						Nnp--;
						NP_RNa2++;
						break;

					//38. NP_RNa_Rr    + H2np --> NP_RNa_Rr_Hr
					case 38:
						NP_RNa_Rr--;
						H2np--;
						NP_RNa_Rr_Hr++;
						break;

					//39. NP_RNa_Rr_Hr + Nnp  --> NP_RNa2_Hr
					case 39:
						NP_RNa_Rr_Hr--;
						Nnp--;
						NP_RNa2_Hr++;
						break;

					//40. NP_RNa2      + H2np --> NP_RNa2_Hr
					case 40:
						NP_RNa2--;
						H2np--;
						NP_RNa2_Hr++;
						break;

					//41. NP_Rr        --> NP_free      + Rnp
					case 41:
						NP_Rr--;
						NP_free++;
						Rnp++;
						break;

					//42. NP_RNa       --> NP_Rr        + Nnp
					case 42:
						NP_RNa--;
						NP_Rr++;
						Nnp++;
						break;

					//43. NP_RNa_Rr    --> NP_RNa       + Rnp
					case 43:
						NP_RNa_Rr--;
						NP_RNa++;
						Rnp++;
						break;

					//44. NP_RNa_Rr    --> NP_Rr2       + Nnp
					case 44:
						NP_RNa_Rr--;
						NP_Rr2++;
						Nnp++;
						break;

					//45. NP_Rr2       --> NP_Rr        + Rnp
					case 45:
						NP_Rr2--;
						NP_Rr++;
						Rnp++;
						break;

					//46. NP_RNa2      --> NP_RNa_Rr    + Nnp
					case 46:
						NP_RNa2--;
						NP_RNa_Rr++;
						Nnp++;
						break;

					//47. NP_Hr        --> NP_free      + H2np
					case 47:
						NP_Hr--;
						NP_free++;
						H2np++;
						break;

					//48. NP_Rr_Hr     --> NP_Rr        + H2np
					case 48:
						NP_Rr_Hr--;
						NP_Rr++;
						H2np++;
						break;

					//49. NP_Rr_Hr     --> NP_Hr        + Rnp
					case 49:
						NP_Rr_Hr--;
						NP_Hr++;
						Rnp++;
						break;

					//50. NP_RNa_Hr    --> NP_Rr_Hr     + Nnp
					case 50:
						NP_RNa_Hr--;
						NP_Rr_Hr++;
						Nnp++;
						break;

					//51. NP_RNa_Hr    --> NP_RNa       + H2np
					case 51:
						NP_RNa_Hr--;
						NP_RNa++;
						H2np++;
						break;

					//52. NP_RNa_Rr_Hr --> NP_RNa_Hr    + Rnp
					case 52:
						NP_RNa_Rr_Hr--;
						NP_RNa_Hr++;
						Rnp++;
						break;

					//53. NP_RNa_Rr_Hr --> NP_Rr2_Hr    + Nnp
					case 53:
						NP_RNa_Rr_Hr--;
						NP_Rr2_Hr++;
						Nnp++;
						break;

					//54. NP_RNa_Rr_Hr --> NP_RNa_Rr    + H2np
					case 54:
						NP_RNa_Rr_Hr--;
						NP_RNa_Rr++;
						H2np++;
						break;

					//55. NP_Rr2_Hr    --> NP_Rr_Hr     + Rnp
					case 55:
						NP_Rr2_Hr--;
						NP_Rr_Hr++;
						Rnp++;
						break;

					//56. NP_Rr2_Hr    --> NP_Rr2       + H2np
					case 56:
						NP_Rr2_Hr--;
						NP_Rr2++;
						H2np++;
						break;

					//57. NP_RNa2_Hr   --> NP_RNa_Rr_Hr + Nnp
					case 57:
						NP_RNa2_Hr--;
						NP_RNa_Rr_Hr++;
						Nnp++;
						break;

					//58. NP_RNa2_Hr   --> NP_RNa2      + H2np
					case 58:
						NP_RNa2_Hr--;
						NP_RNa2++;
						H2np++;
						break;

					//59. HP_free       + Rnp  --> HP_Rr
					case 59:
						HP_free--;
						Rnp--;
						HP_Rr++;
						break;

					//60. HP_free       + H2np --> HP_Hr
					case 60:
						HP_free--;
						H2np--;
						HP_Hr++;
						break;

					//61. HP_Hr         + Rnp  --> HP_Rr_Hr
					case 61:
						HP_Hr--;
						Rnp--;
						HP_Rr_Hr++;
						break;

					//62. HP_Hr         + H2np --> HP_Hr2
					case 62:
						HP_Hr--;
						H2np--;
						HP_Hr2++;
						break;

					//63. HP_Hr2        + Rnp  --> HP_Rr_Hr2
					case 63:
						HP_Hr2--;
						Rnp--;
						HP_Rr_Hr2++;
						break;

					//64. HP_Hr2        + H2np --> HP_Hr3
					case 64:
						HP_Hr2--;
						H2np--;
						HP_Hr3++;
						break;

					//65. HP_Hr3        + Rnp  --> HP_Rr_Hr3
					case 65:
						HP_Hr3--;
						Rnp--;
						HP_Rr_Hr3++;
						break;

					//66. HP_Rr         + Rnp  --> HP_Rr2
					case 66:
						HP_Rr--;
						Rnp--;
						HP_Rr2++;
						break;

					//67. HP_Rr         + Nnp  --> HP_RNa
					case 67:
						HP_Rr--;
						Nnp--;
						HP_RNa++;
						break;

					//68. HP_Rr         + H2np --> HP_Rr_Hr
					case 68:
						HP_Rr--;
						H2np--;
						HP_Rr_Hr++;
						break;

					//69. HP_Rr_Hr      + Rnp  --> HP_Rr2_Hr
					case 69:
						HP_Rr_Hr--;
						Rnp--;
						HP_Rr2_Hr++;
						break;

					//70. HP_Rr_Hr      + Nnp  --> HP_RNa_Hr
					case 70:
						HP_Rr_Hr--;
						Nnp--;
						HP_RNa_Hr++;
						break;

					//71. HP_Rr_Hr      + H2np --> HP_Rr_Hr2
					case 71:
						HP_Rr_Hr--;
						H2np--;
						HP_Rr_Hr2++;
						break;

					//72. HP_Rr_Hr2     + Rnp  --> HP_Rr2_Hr2
					case 72:
						HP_Rr_Hr2--;
						Rnp--;
						HP_Rr2_Hr2++;
						break;

					//73. HP_Rr_Hr2     + Nnp  --> HP_RNa_Hr2
					case 73:
						HP_Rr_Hr2--;
						Nnp--;
						HP_RNa_Hr2++;
						break;

					//74. HP_Rr_Hr2     + H2np --> HP_Rr_Hr3
					case 74:
						HP_Rr_Hr2--;
						H2np--;
						HP_Rr_Hr3++;
						break;

					//75. HP_Rr_Hr3     + Rnp  --> HP_Rr2_Hr3
					case 75:
						HP_Rr_Hr3--;
						Rnp--;
						HP_Rr2_Hr3++;
						break;

					//76. HP_Rr_Hr3     + Nnp  --> HP_RNa_Hr3
					case 76:
						HP_Rr_Hr3--;
						Nnp--;
						HP_RNa_Hr3++;
						break;

					//77. HP_RNa        + Rnp  --> HP_RNa_Rr
					case 77:
						HP_RNa--;
						Rnp--;
						HP_RNa_Rr++;
						break;

					//78. HP_RNa        + H2np --> HP_RNa_Hr
					case 78:
						HP_RNa--;
						H2np--;
						HP_RNa_Hr++;
						break;

					//79. HP_RNa_Hr     + Rnp  --> HP_RNa_Rr_Hr
					case 79:
						HP_RNa_Hr--;
						Rnp--;
						HP_RNa_Rr_Hr++;
						break;

					//80. HP_RNa_Hr     + H2np --> HP_RNa_Hr2
					case 80:
						HP_RNa_Hr--;
						H2np--;
						HP_RNa_Hr2++;
						break;

					//81. HP_RNa_Hr2    + Rnp  --> HP_RNa_Rr_Hr2
					case 81:
						HP_RNa_Hr2--;
						Rnp--;
						HP_RNa_Rr_Hr2++;
						break;

					//82. HP_RNa_Hr2    + H2np --> HP_RNa_Hr3
					case 82:
						HP_RNa_Hr2--;
						H2np--;
						HP_RNa_Hr3++;
						break;

					//83. HP_RNa_Hr3    + Rnp  --> HP_RNa_Rr_Hr3
					case 83:
						HP_RNa_Hr3--;
						Rnp--;
						HP_RNa_Rr_Hr3++;
						break;

					//85. HP_RNa_Rr     + Nnp  --> HP_RNa2
					case 84:
						HP_RNa_Rr--;
						Nnp--;
						HP_RNa2++;
						break;

					//86. HP_RNa_Rr     + H2np --> HP_RNa_Rr_Hr
					case 85:
						HP_RNa_Rr--;
						H2np--;
						HP_RNa_Rr_Hr++;
						break;

					//88. HP_RNa_Rr_HR  + Nnp  --> HP_RNa2_Hr
					case 86:
						HP_RNa_Rr_Hr--;
						Nnp--;
						HP_RNa2_Hr++;
						break;

					//89. HP_RNa_Rr_HR  + H2np --> HP_RNa_Rr_Hr2
					case 87:
						HP_RNa_Rr_Hr--;
						H2np--;
						HP_RNa_Rr_Hr2++;
						break;

					//91. HP_RNa_Rr_HR2 + Nnp  --> HP_RNa2_Hr2
					case 88:
						HP_RNa_Rr_Hr2--;
						Nnp--;
						HP_RNa2_Hr2++;
						break;

					//92. HP_RNa_Rr_HR2 + H2np --> HP_RNa_Rr_Hr3
					case 89:
						HP_RNa_Rr_Hr2--;
						H2np--;
						HP_RNa_Rr_Hr3++;
						break;

					//94. HP_RNa_Rr_HR3 + Nnp  --> HP_RNa2_Hr3
					case 90:
						HP_RNa_Rr_Hr3--;
						Nnp--;
						HP_RNa2_Hr3++;
						break;

					//96. HP_Rr2        + Nnp  --> HP_RNa_Rr
					case 91:
						HP_Rr2--;
						Nnp--;
						HP_RNa_Rr++;
						break;

					//97. HP_Rr2        + H2np --> HP_Rr2_Hr
					case 92:
						HP_Rr2--;
						H2np--;
						HP_Rr2_Hr++;
						break;

					//99. HP_Rr2_Hr     + Nnp  --> HP_RNa_Rr_Hr
					case 93:
						HP_Rr2_Hr--;
						Nnp--;
						HP_RNa_Rr_Hr++;
						break;

					//100. HP_Rr2_Hr     + H2np --> HP_Rr2_Hr2
					case 94:
						HP_Rr2_Hr--;
						H2np--;
						HP_Rr2_Hr2++;
						break;

					//102. HP_Rr2_Hr2    + Nnp  --> HP_RNa_Rr_Hr2
					case 95:
						HP_Rr2_Hr2--;
						Nnp--;
						HP_RNa_Rr_Hr2++;
						break;

					//103. HP_Rr2_Hr2    + H2np --> HP_Rr2_Hr3
					case 96:
						HP_Rr2_Hr2--;
						H2np--;
						HP_Rr2_Hr3++;
						break;

					//105. HP_Rr2_Hr3    + Nnp  --> HP_RNa_Rr_Hr3
					case 97:
						HP_Rr2_Hr3--;
						Nnp--;
						HP_RNa_Rr_Hr3++;
						break;

					//106. HP_RNa2       + H2np --> HP_RNa2_Hr
					case 98:
						HP_RNa2--;
						H2np--;
						HP_RNa2_Hr++;
						break;

					//107. HP_RNa2_Hr    + H2np --> HP_RNa2_Hr2
					case 99:
						HP_RNa2_Hr--;
						H2np--;
						HP_RNa2_Hr2++;
						break;

					//108. HP_RNa2_Hr2   + H2np --> HP_RNa2_Hr3
					case 100:
						HP_RNa2_Hr2--;
						H2np--;
						HP_RNa2_Hr3++;
						break;

					//109. HP_Hr         --> HP_free       + H2np
					case 101:
						HP_Hr--;
						HP_free++;
						H2np++;
						break;

					//110. HP_Hr2        --> HP_Hr         + H2np
					case 102:
						HP_Hr2--;
						HP_Hr++;
						H2np++;
						break;

					//111. HP_Hr3        --> HP_Hr2        + H2np
					case 103:
						HP_Hr3--;
						HP_Hr2++;
						H2np++;
						break;

					//112. HP_Rr         --> HP_free       + Rnp
					case 104:
						HP_Rr--;
						HP_free++;
						Rnp++;
						break;

					//113. HP_Rr_Hr      --> HP_Hr         + Rnp
					case 105:
						HP_Rr_Hr--;
						HP_Hr++;
						Rnp++;
						break;

					//114. HP_Rr_Hr      --> HP_Rr         + H2np
					case 106:
						HP_Rr_Hr--;
						HP_Rr++;
						H2np++;
						break;

					//115. HP_Rr_Hr2     --> HP_Hr2        + Rnp
					case 107:
						HP_Rr_Hr2--;
						HP_Hr2++;
						Rnp++;
						break;

					//116. HP_Rr_Hr2     --> HP_Rr_Hr      + H2np
					case 108:
						HP_Rr_Hr2--;
						HP_Rr_Hr++;
						H2np++;
						break;

					//117. HP_Rr_Hr3     --> HP_Hr3        + Rnp
					case 109:
						HP_Rr_Hr3--;
						HP_Hr3++;
						Rnp++;
						break;

					//118. HP_Rr_Hr3     --> HP_Rr_Hr2     + H2np
					case 110:
						HP_Rr_Hr3--;
						HP_Rr_Hr2++;
						H2np++;
						break;

					//119. HP_RNa        --> HP_Rr         + Nnp
					case 111:
						HP_RNa--;
						HP_Rr++;
						Nnp++;
						break;

					//120. HP_RNa_Hr     --> HP_Rr_Hr      + Nnp
					case 112:
						HP_RNa_Hr--;
						HP_Rr_Hr++;
						Nnp++;
						break;

					//121. HP_RNa_Hr     --> HP_RNa        + H2np
					case 113:
						HP_RNa_Hr--;
						HP_RNa++;
						H2np++;
						break;

					//122. HP_RNa_Hr2    --> HP_Rr_Hr2     + Nnp
					case 114:
						HP_RNa_Hr2--;
						HP_Rr_Hr2++;
						Nnp++;
						break;

					//123. HP_RNa_Hr2    --> HP_RNa_Hr     + H2np
					case 115:
						HP_RNa_Hr2--;
						HP_RNa_Hr++;
						H2np++;
						break;

					//124. HP_RNa_Hr3    --> HP_Rr_Hr3     + Nnp
					case 116:
						HP_RNa_Hr3--;
						HP_Rr_Hr3++;
						Nnp++;
						break;

					//125. HP_RNa_Hr3    --> HP_RNa_Hr2    + H2np
					case 117:
						HP_RNa_Hr3--;
						HP_RNa_Hr2++;
						H2np++;
						break;

					//126. HP_RNa_Rr     --> HP_RNa        + Rnp
					case 118:
						HP_RNa_Rr--;
						HP_RNa++;
						Rnp++;
						break;

					//127. HP_RNa_Rr     --> HP_Rr2        + Nnp
					case 119:
						HP_RNa_Rr--;
						HP_Rr2++;
						Nnp++;
						break;

					//128. HP_RNa_Rr_Hr  --> HP_RNa_Hr     + Rnp
					case 120:
						HP_RNa_Rr_Hr--;
						HP_RNa_Hr++;
						Rnp++;
						break;

					//129. HP_RNa_Rr_Hr  --> HP_Rr2_Hr     + Nnp
					case 121:
						HP_RNa_Rr_Hr--;
						HP_Rr2_Hr++;
						Nnp++;
						break;

					//130. HP_RNa_Rr_Hr  --> HP_RNa_Rr     + H2np
					case 122:
						HP_RNa_Rr_Hr--;
						HP_RNa_Rr++;
						H2np++;
						break;

					//131. HP_RNa_Rr_Hr2 --> HP_RNa_Hr2    + Rnp
					case 123:
						HP_RNa_Rr_Hr2--;
						HP_RNa_Hr2++;
						Rnp++;
						break;

					//132. HP_RNa_Rr_Hr2 --> HP_Rr2_Hr2    + Nnp
					case 124:
						HP_RNa_Rr_Hr2--;
						HP_Rr2_Hr2++;
						Nnp++;
						break;

					//133. HP_RNa_Rr_Hr2 --> HP_RNa_Rr_Hr  + H2np
					case 125:
						HP_RNa_Rr_Hr2--;
						HP_RNa_Rr_Hr++;
						H2np++;
						break;

					//134. HP_RNa_Rr_Hr3 --> HP_RNa_Hr3    + Rnp
					case 126:
						HP_RNa_Rr_Hr3--;
						HP_RNa_Hr3++;
						Rnp++;
						break;

					//135. HP_RNa_Rr_Hr3 --> HP_Rr2_Hr3    + Nnp
					case 127:
						HP_RNa_Rr_Hr3--;
						HP_Rr2_Hr3++;
						Nnp++;
						break;

					//136. HP_RNa_Rr_Hr3 --> HP_RNa_Rr_Hr2 + H2np
					case 128:
						HP_RNa_Rr_Hr3--;
						HP_RNa_Rr_Hr2++;
						H2np++;
						break;

					//137. HP_Rr2        --> HP_Rr         + Rnp
					case 129:
						HP_Rr2--;
						HP_Rr++;
						Rnp++;
						break;

					//138. HP_Rr2_Hr     --> HP_Rr_Hr      + Rnp
					case 130:
						HP_Rr2_Hr--;
						HP_Rr_Hr++;
						Rnp++;
						break;

					//139. HP_Rr2_Hr     --> HP_Rr2        + H2np
					case 131:
						HP_Rr2_Hr--;
						HP_Rr2++;
						H2np++;
						break;

					//140. HP_Rr2_Hr2    --> HP_Rr_Hr2     + Rnp
					case 132:
						HP_Rr2_Hr2--;
						HP_Rr_Hr2++;
						Rnp++;
						break;

					//141. HP_Rr2_Hr2    --> HP_Rr2_Hr     + H2np
					case 133:
						HP_Rr2_Hr2--;
						HP_Rr2_Hr++;
						H2np++;
						break;

					//142. HP_Rr2_Hr3    --> HP_Rr_Hr3     + Rnp
					case 134:
						HP_Rr2_Hr3--;
						HP_Rr_Hr3++;
						Rnp++;
						break;

					//143. HP_Rr2_Hr3    --> HP_Rr2_Hr2    + H2np
					case 135:
						HP_Rr2_Hr3--;
						HP_Rr2_Hr2++;
						H2np++;
						break;

					//144. HP_RNa2       --> HP_RNa_Rr     + Nnp
					case 136:
						HP_RNa2--;
						HP_RNa_Rr++;
						Nnp++;
						break;

					//145. HP_RNa2_Hr    --> HP_RNa_Rr_Hr     + Nnp
					case 137:
						HP_RNa2_Hr--;
						HP_RNa_Rr_Hr++;
						Nnp++;
						break;

					//146. HP_RNa2_Hr    --> HP_RNa2       + H2np
					case 138:
						HP_RNa2_Hr--;
						HP_RNa2++;
						H2np++;
						break;

					//147. HP_RNa2_Hr2   --> HP_RNa_Rr_Hr2 + Nnp
					case 139:
						HP_RNa2_Hr2--;
						HP_RNa_Rr_Hr2++;
						Nnp++;
						break;

					//148. HP_RNa2_Hr2   --> HP_RNa2_Hr    + H2np
					case 140:
						HP_RNa2_Hr2--;
						HP_RNa2_Hr++;
						H2np++;
						break;

					//149. HP_RNa2_Hr3   --> HP_RNa_Rr_Hr3 + Nnp
					case 141:
						HP_RNa2_Hr3--;
						HP_RNa_Rr_Hr3++;
						Nnp++;
						break;

					//150. HP_RNa2_Hr3   --> HP_RNa2_Hr2   + H2np
					case 142:
						HP_RNa2_Hr3--;
						HP_RNa2_Hr2++;
						H2np++;
						break;

					//CBF1 promoter binding
					//151. RP_free + Rnp         --> RP_Rr
					case 143:
						RP_free--;
						Rnp--;
						RP_Rr++;
						break;

					//152. RP_free + H2np        --> RP_Hr
					case 144:
						RP_free--;
						H2np--;
						RP_Hr++;
						break;

					//153. RP_Hr + Rnp           --> RP_Rr_Hr
					case 145:
						RP_Hr--;
						Rnp--;
						RP_Rr_Hr++;
						break;

					//154. RP_Hr + H2np          --> RP_Hr2
					case 146:
						RP_Hr--;
						H2np--;
						RP_Hr2++;
						break;

					//155. RP_Hr2 + Rnp          --> RP_Rr_Hr2
					case 147:
						RP_Hr2--;
						Rnp--;
						RP_Rr_Hr2++;
						break;

					//156. RP_Hr2 + H2np         --> RP_Hr3
					case 148:
						RP_Hr2--;
						H2np--;
						RP_Hr3++;
						break;

					//157. RP_Hr3 + Rnp          --> RP_Rr_Hr3
					case 149:
						RP_Hr3--;
						Rnp--;
						RP_Rr_Hr3++;
						break;

					//158. RP_Rr + Rnp           --> RP_Rr2
					case 150:
						RP_Rr--;
						Rnp--;
						RP_Rr2++;
						break;

					//159. RP_Rr + Nnp           --> RP_RNa
					case 151:
						RP_Rr--;
						Nnp--;
						RP_RNa++;
						break;

					//160. RP_Rr + H2np          --> RP_Rr_Hr
					case 152:
						RP_Rr--;
						H2np--;
						RP_Rr_Hr++;
						break;

					//161. RP_Rr_Hr + Rnp        --> RP_Rr2_Hr
					case 153:
						RP_Rr_Hr--;
						Rnp--;
						RP_Rr2_Hr++;
						break;

					//162. RP_Rr_Hr + Nnp        --> RP_RNa_Hr
					case 154:
						RP_Rr_Hr--;
						Nnp--;
						RP_RNa_Hr++;
						break;

					//163. RP_Rr_Hr + H2np       --> RP_Rr_Hr2
					case 155:
						RP_Rr_Hr--;
						H2np--;
						RP_Rr_Hr2++;
						break;

					//164. RP_Rr_Hr2 + Rnp       --> RP_Rr2_Hr2
					case 156:
						RP_Rr_Hr2--;
						Rnp--;
						RP_Rr2_Hr2++;
						break;

					//165. RP_Rr_Hr2 + Nnp       --> RP_RNa_Hr2
					case 157:
						RP_Rr_Hr2--;
						Nnp--;
						RP_RNa_Hr2++;
						break;

					//166. RP_Rr_Hr2 + H2np      --> RP_Rr_Hr3
					case 158:
						RP_Rr_Hr2--;
						H2np--;
						RP_Rr_Hr3++;
						break;

					//167. RP_Rr_Hr3 + Rnp       --> RP_Rr2_Hr3
					case 159:
						RP_Rr_Hr3--;
						Rnp--;
						RP_Rr2_Hr3++;
						break;

					//168. RP_Rr_Hr3 + Nnp       --> RP_RNa_Hr3
					case 160:
						RP_Rr_Hr3--;
						Nnp--;
						RP_RNa_Hr3++;
						break;

					//169. RP_RNa + Rnp          --> RP_RNa_Rr
					case 161:
						RP_RNa--;
						Rnp--;
						RP_RNa_Rr++;
						break;

					//170. RP_RNa + H2np         --> RP_RNa_Hr
					case 162:
						RP_RNa--;
						H2np--;
						RP_RNa_Hr++;
						break;

					//171. RP_RNa_Hr + Rnp       --> RP_RNa_Rr_Hr
					case 163:
						RP_RNa_Hr--;
						Rnp--;
						RP_RNa_Rr_Hr++;
						break;

					//172. RP_RNa_Hr + H2np      --> RP_RNa_Hr2
					case 164:
						RP_RNa_Hr--;
						H2np--;
						RP_RNa_Hr2++;
						break;

					//173. RP_RNa_Hr2 + Rnp      --> RP_RNa_Rr_Hr2
					case 165:
						RP_RNa_Hr2--;
						Rnp--;
						RP_RNa_Rr_Hr2++;
						break;

					//174. RP_RNa_Hr2 + H2np     --> RP_RNa_Hr3
					case 166:
						RP_RNa_Hr2--;
						H2np--;
						RP_RNa_Hr3++;
						break;

					//175. RP_RNa_Hr3 + Rnp      --> RP_RNa_Rr_Hr3
					case 167:
						RP_RNa_Hr3--;
						Rnp--;
						RP_RNa_Rr_Hr3++;
						break;

					//176. RP_Rr2 + Rnp          --> RP_Rr3
					case 168:
						RP_Rr2--;
						Rnp--;
						RP_Rr3++;
						break;

					//177. RP_Rr2 + Nnp          --> RP_RNa_Rr
					case 169:
						RP_Rr2--;
						Nnp--;
						RP_RNa_Rr++;
						break;

					//178. RP_Rr2 + H2np         --> RP_Rr2_Hr
					case 170:
						RP_Rr2--;
						H2np--;
						RP_Rr2_Hr++;
						break;

					//179. RP_Rr2_Hr + Rnp       --> RP_Rr3_Hr
					case 171:
						RP_Rr2_Hr--;
						Rnp--;
						RP_Rr3_Hr++;
						break;

					//180. RP_Rr2_Hr + Nnp       --> RP_RNa_Rr_Hr
					case 172:
						RP_Rr2_Hr--;
						Nnp--;
						RP_RNa_Rr_Hr++;
						break;

					//181. RP_Rr2_Hr + H2np      --> RP_Rr2_Hr2
					case 173:
						RP_Rr2_Hr--;
						H2np--;
						RP_Rr2_Hr2++;
						break;

					//182. RP_Rr2_Hr2 + Rnp      --> RP_Rr3_Hr2
					case 174:
						RP_Rr2_Hr2--;
						Rnp--;
						RP_Rr3_Hr2++;
						break;

					//183. RP_Rr2_Hr2 + Nnp      --> RP_RNa_Rr_Hr2
					case 175:
						RP_Rr2_Hr2--;
						Nnp--;
						RP_RNa_Rr_Hr2++;
						break;

					//184. RP_Rr2_Hr2 + H2np     --> RP_Rr2_Hr3
					case 176:
						RP_Rr2_Hr2--;
						H2np--;
						RP_Rr2_Hr3++;
						break;

					//185. RP_Rr2_Hr3 + Rnp      --> RP_Rr3_Hr3
					case 177:
						RP_Rr2_Hr3--;
						Rnp--;
						RP_Rr3_Hr3++;
						break;

					//186. RP_Rr2_Hr3 + Nnp      --> RP_RNa_Rr_Hr3
					case 178:
						RP_Rr2_Hr3--;
						Nnp--;
						RP_RNa_Rr_Hr3++;
						break;

					//187. RP_RNa_Rr + Rnp       --> RP_RNa_Rr2
					case 179:
						RP_RNa_Rr--;
						Rnp--;
						RP_RNa_Rr2++;
						break;

					//188. RP_RNa_Rr + Nnp       --> RP_RNa2
					case 180:
						RP_RNa_Rr--;
						Nnp--;
						RP_RNa2++;
						break;

					//189. RP_RNa_Rr + H2np      --> RP_RNa_Rr_Hr
					case 181:
						RP_RNa_Rr--;
						H2np--;
						RP_RNa_Rr_Hr++;
						break;

					//190. RP_RNa_Rr_Hr + Rnp    --> RP_Rr2_Hr
					case 182:
						RP_RNa_Rr_Hr--;
						Rnp--;
						RP_Rr2_Hr++;
						break;

					//191. RP_RNa_Rr_Hr + Nnp    --> RP_RNa2_Hr
					case 183:
						RP_RNa_Rr_Hr--;
						Nnp--;
						RP_RNa2_Hr++;
						break;

					//192. RP_RNa_Rr_Hr + H2np   --> RP_RNa_Rr_Hr2
					case 184:
						RP_RNa_Rr_Hr--;
						H2np--;
						RP_RNa_Rr_Hr2++;
						break;

					//193. RP_RNa_Rr_Hr2 + Rnp   --> RP_RNa_Rr2_Hr2
					case 185:
						RP_RNa_Rr_Hr2--;
						Rnp--;
						RP_RNa_Rr2_Hr2++;
						break;

					//194. RP_RNa_Rr_Hr2 + Nnp   --> RP_RNa2_Hr2
					case 186:
						RP_RNa_Rr_Hr2--;
						Nnp--;
						RP_RNa2_Hr2++;
						break;

					//195. RP_RNa_Rr_Hr2 + H2np  --> RP_RNa_Rr_Hr3
					case 187:
						RP_RNa_Rr_Hr2--;
						H2np--;
						RP_RNa_Rr_Hr3++;
						break;

					//196. RP_RNa_Rr_Hr3 + Rnp   --> RP_RNa_Rr2_Hr3
					case 188:
						RP_RNa_Rr_Hr3--;
						Rnp--;
						RP_RNa_Rr2_Hr3++;
						break;

					//197. RP_RNa_Rr_Hr3 + Nnp   --> RP_RNa2_Hr3
					case 189:
						RP_RNa_Rr_Hr3--;
						Nnp--;
						RP_RNa2_Hr3++;
						break;

					//198. RP_Rr3 + Nnp          --> RP_RNa_Rr2
					case 190:
						RP_Rr3--;
						Nnp--;
						RP_RNa_Rr2++;
						break;

					//199. RP_Rr3 + H2np         --> RP_Rr3_Hr
					case 191:
						RP_Rr3--;
						H2np--;
						RP_Rr3_Hr++;
						break;

					//200. RP_Rr3_Hr + Nnp       --> RP_RNa_Rr2_Hr
					case 192:
						RP_Rr3_Hr--;
						Nnp--;
						RP_RNa_Rr2_Hr++;
						break;

					//201. RP_Rr3_Hr + H2np      --> RP_Rr3_Hr2
					case 193:
						RP_Rr3_Hr--;
						H2np--;
						RP_Rr3_Hr2++;
						break;

					//202. RP_Rr3_Hr2 + Nnp      --> RP_RNa_Rr2_Hr2
					case 194:
						RP_Rr3_Hr2--;
						Nnp--;
						RP_RNa_Rr2_Hr2++;
						break;

					//203. RP_Rr3_Hr2 + H2np     --> RP_Rr3_Hr3
					case 195:
						RP_Rr3_Hr2--;
							H2np--;
							RP_Rr3_Hr3++;
						break;

					//204. RP_Rr3_Hr3 + Nnp      --> RP_RNa_Rr2_Hr3
					case 196:
						RP_Rr3_Hr3--;
						Nnp--;
						RP_RNa_Rr2_Hr3++;
						break;

					//205. RP_RNa_Rr2 + Nnp      --> RP_RNa2_Rr
					case 197:
						RP_RNa_Rr2--;
						Nnp--;
						RP_RNa2_Rr++;
						break;

					//206. RP_RNa_Rr2 + H2np     --> RP_RNa_Rr2_Hr
					case 198:
						RP_RNa_Rr2--;
						H2np--;
						RP_RNa_Rr2_Hr++;
						break;

					//207. RP_RNa_Rr2_Hr + Nnp   --> RP_RNa2_Rr_Hr
					case 199:
						RP_RNa_Rr2_Hr--;
						Nnp--;
						RP_RNa2_Rr_Hr++;
						break;

					//208. RP_RNa_Rr2_Hr + H2np  --> RP_RNa_Rr2_Hr2
					case 200:
						RP_RNa_Rr2_Hr--;
						H2np--;
						RP_RNa_Rr2_Hr2++;
						break;

					//209. RP_RNa_Rr2_Hr2 + Nnp  --> RP_RNa2_Rr_Hr2
					case 201:
						RP_RNa_Rr2_Hr2--;
						Nnp--;
						RP_RNa_Rr_Hr2++;
						break;

					//210. RP_RNa_Rr2_Hr2 + H2np --> RP_RNa_Rr2_Hr3
					case 202:
						RP_RNa_Rr2_Hr2--;
						H2np--;
						RP_RNa_Rr2_Hr3++;
						break;

					//211. RP_RNa_Rr2_Hr3 + Nnp  --> RP_RNa2_Rr_Hr3
					case 203:
						RP_RNa_Rr2_Hr3--;
						Nnp--;
						RP_RNa2_Rr_Hr3++;
						break;

					//212. RP_RNa2 + Rnp         --> RP_RNa2_Rr
					case 204:
						RP_RNa2--;
						Rnp--;
						RP_RNa2_Rr++;
						break;

					//213. RP_RNa2 + H2np        --> RP_RNa2_Hr
					case 205:
						RP_RNa2--;
						H2np--;
						RP_RNa2_Hr++;
						break;

					//214. RP_RNa2_Hr + Rnp      --> RP_RNa2_Rr_Hr
					case 206:
						RP_RNa2_Hr--;
						Rnp--;
						RP_RNa2_Rr_Hr++;
						break;

					//215. RP_RNa2_Hr + H2np     --> RP_RNa2_Hr2
					case 207:
						RP_RNa2_Hr--;
						H2np--;
						RP_RNa2_Hr2++;
						break;

					//216. RP_RNa2_Hr2 + Rnp     --> RP_RNa2_Rr_Hr2
					case 208:
						RP_RNa2_Hr2--;
						Rnp--;
						RP_RNa2_Rr_Hr2++;
						break;

					//217. RP_RNa2_Hr2 + H2np    --> RP_RNa2_Hr3
					case 209:
						RP_RNa2_Hr2--;
						H2np--;
						RP_RNa2_Hr3++;
						break;

					//218. RP_RNa2_Hr3 + Rnp     --> RP_RNa2_Rr_Hr3
					case 210:
						RP_RNa2_Hr3--;
						Rnp--;
						RP_RNa2_Rr_Hr3++;
						break;

					//219. RP_RNa2_Rr + Nnp      --> RP_RNa3
					case 211:
						RP_RNa2_Rr--;
						Nnp--;
						RP_RNa3++;
						break;

					//220. RP_RNa2_Rr + H2np     --> RP_RNa2_Rr_Hr
					case 212:
						RP_RNa2_Rr--;
						H2np--;
						RP_RNa2_Rr_Hr++;
						break;

					//221. RP_RNa2_Rr_Hr + Nnp   --> RP_RNa3_Hr
					case 213:
						RP_RNa2_Rr_Hr--;
						Nnp--;
						RP_RNa3_Hr++;
						break;

					//222. RP_RNa2_Rr_Hr + H2np  --> RP_RNa2_Rr_Hr2
					case 214:
						RP_RNa2_Rr_Hr--;
						H2np--;
						RP_RNa2_Rr_Hr2++;
						break;

					//223. RP_RNa2_Rr_Hr2 + Nnp  --> RP_RNa3_Hr2
					case 215:
						RP_RNa2_Rr_Hr2--;
						Nnp--;
						RP_RNa3_Hr3++;
						break;

					//224. RP_RNa2_Rr_Hr2 + H2np --> RP_RNa2_Rr_Hr3
					case 216:
						RP_RNa2_Rr_Hr2--;
						H2np--;
						RP_RNa2_Rr_Hr3++;
						break;

					//225. RP_RNa2_Rr_Hr3 + Nnp  --> RP_RNa3_Hr3
					case 217:
						RP_RNa2_Rr_Hr3--;
						Nnp--;
						RP_RNa3_Hr3++;
						break;

					//226. RP_RNa3 + H2np        --> RP_RNa3_Hr
					case 218:
						RP_RNa3--;
						H2np--;
						RP_RNa3_Hr++;
						break;

					//227. RP_RNa3_Hr + H2np     --> RP_RNa3_Hr2
					case 219:
						RP_RNa3_Hr--;
						H2np--;
						RP_RNa3_Hr2++;
						break;

					//228. RP_RNa3_Hr2 + H2np    --> RP_RNa3_Hr3
					case 220:
						RP_RNa3_Hr2--;
						H2np--;
						RP_RNa3_Hr3++;
						break;

					//CBF1 promoter disassociation
					//229. RP_Hr			--> RP_free + H2np
					case 221:
						RP_Hr--;
						RP_free++;
						H2np++;
						break;

					//230. RP_Hr2			--> RP_Hr + H2np
					case 222:
						RP_Hr2--;
						RP_Hr++;
						H2np++;
						break;

					//231. RP_Hr3			--> RP_Hr2 + H2np
					case 223:
						RP_Hr3--;
						RP_Hr2++;
						H2np++;
						break;

					//232. RP_Rr			--> RP_free + Rnp
					case 224:
						RP_Rr--;
						RP_free++;
						Rnp++;
						break;

					//233. RP_Rr_Hr			--> RP_Hr + Rnp
					case 225:
						RP_Rr_Hr--;
						RP_Hr++;
						Rnp++;
						break;

					//234. RP_Rr_Hr			--> RP_Rr + H2np
					case 226:
						RP_Rr_Hr--;
						RP_Rr++;
						H2np++;
						break;

					//235. RP_Rr_Hr2		--> RP_Hr2 + Rnp
					case 227:
						RP_Rr_Hr2--;
						RP_Hr2++;
						Rnp++;
						break;

					//236. RP_Rr_Hr2		--> RP_Rr_Hr + H2np
					case 228:
						RP_Rr_Hr2--;
						RP_Rr_Hr++;
						H2np++;
						break;

					//237. RP_Rr_Hr3		--> RP_Hr3 + Rnp
					case 229:
						RP_Rr_Hr3--;
						RP_Hr3++;
						Rnp++;
						break;

					//238. RP_Rr_Hr3		--> RP_Rr_Hr2 + H2np
					case 230:
						RP_Rr_Hr3--;
						RP_Rr_Hr2++;
						H2np++;
						break;

					//239. RP_RNa			--> RP_Rr + Nnp
					case 231:
						RP_RNa--;
						RP_Rr++;
						Nnp++;
						break;

					//240. RP_RNa_Hr		--> RP_Rr_Hr + Nnp
					case 232:
						RP_RNa_Hr--;
						RP_Rr_Hr++;
						Nnp++;
						break;

					//241. RP_RNa_Hr		--> RP_RNa + H2np
					case 233:
						RP_RNa_Hr--;
						RP_RNa++;
						H2np++;
						break;

					//242. RP_RNa_Hr2		--> RP_Rr_Hr2 + Nnp
					case 234:
						RP_RNa_Hr2--;
						RP_Rr_Hr2++;
						Nnp++;
						break;

					//243. RP_RNa_Hr2		--> RP_RNa_Hr + H2np
					case 235:
						RP_RNa_Hr2--;
						RP_RNa_Hr++;
						H2np++;
						break;

					//244. RP_RNa_Hr3		--> RP_Rr_Hr3 + Nnp
					case 236:
						RP_RNa_Hr3--;
						RP_Rr_Hr3++;
						Nnp++;
						break;

					//245. RP_RNa_Hr3		--> RP_RNa_Hr2 + H2np
					case 237:
						RP_RNa_Hr3--;
						RP_RNa_Hr2++;
						H2np++;
						break;

					//246. RP_Rr2			--> RP_Rr + Rnp
					case 238:
						RP_Rr2--;
						RP_Rr++;
						Rnp++;
						break;

					//247. RP_Rr2_Hr		--> RP_Rr_Hr + Rnp
					case 239:
						RP_Rr_Hr--;
						RP_Rr_Hr++;
						Rnp++;
						break;

					//248. RP_Rr2_Hr		--> RP_Rr2 + H2np
					case 240:
						RP_Rr2_Hr--;
						RP_Rr2++;
						H2np++;
						break;

					//249. RP_Rr2_Hr2		--> RP_Rr_Hr2 + Rnp
					case 241:
						RP_Rr2_Hr2--;
						RP_Rr_Hr2++;
						Rnp++;
						break;

					//250. RP_Rr2_Hr2		--> RP_Rr2_Hr + H2np
					case 242:
						RP_Rr2_Hr2--;
						RP_Rr2_Hr++;
						H2np++;
						break;

					//251. RP_Rr2_Hr3		--> RP_Rr_Hr3 + Rnp
					case 243:
						RP_Rr2_Hr3--;
						RP_Rr_Hr3++;
						Rnp++;
						break;

					//252. RP_Rr2_Hr3		--> RP_Rr2_Hr2 + H2np
					case 244:
						RP_Rr2_Hr3--;
						RP_Rr2_Hr2++;
						H2np++;
						break;

					//253. RP_RNa_Rr		--> RP_RNa + Rnp
					case 245:
						RP_RNa_Rr--;
						RP_RNa++;
						Nnp++;
						break;

					//254. RP_RNa_Rr		--> RP_Rr2 + Nnp
					case 246:
						RP_RNa_Rr--;
						RP_Rr2++;
						Nnp++;
						break;

					//255. RP_RNa_Rr_Hr		--> RP_RNa_Hr + Rnp
					case 247:
						RP_RNa_Rr_Hr--;
						RP_RNa_Hr++;
						Rnp++;
						break;

					//256. RP_RNa_Rr_Hr		--> RP_Rr2_Hr + Nnp
					case 248:
						RP_RNa_Rr_Hr--;
						RP_Rr2_Hr++;
						Nnp++;
						break;

					//257. RP_RNa_Rr_Hr		--> RP_RNa_Rr + H2np
					case 249:
						RP_RNa_Rr_Hr--;
						RP_RNa_Rr++;
						H2np++;
						break;

					//258. RP_RNa_Rr_Hr2	--> RP_RNa_Hr2 + Rnp
					case 250:
						RP_RNa_Rr_Hr2--;
						RP_RNa_Hr2++;
						Rnp++;
						break;

					//259. RP_RNa_Rr_Hr2	--> RP_Rr2_Hr2 + Nnp
					case 251:
						RP_RNa_Rr_Hr2--;
						RP_Rr2_Hr2++;
						Nnp++;
						break;

					//260. RP_RNa_Rr_Hr2	--> RP_RNa_Rr_Hr+ H2np
					case 252:
						RP_RNa_Rr_Hr2--;
						RP_RNa_Rr_Hr++;
						H2np++;
						break;

					//261. RP_RNa_Rr_Hr3	--> RP_RNa_Hr3 + Rnp
					case 253:
						RP_RNa_Rr_Hr3--;
						RP_RNa_Hr3++;
						Rnp++;
						break;

					//262. RP_RNa_Rr_Hr3	--> RP_Rr2_Hr3 + Nnp
					case 254:
						RP_RNa_Rr_Hr3--;
						RP_Rr2_Hr3++;
						Nnp++;
						break;

					//263. RP_RNa_Rr_Hr3	--> RP_RNa_Rr_Hr2 + H2np
					case 255:
						RP_RNa_Rr_Hr3--;
						RP_RNa_Rr_Hr2++;
						H2np++;
						break;

					//264. RP_Rr3			--> RP_Rr2 + Rnp
					case 256:
						RP_Rr3--;
						RP_Rr2++;
						Rnp++;
						break;

					//265. RP_Rr3_Hr		--> RP_Rr2_Hr + Rnp
					case 257:
						RP_Rr3_Hr--;
						RP_Rr2_Hr++;
						Rnp++;
						break;

					//266. RP_Rr3_Hr		--> RP_Rr3 + H2np
					case 258:
						RP_Rr3_Hr--;
						RP_Rr3++;
						H2np++;
						break;

					//267. RP_Rr3_Hr2		--> RP_Rr2_Hr2 + Rnp
					case 259:
						RP_Rr3_Hr2--;
						RP_Rr2_Hr2++;
						Rnp++;
						break;

					//268. RP_Rr3_Hr2		--> RP_Rr3_Hr + H2np
					case 260:
						RP_Rr3_Hr2--;
						RP_Rr3_Hr++;
						H2np++;
						break;

					//269. RP_Rr3_Hr3		--> RP_Rr2_Hr3 + Rnp
					case 261:
						RP_Rr3_Hr3--;
						RP_Rr2_Hr3++;
						Rnp++;
						break;

					//270. RP_Rr3_Hr3		--> RP_Rr3_Hr2 + H2np
					case 262:
						RP_Rr3_Hr3--;
						RP_Rr3_Hr2++;
						H2np++;
						break;

					//271. RP_RNa_Rr2		--> RP_RNa_Rr + Rnp
					case 263:
						RP_RNa_Rr2--;
						RP_RNa_Rr++;
						Rnp++;
						break;

					//272. RP_RNa_Rr2		--> RP_Rr3 + Nnp
					case 264:
						RP_RNa_Rr2--;
						RP_Rr3++;
						Nnp++;
						break;

					//273. RP_RNa_Rr2_Hr	--> RP_RNa_Rr_Hr + Rnp
					case 265:
						RP_RNa_Rr2_Hr--;
						RP_RNa_Rr_Hr++;
						Rnp++;
						break;

					//274. RP_RNa_Rr2_Hr	--> RP_Rr3_Hr + Nnp
					case 266:
						RP_RNa_Rr2_Hr--;
						RP_Rr3_Hr++;
						Nnp++;
						break;

					//275. RP_RNa_Rr2_Hr	--> RP_RNa_Rr2 + H2np
					case 267:
						RP_RNa_Rr2_Hr--;
						RP_RNa_Rr2++;
						H2np++;
						break;

					//276. RP_RNa_Rr2_Hr2	--> RP_RNa_Rr_Hr2 + Rnp
					case 268:
						RP_RNa_Rr2_Hr2--;
						RP_RNa_Rr_Hr2++;
						Rnp++;
						break;

					//277. RP_RNa_Rr2_Hr2	--> RP_Rr3_Hr2 + Nnp
					case 269:
						RP_RNa_Rr2_Hr2--;
						RP_Rr3_Hr2++;
						Nnp++;
						break;

					//278. RP_RNa_Rr2_Hr2	--> RP_RNa_Rr2_Hr + H2np
					case 270:
						RP_RNa_Rr2_Hr2--;
						RP_RNa_Rr2_Hr++;
						H2np++;
						break;

					//279. RP_RNa_Rr2_Hr3	--> RP_RNa_Rr_Hr3 + Rnp
					case 271:
						RP_RNa_Rr2_Hr3--;
						RP_RNa_Rr_Hr3++;
						Rnp++;
						break;

					//280. RP_RNa_Rr2_Hr3	--> RP_Rr3_Hr3 + Nnp
					case 272:
						RP_RNa_Rr2_Hr3--;
						RP_Rr3_Hr3++;
						Nnp++;
						break;

					//281. RP_RNa_Rr2_Hr3	--> RP_RNa_Rr2_Hr2 + H2np
					case 273:
						RP_RNa_Rr2_Hr3--;
						RP_RNa_Rr2_Hr2++;
						H2np++;
						break;

					//282. RP_RNa2			--> RP_RNa_Rr + Nnp
					case 274:
						RP_RNa2--;
						RP_RNa_Rr++;
						Nnp++;
						break;

					//283. RP_RNa2_Hr		--> RP_RNa_Rr_Hr + Nnp
					case 275:
						RP_RNa_Hr--;
						RP_RNa_Rr_Hr++;
						Nnp++;
						break;

					//284. RP_RNa2_Hr		--> RP_RNa2 + H2np
					case 276:
						RP_RNa2_Hr--;
						RP_RNa2++;
						H2np++;
						break;

					//285. RP_RNa2_Hr2		--> RP_RNa_Rr_Hr2 + Nnp
					case 277:
						RP_RNa2_Hr2--;
						RP_RNa_Rr_Hr2++;
						Nnp++;
						break;

					//286. RP_RNa2_Hr2		--> RP_RNa2_Hr + H2np
					case 278:
						RP_RNa2_Hr2--;
						RP_RNa2_Hr++;
						H2np++;
						break;

					//287. RP_RNa2_Hr3		--> RP_RNa_Rr_Hr3 + Nnp
					case 279:
						RP_RNa2_Hr3--;
						RP_RNa_Rr_Hr3++;
						Nnp++;
						break;

					//288. RP_RNa2_Hr3		--> RP_RNa2_Hr2 + H2np
					case 280:
						RP_RNa2_Hr3--;
						RP_RNa2_Hr2++;
						H2np++;
						break;

					//289. RP_RNa2_Rr 		--> RP_RNa2 + Rnp
					case 281:
						RP_RNa2_Rr--;
						RP_RNa2++;
						Rnp++;
						break;

					//290. RP_RNa2_Rr 		--> RP_RNa_Rr2 + Nnp
					case 282:
						RP_RNa2_Rr--;
						RP_RNa_Rr2++;
						Nnp++;
						break;

					//291. RP_RNa2_Rr_Hr 	--> RP_RNa2_Hr + Rnp
					case 283:
						RP_RNa2_Rr_Hr--;
						RP_RNa2_Hr++;
						Rnp++;
						break;

					//292. RP_RNa2_Rr_Hr 	--> RP_RNa_Rr2_Hr + Nnp
					case 284:
						RP_RNa2_Rr_Hr--;
						RP_RNa_Rr2_Hr++;
						Nnp++;
						break;

					//293. RP_RNa2_Rr_Hr 	--> RP_RNa2_Rr + H2np
					case 285:
						RP_RNa2_Rr_Hr--;
						RP_RNa2_Rr++;
						H2np++;
						break;

					//294. RP_RNa2_Rr_Hr2	--> RP_RNa2_Hr + Rnp
					case 286:
						RP_RNa2_Rr_Hr2--;
						RP_RNa2_Hr++;
						Rnp++;
						break;

					//295. RP_RNa2_Rr_Hr2	--> RP_RNa_Rr2_Hr2 + Nnp
					case 287:
						RP_RNa2_Rr_Hr2--;
						RP_RNa_Rr2_Hr2++;
						Nnp++;
						break;

					//296. RP_RNa2_Rr_Hr2	--> RP_RNa2_Rr_Hr + H2np
					case 288:
						RP_RNa2_Rr_Hr2--;
						RP_RNa2_Rr_Hr++;
						H2np++;
						break;

					//297. RP_RNa2_Rr_Hr3 	--> RP_RNa2_Hr3 + Rnp
					case 289:
						RP_RNa2_Rr_Hr3--;
						RP_RNa2_Hr3++;
						Rnp++;
						break;

					//298. RP_RNa2_Rr_Hr3 	--> RP_RNa_Rr2_Hr3 + Nnp
					case 290:
						RP_RNa2_Rr_Hr3--;
						RP_RNa_Rr2_Hr3++;
						Nnp++;
						break;

					//299. RP_RNa2_Rr_Hr3 	--> RP_RNa2_Rr_Hr2 + H2np
					case 291:
						RP_RNa2_Rr_Hr3--;
						RP_RNa2_Rr_Hr2++;
						H2np++;
						break;

					//300. RP_RNa3 			--> RP_RNa2_Rr + Nnp
					case 292:
						RP_RNa3--;
						RP_RNa2_Rr++;
						Nnp++;
						break;

					//301. RP_RNa3_Hr 		--> RP_RNa2_Rr_Hr + Nnp
					case 293:
						RP_RNa3_Hr--;
						RP_RNa2_Rr_Hr++;
						Nnp++;
						break;

					//302. RP_RNa3_Hr 		--> RP_RNa3 + H2np
					case 294:
						RP_RNa3_Hr--;
						RP_RNa3++;
						H2np++;
						break;

					//303. RP_RNa3_Hr2 		--> RP_RNa2_Rr_Hr2 + Nnp
					case 295:
						RP_RNa3_Hr2--;
						RP_RNa2_Rr_Hr2++;
						Nnp++;
						break;

					//304. RP_RNa3_Hr2 		--> RP_RNa3_Hr + H2np
					case 296:
						RP_RNa3_Hr2--;
						RP_RNa3_Hr++;
						H2np++;
						break;

					//305. RP_RNa3_Hr3 		--> RP_RNa2_Rr_Hr3 + Nnp
					case 297:
						RP_RNa3_Hr3--;
						RP_RNa2_Rr_Hr3++;
						Nnp++;
						break;

					//306. RP_RNa3_Hr3 		--> RP_RNa3_Hr2 + H2np
					case 298:
						RP_RNa3_Hr3--;
						RP_RNa3_Hr2++;
						H2np++;
						break;

				}
			}

			if(t>tp) {
				cout<<tp<<"\t"<<Rcm<<"\t"<<Hcm<<"\t"<<Nm<<"\t"<<Rnp<<"\t"<<Hnp<<"\t"<<Np<<"\t"<<Hnp<<"\t"<<endl;
				tp++;
			}

			//delta pulse at t=750
			if(t>750 && t<760) {
				delta=10000;
			}
		}




	/*	for(i=0;i<species;i++){
			cout<<"\t"<<X[i];
		}
				cout<<endl;
		*/
		//cout<<n<<endl;

	}

	return 0;

}
