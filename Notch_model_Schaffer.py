from pysb.core import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from pysb.integrate import ScipyOdeSimulator as SOS
from pysb.bng import generate_equations
# "??" denotes issues to be addressed

Model()

## Omitted Species for Schaffer model
#Gamma secretase - bound to ligand-receptor
#Monomer('YSC', ['lr'])
#Notch membrane domain - omitted in Schaffer model
#Monomer('NMD', ['lr'])
#Monomer('ND_endo')
#Monomer('N_endo')
#Monomer('MAM', ['nicd'])
# Monomer('pASCL1')
# Monomer('pNEUROD1')
# Monomer('pSLUG')
# Monomer('pNGN2')
## end section

## Monomer section
Monomer('Delta')
#binding sites for ligand; got rid of gamma secretase site ??
Monomer('mNotch')
Monomer('pNotch') #should the gamma secretase be added back on??
Monomer('NICD', ['rbpj', 'loc'], {'loc' : ['cyt', 'nuc']}) # binds only to RBPJ
Monomer('mRBPJ')
Monomer('pRBPJ', ['nicd', 'Rbox', 'loc'], {'loc' : ['cyt', 'nuc']}) #name I made up for RBPJ-binding site in promoter region
#bad version Monomer('pRBPJ', ['nicd', 'Hprom', 'Nprom', 'Rprom'])
Monomer('mHes1')
Monomer('pHes1', ['pHes1','Hbox', 'loc'], {'loc' : ['cyt', 'nuc']}) # Should it be a dimerizing before binding to Nbox ??
#bad version Monomer('pHES1', ['pHes1','Hprom', 'Nprom', 'Rprom'])
#Promoter Regions (d for DNA)
Monomer('dNotch', ['Hbox1a', 'Hbox1b', 'Rbox1', 'Rbox2'])
Monomer('dHes1', ['Hbox1a', 'Hbox1b', 'Hbox2a', 'Hbox2b', 'Hbox3a', 'Hbox3b', 'Rbox1', 'Rbox2'])
Monomer('dRBPJ', ['Hbox1a', 'Hbox1b', 'Hbox2a', 'Hbox2b', 'Hbox3a', 'Hbox3b', 'Rbox1', 'Rbox2', 'Rbox3'])

Observable('pHes1_cyt', pHes1(loc='cyt'))

# a[62] = c[24] * 2 * HP_Hr * H2np; // 63.  HP_Hr + H2np --> HP_Hr2
# versus
#a[64]=c[24] * HP_Hr2 * H2np;			//65. HP_Hr2        + H2np --> HP_Hr3
#

# print model.monomers

fill = 1

#Default forward, reverse, and catalytic rates
KF = 1e-6
KR = 1e-3
KC = 1

# cell_size = ??
# nucleus_size = ??
# membrane_size = ??

Parameter('Delta_0', 0)
# Parameter('NICD_cyt_0', 0)
# Parameter('NICD_nuc_0', 0)
# get rid of ?? Parameter('NMD_0', fill)
Parameter('mHes1_0', 1) 			# Hcm=1;			//3.64x10^(-12)
Parameter('pHes1_cyt_0', 35) 	# Hcp=35;		//1.25*10^(-10)
Parameter('pHes1_nuc_0', 109) 	# Hnp=109;		//2.77*10^(-9)
Parameter('dimHes1_nuc_0', 101) 	# H2np=101;		//2.5662*10^(-9)
Parameter('mRBPJ_0', 1) 			# Rcm=1;			//1.2*10^(-12)
Parameter('pRBPJ_cyt_0', 11) 	# Rcp=11;		//3.74*10^(-11)
Parameter('pRBPJ_nuc_0', 447) 	# Rnp=447;		//1.13*10^(-8)
Parameter('mNotch_0', 70) 		# Nm=70;			//2.51*10^(-10)	//Notch mRNA
Parameter('pNotch_0', 8130) 		#??# Np=8130;	//2.95*10^(-8)	//Notch protein in the cytoplasm
Parameter('dNotch_0', 2) 		# NP_free=2;		//Notch promoter (diploid so start with 2 free promoter regions)
Parameter('dHes1_0', 2) 			# HP_free=2;		//Hes1 promoter (diploid so start with 2 free promoter regions)
Parameter('dRBPJ_0', 2) 			# RP_free=2;		//RPJK promoter (diploid so start with 2 free promoter regions)

#How much should I specify in the initials ??
#do i put location or no ??
Initial(Delta, Delta_0)
# Initial(NICD(rbpj=None, loc=cyt), NICD_cyt_0)
# Initial(NICD(rbpj=None, loc=nuc), NICD_nuc_0)
#Initial(NMD(), NMD_0)
Initial(mHes1(), mHes1_0)
Initial(pHes1(pHes1=None, Hbox=None, loc='cyt'), pHes1_cyt_0)
Initial(pHes1(pHes1=None, Hbox=None, loc='nuc'), pHes1_nuc_0)
Initial(pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox= None, loc='nuc'), dimHes1_nuc_0)
Initial(mRBPJ(), mRBPJ_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc='cyt'), pRBPJ_cyt_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc='nuc'), pRBPJ_nuc_0)
Initial(mNotch(), mNotch_0)
Initial(pNotch(), pNotch_0)
Initial(dNotch(Hbox1a=None, Hbox1b=None, Rbox1=None, Rbox2=None), dNotch_0)
Initial(dHes1(Hbox1a=None, Hbox1b=None, Hbox2a=None, Hbox2b=None, Hbox3a=None, Hbox3b=None, Rbox1=None, Rbox2=None), dHes1_0)
Initial(dRBPJ(Hbox1a=None, Hbox1b=None, Hbox2a=None, Hbox2b=None, Hbox3a=None, Hbox3b=None, Rbox1=None, Rbox2=None, Rbox3=None), dRBPJ_0)

##
#Add section for times delays ??
###

delay=True	#//set true here if you want to implement delays

#//Delay times for the 6 delayed reactions
#//units here are minutes delayed
TpNc = 21
TpRc = 4.3
TpHc = 2.35
TmNc = 70
TmRc = 20
TmHc = 10

#Degradation parameters
#	//Units here are per minute
Parameter('KdNm', 0.00035) 		#degradation Notch mRNA
Parameter('KdNp', 0.017) 		#degradation Notch protein

Parameter('KdHcm', 0.029) 		#// deg Hes1 mRNA
# Parameter('KdHcp', 0.0315) 	#//degradation Hes1 protein
# Parameter('KdHnp', 0.0315) 	#//
Parameter('KdHes1', 0.0315)		#//degradation Hes1 protein
Parameter('KdH2np', 0) 			#//disable Hes1 dimer degradation <-- I didn't write that -MS

# Parameter('KdNcp', 0.00385) 	#degradation NICD protein
# Parameter('KdNnp', 0.00385) 	#""
# Parameter('KdNICD', 0.00385) 	#degradation NICD protein
Expression('KdNICD', 0.00385 + (0.0014-0.00385)*(pHes1_cyt > 100)) #degradation NICD protein

Parameter('KdRcm', 0.0075) 		#degradation RBPJ mRNA
# Parameter('KdRcp', 0.00231) 	#degradation RBPJ protein
# Parameter('KdRnp', 0.00231) 	# ""
Parameter('KdRBPJ', 0.00231) 	#degradation RBPJ protein

#	//units here are per minute
#protein translation rates
Parameter('KtrHc', 4.5)
Parameter('KtrN', 2)
Parameter('KtrRc', 3.2)

# units here are per minute
#nuclear rate of import
Parameter('KniRcp', 0.1)
Parameter('KniHcp', 0.1)
Parameter('KniNcp', 0.1)

#//7.6*10^8 per molar per minute, forward rate of NICD binding to Delta
Parameter('KfNcp', 0.000276)
# Hes1 dimer association constant
Parameter('KaHp', 0.02535)
# // disable disassociation of Hes1 dimer
Parameter('KrHp', 0)
#RBPJ-DNA association constant; //Keqr=8.19e-3; //float Kr=0.001628;
Parameter('Kr', 0.0410)
Parameter('Kr_r', 5)
#Kn = Hes1-DNA association constant; // Keqn = 5.07e-3; // float Kn = 0.00072;
Parameter('Kn', 0.0254)
Parameter('Kn_r', 5)
#Ka = RBPJ-NICD association constant; // Keqa = 2.54e-3;  // float Ka = 0.00036;
Parameter('Ka', 0.0127)
Parameter('Ka_r', 5)

#//These generation rates are based on molecules generated per cell per minute <= from original code
Parameter('Vmaxh', 197) #//9.98;
Parameter('Vmaxr', 79) #//4.99
Parameter('Vmaxn', 21.6) #//1.43
Parameter('Vbh', 4.5) #//2.53*10^(-12)M/min
Parameter('Vbr', 1.7) #//1.27*10^(-12)M/min
Parameter('Vbn', 0.5) #//3.6*10^(-13)M/min
Parameter('rNbox', 0.3) #
Parameter('rRbox', 0.2)
Parameter('tc', 0.5)

# //Generation reactions
# a[0]=c[0];		//1. 0 --> Nm		0.RfNm
# a[1]=c[1];		//2. 0 --> Rcm		1.RfRm
# a[2]=c[2];		//3. 0 --> Hcm		2.RfHm

# .... Add these in last --LAH

# //Translation channels (delay these if they occur in delay model)
# a[3]=c[3]*Nm;		//4. Nm --> Np		3.KtrN
# a[4]=c[4]*Rcm;		//5. Rcm --> Rcp		4.KtrRc
# a[5]=c[5]*Hcm;		//6. Hcm --> Hcp		5.KtrHc

Rule('mNotch_translation', mNotch() >> mNotch() + pNotch() , KtrN)
Rule('mRBPJ_translation', mRBPJ() >> mRBPJ() + pRBPJ(nicd=None, Rbox=None, loc='cyt'), KtrRc)
Rule('mHes1_translation', mHes1() >> mHes1() + pHes1(pHes1=None, Hbox=None, loc='cyt'), KtrHc)
 
# //Other reactions in the cytoplasm
# //a[6]=c[6]*Hnp*Hnp;	//7. Hnp + Hnp --> H2np		6.KaHp
# //a[7]=c[7]*H2np;		//8. H2np --> Hnp + Hnp		7.KrHp
# a[8]=c[8]*delta*Np;	//9. delta + Np --> Ncp		8.KfNcp

Rule('Delta_Notch_to_NICD', Delta() + pNotch() >> NICD(rbpj=None, loc='cyt'), KfNcp)

# // Effect of GSK3b on NICD half life
# 
# if(Hcp>100)
# {
#   c[12]=c[15]=0.0014;   // KdNnp = KdNcp = 0.0014 - 3hr half life
# }
# else
# {
#   c[12]=c[15]=0.00385;   // KdNnp = KdNcp = 0.00385 - 8hr half life
# }
# 
# //Degeneration channels
# a[9]=c[9]*Nm;			//10. Nm --> 0		9.KdNm
# a[10]=c[10]*Rcm;		//11. Rcm --> 0		10.KdRcm
# a[11]=c[11]*Hcm;		//12. Hcm --> 0		11.KdHcm
# a[12]=c[12]*Nnp;		//13. Nnp --> 0		12.KdNnp
# a[13]=c[13]*Rnp;		//14. Rnp --> 0		13.KdRnp
# a[14]=c[14]*Np;		//15. Np --> 0		14.KdNp
# a[15]=c[15]*Ncp;		//16. Ncp --> 0		15.KdNcp = KdNnp
# a[16]=c[16]*Rcp;		//17. Rcp --> 0		16.KdRcp = KdRnp
# a[17]=c[17]*Hcp;		//18. Hcp --> 0		17.KdHcp = KdHnp
# a[18]=KdHnp*Hnp;		//19. Hnp --> 0		xx.KdHnp
# //a[19]=0;				//20. H2np --> 0		18.KdH2np

Rule('mNotch_degradation', mNotch() >> None, KdNm) # //10. Nm --> 0  KdNm
Rule('pNotch_degradation', pNotch() >> None, KdNp) # //15. Np --> 0  KdNp

Rule('mHes1_degradation', mHes1() >> None, KdHcm) # //12. Hcm --> 0  KdHcm
Rule('pHes1_monomer_degradation', pHes1(pHes1=None, Hbox=None) >> None, KdHes1) # //18. Hcp --> 0  KdHcp; //19. Hnp --> 0  KdHnp
Rule('pHes1_dimer_degradation', pHes1(pHes1=ANY, Hbox=None, loc='nuc') >> None, KdH2np) # //20. H2np --> 0  KdH2np

Rule('NICD_degradation', NICD(rbpj=None) >> None, KdNICD) # //13. Nnp --> 0  KdNnp; //16. Ncp --> 0  KdNcp 

Rule('mRBPJ_degradation', mRBPJ() >> None, KdRcm) # //11. Rcm --> 0  KdRcm
Rule('pRBPJ_degradation', pRBPJ(nicd=None, Rbox=None) >> None, KdRBPJ) # //14. Rnp --> 0  KdRnp; //17. Rcp --> 0  KdRcp

print
print model.rules
	