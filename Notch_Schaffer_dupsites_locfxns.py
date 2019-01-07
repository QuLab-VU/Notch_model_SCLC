from pysb import *
from pysb.simulator.scipyode import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

# delta: Delta protein (stimulus)

# HP_free: Hes1 promoter
# Hcm: Hes1 mRNA
# Hcp: Hes1 cytoplasmic protein
# Hnp: Hes1 nuclear protein
# H2np: Hes1 nuclear dimer

# RP_free: RBPJ promoter
# Rcm: RBPJ mRNA
# Rcp: RBPJ cytoplasmic protein
# Rnp: RBPJ nuclear protein

# NP_free: Notch promoter
# Nm: Notch1 mRNA
# Np: Notch1 protein

# Ncp: NICD cytoplasmic protein
# Nnp: NICD nuclear protein

Vcyt = 4.581e-13 # liter
Vnuc = 6.545e-14 # liter

Model()

Monomer('Delta')
Monomer('mNotch')
Monomer('pNotch') 
Monomer('NICD', ['rbpj', 'loc'], {'loc' : ['cyt', 'nuc']}) 
Monomer('mRBPJ')
Monomer('pRBPJ', ['nicd', 'Rbox', 'loc'], {'loc' : ['cyt', 'nuc']})
Monomer('mHes1')
Monomer('pHes1', ['pHes1', 'Hbox', 'loc'], {'loc' : ['cyt', 'nuc']}) 
Monomer('dNotch', ['Hbox', 'Rbox', 'Rbox'])
Monomer('dHes1', ['Hbox', 'Hbox', 'Hbox', 'Rbox', 'Rbox'])
Monomer('dRBPJ', ['Hbox', 'Hbox', 'Hbox', 'Rbox', 'Rbox', 'Rbox'])

Monomer('mNotch_delayed')
Monomer('pNotch_delayed')
Monomer('mRBPJ_delayed')
Monomer('pRBPJ_delayed')
Monomer('mHes1_delayed')
Monomer('pHes1_delayed')

Parameter('Delta_0', 0)
Parameter('mHes1_0', 1) 		# Hcm=1;		//3.64x10^(-12)
Parameter('pHes1_cyt_0', 35) 	# Hcp=35;		//1.25*10^(-10)
Parameter('pHes1_nuc_0', 109) 	# Hnp=109;		//2.77*10^(-9)
Parameter('dimHes1_nuc_0', int(0.02535*pHes1_nuc_0.value**2)) # 101) # H2np=101;		//2.5662*10^(-9)
Parameter('mRBPJ_0', 1) 		# Rcm=1;		//1.2*10^(-12)
Parameter('pRBPJ_cyt_0', 11) 	# Rcp=11;		//3.74*10^(-11)
Parameter('pRBPJ_nuc_0', 447) 	# Rnp=447;		//1.13*10^(-8)
Parameter('mNotch_0', 70) 		# Nm=70;		//2.51*10^(-10)	//Notch mRNA
Parameter('pNotch_0', 8130) 	#??# Np=8130;	//2.95*10^(-8)	//Notch protein in the cytoplasm
Parameter('dNotch_0', 2) 		# NP_free=2;	//Notch promoter (diploid so start with 2 free promoter regions)
Parameter('dHes1_0', 2) 		# HP_free=2;	//Hes1 promoter (diploid so start with 2 free promoter regions)
Parameter('dRBPJ_0', 2) 		# RP_free=2;	//RPJK promoter (diploid so start with 2 free promoter regions)

Initial(Delta(), Delta_0)
Initial(mHes1(), mHes1_0)
Initial(pHes1(pHes1=None, Hbox=None, loc='cyt'), pHes1_cyt_0)
Initial(pHes1(pHes1=None, Hbox=None, loc='nuc'), pHes1_nuc_0)
Initial(pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox= None, loc='nuc'), dimHes1_nuc_0)
Initial(mRBPJ(), mRBPJ_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc='cyt'), pRBPJ_cyt_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc='nuc'), pRBPJ_nuc_0)
Initial(mNotch(), mNotch_0)
Initial(pNotch(), pNotch_0)
Initial(dNotch(Hbox=None, Rbox=(None,None)), dNotch_0)
Initial(dHes1(Hbox=(None,None,None), Rbox=(None,None)), dHes1_0)
Initial(dRBPJ(Hbox=(None,None,None), Rbox=(None,None,None)), dRBPJ_0)

delay=True	#//set true here if you want to implement delays

#//Delay times for the 6 delayed reactions
#//units here are minutes delayed

# Degradation constant of Hes1 protein (min^-1) kdHcp, kdHnp 0.0315 [30]
Parameter('KdHes1', 0.0315)		#//degradation Hes1 protein
# Degradation constant of Hes1 mRNA (min^-1) kdHcm 0.029 [30]
Parameter('KdHcm', 0.029) 		#// deg Hes1 mRNA
# Degradation constant of RBP-Jk protein (min^-1) kdRcp 0.00231 [82]
Parameter('KdRBPJ', 0.00231) 	#degradation RBPJ protein
# Degradation constant of RBP-JK mRNA (min^-1) kdRcm 0.0075 [82]
Parameter('KdRcm', 0.0075) 		#degradation RBPJ mRNA
# Degradation constant of full-length Notch1 protein (min^-1) kdNp 0.017 [54], Text
Parameter('KdNp', 0.017) 		#degradation Notch protein
# Degradation constant of NICD protein (min^-1) kdNcp,kdNnp 0.0014 or 0.00385 [55], Text
Observable('pHes1_cyt', pHes1(loc='cyt'))
Expression('KdNICD', 0.00385 + (0.0014-0.00385)*(pHes1_cyt > 100)) #degradation NICD protein
# Degradation of Notch mRNA (min^-1) kdNm 0.0058 [83]
Parameter('KdNm', 0.0058) #0.00035) 		#degradation Notch mRNA
# Cooperativity factor for Hes1-DNA binding Cn 1 Text
# Cooperativity factor for RBP-Jk DNA binding Cr 1 Text
# Cooperativity factor for RBP-Jk Hes1 DNA binding Cnr 1 Text
# Rate of protein translation from Hes1 mRNA (min^-1) KtrHc 4.5 [45], Text
Parameter('KtrHc', 4.5)
# Rate of protein translation from RBP-Jk mRNA (min^-1) KtrRc 2.5 Text
Parameter('KtrRc', 2.5) #3.2)
# Rate of protein translation from Notch1 mRNA (min^-1) KtrN 1 Text
Parameter('KtrN', 1) #2)
# RBP-Jk DNA association constant (M^-1) Kr 3.23x10^8 [84]
Parameter('Kr', 0.0410)
Parameter('Kr_r', 5)
# Hes1 DNA association constant (M^-1) Kn 2x10^8 Text
Parameter('Kn', 0.0254)
Parameter('Kn_r', 5)
# RBP-Jk NICD association constant (M^-1) Ka 1x10^8 Text
Parameter('Ka', 0.0127)
Parameter('Ka_r', 5)
# Hes1 dimer association constant (M^-1) KaHp 1x10^9 Text
Parameter('KaHp', 2*0.02535)
Parameter('KrHp', 1) #0)
# Transcriptional time delay for Hes1 (min) TmHc 10 [52], Text
TmHc = 10 
# Translational time delay for Hes1 (min) TpHc 2.35 [52], Text
TpHc = 2.35 
# Transcriptional time delay for RBP-Jk (min) TmRc 20 [52], Text
TmRc = 20 
# Translational time delay for RBP-Jk (min) TpRc 4.3 [52], Text
TpRc = 4.3 
# Transcriptional time delay for Notch1 (min) TmNc 70 [52], Text
TmNc = 70 
# Translational time delay for Notch1 (min) TpNc 21 [52], Text
TpNc = 21 
# Basal transcriptional rate for Hes1 (M/min) Vbh 1.14x10^-10 [45], Text
Parameter('Vbh', 4.5) 		#//2.53*10^(-12)M/min
# Basal transcriptional rate for RBP-Jk (M/min) Vbr 4.3x10^-11 Text
Parameter('Vbr', 1.7) 		#//1.27*10^(-12)M/min
# Basal transcriptional rate for Notch1 (M/min) Vbn 1.23x10^-11 Text
Parameter('Vbn', 0.5) 		#//3.6*10^(-13)M/min
# Maximal transcriptional rate for Hes1 (M/min) Vmaxh 5x10^-10 [53], Text
Parameter('Vmaxh', 19.7) #197) 	#//9.98;
# Maximal transcriptional rate for RBP-Jk (M/min) Vmaxr 2x10^-10 Text
Parameter('Vmaxr', 7.9) #79) 		#//4.99
# Maximal transcriptional rate for Notch1 (M/min) Vmaxn 5.5x10^-11 Text
Parameter('Vmaxn', 2.16) #21.6) 	#//1.43
# Nuclear import rate of Hes1 protein (min^-1) kniHcp 0.1 [85], Text
Parameter('KniHcp', 0.1)
# Nuclear import rate of RBP-Jk protein (min^-1) kniRcp 0.1 [85], Text
Parameter('KniRcp', 0.1)
# Nuclear import rate of NICD protein (min^-1) kniNcp 0.1 [85], Text
Parameter('KniNcp', 0.1)
# NICD generation constant upon Delta binding (M^-1 min^-1) KfNcp 7.6x10^7 Text
Parameter('KfNcp', 0.000276)
# Repression constant of Hes1 bound to N-box rNbox 0.3 [30], Text
Parameter('rNbox', 0.03)
# Repression constant of RBP-Jk alone bound to promoter rR 0.2 Text
Parameter('rRbox', 0.2)

#Degradation parameters
#	//Units here are per minute
Parameter('KdH2np', 0) 			#//disable Hes1 dimer degradation <-- I didn't write that -MS

#	//units here are per minute
#protein translation rates

# units here are per minute
#nuclear rate of import

#//7.6*10^8 per molar per minute, forward rate of NICD binding to Delta

# Hes1 dimer association constant

#RBPJ-DNA association constant; //Keqr=8.19e-3; //float Kr=0.001628;

#Kn = Hes1-DNA association constant; // Keqn = 5.07e-3; // float Kn = 0.00072;

#Ka = RBPJ-NICD association constant; // Keqa = 2.54e-3;  // float Ka = 0.00036;

#//These generation rates are based on molecules generated per cell per minute <= from original code

Parameter('tc', 0.5)

# Hes1 dimerization
# H2np = KaHp*Hnp*Hnp;  // THIS MIGHT NEED ADJUSTMENT
Rule('Hes1_dimerization', pHes1(pHes1=None, Hbox=None, loc='nuc') + pHes1(pHes1=None, Hbox=None, loc='nuc') | 
	pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox=None, loc='nuc'), KaHp, KrHp)

# //Generation reactions
# a[0]=c[0];		//1. 0 --> Nm		0.RfNm
# a[1]=c[1];		//2. 0 --> Rcm		1.RfRm
# a[2]=c[2];		//3. 0 --> Hcm		2.RfHm

Tag('x')

Observable('dNotch_pHes1', dNotch(Hbox=ANY))
Observable('dNotch_pRBPJ', dNotch(Rbox=1) % pRBPJ(Rbox=1, nicd=None))
Observable('dNotch_pNICD', dNotch(Rbox=1) % pRBPJ(Rbox=1, nicd=ANY))
Expression('ktrN', rNbox**dNotch_pHes1(x) * rRbox**dNotch_pRBPJ(x) * ((dNotch_pNICD(x) > 0)*tc**(2-dNotch_pNICD(x))*Vmaxn + Vbn))
Rule('mNotch_transcription', dNotch() @x >> dNotch() @x + mNotch_delayed(), ktrN)
Rule('mNotch_appears', mNotch_delayed() >> mNotch(), Parameter('kmN_delay', 1./TmNc))

Observable('dRBPJ_pHes1', dRBPJ(Hbox=ANY))
Observable('dRBPJ_pRBPJ', dRBPJ(Rbox=1) % pRBPJ(Rbox=1, nicd=None))
Observable('dRBPJ_pNICD', dRBPJ(Rbox=1) % pRBPJ(Rbox=1, nicd=ANY))
Expression('ktrR', rNbox**dRBPJ_pHes1(x) * rRbox**dRBPJ_pRBPJ(x) * ((dRBPJ_pNICD(x) > 0)*tc**(3-dRBPJ_pNICD(x))*Vmaxr + Vbr))
Rule('mRBPJ_transcription',  dRBPJ() @x  >> dRBPJ() @x + mRBPJ_delayed(), ktrR)
Rule('mRBPJ_appears', mRBPJ_delayed() >> mRBPJ(), Parameter('kmR_delay', 1./TmRc))

Observable('dHes1_pHes1', dHes1(Hbox=ANY))
Observable('dHes1_pRBPJ', dHes1(Rbox=1) % pRBPJ(Rbox=1, nicd=None))
Observable('dHes1_pNICD', dHes1(Rbox=1) % pRBPJ(Rbox=1, nicd=ANY))
Expression('ktrH', rNbox**dHes1_pHes1(x) * rRbox**dHes1_pRBPJ(x) * ((dHes1_pNICD(x) > 0)*tc**(2-dHes1_pNICD(x))*Vmaxh + Vbh))
Rule('mHes1_transcription', dHes1() @x  >> dHes1() @x + mHes1_delayed(), ktrH)
Rule('mHes1_appears', mHes1_delayed() >> mHes1(), Parameter('kmH_delay', 1./TmHc))

# //Translation channels (delay these if they occur in delay model)
# a[3]=c[3]*Nm;			//4. Nm --> Np		3.KtrN
# a[4]=c[4]*Rcm;		//5. Rcm --> Rcp	4.KtrRc
# a[5]=c[5]*Hcm;		//6. Hcm --> Hcp	5.KtrHc

Rule('mNotch_translation', mNotch() >> mNotch() + pNotch_delayed() , KtrN)
Rule('mRBPJ_translation', mRBPJ() >> mRBPJ() + pRBPJ_delayed(), KtrRc)
Rule('mHes1_translation', mHes1() >> mHes1() + pHes1_delayed(), KtrHc)
 
Rule('pNotch_appears', pNotch_delayed() >> pNotch(), Parameter('kpN_delay', 1./TpNc))
Rule('pRBPJ_appears',  pRBPJ_delayed()  >> pRBPJ(nicd=None, Rbox=None, loc='cyt'), Parameter('kpR_delay', 1./TpRc))
Rule('pHes1_appears',  pHes1_delayed()  >> pHes1(pHes1=None, Hbox=None, loc='cyt'), Parameter('kpH_delay', 1./TpHc))
 
# //Other reactions in the cytoplasm
# //a[6]=c[6]*Hnp*Hnp;	//7. Hnp + Hnp --> H2np		6.KaHp
# //a[7]=c[7]*H2np;		//8. H2np --> Hnp + Hnp		7.KrHp
# a[8]=c[8]*delta*Np;	//9. delta + Np --> Ncp		8.KfNcp

Rule('Delta_Notch_to_NICD', Delta() + pNotch() >> Delta() + NICD(rbpj=None, loc='cyt'), KfNcp)

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
# //a[19]=0;			//20. H2np --> 0	18.KdH2np

Rule('mNotch_degradation', mNotch() >> None, KdNm) # //10. Nm --> 0  KdNm
Rule('pNotch_degradation', pNotch() >> None, KdNp) # //15. Np --> 0  KdNp

Rule('mHes1_degradation', mHes1() >> None, KdHcm) # //12. Hcm --> 0  KdHcm
Rule('pHes1_monomer_degradation', pHes1(pHes1=None, Hbox=None) >> None, KdHes1) # //18. Hcp --> 0  KdHcp; //19. Hnp --> 0  KdHnp
Rule('pHes1_dimer_degradation', pHes1(pHes1=ANY, Hbox=None, loc='nuc') >> None, KdH2np) # //20. H2np --> 0  KdH2np

Rule('NICD_degradation', NICD(rbpj=None) >> None, KdNICD) # //13. Nnp --> 0  KdNnp; //16. Ncp --> 0  KdNcp 

Rule('mRBPJ_degradation', mRBPJ() >> None, KdRcm) # //11. Rcm --> 0  KdRcm
Rule('pRBPJ_degradation', pRBPJ(nicd=None, Rbox=None) >> None, KdRBPJ) # //14. Rnp --> 0  KdRnp; //17. Rcp --> 0  KdRcp

# //Protein transport channels cyt --> "nuc only transport into nucleus"
# a[20]=c[19]*Ncp;		//21. Ncp --> Nnp			19.KniNcp
# a[21]=c[20]*Rcp;		//22. Rcp --> Rnp			20.KniRcp
# a[22]=c[21]*Hcp;		//23. Hcp --> Hnp			21.KniHcp

Rule('NICD_nuclear_transport', NICD(rbpj=None, loc='cyt') >> NICD(rbpj=None, loc='nuc'), KniNcp)
Rule('RBPJ_nuclear_transport', pRBPJ(nicd=None, Rbox=None, loc='cyt') >> pRBPJ(nicd=None, Rbox=None, loc='nuc'), KniRcp)
Rule('Hes1_nuclear_transport', pHes1(pHes1=None, Hbox=None, loc='cyt') >> pHes1(pHes1=None, Hbox=None, loc='nuc'), KniHcp)

# //Notch Promoter binding / unbinding (18*2 reactions)
# a[23]=c[22]*2*NP_free*Rnp;		//24. NP_free      + Rnp  --> NP_Rr
# a[27]=c[22]*NP_Rr*Rnp;			//28. NP_Rr        + Rnp  --> NP_Rr2
# a[28]=c[22]*2*NP_Hr*Rnp;			//29. NP_Hr        + Rnp  --> NP_Rr_Hr
# a[29]=c[22]*NP_Rr_Hr*Rnp;			//30. NP_Rr_Hr     + Rnp  --> NP_Rr2_Hr
# a[31]=c[22]*NP_RNa*Rnp;			//32. NP_RNa       + Rnp  --> NP_RNa_Rr
# a[36]=c[22]*NP_RNa_Hr*Rnp;		//37. NP_RNa_Hr    + Rnp  --> NP_RNa_Rr_Hr
#
# a[41]=c[23]*NP_Rr;				//42. NP_Rr        --> NP_free      + Rnp
# a[45]=c[25]*2*NP_Rr2;				//46. NP_Rr2       --> NP_Rr        + Rnp
# a[49]=c[25]*NP_Rr_Hr;				//50. NP_Rr_Hr     --> NP_Hr        + Rnp
# a[55]=c[25]*2*NP_Rr2_Hr;			//56. NP_Rr2_Hr    --> NP_Rr_Hr     + Rnp
# a[43]=c[25]*NP_RNa_Rr;			//44. NP_RNa_Rr    --> NP_RNa       + Rnp
# a[52]=c[25]*NP_RNa_Rr_Hr;			//53. NP_RNa_Rr_Hr --> NP_RNa_Hr    + Rnp

Rule('dNotch_bind_Rbox', dNotch(Rbox=None) + pRBPJ(nicd=None, Rbox=None, loc='nuc') |
	dNotch(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc'), Kr,Kr_r)

# a[24]=c[24]*NP_free*H2np;			//25. NP_free      + H2np --> NP_Hr
# a[25]=c[24]*NP_Rr*H2np;			//26. NP_Rr        + H2np --> NP_Rr_Hr
# a[32]=c[24]*NP_RNa*H2np;			//33. NP_RNa       + H2np --> NP_RNa_Hr
# a[34]=c[24]*NP_Rr2*H2np;			//35. NP_Rr2       + H2np --> NP_Rr2_Hr
# a[38]=c[24]*NP_RNa_Rr*H2np;		//39. NP_RNa_Rr    + H2np --> NP_RNa_Rr_Hr
# a[40]=c[24]*NP_RNa2*H2np;			//41. NP_RNa2      + H2np --> NP_RNa2_Hr
#
# a[47]=c[23]*NP_Hr;				//48. NP_Hr        --> NP_free      + H2np
# a[48]=c[25]*NP_Rr_Hr;				//49. NP_Rr_Hr     --> NP_Rr        + H2np
# a[51]=c[25]*NP_RNa_Hr;			//52. NP_RNa_Hr    --> NP_RNa       + H2np
# a[56]=c[25]*NP_Rr2_Hr;			//57. NP_Rr2_Hr    --> NP_Rr2       + H2np
# a[54]=c[25]*NP_RNa_Rr_Hr;			//55. NP_RNa_Rr_Hr --> NP_RNa_Rr    + H2np
# a[58]=c[25]*NP_RNa2_Hr;			//59. NP_RNa2_Hr   --> NP_RNa2      + H2np

Rule('dNotch_bind_Hes1_dimer', dNotch(Hbox=None) + pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox=None, loc='nuc') |
	dNotch(Hbox=[2,3]) % pHes1(pHes1=1, Hbox=2, loc='nuc') % pHes1(pHes1=1, Hbox=3, loc='nuc'), Kn,Kn_r)

# a[26]=c[26]*NP_Rr*Nnp;			//27. NP_Rr        + Nnp  --> NP_RNa
# a[30]=c[26]*NP_Rr_Hr*Nnp;			//31. NP_Rr_Hr     + Nnp  --> NP_RNa_Hr
# a[33]=c[26]*2*NP_Rr2*Nnp;			//34. NP_Rr2       + Nnp  --> NP_RNa_Rr
# a[35]=c[26]*2*NP_Rr2_Hr*Nnp;		//36. NP_Rr2_Hr    + Nnp  --> NP_RNA_Rr_Hr
# a[37]=c[26]*NP_RNa_Rr*Nnp;		//38. NP_RNa_Rr    + Nnp  --> NP_RNa2
# a[39]=c[26]*NP_RNa_Rr_Hr*Nnp;		//40. NP_RNa_Rr_Hr + Nnp  --> NP_RNa2_Hr
#
# a[42]=c[25]*NP_RNa;				//43. NP_RNa       --> NP_Rr        + Nnp
# a[50]=c[25]*NP_RNa_Hr;			//51. NP_RNa_Hr    --> NP_Rr_Hr     + Nnp
# a[44]=c[25]*NP_RNa_Rr;			//45. NP_RNa_Rr    --> NP_Rr2       + Nnp
# a[53]=c[25]*NP_RNa_Rr_Hr;			//54. NP_RNa_Rr_Hr --> NP_Rr2_Hr    + Nnp
# a[46]=c[25]*NP_RNa2;				//47. NP_RNa2      --> NP_RNa_Rr    + Nnp
# a[57]=c[25]*2*NP_RNa2_Hr;			//58. NP_RNa2_Hr   --> NP_RNa_Rr_Hr + Nnp

Rule('dNotch_pRBPJ_bind_NICD', dNotch(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc') + NICD(rbpj=None, loc='nuc') |
	dNotch(Rbox=1) % pRBPJ(nicd=2, Rbox=1, loc='nuc') % NICD(rbpj=2, loc='nuc'), Ka,Ka_r)

# // Hes1 promoter binding / unbinding
# a[59]=c[22]*2*HP_free*Rnp;		//60. HP_free       + Rnp  --> HP_Rr
# a[61]=c[22]*2*HP_Hr*Rnp;			//62. HP_Hr         + Rnp  --> HP_Rr_Hr
# a[63]=c[22]*2*HP_Hr2*Rnp;			//64. HP_Hr2        + Rnp  --> HP_Rr_Hr2
# a[65]=c[22]*2*HP_Hr3*Rnp;			//66. HP_Hr3        + Rnp  --> HP_Rr_Hr3
# a[66]=c[22]*HP_Rr*Rnp;			//67. HP_Rr         + Rnp  --> HP_Rr2
# a[69]=c[22]*HP_Rr_Hr*Rnp;			//70. HP_Rr_Hr      + Rnp  --> HP_Rr2_Hr
# a[72]=c[22]*HP_Rr_Hr2*Rnp;		//73. HP_Rr_Hr2     + Rnp  --> HP_Rr2_Hr2
# a[75]=c[22]*HP_Rr_Hr3*Rnp;		//76. HP_Rr_Hr3     + Rnp  --> HP_Rr2_Hr3
# a[77]=c[22]*HP_RNa*Rnp;			//78. HP_RNa        + Rnp  --> HP_RNa_Rr
# a[79]=c[22]*HP_RNa_Hr*Rnp;		//80. HP_RNa_Hr     + Rnp  --> HP_RNa_Rr_Hr
# a[81]=c[22]*HP_RNa_Hr2*Rnp;		//82. HP_RNa_Hr2    + Rnp  --> HP_RNa_Rr_Hr2
# a[83]=c[22]*HP_RNa_Hr3*Rnp;		//84. HP_RNa_Hr3    + Rnp  --> HP_RNa_Rr_Hr3
#
# a[104]=c[23]*HP_Rr;				//112. HP_Rr         --> HP_free       + Rnp
# a[105]=c[25]*HP_Rr_Hr;			//113. HP_Rr_Hr      --> HP_Hr         + Rnp
# a[107]=c[25]*HP_Rr_Hr2;			//115. HP_Rr_Hr2     --> HP_Hr2        + Rnp
# a[109]=c[25]*HP_Rr_Hr3;			//117. HP_Rr_Hr3     --> HP_Hr3        + Rnp
# a[118]=c[25]*HP_RNa_Rr;			//126. HP_RNa_Rr     --> HP_RNa        + Rnp
# a[120]=c[25]*HP_RNa_Rr_Hr;		//128. HP_RNa_Rr_Hr  --> HP_RNa_Hr     + Rnp
# a[123]=c[25]*HP_RNa_Rr_Hr2;		//131. HP_RNa_Rr_Hr2 --> HP_RNa_Hr2    + Rnp
# a[126]=c[25]*HP_RNa_Rr_Hr3;		//134. HP_RNa_Rr_Hr3 --> HP_RNa_Hr3    + Rnp
# a[129]=c[25]*2*HP_Rr2;			//137. HP_Rr2        --> HP_Rr         + Rnp
# a[130]=c[25]*HP_Rr2_Hr;			//138. HP_Rr2_Hr     --> HP_Rr_Hr      + Rnp
# a[132]=c[25]*HP_Rr2_Hr2;			//140. HP_Rr2_Hr2    --> HP_Rr_Hr2     + Rnp
# a[134]=c[25]*HP_Rr2_Hr3;			//142. HP_Rr2_Hr3    --> HP_Rr_Hr3     + Rnp

Rule('dHes1_bind_Rbox', dHes1(Rbox=None) + pRBPJ(nicd=None, Rbox=None, loc='nuc') |
	dHes1(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc'), Kr,Kr_r)

# a[60]=c[24]*3*HP_free*H2np;		//61. HP_free       + H2np --> HP_Hr
# a[62]=c[24]*2*HP_Hr*H2np;			//63. HP_Hr         + H2np --> HP_Hr2
# a[64]=c[24]*HP_Hr2*H2np;			//65. HP_Hr2        + H2np --> HP_Hr3
# a[68]=c[24]*3*HP_Rr*H2np;			//69. HP_Rr         + H2np --> HP_Rr_Hr
# a[71]=c[24]*2*HP_Rr_Hr*H2np;		//72. HP_Rr_Hr      + H2np --> HP_Rr_Hr2
# a[74]=c[24]*HP_Rr_Hr2*H2np;		//75. HP_Rr_Hr2     + H2np --> HP_Rr_Hr3
# a[78]=c[24]*3*HP_RNa*H2np;		//79. HP_RNa        + H2np --> HP_RNa_Hr
# a[80]=c[24]*2*HP_RNa_Hr*H2np;		//81. HP_RNa_Hr     + H2np --> HP_RNa_Hr2
# a[82]=c[24]*HP_RNa_Hr2*H2np;		//83. HP_RNa_Hr2    + H2np --> HP_RNa_Hr3
# a[85]=c[24]*3*HP_RNa_Rr*H2np;		//86. HP_RNa_Rr     + H2np --> HP_RNa_Rr_Hr
# a[87]=c[24]*2*HP_RNa_Rr_Hr*H2np;	//89. HP_RNa_Rr_Hr  + H2np --> HP_RNa_Rr_Hr2
# a[89]=c[24]*HP_RNa_Rr_Hr2*H2np;	//92. HP_RNa_Rr_Hr2 + H2np --> HP_RNa_Rr_Hr3
# a[92]=c[24]*3*HP_Rr2*H2np;		//97. HP_Rr2        + H2np --> HP_Rr2_Hr
# a[94]=c[24]*2*HP_Rr2_Hr*H2np;		//100. HP_Rr2_Hr     + H2np --> HP_Rr2_Hr2
# a[96]=c[24]*HP_Rr2_Hr2*H2np;		//103. HP_Rr2_Hr2    + H2np --> HP_Rr2_Hr3
# a[98]=c[24]*3*HP_RNa2*H2np;		//106. HP_RNa2       + H2np --> HP_RNa2_Hr
# a[99]=c[24]*2*HP_RNa2_Hr*H2np;	//107. HP_RNa2_Hr    + H2np --> HP_RNa2_Hr2
# a[100]=c[24]*HP_RNa2_Hr2*H2np;	//108. HP_RNa2_Hr2   + H2np --> HP_RNa2_Hr3
#
# a[101]=c[23]*HP_Hr;				//109. HP_Hr         --> HP_free       + H2np
# a[102]=c[25]*2*HP_Hr2;			//110. HP_Hr2        --> HP_Hr         + H2np
# a[103]=c[25]*3*HP_Hr3;			//111. HP_Hr3        --> HP_Hr2        + H2np
# a[106]=c[25]*HP_Rr_Hr;			//114. HP_Rr_Hr      --> HP_Rr         + H2np
# a[108]=c[25]*2*HP_Rr_Hr2;			//116. HP_Rr_Hr2     --> HP_Rr_Hr      + H2np
# a[110]=c[25]*3*HP_Rr_Hr3;			//118. HP_Rr_Hr3     --> HP_Rr_Hr2     + H2np
# a[113]=c[25]*HP_RNa_Hr;			//121. HP_RNa_Hr     --> HP_RNa        + H2np
# a[115]=c[25]*2*HP_RNa_Hr2;		//123. HP_RNa_Hr2    --> HP_RNa_Hr     + H2np
# a[117]=c[25]*3*HP_RNa_Hr3;		//125. HP_RNa_Hr3    --> HP_RNa_Hr2    + H2np
# a[122]=c[25]*HP_RNa_Rr_Hr;		//130. HP_RNa_Rr_Hr  --> HP_RNa_Rr     + H2np
# a[125]=c[25]*2*HP_RNa_Rr_Hr2;		//133. HP_RNa_Rr_Hr2 --> HP_RNa_Rr_Hr  + H2np
# a[128]=c[25]*3*HP_RNa_Rr_Hr3;		//136. HP_RNa_Rr_Hr3 --> HP_RNa_Rr_Hr2 + H2np
# a[131]=c[25]*HP_Rr2_Hr;			//139. HP_Rr2_Hr     --> HP_Rr2        + H2np
# a[133]=c[25]*2*HP_Rr2_Hr2;		//141. HP_Rr2_Hr2    --> HP_Rr2_Hr     + H2np
# a[135]=c[25]*3*HP_Rr2_Hr3;		//143. HP_Rr2_Hr3    --> HP_Rr2_Hr2    + H2np
# a[138]=c[25]*HP_RNa2_Hr;			//146. HP_RNa2_Hr    --> HP_RNa2       + H2np
# a[140]=c[25]*2*HP_RNa2_Hr2;		//148. HP_RNa2_Hr2   --> HP_RNa2_Hr    + H2np
# a[142]=c[25]*3*HP_RNa2_Hr3;		//150. HP_RNa2_Hr3   --> HP_RNa2_Hr2   + H2np

Rule('dHes1_bind_Hes1_dimer', dHes1(Hbox=None) + pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox=None, loc='nuc') |
	dHes1(Hbox=[2,3]) % pHes1(pHes1=1, Hbox=2, loc='nuc') % pHes1(pHes1=1, Hbox=3, loc='nuc'), Kn,Kn_r)

# a[67]=c[22]*HP_Rr*Nnp;			//68. HP_Rr         + Nnp  --> HP_RNa
# a[70]=c[26]*HP_Rr_Hr*Nnp;			//71. HP_Rr_Hr      + Nnp  --> HP_RNa_Hr
# a[73]=c[26]*HP_Rr_Hr2*Nnp;		//74. HP_Rr_Hr2     + Nnp  --> HP_RNa_Hr2
# a[76]=c[26]*HP_Rr_Hr3*Nnp;		//77. HP_Rr_Hr3     + Nnp  --> HP_RNa_Hr3
# a[84]=c[26]*HP_RNa_Rr*Nnp;		//85. HP_RNa_Rr     + Nnp  --> HP_RNa2
# a[86]=c[26]*HP_RNa_Rr_Hr*Nnp;		//88. HP_RNa_Rr_Hr  + Nnp  --> HP_RNa2_Hr
# a[88]=c[26]*HP_RNa_Rr_Hr2*Nnp;	//91. HP_RNa_Rr_Hr2 + Nnp  --> HP_RNa2_Hr2
# a[90]=c[26]*HP_RNa_Rr_Hr3*Nnp;	//94. HP_RNa_Rr_Hr3 + Nnp  --> HP_RNa2_Hr3
# a[91]=c[26]*2*HP_Rr2*Nnp;			//96. HP_Rr2        + Nnp  --> HP_RNa_Rr
# a[93]=c[26]*2*HP_Rr2_Hr*Nnp;		//99. HP_Rr2_Hr     + Nnp  --> HP_RNa_Rr_Hr
# a[95]=c[26]*2*HP_Rr2_Hr2*Nnp;		//102. HP_Rr2_Hr2    + Nnp  --> HP_RNa_Rr_Hr2
# a[97]=c[26]*2*HP_Rr2_Hr3*Nnp;		//105. HP_Rr2_Hr3    + Nnp  --> HP_RNa_Rr_Hr3
#
# a[111]=c[25]*HP_RNa;				//119. HP_RNa        --> HP_Rr         + Nnp
# a[112]=c[25]*HP_RNa_Hr;			//120. HP_RNa_Hr     --> HP_Rr_Hr      + Nnp
# a[114]=c[25]*HP_RNa_Hr2;			//122. HP_RNa_Hr2    --> HP_Rr_Hr2     + Nnp
# a[116]=c[25]*HP_RNa_Hr3;			//124. HP_RNa_Hr3    --> HP_Rr_Hr3     + Nnp
# a[119]=c[25]*HP_RNa_Rr;			//127. HP_RNa_Rr     --> HP_Rr2        + Nnp
# a[121]=c[25]*HP_RNa_Rr_Hr;		//129. HP_RNa_Rr_Hr  --> HP_Rr2_Hr     + Nnp
# a[124]=c[25]*HP_RNa_Rr_Hr2;		//132. HP_RNa_Rr_Hr2 --> HP_Rr2_Hr2    + Nnp
# a[127]=c[25]*HP_RNa_Rr_Hr3;		//135. HP_RNa_Rr_Hr3 --> HP_Rr2_Hr3    + Nnp
# a[136]=c[25]*2*HP_RNa2;			//144. HP_RNa2       --> HP_RNa_Rr     + Nnp
# a[137]=c[25]*2*HP_RNa2_Hr;		//145. HP_RNa2_Hr    --> HP_RNa_Rr_Hr  + Nnp
# a[139]=c[25]*2*HP_RNa2_Hr2;		//147. HP_RNa2_Hr2   --> HP_RNa_Rr_Hr2 + Nnp
# a[141]=c[25]*2*HP_RNa2_Hr3;		//149. HP_RNa2_Hr3   --> HP_RNa_Rr_Hr3 + Nnp

Rule('dHes1_pRBPJ_bind_NICD', dHes1(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc') + NICD(rbpj=None, loc='nuc') |
	dHes1(Rbox=1) % pRBPJ(nicd=2, Rbox=1, loc='nuc') % NICD(rbpj=2, loc='nuc'), Ka,Ka_r)

# //CBF1 promoter binding / unbinding
# a[143]=c[22]*3*RP_free*Rnp;		//151. RP_free + Rnp         --> RP_Rr
# a[145]=c[22]*3*RP_Hr*Rnp;			//153. RP_Hr + Rnp           --> RP_Rr_Hr
# a[147]=c[22]*3*RP_Hr2*Rnp;		//155. RP_Hr2 + Rnp          --> RP_Rr_Hr2
# a[149]=c[22]*3*RP_Hr3*Rnp;		//157. RP_Hr3 + Rnp          --> RP_Rr_Hr3
# a[150]=c[22]*2*RP_Rr*Rnp;			//158. RP_Rr + Rnp           --> RP_Rr2
# a[153]=c[22]*2*RP_Rr_Hr*Rnp;		//161. RP_Rr_Hr + Rnp        --> RP_Rr2_Hr
# a[156]=c[22]*2*RP_Rr_Hr2*Rnp;		//164. RP_Rr_Hr2 + Rnp       --> RP_Rr2_Hr2
# a[159]=c[22]*2*RP_Rr_Hr3*Rnp;		//167. RP_Rr_Hr3 + Rnp       --> RP_Rr2_Hr3
# a[161]=c[22]*2*RP_RNa*Rnp;		//169. RP_RNa + Rnp          --> RP_RNa_Rr
# a[163]=c[22]*2*RP_RNa_Hr*Rnp;		//171. RP_RNa_Hr + Rnp       --> RP_RNa_Rr_Hr
# a[165]=c[22]*2*RP_RNa_Hr2*Rnp;	//173. RP_RNa_Hr2 + Rnp      --> RP_RNa_Rr_Hr2
# a[167]=c[22]*2*RP_RNa_Hr3*Rnp;	//175. RP_RNa_Hr3 + Rnp      --> RP_RNa_Rr_Hr3
# a[168]=c[22]*RP_Rr2*Rnp;			//176. RP_Rr2 + Rnp          --> RP_Rr3
# a[171]=c[24]*RP_Rr2_Hr*Rnp;		//179. RP_Rr2_Hr + Rnp       --> RP_Rr3_Hr
# a[174]=c[22]*RP_Rr2_Hr2*Rnp;		//182. RP_Rr2_Hr2 + Rnp      --> RP_Rr3_Hr2
# a[177]=c[22]*RP_Rr2_Hr3*Rnp;		//185. RP_Rr2_Hr3 + Rnp      --> RP_Rr3_Hr3
# a[179]=c[22]*RP_RNa_Rr*Rnp;		//187. RP_RNa_Rr + Rnp       --> RP_RNa_Rr2
# a[182]=c[22]*RP_RNa_Rr_Hr*Rnp;	//190. RP_RNa_Rr_Hr + Rnp    --> RP_Rr2_Hr
# a[185]=c[22]*RP_RNa_Rr_Hr2*Rnp;	//193. RP_RNa_Rr_Hr2 + Rnp   --> RP_RNa_Rr2_Hr2
# a[188]=c[22]*RP_RNa_Rr_Hr3*Rnp;	//196. RP_RNa_Rr_Hr3 + Rnp   --> RP_RNa_Rr2_Hr3
# a[204]=c[22]*RP_RNa2*Rnp;			//212. RP_RNa2 + Rnp         --> RP_RNa2_Rr
# a[206]=c[22]*RP_RNa2_Hr*Rnp;		//214. RP_RNa2_Hr + Rnp      --> RP_RNa2_Rr_Hr
# a[208]=c[22]*RP_RNa2_Hr2*Rnp;		//216. RP_RNa2_Hr2 + Rnp     --> RP_RNa2_Rr_Hr2
# a[210]=c[22]*RP_RNa2_Hr3*Rnp;		//218. RP_RNa2_Hr3 + Rnp     --> RP_RNa2_Rr_Hr3
#
# a[224]=c[23]*RP_Rr;				//232. RP_Rr				--> RP_free + Rnp
# a[225]=c[25]*RP_Rr_Hr;			//233. RP_Rr_Hr				--> RP_Hr + Rnp
# a[227]=c[25]*RP_Rr_Hr2;			//235. RP_Rr_Hr2			--> RP_Hr2 + Rnp
# a[229]=c[25]*RP_Rr_Hr3;			//237. RP_Rr_Hr3			--> RP_Hr3 + Rnp
# a[238]=c[25]*2*RP_Rr2;			//246. RP_Rr2				--> RP_Rr + Rnp
# a[239]=c[25]*2*RP_Rr2_Hr;			//247. RP_Rr2_Hr			--> RP_Rr_Hr + Rnp
# a[241]=c[25]*2*RP_Rr2_Hr2;		//249. RP_Rr2_Hr2			--> RP_Rr_Hr2 + Rnp
# a[243]=c[25]*2*RP_Rr2_Hr3;		//251. RP_Rr2_Hr3			--> RP_Rr_Hr3 + Rnp
# a[245]=c[25]*RP_RNa_Rr;			//253. RP_RNa_Rr			--> RP_RNa + Rnp
# a[247]=c[25]*RP_RNa_Rr_Hr;		//255. RP_RNa_Rr_Hr			--> RP_RNa_Hr + Rnp
# a[250]=c[25]*RP_RNa_Rr_Hr2;		//258. RP_RNa_Rr_Hr2		--> RP_RNa_Hr2 + Rnp
# a[253]=c[25]*RP_RNa_Rr_Hr3;		//261. RP_RNa_Rr_Hr3		--> RP_RNa_Hr3 + Rnp
# a[256]=c[25]*3*RP_Rr3;			//264. RP_Rr3				--> RP_Rr2 + Rnp
# a[257]=c[25]*3*RP_Rr3_Hr;			//265. RP_Rr3_Hr			--> RP_Rr2_Hr + Rnp
# a[259]=c[25]*3*RP_Rr3_Hr2;		//267. RP_Rr3_Hr2			--> RP_Rr2_Hr2 + Rnp
# a[261]=c[25]*3*RP_Rr3_Hr3;		//269. RP_Rr3_Hr3			--> RP_Rr2_Hr3 + Rnp
# a[263]=c[25]*2*RP_RNa_Rr2;		//271. RP_RNa_Rr2			--> RP_RNa_Rr + Rnp
# a[265]=c[25]*2*RP_RNa_Rr2_Hr;		//273. RP_RNa_Rr2_Hr		--> RP_RNa_Rr_Hr + Rnp
# a[268]=c[25]*2*RP_RNa_Rr2_Hr2;	//276. RP_RNa_Rr2_Hr2		--> RP_RNa_Rr_Hr2 + Rnp
# a[271]=c[25]*2*RP_RNa_Rr2_Hr3;	//279. RP_RNa_Rr2_Hr3		--> RP_RNa_Rr_Hr3 + Rnp
# a[281]=c[25]*RP_RNa2_Rr;			//289. RP_RNa2_Rr 			--> RP_RNa2 + Rnp
# a[283]=c[25]*RP_RNa2_Rr_Hr;		//291. RP_RNa2_Rr_Hr 		--> RP_RNa2_Hr + Rnp
# a[286]=c[25]*RP_RNa2_Rr_Hr2;		//294. RP_RNa2_Rr_Hr2		--> RP_RNa2_Hr + Rnp
# a[289]=c[25]*RP_RNa2_Rr_Hr3;		//297. RP_RNa2_Rr_Hr3 		--> RP_RNa2_Hr3 + Rnp

Rule('dRBPJ_bind_Rbox', dRBPJ(Rbox=None) + pRBPJ(nicd=None, Rbox=None, loc='nuc') |
	dRBPJ(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc'), Kr,Kr_r)

# a[144]=c[24]*3*RP_free*H2np;		//152. RP_free + H2np        --> RP_Hr
# a[146]=c[24]*2*RP_Hr*H2np;		//154. RP_Hr + H2np          --> RP_Hr2
# a[148]=c[24]*RP_Hr2*H2np;			//156. RP_Hr2 + H2np         --> RP_Hr3
# a[152]=c[24]*3*RP_Rr*H2np;		//160. RP_Rr + H2np          --> RP_Rr_Hr
# a[155]=c[24]*2*RP_Rr_Hr*H2np;		//163. RP_Rr_Hr + H2np       --> RP_Rr_Hr2
# a[158]=c[24]*RP_Rr_Hr2*H2np;		//166. RP_Rr_Hr2 + H2np      --> RP_Rr_Hr3
# a[162]=c[24]*3*RP_RNa*H2np;		//170. RP_RNa + H2np         --> RP_RNa_Hr
# a[164]=c[24]*2*RP_RNa_Hr*H2np;	//172. RP_RNa_Hr + H2np      --> RP_RNa_Hr2
# a[166]=c[24]*RP_RNa_Hr2*H2np;		//174. RP_RNa_Hr2 + H2np     --> RP_RNa_Hr3
# a[170]=c[24]*3*RP_Rr2*H2np;		//178. RP_Rr2 + H2np         --> RP_Rr2_Hr
# a[173]=c[24]*2*RP_Rr2_Hr*H2np;	//181. RP_Rr2_Hr + H2np      --> RP_Rr2_Hr2
# a[176]=c[24]*RP_Rr2_Hr2*H2np;		//184. RP_Rr2_Hr2 + H2np     --> RP_Rr2_Hr3
# a[181]=c[24]*3*RP_RNa_Rr*H2np;	//189. RP_RNa_Rr + H2np      --> RP_RNa_Rr_Hr
# a[184]=c[24]*2*RP_RNa_Rr_Hr*H2np;	//192. RP_RNa_Rr_Hr + H2np   --> RP_RNa_Rr_Hr2
# a[187]=c[24]*RP_RNa_Rr_Hr2*H2np;	//195. RP_RNa_Rr_Hr2 + H2np  --> RP_RNa_Rr_Hr3
# a[191]=c[24]*3*RP_Rr3*H2np;		//199. RP_Rr3 + H2np         --> RP_Rr3_Hr
# a[193]=c[24]*2*RP_Rr3_Hr*H2np;	//201. RP_Rr3_Hr + H2np      --> RP_Rr3_Hr2
# a[195]=c[24]*RP_Rr3_Hr2*H2np;		//203. RP_Rr3_Hr2 + H2np     --> RP_Rr3_Hr3
# a[198]=c[24]*3*RP_RNa_Rr2*H2np;	//206. RP_RNa_Rr2 + H2np     --> RP_RNa_Rr2_Hr
# a[200]=c[24]*2*RP_RNa_Rr2_Hr*H2np;//208. RP_RNa_Rr2_Hr + H2np  --> RP_RNa_Rr2_Hr2
# a[202]=c[24]*RP_RNa_Rr2_Hr2*H2np;	//210. RP_RNa_Rr2_Hr2 + H2np --> RP_RNa_Rr2_Hr3
# a[205]=c[24]*3*RP_RNa2*H2np;		//213. RP_RNa2 + H2np        --> RP_RNa2_Hr
# a[207]=c[24]*2*RP_RNa2_Hr*H2np;	//215. RP_RNa2_Hr + H2np     --> RP_RNa2_Hr2
# a[209]=c[24]*RP_RNa2_Hr2*H2np;	//217. RP_RNa2_Hr2 + H2np    --> RP_RNa2_Hr3
# a[212]=c[24]*3*RP_RNa2_Rr*H2np;	//220. RP_RNa2_Rr + H2np     --> RP_RNa2_Rr_Hr
# a[214]=c[24]*2*RP_RNa2_Rr_Hr*H2np;//222. RP_RNa2_Rr_Hr + H2np  --> RP_RNa2_Rr_Hr2
# a[216]=c[24]*RP_RNa2_Rr_Hr2*H2np;	//224. RP_RNa2_Rr_Hr2 + H2np --> RP_RNa2_Rr_Hr3
# a[218]=c[24]*3*RP_RNa3*H2np;		//226. RP_RNa3 + H2np        --> RP_RNa3_Hr
# a[219]=c[24]*2*RP_RNa3_Hr*H2np;	//227. RP_RNa3_Hr + H2np     --> RP_RNa3_Hr2
# a[220]=c[24]*RP_RNa3_Hr2*H2np;	//228. RP_RNa3_Hr2 + H2np    --> RP_RNa3_Hr3
#
# a[221]=c[23]*RP_Hr;				//229. RP_Hr				--> RP_free + H2np
# a[222]=c[25]*2*RP_Hr2;			//230. RP_Hr2				--> RP_Hr + H2np
# a[223]=c[25]*3*RP_Hr3;			//231. RP_Hr3				--> RP_Hr2 + H2np
# a[226]=c[25]*RP_Rr_Hr;			//234. RP_Rr_Hr				--> RP_Rr + H2np
# a[228]=c[25]*2*RP_Rr_Hr2;			//236. RP_Rr_Hr2			--> RP_Rr_Hr + H2np
# a[230]=c[25]*3*RP_Rr_Hr3;			//238. RP_Rr_Hr3			--> RP_Rr_Hr2 + H2np
# a[233]=c[25]*RP_RNa_Hr;			//241. RP_RNa_Hr			--> RP_RNa + H2np
# a[235]=c[25]*2*RP_RNa_Hr2;		//243. RP_RNa_Hr2			--> RP_RNa_Hr + H2np
# a[237]=c[25]*3*RP_RNa_Hr3;		//245. RP_RNa_Hr3			--> RP_RNa_Hr2 + H2np
# a[240]=c[25]*RP_Rr2_Hr;			//248. RP_Rr2_Hr			--> RP_Rr2 + H2np
# a[242]=c[25]*2*RP_Rr2_Hr2;		//250. RP_Rr2_Hr2			--> RP_Rr2_Hr + H2np
# a[244]=c[25]*3*RP_Rr2_Hr3;		//252. RP_Rr2_Hr3			--> RP_Rr2_Hr2 + H2np
# a[249]=c[25]*RP_RNa_Rr_Hr;		//257. RP_RNa_Rr_Hr			--> RP_RNa_Rr + H2np
# a[252]=c[25]*2*RP_RNa_Rr_Hr2;		//260. RP_RNa_Rr_Hr2		--> RP_RNa_Rr_Hr+ H2np
# a[255]=c[25]*3*RP_RNa_Rr_Hr3;		//263. RP_RNa_Rr_Hr3		--> RP_RNa_Rr_Hr2 + H2np
# a[258]=c[25]*RP_Rr3_Hr;			//266. RP_Rr3_Hr			--> RP_Rr3 + H2np
# a[260]=c[25]*2*RP_Rr3_Hr2;		//268. RP_Rr3_Hr2			--> RP_Rr3_Hr + H2np
# a[262]=c[25]*3*RP_Rr3_Hr3;		//270. RP_Rr3_Hr3			--> RP_Rr3_Hr2 + H2np
# a[267]=c[25]*RP_RNa_Rr2_Hr;		//275. RP_RNa_Rr2_Hr		--> RP_RNa_Rr2 + H2np
# a[270]=c[25]*2*RP_RNa_Rr2_Hr2;	//278. RP_RNa_Rr2_Hr2		--> RP_RNa_Rr2_Hr + H2np
# a[273]=c[25]*3*RP_RNa_Rr2_Hr3;	//281. RP_RNa_Rr2_Hr3		--> RP_RNa_Rr2_Hr2 + H2np
# a[276]=c[25]*RP_RNa2_Hr;			//284. RP_RNa2_Hr			--> RP_RNa2 + H2np
# a[278]=c[25]*2*RP_RNa2_Hr2;		//286. RP_RNa2_Hr2			--> RP_RNa2_Hr + H2np
# a[280]=c[25]*3*RP_RNa2_Hr3;		//288. RP_RNa2_Hr3			--> RP_RNa2_Hr2 + H2np
# a[285]=c[25]*RP_RNa2_Rr_Hr;		//293. RP_RNa2_Rr_Hr 		--> RP_RNa2_Rr + H2np
# a[288]=c[25]*2*RP_RNa2_Rr_Hr2;	//296. RP_RNa2_Rr_Hr2		--> RP_RNa2_Rr_Hr + H2np
# a[291]=c[25]*3*RP_RNa2_Rr_Hr3;	//299. RP_RNa2_Rr_Hr3 		--> RP_RNa2_Rr_Hr2 + H2np
# a[294]=c[25]*RP_RNa3_Hr;			//302. RP_RNa3_Hr 			--> RP_RNa3 + H2np
# a[296]=c[25]*2*RP_RNa3_Hr2;		//304. RP_RNa3_Hr2 			--> RP_RNa3_Hr + H2np
# a[298]=c[25]*3*RP_RNa3_Hr3;		//306. RP_RNa3_Hr3 			--> RP_RNa3_Hr2 + H2np

Rule('dRBPJ_bind_Hes1_dimer', dRBPJ(Hbox=None) + pHes1(pHes1=1, Hbox=None, loc='nuc') % pHes1(pHes1=1, Hbox=None, loc='nuc') |
	dRBPJ(Hbox=[2,3]) % pHes1(pHes1=1, Hbox=2, loc='nuc') % pHes1(pHes1=1, Hbox=3, loc='nuc'), Kn,Kn_r)

# a[151]=c[26]*RP_Rr*Nnp;			//159. RP_Rr + Nnp           --> RP_RNa
# a[154]=c[26]*RP_Rr_Hr*Nnp;		//162. RP_Rr_Hr + Nnp        --> RP_RNa_Hr
# a[157]=c[26]*RP_Rr_Hr2*Nnp;		//165. RP_Rr_Hr2 + Nnp       --> RP_RNa_Hr2
# a[160]=c[26]*RP_Rr_Hr3*Nnp;		//168. RP_Rr_Hr3 + Nnp       --> RP_RNa_Hr3
# a[169]=c[26]*2*RP_Rr2*Nnp;		//177. RP_Rr2 + Nnp          --> RP_RNa_Rr
# a[172]=c[26]*2*RP_Rr2_Hr*Nnp;		//180. RP_Rr2_Hr + Nnp       --> RP_RNa_Rr_Hr
# a[175]=c[26]*2*RP_Rr2_Hr2*Nnp;	//183. RP_Rr2_Hr2 + Nnp      --> RP_RNa_Rr_Hr2
# a[178]=c[26]*2*RP_Rr2_Hr3*Nnp;	//186. RP_Rr2_Hr3 + Nnp      --> RP_RNa_Rr_Hr3 
# a[180]=c[26]*RP_RNa_Rr*Nnp;		//188. RP_RNa_Rr + Nnp       --> RP_RNa2
# a[183]=c[26]*RP_RNa_Rr_Hr*Nnp;	//191. RP_RNa_Rr_Hr + Nnp    --> RP_RNa2_Hr
# a[186]=c[26]*RP_RNa_Rr_Hr2*Nnp;	//194. RP_RNa_Rr_Hr2 + Nnp   --> RP_RNa2_Hr2
# a[189]=c[26]*RP_RNa_Rr_Hr3*Nnp;	//197. RP_RNa_Rr_Hr3 + Nnp   --> RP_RNa2_Hr3
# a[190]=c[26]*3*RP_Rr3*Nnp;		//198. RP_Rr3 + Nnp          --> RP_RNa_Rr2
# a[192]=c[26]*3*RP_Rr3_Hr*Nnp;		//200. RP_Rr3_Hr + Nnp       --> RP_RNa_Rr2_Hr
# a[194]=c[26]*3*RP_Rr3_Hr2*Nnp;	//202. RP_Rr3_Hr2 + Nnp      --> RP_RNa_Rr2_Hr2
# a[196]=c[26]*3*RP_Rr3_Hr3*Nnp;	//204. RP_Rr3_Hr3 + Nnp      --> RP_RNa_Rr2_Hr3
# a[197]=c[26]*2*RP_RNa_Rr2*Nnp;	//205. RP_RNa_Rr2 + Nnp      --> RP_RNa2_Rr
# a[199]=c[26]*2*RP_RNa_Rr2_Hr*Nnp;	//207. RP_RNa_Rr2_Hr + Nnp   --> RP_RNa2_Rr_Hr
# a[201]=c[26]*2*RP_RNa_Rr2_Hr2*Nnp;//209. RP_RNa_Rr2_Hr2 + Nnp  --> RP_RNa2_Rr_Hr2
# a[203]=c[26]*2*RP_RNa_Rr2_Hr3*Nnp;//211. RP_RNa_Rr2_Hr3 + Nnp  --> RP_RNa2_Rr_Hr3
# a[211]=c[26]*RP_RNa2_Rr*Nnp;		//219. RP_RNa2_Rr + Nnp      --> RP_RNa3
# a[213]=c[26]*RP_RNa2_Rr_Hr*Nnp;	//221. RP_RNa2_Rr_Hr + Nnp   --> RP_RNa3_Hr
# a[215]=c[26]*RP_RNa2_Rr_Hr2*Nnp;	//223. RP_RNa2_Rr_Hr2 + Nnp  --> RP_RNa3_Hr2
# a[217]=c[26]*RP_RNa2_Rr_Hr3*Nnp;	//225. RP_RNa2_Rr_Hr3 + Nnp  --> RP_RNa3_Hr3
#
# a[231]=c[25]*RP_RNa;				//239. RP_RNa				--> RP_Rr + Nnp
# a[232]=c[25]*RP_RNa_Hr;			//240. RP_RNa_Hr			--> RP_Rr_Hr + Nnp
# a[234]=c[25]*RP_RNa_Hr2;			//242. RP_RNa_Hr2			--> RP_Rr_Hr2 + Nnp
# a[236]=c[25]*RP_RNa_Hr3;			//244. RP_RNa_Hr3			--> RP_Rr_Hr3 + Nnp
# a[246]=c[25]*RP_RNa_Rr;			//254. RP_RNa_Rr			--> RP_Rr2 + Nnp
# a[248]=c[25]*RP_RNa_Rr_Hr;		//256. RP_RNa_Rr_Hr			--> RP_Rr2_Hr + Nnp
# a[251]=c[25]*RP_RNa_Rr_Hr2;		//259. RP_RNa_Rr_Hr2		--> RP_Rr2_Hr2 + Nnp
# a[254]=c[25]*RP_RNa_Rr_Hr3;		//262. RP_RNa_Rr_Hr3		--> RP_Rr2_Hr3 + Nnp
# a[264]=c[25]*RP_RNa_Rr2;			//272. RP_RNa_Rr2			--> RP_Rr3 + Nnp
# a[266]=c[25]*RP_RNa_Rr2_Hr;		//274. RP_RNa_Rr2_Hr		--> RP_Rr3_Hr + Nnp
# a[269]=c[25]*RP_RNa_Rr2_Hr2;		//277. RP_RNa_Rr2_Hr2		--> RP_Rr3_Hr2 + Nnp
# a[272]=c[25]*RP_RNa_Rr2_Hr3;		//280. RP_RNa_Rr2_Hr3		--> RP_Rr3_Hr3 + Nnp
# a[274]=c[25]*2*RP_RNa2;			//282. RP_RNa2				--> RP_RNa_Rr + Nnp
# a[275]=c[25]*2*RP_RNa2_Hr;		//283. RP_RNa2_Hr			--> RP_RNa_Rr_Hr + Nnp
# a[277]=c[25]*2*RP_RNa2_Hr2;		//285. RP_RNa2_Hr2			--> RP_RNa_Rr_Hr2 + Nnp
# a[279]=c[25]*2*RP_RNa2_Hr3;		//287. RP_RNa2_Hr3			--> RP_RNa_Rr_Hr3 + Nnp
# a[282]=c[25]*2*RP_RNa2_Rr;		//290. RP_RNa2_Rr 			--> RP_RNa_Rr2 + Nnp
# a[284]=c[25]*2*RP_RNa2_Rr_Hr;		//292. RP_RNa2_Rr_Hr 		--> RP_RNa_Rr2_Hr + Nnp
# a[287]=c[25]*2*RP_RNa2_Rr_Hr2;	//295. RP_RNa2_Rr_Hr2		--> RP_RNa_Rr2_Hr2 + Nnp
# a[290]=c[25]*2*RP_RNa2_Rr_Hr3;	//298. RP_RNa2_Rr_Hr3 		--> RP_RNa_Rr2_Hr3 + Nnp
# a[292]=c[25]*3*RP_RNa3;			//300. RP_RNa3 				--> RP_RNa2_Rr + Nnp
# a[293]=c[25]*3*RP_RNa3_Hr;		//301. RP_RNa3_Hr 			--> RP_RNa2_Rr_Hr + Nnp
# a[295]=c[25]*3*RP_RNa3_Hr2;		//303. RP_RNa3_Hr2 			--> RP_RNa2_Rr_Hr2 + Nnp
# a[297]=c[25]*3*RP_RNa3_Hr3;		//305. RP_RNa3_Hr3 			--> RP_RNa2_Rr_Hr3 + Nnp

Rule('dRBPJ_pRBPJ_bind_NICD', dRBPJ(Rbox=1) % pRBPJ(nicd=None, Rbox=1, loc='nuc') + NICD(rbpj=None, loc='nuc') |
	dRBPJ(Rbox=1) % pRBPJ(nicd=2, Rbox=1, loc='nuc') % NICD(rbpj=2, loc='nuc'), Ka,Ka_r)

Observable('delta',Delta())
Observable('Rcm',  mRBPJ())
Observable('Hcm',  mHes1())
Observable('Nm',   mNotch())
Observable('Rcp',  pRBPJ(nicd=None,Rbox=None,loc='cyt'))
Observable('Hcp',  pHes1(pHes1=None,Hbox=None,loc='cyt'))
Observable('Np',   pNotch())
Observable('Hnp',  pHes1(pHes1=None,Hbox=None,loc='nuc'))
Observable('H2np', pHes1(pHes1=ANY, Hbox=None, loc='nuc'), match='species')

### SIMULATION ###

sim = ScipyOdeSimulator(model, verbose=True, integrator_options={'atol' : 1e-8, 'rtol' : 1e-6})
# sim = BngSimulator(model, verbose=5)

tspan1 = np.linspace(0, 750, 751)
x = sim.run(tspan=tspan1)

for obs in ['delta', 'Rcm', 'Hcm', 'Nm', 'Rcp', 'Hcp', 'Np', 'Hnp', 'H2np']:
	plt.figure(obs)
	vol = Vnuc if obs in ['Hnp', 'H2np'] else Vcyt
	plt.plot(tspan1, x.observables[obs]/6.022e23/vol, lw=2, label=obs)
# 	if obs == 'H2np':
# 		plt.plot(tspan1, 0.5*KaHp.value/KrHp.value*x.observables['Hnp']**2, '--k')
	plt.legend()

tspan2 = np.linspace(750, 100000, (100000-750)+1)
initials = x.species[-1]
delta_index = [str(sp) for sp in model.species].index('Delta()')
initials[delta_index] = 1e-3/KfNcp.value # 0.36 #10000
x = sim.run(tspan=tspan2, initials=initials)

for obs in ['delta', 'Rcm', 'Hcm', 'Nm', 'Rcp', 'Hcp', 'Np', 'Hnp', 'H2np']:
	plt.figure(obs)
	vol = Vnuc if obs in ['Hnp', 'H2np'] else Vcyt
	plt.plot(tspan2, x.observables[obs]/6.022e23/vol, lw=2)
# 	if obs == 'H2np':
# 		plt.plot(tspan2, 0.5*KaHp.value/KrHp.value*x.observables['Hnp']**2, '--k')

plt.show()


