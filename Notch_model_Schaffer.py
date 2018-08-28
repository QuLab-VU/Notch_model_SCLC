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
Monomer('Delta', ['rec']
#binding sites for ligand; got rid of gamme secretase site ??
Monomer('mNotch')
Monomer('pNotch', ['lig']) #should the gamma secretase be added back on??
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
Parameter('mNotch_0', 70)
Parameter('pNotch_0', 8130) #??
Parameter('NICD_cyt_0', 0)
Parameter('NICD_nuc_0', 0)
# get rid of ?? Parameter('NMD_0', fill)
Parameter('mRBPJ_0', 1)
Parameter('pRBPJ_cyt_0', 11)
Parameter('pRBPJ_nuc_0', 447)
Parameter('mHes1_0', 1)
Parameter('pHes1_cyt_0', 35)
Parameter('pHes1_nuc_0', 109)
Parameter('dNotch_0', 2)
Parameter('dHes1_0', 2)
Parameter('dRBPJ_0', 2)
Parameter('dimHes1_nuc_0', 101)


#HOw much should I specify in the initials ??
#do i put location or no ??
Initial(Delta(rec=None), Delta_0)
Initial(mNotch(), mNotch_0)
Initial(pNotch(lig=None), pNotch_0)
Initial(NICD(rbpj=None, loc=cyt), NICD_cyt_0)
Initial(NICD(rbpj=None, loc=nuc), NICD_nuc_0)
#Initial(NMD(), NMD_0)
Initial(mRBPJ(), mRBPJ_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc=cyt), pRBPJ_cyt_0)
Initial(pRBPJ(nicd=None, Rbox=None, loc=nuc), pRBPJ_nuc_0)
Initial(mHes1(), mHES1_0)
Initial(pHes1(pHes1=None, Hbox=None, loc=cyt), pHES1_cyt_0)
Initial(pHes1(pHes1=None, Hbox=None, loc=nuc), pHES1_nuc_0)
Initial(dNotch(Hbox1a=None, Hbox1b=None, Rbox1=None, Rbox2=None), dNotch_0)
Initial(dHes1(Hbox1a=None, Hbox1b=None, Hbox2a=None, Hbox2b=None, Hbox3a=None, Hbox3b=None, Rbox1=None, Rbox2=None), mHes1_0)
Initial(dRBPJ(Hbox1a=None, Hbox1b=None, Hbox2a=None, Hbox2b=None, Hbox3a=None, Hbox3b=None, Rbox1=None, Rbox2=None, Rbox3=None), pRBPJ_0)
Initial(pHes1(pHes1=1, Hbox=None, loc=nuc) % pHes1(pHes1=1, Hbox= None, loc=nuc), dimHes1_nuc_0)


##
#Add section for times delays ??
###


#Degradation parameters

#	//Units here are per minute

#degradation constant of full length Notch1
Parameter('KdNp', 0.017)
#//degradation HCP
Parameter('KdHcp', 0.0315)
#//deg Hnp
Parameter('KdHnp', 0.0315)
#// deg mHes1
Parameter('KdHcm', 0.029)
#//disable Hes1 dimer degradation <-- I didn't write that -MS
Parameter('KdH2np' 0)
Parameter('KdNcp', 0.00385)
Parameter('KdNnp', 0.00385)
Parameter('KdNm', 0.00035)
Parameter('KdRcp', 0.00231)
Parameter('KdRnp', 0.00231)
Parameter('KdRcm', 0.0075)

#	//units here are per minute
#protein translation rates
Parameter('KtrHc', 4.5)
Parameter('KtrN' 2)
Parameter('KtrRc' 3.2)

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
	float rRbox=0.2;
	float tc=0.5;