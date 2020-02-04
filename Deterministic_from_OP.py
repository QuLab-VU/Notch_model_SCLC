from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, N_A
import sys

Model()

r_cell= (10./2.)/1e6 #m
r_nucleus = (5./2.)/1e6 #m

Vol_cell = (4/3)*pi*(r_cell)**3*1000
Vol_nucleus = (4/3)*pi*(r_nucleus)**3*1000
Vol_cyto = Vol_cell-Vol_nucleus

print('Volumes')
print(Vol_cell)
print(Vol_nucleus)
print(Vol_cyto)
print(pi)
print(N_A)




######MONOMERS#######

#Hes1
Monomer('HcmD')
Monomer('Hcm')
Monomer('HcpD')
Monomer('Hcp')
Monomer('Hnp')#MS:add binding sites or no? As of now, no.
#Monomer('Hdnp')#MS: rewrite as dimer or no?
#RBPJ
Monomer('RcmD')
Monomer('Rcm')
Monomer('RcpD')
Monomer('Rcp')
Monomer('Rnp')
#Notch/NICD
Monomer('NmD')
Monomer('Nm')
Monomer('NpD')
Monomer('Np')
Monomer('Ncp')
Monomer('Nnp')

######INITIAL PARAMETERS#######
#initial params left open are those not in Table 2
#Hes1
Parameter('HcmD_0', 0.0)#TBD
Parameter('Hcm_0', 1)
Parameter('HcpD_0', 0.0)
Parameter('Hcp_0', 41)
Parameter('Hnp_0', 130)#MS:add binding sites or no? As of now, no.
#Parameter('Hdnp_0', )#MS: rewrite as dimer or no?
#RBPJ
Parameter('RcmD_0', 0.0)
Parameter('Rcm_0', 1.0)
Parameter('RcpD_0', 0.0)
Parameter('Rcp_0', 10)
Parameter('Rnp_0', 422)
#Notch/NICD
Parameter('NmD_0', 0.0)
Parameter('Nm_0', 5)
Parameter('NpD_0', 0.0)
Parameter('Np_0', 269)
Parameter('Ncp_0', 0.0)
Parameter('Nnp_0', 0.0)

######INITIALS######
#Hes1
Initial(HcmD(), HcmD_0)
Initial(Hcm(), Hcm_0)
Initial(HcpD(), HcpD_0)
Initial(Hcp(), Hcp_0)
Initial(Hnp(), Hnp_0)#MS:add binding sites or no? As of now, no.
#Initial(Hdnp, Hdnp_0)
#RBPJ
Initial(RcmD(), RcmD_0)
Initial(Rcm(), Rcm_0)
Initial(RcpD(), RcpD_0)
Initial(Rcp(), Rcp_0)
Initial(Rnp(), Rnp_0)
#Notch/NICD
Initial(NmD(), NmD_0)
Initial(Nm(), Nm_0)
Initial(NpD(), NpD_0)
Initial(Np(), Np_0)
Initial(Ncp(), Ncp_0)
Initial(Nnp(), Nnp_0)

#####OBSERVABLES######
Observable('Hnp_obs', Hnp())

######RATE CONSTANT PARAMETERS#######
#Numbers are placeholders
#Hes1 params
Parameter('dummy', 1.0)
Expression('RfHcm', dummy)#M/min are original units, crazy equation, MS
# [I will be converting ALL the other parameters first, so no need to convert within RfHcm
Parameter('TmHc', 1.0)#min
Parameter('kdHcm', 1.0)#/min
Parameter('KtrHc', 1.0)#/min
Parameter('TpHc', 1.0)#min
Parameter('kniHcp', 1.0)#/min
Parameter('kdHcp', 1.0)#/min
Parameter('kdHnp', 1.0)#/min
Parameter('KaHp', 1.0e9/(N_A*Vol_nucleus)) #has to be large value in order #/M
#Parameter('kf', KaHp.value) #MS: not listed on any table
#Parameter('kr', 1.0) #MS:not listed on any table #unsure of unit
#Parameter('kdH2np', 0.0) #not in code

#ptxa = pre-transcriptional assembly
#actual_tx = actual transcription
#mrna_deg = mRNA degradation
#ptla = pre-translational assembly
#actual_tl = actual translation
#nuclear_import
#cyto_degradation = protein degradation of in cytosol

Rule('Hes1_ptxa', None >> HcmD(), RfHcm)
Rule('Hes1_actual_tx', HcmD() >> Hcm(), Parameter('k_Hcm_tx', 1.0/TmHc.value))
Rule('Hes1_mrna_deg', Hcm() >> None, kdHcm)
Rule('Hes1_ptla', HcmD() >> HcmD() + HcpD(), KtrHc)
Rule('Hes1_actual_tl', HcpD() >> Hcp(), Parameter('k_Hcp_tl', 1.0/TpHc.value))
#Rule('Hes1_nuclear_import', Hcp() >> Hnp(), Parameter('k_nimport_vol', kniHcp.value*7.0))
Rule('Hes1_nuclear_import_pt1', Hcp() >> None, kniHcp)
Rule('Hes1_nuclear_import_pt2', Hcp() >> Hcp() + Hnp(), Parameter('k_nimport_vol', kniHcp.value))
Rule('Hes1_cyto_degradation', Hcp() >> None, kdHcp)
Rule('Hes1_nuc_degradation', Hnp() >> None, kdHnp)
#Rule('Hes1_dimerization', Hnp() + Hnp() | Hdnp(), kf, kr)
Expression('Hdnp_exp', KaHp*Hnp_obs**2)
#Rule('Hes1_dim_degradation', Hdnp() >> None, kdH2np)#not in code


#RBPJ params
Expression('RfRcm', dummy)#M/min are original units, crazy equation, MS
Parameter('TmRc', 1.0)#min
Parameter('kdRcm', 1.0)#/min
Parameter('KtrRc', 1.0)#/min
Parameter('TpRc', 1.0)#min
Parameter('kniRcp', 1.0)#/min
Parameter('kdRcp', 1.0)#/min
Parameter('kdRnp', 1.0)#/min

Rule('RBPJ_ptxa', None >> RcmD(), RfRcm)
Rule('RBPJ_actual_tx', RcmD() >> Rcm(), Parameter('k_Rcm_tx', 1/TmRc.value))
Rule('RBPJ_mrna_deg', Rcm() >> None, kdRcm)
Rule('RBPJ_ptla', RcmD() >> RcmD() + RcpD(), KtrRc)
Rule('RBPJ_actual_tl', RcpD() >> Rcp(), Parameter('k_Rcp_tl', 1/TpRc.value))
Rule('RBPJ_nuclear_import', Rcp() >> Rnp(), kniRcp)
Rule('RBPJ_cyto_degradation', Rcp() >> None, kdRcp)
Rule('RBPJ_nuc_degradation', Rnp() >> None, kdRnp)

#Notch/NICD params
Expression('RfNm', dummy)#M/min are original units, crazy equation, MS
Parameter('TmNc', 1.0)#min
Parameter('kdNm', 1.0)#/min
Parameter('KtrN', 1.0)#/min
Parameter('TpNc', 1.0)#min
Parameter('kdNp', 1.0)#/min
#Parameter('kfNcp', 7.6e7/(N_A*Vol_cyto))#/(M*min)
Parameter('kDelp', 1e-3)#/min
Parameter('kniNcp', 1.0)#/min
Parameter('kdNcp', 1.0)#/min
Parameter('kdNnp', 1.0)#/min

Rule('Notch1_ptxa', None >> NmD(), RfNm)
Rule('Notch1_actual_tx', RcmD() >> Rcm(), Parameter('k_Ncm_tx', 1/TmNc.value))
Rule('Notch1_mrna_deg', Nm() >> None, kdNm)
Rule('Notch1_ptla', NmD() >> NmD() + NpD(), KtrN)
Rule('Notch1_actual_tl', NpD() >> Np(), Parameter('k_Ncp_tx', 1/TpNc.value))
Rule('Notch1_membrane_degradation', Np() >> None, kdNp)
Rule('Notch1_cleavage', Np() >> Ncp(), kDelp)
Rule('Notch1_nuclear_import', Ncp() >> Nnp(), kniNcp)
Rule('Notch1_cyto_degradation', Ncp() >> None, kdNcp)
Rule('Notch1_nuc_degradation', Nnp() >> None, kdNnp)
sys.exit()
######OBSERVABLE#######
Observable('name_of_observable', species_of_interest)#MS

###########################################################
### Equations from "MATLAB" code (Berkeley Madonna code)###
###########################################################

#Parameters within RfHcm that are implemented above but had to be changed
#
#Parameters within RfHcm that are implemented above/from the old code
#

#New Parameters

Parameter('Cr', 1.0)#
Parameter('Ka', 1e8) #
Parameter('Kr', 3.23e8)#
Parameter('rR', 0.2)#
Parameter('Kn', 2e8)#
Parameter('rNbox', 0.3)#
Parameter('Cnr', 1.0)#
Parameter('Cn', 1.0)#
Parameter('', )#
Parameter('', )#



#New Observables

Observable('NnpD_obs', )
Observable('RnpD_obs', )
Observable('HdnpD_obs', )

#######EQUATION######
RfHcm = 1/D(alpha*Vb+2*Ka*Kr*NnpD*(1+3*Cnr*HdnpD*Kn*rNbox*(1+2*Cn*HdnpD*Kn*rNbox+2*Cn2*HdnpD2*Kn2*rNbox2))*RnpD*(tc+2*Cr*Kr*RnpD*(Ka*NnpD+rR*tc))*Vmax)
# alpha = alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7 + alpha8
alpha1 = 1
alpha2 = 4*Cr*Ka**2*Kr**2*NnpD**2*RnpD**2 #MS: check these species/observables
alpha3 = 2*Kr*RnpD*rR
alpha4 = 2*Cr*Kr**2*RnpD**2*rR**2
alpha5 = 2*Ka*Kr*NnpD*RnpD*(1+2*Cr*Kr*RnpD*rR)
alpha6 = 3*HdnpD*Kn*rNbox*(1+2*Cnr*Kr*RnpD*(Ka*NnpD+2*Cr*Ka**2*Kr*NnpD**2*RnpD+rR+2*Cr*Ka*Kr*NnpD*RnpD*rR+Cr*Kr*RnpD*rR**2))
alpha7 = 6*Cn*HdnpD**2*Kn**2*rNbox**2*(1+2*Cnr*Kr*RnpD*(Ka*NnpD+2*Cr*Ka**2*Kr*NnpD**2*RnpD+rR+2*Cr*Ka*Kr*NnpD*RnpD*rR+Cr*Kr*RnpD*rR**2))
#alpha8 = 6*Cn**2*HdnpD**3*Kn**3*rNbox**3*(1+2*Cnr*Kr*RnpD*(Ka*NnpD+2*Cr*Ka**2*Kr*NnpD**2*RnpD+rR+2*Cr*Ka*Kr*NnpD*RnpD*rR+Cr*Kr*RnpD*rR**2)))









#####Copied from LH version#######


Parameter('KdHes1', 0.0315)  # //degradation Hes1 protein

#new parameters -q
# Degradation constant of Hes1 protein (min^-1) kdHcp, kdHnp 0.0315 [30]
Parameter('KdHcp', 0.0315)  # //degradation Hes1 protein
# Degradation constant of Hes1 protein (min^-1) kdHcp, kdHnp 0.0315 [30]
Parameter('KdHnp', 0.0315)  # //degradation Hes1 protein
# Degradation constant of Hes1 mRNA (min^-1) kdHcm 0.029 [30]
Parameter('KdHcm', 0.029)  # // deg Hes1 mRNA
# Degradation constant of RBP-Jk protein (min^-1) kdRcp 0.00231 [82]
Parameter('KdRcp', 0.00231)  # degradation RBPJ protein
# Degradation constant of RBP-JK mRNA (min^-1) kdRcm 0.0075 [82]
Parameter('KdRcm', 0.0075)  # degradation RBPJ mRNA
# Degradation constant of full-length Notch1 protein (min^-1) kdNp 0.017 [54], Text
Parameter('KdNp', 0.017)  # degradation Notch protein
# Degradation constant of NICD protein (min^-1) kdNcp,kdNnp 0.0014 or 0.00385 [55], Text
Parameter('KdNcp', 0.0014)#degradation constant of NICD protein (min^-1)    #MS:this was not in the code from LH
Parameter('KdNnp', 0.0014)#degradation constant of NICD protein (min^-1)    #MS:this was not in the code from LH
Parameter('KdNm', 0.0058)  # 0.00035) 		#degradation Notch mRNA
# Cooperativity factor for Hes1-DNA binding Cn 1 Text
# Cooperativity factor for RBP-Jk DNA binding Cr 1 Text
# Cooperativity factor for RBP-Jk Hes1 DNA binding Cnr 1 Text
# Rate of protein translation from Hes1 mRNA (min^-1) KtrHc 4.5 [45], Text
Parameter('KtrHc', 4.5)
# Rate of protein translation from RBP-Jk mRNA (min^-1) KtrRc 2.5 Text
Parameter('KtrRc', 2.5)  # 3.2)
# Rate of protein translation from Notch1 mRNA (min^-1) KtrN 1 Text
Parameter('KtrN', 1)  # 2)








