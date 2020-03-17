from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, N_A
import sys

Model()

# <editor-fold desc="Constants">
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
# </editor-fold>

# <editor-fold desc="Monomers">

# <editor-fold desc="True species">
#True species, '_T' represents "true" or current concentration
Monomer('Hcm_t')
Monomer('Hcp_t')
Monomer('Hnp_t')#MS:add binding sites or no? As of now, no.
#Monomer('Hdnp_r')#MS: rewrite as dimer or no?
#RBPJ
Monomer('Rcm_t')
Monomer('Rcp_t')
Monomer('Rnp_t')
#Notch/NICD
Monomer('Nm_t')
Monomer('Np_t')
Monomer('Ncp_t')
Monomer('Nnp_t')
# </editor-fold>

get rid:
#Monomer('X_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_10', '_20', '_70']})

# <editor-fold desc="Delay/Sham species">
#Hes1
Fix the time/delay options
Monomer('Hcm_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Hcp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time':['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Hnp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time':['_235', '_43', '_10', '_20', '_21', '_70']})#MS:add binding sites or no? As of now, no.
#Monomer('Hdnp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time':['_235', '_43', '_10', '_20', '_21', '_70']})#MS: rewrite as dimer or no?
#RBPJ
Monomer('Rcm_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Rcp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Rnp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
#Notch/NICD
Monomer('Nm_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Np_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Ncp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
Monomer('Nnp_d', ['delay', 'time'], {'delay': ['y', 'n'], 'time': ['_235', '_43', '_10', '_20', '_21', '_70']})
# </editor-fold>

# </editor-fold>

# <editor-fold desc="Parameters">
###### PARAMETERS#######
Parameter('dummy', 1.0)
#Table 1 Parameters
Parameter('kdHcp', 0.0315)#/min
Parameter('kdHnp', 0.0315)#/min
Parameter('kdHcm', 0.029)#/min
Parameter('kdRcp', 0.00231)#/min
Parameter('kdRcm', 0.0075)#/min
Parameter('kdNp', 0.017)#/min
Parameter('kdNcp', 0.0014)#/min MS: why are there two different values
Parameter('kdNnp', 0.0014)#/min MS: why are there two different values
Parameter('kdNm', 0.0058)#/min
Parameter('Cn', 1.0)#no unit
Parameter('Cr', 1.0)#no unit
Parameter('Cnr', 1.0)#no unit
Parameter('KtrHc', 4.5) #/min
Parameter('KtrRc', 2.5) #/min
Parameter('KtrN', 1) #/min
Parameter('Kr', 3.23e8/(Vol_nucleus*N_A)) #/M, but converted to 1/molecule (in the nucleus)
Parameter('Kn', 2.0e8/(Vol_nucleus*N_A)) #/M, but converted to 1/molecule (in the nucleus)
Parameter('Ka', 1.0e8/(Vol_nucleus*N_A)) #/M, but converted to 1/molecule (in the nucleus)
Parameter('KaHp', 1.0e9/(Vol_nucleus*N_A)) #/M, but converted to 1/molecule (in the nucleus)
Parameter('TmHc', 10.0)#min
Parameter('TpHc', 2.35)#min
Parameter('TmRc', 20) #min
Parameter('TpRc', 2.3)#min
Parameter('TmNc', 70)#min
Parameter('TpNc', 21)#min
Parameter('Vbh', 1.14e-10*Vol_nucleus*N_A)#M/min, but converted to molecules/min (in the nucleus)
Parameter('Vbr', 4.3e-11*Vol_nucleus*N_A) #M/min, but converted to molecules/min (in the nucleus)
Parameter('Vbn', 1.23e-11*Vol_nucleus*N_A)#M/min, but converted to molecules/min (in the nucleus)
Parameter('Vmaxh', 5e-10*Vol_nucleus*N_A) #M/min, but converted to molecules/min (in the nucleus)
Parameter('Vmaxr', 2e-10*Vol_nucleus*N_A)#M/min, but converted to molecules/min (in the nucleus)
Parameter('Vmaxn', 5.5e-11*Vol_nucleus*N_A)#M/min, but converted to molecules/min (in the nucleus)
Parameter('kniHcp', 0.1)#/min
Parameter('kniRcp', 0.1)#/min
Parameter('kniNcp', 0.1)#/min
Parameter('KfNcp', 7.6e7/(Vol_cyto*N_A))#/M/min, but converted to 1/molecules/min (in the cytosol)
#MS: why is KfNCP also listed below
Parameter('rNbox', 0.3)#no unit
Parameter('rR', 0.2)#no unit


#More Parameters, but not in Table 1
look into this:
Parameter('kdRnp', 1.)# MS is this needed? it's not in table 1 !!!
#Parameter('kfNcp', 7.6e7/(N_A*Vol_cyto))#/(M*min)
Parameter('kDelp', 1e-3)#/min
Parameter('tc', 0.5)
# </editor-fold>

 # <editor-fold desc="Initial Parameters">
# <editor-fold desc="Initial Parameters from OP">
###### INITIAL PARAMETERS #######
#initial params left open are those not in Table 2
#these are from Table 2, column 2 (also see Excel sheet)
#They are all in molecules/cell
#Hes1

Parameter('Hcm_t_0', 1)
Parameter('Hcp_t_0', 41)
Parameter('Hnp_t_0', 130)#MS:add binding sites or no? As of now, no.
#Parameter('Hdnp_t_0', )#MS: rewrite as dimer or no?
#RBPJ
Parameter('Rcm_t_0', 1.0)
Parameter('Rcp_t_0', 10)
Parameter('Rnp_t_0', 422)
#Notch/NICD
Parameter('Nm_t_0', 5)
Parameter('Np_t_0', 269)
Parameter('Ncp_t_0', 0.0)
Parameter('Nnp_t_0', 0.0)
# </editor-fold>

??all of them need an initial parameter??

###DELAY INITIAL PARAMETERS###
#only these six delay species were mentioned in the ODE's
#Hes1
Parameter('HcmD_0', 0.0)#TBD
Parameter('HcpD_0', 0.0)
#RBPJ
Parameter('RcmD_0', 0.0)
Parameter('RcpD_0', 0.0)
#Notch
Parameter('NmD_0', 0.0)
Parameter('NpD_0', 0.0)





# </editor-fold>



# <editor-fold desc="Actual Initials">
# , ['delay'], {'delay': ['_now', '_235', '_43', '_10', '_20', '_21', '_70']}
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
# </editor-fold>

#New Observables/Parameters for equations
# Observable('NnpD_obs', X1(delay = '_0'))
# Observable('RnpD_obs', )
# Observable('HdnpD_obs', )
#MS: these shouldn't really be parameters
Parameter('NnpD_obs_10', 1.)
Parameter('RnpD_obs_10', 1.)
Parameter('HdnpD_obs_10', 1.)
#####OBSERVABLES######
Observable('Hnp_obs', Hnp())

######RATE CONSTANT PARAMETERS#######
#Numbers are placeholders
#Hes1 params
Expression('D', 1+ 2*Cr*Kr**2*(1+2*Ka*NnpD_obs+2*Ka**2*NnpD_obs**2)
           *RnpD_obs**2+2*Kr*(RnpD_obs+Ka*NnpD_obs*RnpD_obs)
           +3*HdnpD_obs*Kn
           *(1+2*Cnr*Kr*RnpD_obs
             *(1+Cr*Kr*RnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+Ka
               *(NnpD_obs+2*Cr*Kr*NnpD_obs*RnpD_obs)))
           +6*Cn*HdnpD_obs**2*Kn**2
           *(1+2*Cnr*Kr*RnpD_obs
             *(1+Cr*Kr*RnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+Ka
               *(NnpD_obs+2*Cr*Kr*NnpD_obs*RnpD_obs)))
           +6*Cn**2*HdnpD_obs**3*Kn**3
           *(1+2*Cnr*Kr*RnpD_obs
             *(1+Cr*Kr*RnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+Ka
               *(NnpD_obs+2*Cr*Kr*NnpD_obs*RnpD_obs))))
Expression('RfHcm',
           1/D
           *((1+4*Cr*Ka**2*Kr**2*NnpD_obs**2*RnpD_obs**2
              +2*Kr*RnpD_obs*rR+ 2*Cr*Kr**2*RnpD_obs**2*rR**2
              +2*Ka*Kr*NnpD_obs*RnpD_obs*(1+2*Cr*Kr*RnpD_obs*rR)
              +3*HdnpD_obs*Kn*rNbox
              *(1+2*Cnr*Kr*RnpD_obs
                *(Ka*NnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+rR
                  +2*Cr*Ka*Kr*NnpD_obs*RnpD_obs*rR+Cr*Kr*RnpD_obs*rR**2))
              + 6*Cn*HdnpD_obs**2*Kn**2*rNbox**2
              *(1+2*Cnr*Kr*RnpD_obs
                *(Ka*NnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+rR
                  +2*Cr*Ka*Kr*NnpD_obs*RnpD_obs*rR+Cr*Kr*RnpD_obs*rR**2))
              + 6*Cn**2*HdnpD_obs**3*Kn**3*rNbox**3
              *(1+2*Cnr*Kr*RnpD_obs
                *(Ka*NnpD_obs+2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+rR
                  +2*Cr*Ka*Kr*NnpD_obs*RnpD_obs*rR+Cr*Kr*RnpD_obs*rR**2)))
             *Vbh+2*Ka*Kr*NnpD_obs
             *(1+3*Cnr*HdnpD_obs*Kn*rNbox
               *(1+2*Cn*HdnpD_obs*Kn*rNbox+2*Cn**2*HdnpD_obs**2
                 *Kn**2*rNbox**2))*RnpD_obs
             *(tc+2*Cr*Kr*RnpD_obs*(Ka*NnpD_obs+rR*tc))*Vmaxh))
#M/min are original units, crazy equation, MS

# [I will be converting ALL the other parameters first, so no need to convert within RfHcm]
# Parameter('kf', KaHp.value) #MS: not listed on any table
# Parameter('kr', 1.0) #MS:not listed on any table #unsure of unit
# Parameter('kdH2np', 0.0) #not in code

# Have converted units up to here

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
Expression('C', 1+Kr*(3+2*Ka*NnpD_obs)*RnpD_obs+6*Cr*Kr**2
           *(1+2*Ka*NnpD_obs+2*Ka**2*NnpD_obs**2)*RnpD_obs**2
           +6*Cr**2*Kr**3*(1+3*Ka*NnpD_obs+6*Ka**2*NnpD_obs**2
            +6*Ka**3*NnpD_obs**3)*RnpD_obs**3+3*HdnpD_obs*Kn
           *(1+3*Cnr*Kr*RnpD_obs*(1+2*Cr*Kr*RnpD_obs
            +2*Cr**2*Kr**2*RnpD_obs**2
                +12*Cr**2*Ka**3*Kr**2*NnpD_obs**3*RnpD_obs**2
                +4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs*(1+3*Cr*Kr*RnpD_obs)
                +Ka*NnpD_obs*(1+4*Cr*Kr*RnpD_obs+6*Cr**2*Kr**2*RnpD_obs**2)))
           +6*Cn*HdnpD_obs**2*Kn**2*(1+3*Cnr*Kr*RnpD_obs
                *(1+2*Cr*Kr*RnpD_obs+2*Cr**2*Kr**2*RnpD_obs**2
                  +12*Cr**2*Ka**3*Kr**2*NnpD_obs**3*RnpD_obs**2
                  +4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs*(1+3*Cr*Kr*RnpD_obs)
                  +Ka*NnpD_obs*(1+4*Cr*Kr*RnpD_obs+6*Cr**2*Kr**2*RnpD_obs**2)))
           +6*Cn**2*HdnpD_obs**3*Kn**3*(1+3*Cnr*Kr*RnpD_obs
                *(1+2*Cr*Kr*RnpD_obs+2*Cr**2*Kr**2*RnpD_obs**2
                  +12*Cr**2*Ka**3*Kr**2*NnpD_obs**3*RnpD_obs**2
                  +4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs*(1+3*Cr*Kr*RnpD_obs)
                  +Ka*NnpD_obs
                  *(1+4*Cr*Kr*RnpD_obs+6*Cr**2*Kr**2*RnpD_obs**2))))
Expression('RfRcm',
           1/C
           *((1+36*Cr**2*Ka**3*Kr**3*NnpD_obs**3*RnpD_obs**3)
             +3*Kr*RnpD_obs*rR+6*Cr*Kr**2*RnpD_obs**2*rR**2
             +6*Cr**2*Kr**3*RnpD_obs**3*rR**3
             +12*Cr*Ka**2*Kr**2*NnpD_obs**2*RnpD_obs**2
             *(1+3*Cr*Kr*RnpD_obs*rR)+2*Ka*Kr*NnpD_obs*RnpD_obs
             *(1+3*Cr*Kr*RnpD_obs*rR)**2+3*HdnpD_obs*Kn*rNbox
             *(1+3*Cnr*Kr*RnpD_obs*(12*Cr**2*Ka**3*Kr**2*NnpD_obs**3
           *RnpD_obs**2+4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs
           *(1+3*Cr*Kr*RnpD_obs*rR)
           +rR*(1+2*Cr*Kr*RnpD_obs*rR+2*Cr**2*Kr**2*RnpD_obs**2*rR**2)
            +Ka*NnpD_obs*(1+4*Cr*Kr*RnpD_obs*rR
            +6*Cr**2*Kr**2*RnpD_obs**2*rR**2)))
             +6*Cn*HdnpD_obs**2*Kn**2*rNbox**2
             *(1+3*Cnr*Kr*RnpD_obs
               *(12*Cr**2*Ka**3*Kr**2*NnpD_obs**3*RnpD_obs**2
                 +4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs
                 *(1+3*Cr*Kr*RnpD_obs*rR)+rR
                 *(1+2*Cr*Kr*RnpD_obs*rR+2*Cr**2*Kr**2*RnpD_obs**2*rR**2)
                 +Ka*NnpD_obs
                 *(1+4*Cr*Kr*RnpD_obs*rR+6*Cr**2*Kr**2*RnpD_obs**2*rR**2)))
             +6*Cn**2*HdnpD_obs**3*Kn**3*rNbox**3
             *(1+3*Cnr*Kr*RnpD_obs
               *(12*Cr**2*Ka**3*Kr**2*NnpD_obs**3*RnpD_obs**2
                 +4*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs
                 *(1+3*Cr*Kr*RnpD_obs*rR)+rR
                 *(1+2*Cr*Kr*RnpD_obs*rR
                   +2*Cr**2*Kr**2*RnpD_obs**2*rR**2)+Ka*NnpD_obs
                 *(1+4*Cr*Kr*RnpD_obs*rR+6*Cr**2*Kr**2*RnpD_obs**2*rR**2))))
           *Vbr+Ka*Kr*NnpD_obs*RnpD_obs
           *((2+9*Cnr*HdnpD_obs*Kn*rNbox*(1+2*Cn*HdnpD_obs*Kn*rNbox
                +2*Cn**2*HdnpD_obs**2*Kn**2*rNbox**2))*tc**2+12*Cr*Kr
             *(1+3*Cnr*HdnpD_obs*Kn*rNbox*(1+2*Cn*HdnpD_obs*Kn*rNbox
                 +2*Cn**2*HdnpD_obs**2*Kn**2*rNbox**2))*RnpD_obs*tc
             *(Ka*NnpD_obs+rR*tc)+18*Cr**2*Kr**2
             *(1+3*Cnr*HdnpD_obs*Kn*rNbox*(1+2*Cn*HdnpD_obs*Kn*rNbox
                 +2*Cn**2*HdnpD_obs**2*Kn**2*rNbox**2))*RnpD_obs**2
             *(2*Ka**2*NnpD_obs**2+2*Ka*NnpD_obs*rR*tc+rR**2*tc**2))
           *Vmaxr)#M/min are original units, crazy equation, MS


Rule('RBPJ_ptxa', None >> RcmD(), RfRcm)
Rule('RBPJ_actual_tx', RcmD() >> Rcm(), Parameter('k_Rcm_tx', 1/TmRc.value))
Rule('RBPJ_mrna_deg', Rcm() >> None, kdRcm)
Rule('RBPJ_ptla', RcmD() >> RcmD() + RcpD(), KtrRc)
Rule('RBPJ_actual_tl', RcpD() >> Rcp(), Parameter('k_Rcp_tl', 1/TpRc.value))
Rule('RBPJ_nuclear_import', Rcp() >> Rnp(), kniRcp)
Rule('RBPJ_cyto_degradation', Rcp() >> None, kdRcp)
Rule('RBPJ_nuc_degradation', Rnp() >> None, kdRnp)

#Notch/NICD params
Expression('E', 1+2*Cr*Kr**2*(1+2*Ka*NnpD_obs+2*Ka**2*NnpD_obs**2)
           *RnpD_obs**2+2*Kr*(RnpD_obs+Ka*NnpD_obs*RnpD_obs)
           +HdnpD_obs*Kn*(1+2*Cnr*Kr*RnpD_obs
                          *(1+Ka*NnpD_obs+Cr*Kr*RnpD_obs
                            +2*Cr*Ka*Kr*NnpD_obs*RnpD_obs
                            +2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs)))
Expression('RfNm',
           1/E
           *((1+4*Cr*Ka**2*Kr**2*NnpD_obs**2*RnpD_obs**2
             +2*Kr*RnpD_obs*rR+2*Cr*Kr**2*RnpD_obs**2*rR**2
             +2*Ka*Kr*NnpD_obs*RnpD_obs*(1+2*Cr*Kr*RnpD_obs*rR)
             +HdnpD_obs*Kn*rNbox*(1+2*Cnr*Kr*RnpD_obs*(Ka*NnpD_obs
                    +2*Cr*Ka**2*Kr*NnpD_obs**2*RnpD_obs+rR
                    +2*Cr*Ka*Kr*NnpD_obs*RnpD_obs*rR+Cr*Kr*RnpD_obs*rR**2)))
           *Vbn
           +2*Ka*Kr*NnpD_obs*(1+Cnr*HdnpD_obs*Kn*rNbox)*RnpD_obs
           *(tc+2*Cr*Kr*RnpD_obs*(Ka*NnpD_obs+rR*tc))
           *Vmaxn))
            #M/min are original units, crazy equation, MS


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

######OBSERVABLE#######
#Observable('name_of_observable', species_of_interest)#MS
print('Exit here')
sys.exit()

##set up the delays



###########################################################
### Equations from "MATLAB" code (Berkeley Madonna code)###
###########################################################

#Parameters within RfHcm that are implemented above but had to be changed
#
#Parameters within RfHcm that are implemented above/from the old code
#


