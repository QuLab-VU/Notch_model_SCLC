from pysb.core import *

from pysb.integrate import *

import matplotlib.pyplot as plt

import numpy as np

from pysb.integrate import ScipyOdeSimulator as SOS

from pysb.bng import generate_equations

Model()

Monomer('Delta', ['r', 'loc'], {'loc': ['cell_1', 'cell_2']})
Monomer('Notch', ['l', 'ysc'])
#Gamma secretase
Monomer('YSC', ['lr'])
Monomer('NICD', ['csl'])
#Notch membrane domain
Monomer('NMD')
Monomer('NMD', ['lr'])
Monomer('ND_endo')
Monomer('N_endo')
Monomer('CSL', ['nicd'])
Monomer('MAM', ['nicd'])
Monomer('mHes1')
Monomer('pHES1')
Monomer('pASCL1')
Monomer('pNEUROD1')
Monomer('pSLUG')
Monomer('pNGN2')

# print model.monomers

fill = 1

#Default forward, reverse, and catalytic rates
KF = 1e-6
KR = 1e-3
KC = 1

Parameter('Delta_0', fill)
Parameter('Notch_0', fill)
Parameter('YSC_0', fill)
Parameter('NICD_0', fill)
Parameter('NMD_0', fill)
Parameter('ND_endo_0', fill)
Parameter('N_endo_0', fill)
Parameter('CSL_0', fill)
Parameter('MAM_0', fill)
Parameter('mHES1_0', fill)
Parameter('pHES1_0', fill)
Parameter('pASCL1_0', fill)
Parameter('pNEUROD1_0', fill)
Parameter('pSLUG_0', fill)
Parameter('pNGN2_0', fill)

# Compartments

Parameter('VolECM', 10000)
Parameter('VappCM', 100)
Parameter('VCyto1', 4000)
Parameter('VCyto2', 4000)
Parameter('VappNM', 40)
Parameter('VNuc1', 500)
Parameter('VNuc2', 500)

#Two different cell membranes, cytosol volumes, nuclear membranes, and nuclear volumes
Compartment('ECM', None, 3, VolECM)
Compartment('CellMem1', ECM, 2, VappCM)
Compartment('CellMem2', ECM, 2, VappCM)
Compartment('Cyto1', CellMem1, 3, VCyto1)
Compartment('Cyto2', CellMem2, 3, VCyto2)
Compartment('NM1', Cyto1, 2, VappNM)
Compartment('NM2', Cyto2, 2, VappNM)
Compartment('Nuc1', NM1, 3, VNuc1)
Compartment('Nuc2', NM2, 3, VNuc2)



Initial(Delta(r=None, loc='cell_1') ** ECM, Delta_0)
Initial(Delta(r=None, loc='cell_2') ** ECM, Delta_0)

Initial(Notch(), Notch_0)
Initial(Notch(l=None, ysc=None) ** CellMem1, Notch_0)
Initial(R(l=None, ysc=None) ** CellMem2, Notch_0)
Initial(YSC(lr=None) ** CellMem1, YSC_0)
Initial(YSC(lr=None) ** CellMem2, YSC_0)

Initial(NICD(), NICD_0)
Initial(NMD(), NMD_0)
Initial(ND_endo(), ND_endo_0)
Initial(N_endo(), N_endo_0)
Initial(CSL(nicd=None), CSL_0)
Initial(MAM(nicd=None), MAM_0)
Initial(mHES1(), mHES1_0)
Initial(pHES1(), pHES1_0)
Initial(pASCL1, pASCL1_0)
Initial(pNEUROD1(), pNEUROD1_0)
Initial(pSLUG(), pSLUG_0)
Initial(pNGN2(), pNGN2_0)

#Part 1
#rate Parameters
Parameter('kNDf', fill)
Parameter('kNDr', fill)
Parameter('kNDYf', fill)
Parameter('kNDYr', fill)
Parameter('kYSCICD', fill)

#Rules
Rule()

# #the delta in one cell binds to the notch in another cell, forming a complex
Rule('N_1_binds_D_2', Notch(r=None, loc='cell_1') ** ECM + Delta(l=None, ysc=None) ** CellMem2 |
     Notch(r=1, loc='cell_1') ** ECM % Delta(l=1, ysc=None) ** CellMem2, kNDf, kNDr)

#This is the reverse polarity binding. Should we only do a one-way ligand binding for the first run???
Rule('N_2_binds_D_1', Notch(r=None, loc='cell_2') ** ECM + R(l=None, ysc=None) ** CellMem1 |
     Notch(r=1, loc='cell_2') ** ECM % Delta(l=1, ysc=None) ** CellMem1, kNDf, kNDr)

#the Gamma secretase then binds the receptor complex
#Should l really be equal to ANY??? should the ND complex be indicated???
Rule('YSC_binding', YSC(lr=None) + Notch(l=ANY, ysc=None) |
     YSC(lr=1) % Notch(l=ANY, ysc=1), kNDYf, kNDYr)

#the complex then gets cleaved into NICD and NMD (Simplifying assumption, combining the ADAM and gamma secretase cleavage steps)
Rule('NICD_formation_1', Delta(r=1, loc='cell_1') ** ECM % Notch(l=1, ysc=2) ** CellMem2 % YSC(lr=2) ** CellMem2
     >> NICD(csl=None) ** Cyto2 + YSC(lr=None) ** CellMem2 + ND_endo() ** Cyto1 + N_endo() ** Cyto2, kYSCICD)

Rule('NICD_formation_2', Delta(r=1, loc='cell_2') ** ECM % Notch(l=1, ysc=2) ** CellMem1 % YSC(lr=2) ** CellMem1
     >> NICD(csl=None) ** Cyto1 + YSC(lr=None) ** CellMem1 + ND_endo() ** Cyto2 + N_endo() ** Cyto1, kYSCICD)

#Part 2
Parameter('kCytNucf', fill)
Parameter('kCytNucr', fill)
Parameter('kNICDCSLf', fill)
Parameter('kNICDCSLr', fill)

#the NICD translocates to the nucleus

Rule('NICD_cyto_to_nuc_1', NICD(csl=None) ** Cyto1 | NICD(csl=None) ** Nuc1, kCytNucf, kCytNucr)
Rule('NICD_cyto_to_nuc_2', NICD(csl=None) ** Cyto2 | NICD(csl=None) ** Nuc2, kCytNucf, kCytNucr)

#the NICD binds with CSL
Rule('ICD_CSL_bind', NICD(csl=None) + CSL(nicd=None) | NICD(csl=1) % CSL(nicd=1), kNICDCSLf, kNICDCSLr)


#the NICD*CSL complex binds with MAM1

#can the reverse order happen? Or does that even matter?


##synthesis of Hes1

Parameter('ktrcxn', 1e-1)
Parameter('ktrlxn', 1e-1)
Parameter('kmrnadeg', 2.3e-2)
Parameter('km_n_to_c', 0.5)
Parameter('ktf_c_to_n', 0.3)
#Kprotsynth = rate constant for protein production
Parameter('Kprotsynth', ktrcxn.value * ktrlxn.value * km_n_to_c.value / kmrnadeg.value)
Parameter('kgICf', 3e-1)
Parameter('kgICr', 2e-4)
Parameter('kgpHf', 3e-1)
Parameter('kgpHr', 2e-4)
Parameter('kgICpHf', 1e-1)
Parameter('kgICpHr', 2e-4)
Parameter('kgpHICf', 1e-1)
Parameter('kgpHICr', 2e-4)
Parameter('kphes1_deg', 1e-2)


Expression('kHes1_1',
           ktrcxn * (kgICpHr / kgICpHf) * CSL_on_1 /
           (

               (kgICr / kgICf) * (kgICpHr / kgICpHf) +

               (kgICpHr / kgICpHf) * CSL_on_1 +

               (kgpHICr / kgpHICf) * prot_Hes1_1 +

               prot_Hes1_1 * CSL_on_1

           ))



Expression('kHes1_2',
           ktrcxn * (kgICpHr / kgICpHf) * CSL_on_2 /
           (

               (kgICr / kgICf) * (kgICpHr / kgICpHf) +

               (kgICpHr / kgICpHf) * CSL_on_2 +

               (kgpHICr / kgpHICf) * prot_Hes1_2 +

               prot_Hes1_2 * CSL_on_2

           ))

Rule('mHes1_synth_1', None >> mHes1() ** Cyto1, kHes1_1)
Rule('mHes1_synth_2', None >> mHes1() ** Cyto2, kHes1_2)

Parameter('km_to_p', ktrlxn.value * km_n_to_c.value * ktf_c_to_n.value)

Rule('pHes1_synth_1', mHes1() ** Cyto1 >> pHES1() ** Nuc1, km_to_p)
Rule('pHes1_synth_2', mHes1() ** Cyto2 >> pHES1() ** Nuc2, km_to_p)

Rule('mHes1_decay', mHes1() >> None, kmrnadeg)
Rule('pHES1_decay', pHES1() >> None, kphes1_deg)

## SLUG synthesis

Parameter('kgSlICf', 1e-2)
Parameter('kgSlICr', 2e-4)

Parameter('n', 1)
Expression('kSlug_1', Kprotsynth * CSL_on_1 ** n / ((kgSlICr / kgSlICf) ** n + CSL_on_1 ** n))
Expression('kSlug_2', Kprotsynth * CSL_on_2 ** n / ((kgSlICr / kgSlICf) ** n + CSL_on_2 ** n))
#

Rule('synth_Slug_1', None >> pSLUG() ** Cyto1, kSlug_1)
Rule('synth_Slug_2', None >> pSLUG() ** Cyto2, kSlug_2)

# ASCL1 synthesis

Parameter('kgAs1pHf', 1e-2)
Parameter('kgAs1pHr', 2e-5)

Expression('kAscl1_1', Kprotsynth * (kgAs1pHr / kgAs1pHf) ** n / ((kgAs1pHr / kgAs1pHf) ** n + prot_Hes1_1 ** n))
Expression('kAscl1_2', Kprotsynth * (kgAs1pHr / kgAs1pHf) ** n / ((kgAs1pHr / kgAs1pHf) ** n + prot_Hes1_2 ** n))

#

Rule('synth_Ascl1_1', None >> pASCL1() ** Cyto1, kAscl1_1)
Rule('synth_Ascl1_2', None >> pASCL1() ** Cyto2, kAscl1_2)

# NGN2 Synthesis

Parameter('kgNg2pHf', 1e-2)
Parameter('kgNg2pHr', 2e-5)

#
Expression('kNgn2_1', Kprotsynth * ktf_c_to_n / ((kgNg2pHr / kgNg2pHf) ** n + prot_Hes1_1 ** n))
Expression('kNgn2_2', Kprotsynth * ktf_c_to_n / ((kgNg2pHr / kgNg2pHf) ** n + prot_Hes1_2 ** n))
#

Rule('synth_Ngn2_1', None >> pNGN2() ** Nuc1, kNgn2_1)
Rule('synth_Ngn2_2', None >> pNGN2() ** Nuc2, kNgn2_2)

#----------------------------
#----------------------------
#-------Delete this section when done


Rule('mHes1_synth_1', None >> mHes1() ** Cyto1, kHes1_1)
Rule('mHes1_decay', mHes1() >> None, kmrnadeg)

Parameter('km_to_p', ktrlxn.value * km_n_to_c.value * ktf_c_to_n.value)

Rule('pHes1_synth_1', mHes1() ** Cyto1 >> pHES1() ** Cyto1, km_to_p)
Rule('pHes1_nuc_import_1', pHES1() ** Cyto1 >> pHES1() ** Nuc1, kn_to_c)
Rule('pHES1_decay', pHES1() >> None, kphes1_deg)


#----------------------------
#----------------------------
#-------