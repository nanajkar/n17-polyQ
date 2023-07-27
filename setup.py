import scipy.spatial
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns;
from MDAnalysis import *
import MDAnalysis as mda
import argparse
import sys


## Set the baseline simulation temperature
temp = 380

## Change the mdp files to reflect system components

with open('em.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein \n'
content[125] = 'tau_t\t\t\t= 1.0 \n'
content[126] = 'ref_t\t\t\t= %d  \n' %(temp)

f = open('em.mdp','w')
for s in content:
    f.write(s)
f.close()

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein \n'
content[125] = 'tau_t\t\t\t= 1.0 \n'
content[126] = 'ref_t\t\t\t= %d  \n' %(temp)

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()

## Minimize and equilibriate the protein in vaccum

os.system("gmx grompp -f em.mdp -p newprotein.top -c protein.gro -o em_md/proem1.tpr -maxwarn 1")
os.system("gmx mdrun -s em_md/proem1.tpr -c em_md/proem1.gro -tableb ../tables/table_d0.xvg ../tables/table_a*.xvg")


os.system("gmx grompp -f npt_eq.mdp -p newprotein.top -c em_md/proem1.gro -o em_md/proeq1.tpr -maxwarn 2")
os.system("gmx mdrun -s em_md/proeq1.tpr -c proeq1.gro -tableb ../tables/table_d0.xvg ../tables/table_a*.xvg")

## proeq1.gro is the equilibriated protein structure that can be added to the waterbox 



seq_to_resname = { \
	'A': 'ALA',
	'G': 'GLY',
	'S': 'SER',
	'I': 'ILE',
	'L': 'LEU',
	'V': 'VAL',
	'M': 'MET',
	'D': 'ASP',
	'Q': 'GLN',
	'N': 'ASN',
	'P': 'PRO',
	'T': 'THR',
	'K': 'LYS',
	'R': 'ARG',
	'E': 'GLU',
	'W': 'TRP',
	'F': 'PHE',
	'Y': 'TYR',
	'H': 'HIS',
	'Rsd': 'RSD'
}
seq_to_resno = { \
	'A': 3,
	'G': 3,
	'S': 6,
	'I': 4,
	'L': 4,
	'V': 4,
	'M': 4,
	'D': 5,
	'Q': 6,
	'N': 6,
	'P': 4,
	'T': 6,
	'K': 5,
	'R': 5,
	'E': 5,
	'W': 9,
	'F': 6,
	'Y': 8,
	'H': 7,
	'Rsd': 6
}
f = open('./gen_conf/seq.txt', 'r')
A = f.readlines()
f.close()
string = A[0].strip()
Seq = [seq_to_resname[i] for i in string]
Seqno = [seq_to_resno[i] for i in string]


### load in the water box
u = Universe('water_12.gro')

#Store initial dimensions
boxz = u.dimensions[2]
boxx = u.dimensions[0]
boxy = u.dimensions[1]


################### Add 6 peptides to water_12.gro

nx = 2
ny = 2
x = np.linspace(0.2*boxx/10, 0.8*(boxx)/10.0, nx)
y = np.linspace(0.2*boxx/10, 0.8*(boxx)/10.0, ny)
xv, yv = np.meshgrid(x, y)
Peploc = np.array([zip(xv[i], yv[i]) for i in range(nx)])
PROT=[]
a=0
for i in Peploc:
	b=0
	a+=1
	for j in i:
		b+=1
		os.system("gmx editconf -f proeq1.gro -o Nsvs_mod_" + str(a)+str(b)+".gro -center "+str(j[0])+" "+str(j[1])+" "+str(0.5*boxz/10.0)+" -box "+str(boxx/10.0)+" "+str(boxy/10.0)+" "+str(boxz/10.0))
		newgro = "Nsvs_mod_" + str(a)+str(b)+".gro"
		PROT.append(Universe(newgro).atoms)


os.system("gmx editconf -f proeq1.gro -o 2prot_mod1.gro -center "+str(0.5*boxx/10.0)+" "+str(0.5*boxy/10.0)+" "+str(0.25*boxz/10.0)+" -box "+str(boxx/10.0)+" "+str(boxy/10.0)+" "+str(boxz/10.0))
PROT.append(Universe("2prot_mod1.gro").atoms)

os.system("gmx editconf -f proeq1.gro -o 2prot_mod2.gro -center "+str(0.5*boxx/10.0)+" "+str(0.5*boxy/10.0)+" "+str(0.75*boxz/10.0)+" -box "+str(boxx/10.0)+" "+str(boxy/10.0)+" "+str(boxz/10.0))
PROT.append(Universe("2prot_mod2.gro").atoms)

proteins = mda.Merge(PROT[0])
for p in PROT[1:]:    
    proteins = mda.Merge(proteins.atoms, p)
wat = u.select_atoms("resname PW")

combined = mda.Merge( proteins.atoms, wat.atoms)   
no_overlap = combined.select_atoms("same resid as (not around 6 name BB BBp BBm S1 S1p S1m S2 S2p S2m)")
no_overlap.write("overlap_fix.gro")
with open("overlap_fix.gro") as f:
	content = f.readlines()
content[-1] = '   ' + str(boxx/10.0) + '   ' + str(boxy/10.0) + '   '+ str(boxz/10.0)
f = open("pepappend.gro", "w")
for s in content:
    f.write(s)
f.close()
os.system("gmx genconf -f pepappend.gro -o pepappend.gro")

### pepappend.gro contains solvated protein. 


#############################
### Now update mdp files


with open('em.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein PW \n'
content[125] = 'tau_t\t\t\t= 1.0 1.0 \n'
content[126] = 'ref_t\t\t\t= %d %d  \n' %(temp,temp)

f = open('em.mdp','w')
for s in content:
    f.write(s)
f.close()

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein PW \n'
content[125] = 'tau_t\t\t\t= 1.0 1.0 \n'
content[126] = 'ref_t\t\t\t= %d %d  \n' %(temp,temp)

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()

#############################
### Update topology files

u = Universe("pepappend.gro")

pep = u.select_atoms("name BB")
wat = u.select_atoms("resname PW")
pepno = int(pep.n_atoms/len(Seq))
watno = int(wat.n_atoms/3)
f_top = open("em.top", "r")
with f_top as f:
	content = f.readlines()

content_mod = content[0:14]
content_mod.append("Protein"+"\t"+str(pepno)+"\n")
content_mod.append("PW"+"\t"+str(watno)+"\n")
f = open("em.top", "w")
for s in content_mod:
    f.write(s)
f.close()
content_mod[2] = '#include "itps/water.md.itp" \n'
f = open("md.top", "w")
for s in content_mod:
    f.write(s)
f.close()

#############################
## A minimization and addition of ions


os.system("gmx grompp -f em.mdp -p em.top -c pepappend.gro -o em_md/em1.tpr -maxwarn 1")
os.system("gmx mdrun -s em_md/em1.tpr -c em_md/em.gro -tableb ../tables/table_d0.xvg ../tables/table_a*.xvg")

os.system("gmx grompp -f npt_eq.mdp -p md.top -c em_md/em.gro -r em_md/em.gro -o em_md/ion.tpr -maxwarn 1")
os.system("gmx genion -s em_md/ion.tpr -nn 119 -nname CL -np 113 -pname ION -o mem_w_ion.gro") 

## mem_w_ion.gro contains Protein, PW, ION, CL



#############################
### Now update mdp files


with open('em.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein PW ION CL\n'
content[125] = 'tau_t\t\t\t= 1.0 1.0 1.0 1.0 \n'
content[126] = 'ref_t\t\t\t= %d %d %d %d \n' %(temp,temp,temp,temp)

f = open('em.mdp','w')
for s in content:
    f.write(s)
f.close()

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[123] = 'tc_grps\t\t\t= Protein PW ION CL\n'
content[125] = 'tau_t\t\t\t= 1.0 1.0 1.0 1.0 \n'
content[126] = 'ref_t\t\t\t= %d %d %d %d \n' %(temp,temp,temp,temp)

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()



#############################
### Update topology files

u = Universe("mem_w_ion.gro")

pep = u.select_atoms("name BB")
wat = u.select_atoms("resname PW")
ion = u.select_atoms("name ION")
cl = u.select_atoms("name CL")
pepno = int(pep.n_atoms/len(Seq))
watno = int(wat.n_atoms/3)
ionno = int(ion.n_atoms)
clno = int(cl.n_atoms)

with open("em.top") as f:
	content = f.readlines()

content_mod = content[0:14]
content_mod.append("Protein"+"\t"+str(pepno)+"\n")
content_mod.append("PW"+"\t"+str(watno)+"\n")
content_mod.append("ION"+"\t"+str(ionno)+"\n")
content_mod.append("CL"+"\t"+str(clno)+"\n")

f = open("em.top", "w")
for s in content_mod:
    f.write(s)
f.close()
content_mod[2] = '#include "itps/water.md.itp" \n'
f = open("md.top", "w")
for s in content_mod:
    f.write(s)
f.close()

#############################

### Multiple rounds of NPT equilibriation. 

### 10000 step equilibriation with a timestep of 1 femtosecond

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[11] = 'dt\t\t\t\t\t= 0.001\n'
content[12] = 'nsteps\t\t\t\t\t= 10000\n'

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()

os.system("gmx grompp -f npt_eq.mdp -p md.top -c mem_w_ion.gro -o em_md/mem_w_ion.md.tpr -r mem_w_ion.gro -maxwarn 1")
os.system("gmx mdrun -s em_md/mem_w_ion.md.tpr -c em_md/mem_w_ion.md.gro -tableb  ../tables/table_d0.xvg ../tables/table_a*.xvg")

### 10000 step equilibriation with a timestep of 5 femtosecond

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[11] = 'dt\t\t\t\t\t= 0.005\n'
content[12] = 'nsteps\t\t\t\t\t= 10000\n'

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()

os.system("gmx grompp -f npt_eq.mdp -p md.top -c em_md/mem_w_ion.md.gro -o em_md/mem_w_ion.md2.tpr -r em_md/mem_w_ion.md.gro -maxwarn 1")
os.system("gmx mdrun -s em_md/mem_w_ion.md2.tpr -c em_md/mem_w_ion.md2.gro -tableb  ../tables/table_d0.xvg ../tables/table_a*.xvg")


### 5000 step equilibriation with a timestep of 10 femtosecond

with open('npt_eq.mdp') as f:
    content = f.readlines()

content[11] = 'dt\t\t\t\t\t = 0.01\n'
content[12] = 'nsteps\t\t\t\t\t= 5000\n'

f = open('npt_eq.mdp','w')
for s in content:
    f.write(s)
f.close()

os.system("gmx grompp -f npt_eq.mdp -p md.top -c em_md/mem_w_ion.md2.gro -o em_md/mem_w_ion.md3.tpr -r em_md/mem_w_ion.md2.gro -maxwarn 1")
os.system("gmx mdrun -s em_md/mem_w_ion.md2.tpr -c em_md/mem_w_ion.md3.gro -tableb  ../tables/table_d0.xvg ../tables/table_a*.xvg")

### Remove unnecessary files
os.system("rm -rf *.xtc *.trr *.edr *.log *.pdb *.cpt") 
os.system("rm Nsvs*.gro overlap* 2prot_mod*.gro \#*")



## Create the final tpr
## Change the npt_eq.mdp 


os.system("gmx grompp -f npt_eq.mdp -p md.top -c em_md/mem_w_ion.md3.gro -o FINAL.tpr -r em_md/mem_w_ion.md3.gro -maxwarn 1")




