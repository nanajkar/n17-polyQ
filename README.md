# n17-polyQ


This repository contains all files required to set up protein aggregation simulations similar to those carried out in  **Conformations of the Huntingtin protein and variations in its aggregate morphology**. This tutorial is used to set up N17+7Q aggregation system, but other lengths of polyQ can be created with minor tweaks to the protocol. 


## Files
- ```itps/```
  + ```fixedff_noCP.itp``` - A version of the PRompt forcefield with no cation pi interactions
  + ```martini_v2.0_ions.itp``` - Topology file for monovalent cations
  + ```protein.itp``` - Protein topology file, generated by ```gen_prot_em.itp.py```
  + ```posre_protein.itp``` - Restrains the peptide backbone at positions K6 and L14
  + ```water.em.itp``` - Topology for Martini Polarizable Water. Interaction between the backbone and dummy charges are treated as contraints.
  + ```water.md.itp``` - Topology for Martini Polarizable Water. Interaction between the backbone and dummy charges are modelled as harmonic bonds.
- ```gen_conf/```
  + ```create_prot_gro.py``` - Script to construct the initial gro file of the protein from sequence 
  + ```gen_prot_em.itp.py``` - Script to contruct the topology file for the protein
  + ```seq.txt``` - Input for both ```create_prot_gro.py``` and ```+gen_prot_em.itp.py``` 
- ```tables/``` - Contains angular and dihedral potentials

### Generate the protein gro file 
Navigate to ```gen_conf```. Create a file called seq.txt containing the sequence of the protein to simulate. For this tutorial, we use the sequence for N17+7Q, MATLEKLMKAFESLKSFQQQQQQQ.

The following command will create a gro file entitled protein.gro, corresponding to the sequence supplied in seq.txt
> python create_prot_gro.py -f seq.txt -o protein.gro

### Generate the protein itp file

Use the following command to create an itp file for the sequence provided in seq.txt
> python gen_prot_em.itp.py -f seq.txt -o protein.itp

In our model, a dihedral potential is applied to the N17 domain, but not to the polyQ domain. This needs to be changed manually in the protein.itp file. In the **[dihedrals]** section, delete the lines corresponding to dihedral angles in the polyQ domain. Thus, the last line in the **[dihedrals]** section should be as follows-
> 61 &nbsp; 65 &nbsp; 70 &nbsp; 76 &nbsp; 8 &nbsp; 0 &nbsp; 10

Navigate back to the outer directory using ```cd ..```

### Set up script

setup.py is a comprehensive script that sets up simulation with minor ( with minor tweaks). The script is set up as follows

1. Minimizing and equilibrating the extended protein structure in vaccum
2. Addition of protein to the water box
3. Updating the mdp and topology files to reflect the particles in the system 
4. Addition of counterions to neutralize the system and a salt concentration amounting to 125 mM. The number of ions to be added will change depending on the size  of the simulation box.  
5. Multiple rounds of equilibration
6. Creation of the final tpr file


To run it, simply run the following command.
> python setup.py

To calculate the number of ions to be added, use the following formula
$$ \frac{4*number of waters*ion concentration}{55.5+ion concentration} $$



