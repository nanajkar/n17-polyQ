# n17-polyQ


Thhis repository contains all files required to set up protein aggregation simulations similar to those carried out in  ** Conformations of the Huntingtin protein and variations in its aggregate morphology**. This tutorial is used to set up N17+7Q aggregation system, but other lengths of polyQ can be created by following the same protocol with minor tweaks. 

**Need to include a file tree so people know contents of each folder**

### Generate the protein gro file 
Create a file called seq.txt containing the sequence of the protein to simulate. For N17+7Q system the sequence is MATLEKLMKAFESLKSFQQQQQQQ.

The following command will create a gro file entitled protein.gro, corresponding to the sequence supplied in seq.txt
> python create_prot_gro.py -f seq.txt -o protein.gro

### Generate the protein itp file

Use the following command to create an itp file for the sequence provided in seq.txt
> python gen_prot_em.itp.py -f seq.txt -o protein.itp

In our model, a dihedral potential is applied to the N17 domain, but not to the polyQ domain. This needs to be changed manually in the protein.itp file. In the [dihedrals] section, delete the lines corresponding to dihedral angles in the polyQ domain. Thus, the last line in the [dihedrals] section should be as follows-
> 61              65              70              76              8                 0      10

Move this file to the itps folder.

### Set up script

protocol.py is a comprehensive script that can be used to setup all simulations(with minor tweaks). The script is set up as follows

1. Minimizing and equilibrating the extended protein structure in vaccum
2. Addition of proteins to the water box
3. Updating the mdp and topology files to reflect the particles in the system 
4. Addition of counterions and a salt concentration that amounts to 125 mM
5. Multiple rounds of equilibration
6. Creation of the final tpr file



