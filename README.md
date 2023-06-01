# n17-polyQ

This tutorial explains how to set up the N17+7Q system discussed in ** Conformations of the Huntingtin protein and variations in its aggregate morphology**. Other lengths of polyQ can be created by following the same protocol. 

### Generate the protein gro file 
Create a file called seq.txt containing the sequence of the protein to simulate. For N17+7Q system the sequence is MATLEKLMKAFESLKSFQQQQQQQ.

The following command will create a gro file entitled protein.gro that corresponds to the sequence supplied in seq.txt
> python create_prot_gro.py -f seq.txt -o protein.gro

### Generate the protein itp file

Use the following command to create an itp file for the sequence provided in seq.txt
> python gen_prot_em.itp.py -f seq.txt -o protein.itp

In our model, a dihedral potential is applied to the N17 domain, but not to the polyQ domain. This needs to be changed manually in the protein.itp file. In the [dihedrals] section, delete the lines corresponding to dihedral angles in the polyQ domain. Thus, the last line in the [dihedrals] section should be as follows-
> 61              65              70              76              8                 0      10

Move this file to the itps folder.

### Set up script

The simulation will be set up using protocol.py. This script adds proteins and ions to the water box and performs many rounds of equilibriation. 


