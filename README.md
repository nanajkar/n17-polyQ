# n17-polyQ

This tutorial explains how to set up the N17+7Q system discussed in ** Conformations of the Huntingtin protein and variations in its aggregate morphology**. Other lengths of polyQ can be created with the same protocol. 

## Generate the protein gro file 
Create a file called seq.txt containing the sequence of the protein to simulate. For N17+7Q system the sequence is MATLEKLMKAFESLKSFQQQQQQQ.

> python create_prot_gro.py -f seq.txt -o protein.gro

## Generate the itp files 

Use the following command to create an itp file for the sequence provided in seq.txt
> python gen_prot_em.itp.py -f seq.txt -o protein.itp

In our model, a dihedral potential is applied to the N17 domain, but not to the polyQ domain. This needs to be changed manually in the protein.itp file. In the [dihedrals] section, delete the lines corresponding to dihedral angles in the polyQ domain. 


3. Use gen_prot.py to create an itp file for the protein. Move this file to the itp folder
4. Create the gro file using the xx script. This will generate a gro file from the supplied seq.txt file. 
5. M
## Add proteins to water box
## Run the steup.py script
