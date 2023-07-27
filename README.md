# n17-polyQ


This repository contains all files required to set up protein aggregation simulations similar to those carried out in  **Conformations of the Huntingtin protein and variations in its aggregate morphology**. This tutorial is used to set up N17+7Q aggregate system, but other lengths of polyQ can be created with minimal changes. 


## Software Requirements
- Python 2 or Python 3
- Gromacs 2019.4 (Will also work with Gromacs 2018.3). This tutorial assumes a non-MPI installation of Gromacs, however, by changing the gmx binary to `gmx_mpi`, the tutorial should work with MPI installations.

 ## Contributions and  Acknowledegments
 This tutorial was created by Neha Nanajkar. Simulations referenced in the paper were carried out by Neha Nanajkar. Subsequent analysis was jointly carried out by Neha Nanajkar and Abhilash Sahoo, under the guidance of Silvina Matysiak. 

 For questions or clarifications please reach out to Silvina Matysiak(matysiak@umd.edu) or Neha Nanajkar(nanajkar@terpmail.umd.edu)


## Files
- ```itps/```
  + ```fixedff_noCP.itp``` - A version of the ProMPT forcefield with no cation pi interactions. See [Sahoo et al. Transferable and Polarizable Coarse Grained Model for Proteins─ProMPT J. Chem. Theory Comput (2022)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00269) for more details.
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

### 1. Generate the protein gro file 
Navigate to ```gen_conf```. Create a file called seq.txt containing the sequence of the protein to simulate. For this tutorial, we use the sequence for N17+7Q, MATLEKLMKAFESLKSFQQQQQQQ.

The following command will create a gro file entitled protein.gro, corresponding to the sequence supplied in seq.txt
> python create_prot_gro.py -f seq.txt -o protein.gro

### 2. Generate the protein itp file

Use the following command to create an itp file for the sequence provided in seq.txt
> python gen_prot_em.itp.py -f seq.txt -o protein.itp

In our model, a dihedral potential is applied to the backbone of the N17 domain, but not to the polyQ domain. This needs to be changed manually in the protein.itp file. In the **[dihedrals]** section, delete the lines corresponding to dihedral angles in the polyQ domain. Thus, the last line in the `tabulated dihedrals between BB and BB beads` section should be between beads 61, 65, 70, and 76. 

Navigate back to the outer directory using ```cd ..```

### 3. Set up script

setup.py is a comprehensive script that sets up simulations. The script is set up as follows - 

1. Minimizing and equilibrating the extended protein structure in vaccum
2. Addition of protein to the water box
3. Updating the mdp and topology files to reflect the particles in the system 
4. Addition of counterions to neutralize the system and a salt concentration amounting to 125 mM using the ```genion``` command. The number of ions to be added will change depending on the size  of the simulation box. To calculate the number of ions to be added, use the following formula :  $` \frac{4*N_{PW}*M_{NaCl}}{55.5+M_{NaCl}} `$ , where $` N_{PW}`$ represents the number of coarse-grain water beads in the system and $` M_{NaCl}`$ denotes the desired salt molar concentration, in this case 125mM. 
5. Multiple rounds of NPT equilibration
6. Creation of the final tpr file


To run it, simply run the following command.
> python setup.py







