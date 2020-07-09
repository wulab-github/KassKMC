# KassKMC
This program (KassKMC), uses kinetic Monte-Carlo simulation algorithm to calculate the rates of protein-protein association

Following files are included in this package:

1---------------------------------------------------------------------------------------->
The executable file of the KassKMC simulation: 

KassPBato

To run the executable file, type in ./kassPBato para1 para2 listname.
para1 indicates the value of the weight constant w-alpha in the equation (3) of the paper between the electrostatic and hydrophobic interactions.
para2 indicates the value of the weight constant w in the equation (2) of the paper between the physics-based and statistics-based potentials.
listname indicates the list of protein complexes as inputs for simulation, the format of the listname file can be found below.
For an example: ./kassPBato 04 06 pdb62BCHMKS, in which the value of w-alpha equals 0.4, the value of w equals 0.6 and pdb62BCHMKS.txt contains the list for all the protein complexes as inputs.


2---------------------------------------------------------------------------------------->
Fortran source code for the KassKMC simulation algorithm:

KassCal_PB_StatFF_ParaAuto_main.f
rotatefit.f
rs.f
diag.f



3--------------------------------------------------------------------------------------->
The energy parameters used in the simulation:

3didPotential_ResInterfaceContact_04232018.dat
res_index.dat



4--------------------------------------------------------------------------------------->
The list of the 62 protein complexes in the benchmark:

pdb62BCHMKS.txt



5--------------------------------------------------------------------------------------->
The sample PDB input for the test system:

2VLN.pdb



6--------------------------------------------------------------------------------------->
The sample output from the simulation for the test system:

Result_ResPB_rec_2VLN_A00_B00_WS04WK06.dat


