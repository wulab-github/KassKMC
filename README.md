# KassKMC
This program (KassKMC), uses kinetic Monte-Carlo simulation algorithm to calculate the rates of protein-protein association

Following files are included in this package:

1----------------------------------------------------->
The executable file of the KassKMC simulation: 

KassPBato

To run the executable file, type in:

	./kassPBato para1 para2 listname
	
para1 indicates the value of the constant w-alpha that balances the weights between the electrostatic and hydrophobic interactions.
para2 indicates the value of the constant w that balances the wwights between the physics-based and statistics-based potentials.
listname indicates the list of protein complexes as inputs for simulation, the format of the listname file can be found below.
For an example: ./kassPBato 04 06 pdb62BCHMKS, in which the value of w-alpha equals 0.4, the value of w equals 0.6 and pdb62BCHMKS.txt contains the list for all the protein complexes as inputs.


2--------------------------------------------------->
Fortran source code for the KassKMC simulation algorithm:

	KassCal_PB_StatFF_ParaAuto_main.f
	rotatefit.f
	rs.f
	diag.f

The first file is the FORTRAN77 codes for the main simulation program of the Monte-Carlo algorithm. The second file contains the FORTRAN77 subroutine deals with the rigid body superposition between two odjects. The last two files contain the FORTRAN77 subrountines deal with the matrix diagonalization.

3-------------------------------------------------->
The energy parameters used in the simulation:

	3didPotential_ResInterfaceContact_04232018.dat
	res_index.dat

The first file contains all the energy parameters derived from the statistical potential
It has the following format:

	GLY GLY   3.000      3.000
	ALA GLY   3.000      3.000
	ALA ALA  -0.076      4.392
	VAL GLY   3.000      3.000
	......
	......

The first two columns are the names of the two amino acids that form a contact. The third column is the strength of the interaction between the specific two types of amino acids indicating in the first two columns in the unit of kT, and the fourth column is the distance cutoff of the interaction.


4-------------------------------------------------->
The list of the 62 protein complexes in the benchmark:

pdb62BCHMKS.txt

The list file has the following format:

	2VLN 1 A 1 B 0.1 100000000.0
	1BRS 1 A 1 D 0.103 250000000.0 
	1UEA 1 B 1 A 0.23 200000.0	 
	1QA9 1 A 1 B 0.166 400000.0	 
	2B4J 2 A B 1 C 0.181 480000.0	
	1SGN 1 E 1 I 0.26 1200000.0	 
	1VFB 2 A B 1 C 0.15 1400000.0 
	1FLE 1 E 1 I 0.25 3600000.0	 
	1TLU 1 A 1 B 0.01 5600000.0  
	1FFW 1 A 1 B 0.15 62000000.0 
		......
		......

The first column indicates the PDB id of the protein complex. The number next indicates how many subunits in the first binding partner of the complex. 
For instance, the first binding partner of 2VLN in the list above has one subunit, while the first binding partner of 2B4J has two subunits.
The letters after the number indicate the chain ids of the subunits in the first binding partners.
Similarly, tte number after these letters indicates how many subunits in the second binding partner of the complex, while the letters after this number indicate the chain ids of the subunits in the second binding partners. Afterwards, the real number after the chain ids of the second binding partner indicates the ionic strength in which the experimental association rate was measured. The unit of the ionic strength is M. Finally, the real numer in the last column is the experimentally measured value of association rate for the corresponding protein complex between the first and second binding partners which chain ids are indicated in the list. 


5--------------------------------------------------->
The sample PDB input for the test system:

2VLN.pdb

The input files have the standard PDB format as following:

	ATOM      1  N   SER A   6      -9.008  52.410  14.065  1.00 22.23           N  
	ATOM      2  CA  SER A   6      -8.896  52.304  12.583  1.00 20.74           C  
	ATOM      3  C   SER A   6      -7.955  53.416  12.154  1.00 18.44           C  
	ATOM      4  O   SER A   6      -7.592  54.266  12.955  1.00 18.02           O  
	ATOM      5  CB  SER A   6     -10.289  52.471  11.966  1.00 22.99           C  
	ATOM      6  OG  SER A   6     -10.554  53.797  11.560  1.00 26.51           O  
	ATOM      7  N   ILE A   7      -7.603  53.435  10.874  1.00 14.65           N  
	ATOM      8  CA  ILE A   7      -6.620  54.443  10.387  1.00 13.96           C  
	ATOM      9  C   ILE A   7      -7.187  55.856  10.588  1.00 14.95           C  
	ATOM     10  O   ILE A   7      -6.438  56.802  10.792  1.00 14.72           O  
	ATOM     11  CB  ILE A   7      -6.182  54.168   8.945  1.00 15.26           C  
	ATOM     12  CG1 ILE A   7      -4.756  54.727   8.682  1.00 11.56           C  
	ATOM     13  CG2 ILE A   7      -7.217  54.705   7.959  1.00 14.47           C  
	ATOM     14  CD1 ILE A   7      -4.065  54.315   7.300  1.00 13.73           C  
				......
				......

They can be either downloaded from the Protein Data Bank, or computationally modeled. The input PDB structures require backbone heavy atoms, as well as sidechain heavy atoms. HOwever, the files might but not necessarily contain hydrogen atoms. If the original PDB files contain hydrogen atoms, you don't have to remove them. On the other hand, if the original PDB files don't contain hydrogen atoms, you don't have to add extra hydrogen atoms.


6--------------------------------------------------->
The sample output from the simulation for the test system:

Result_ResPB_rec_2VLN_A00_B00_WS04WK06.dat



    1  0      1000.00000        0.00000       11.47618       41.55275        8.65784    0       -2.40446
    2  1       335.00000        0.07895       12.04541       32.26528        5.03485    3      -12.18327
    3  1       642.00000        0.07895       13.43522       28.28039        6.42724    3      -13.69676
    4  0      1000.00000        0.00000       24.31463       66.61725        9.88474    0       -0.09261
    5  0      1000.00000        0.00000       39.63525      107.46005       16.95219    0       -0.00072
	......
	......
