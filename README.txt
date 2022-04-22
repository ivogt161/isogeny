This is an outline of computations performed by various files in this folder.  

The two main files are:

- The file verify_Proposition_5-1.py contains sage functions written to check every subgroup of GL_2(Z/2^nZ) and determine if it is exceptional, as in Definition 3.5. It will take approximately 3 hours to run this file, and all of the data files in the data subfolder must be present.

- The file verify_Theorem_1.py contains sage functions necessary to compute the genera of fiber products of modular curves X \times_{X(1)} Y used in the proof of Theorem 1 part (2).  The method is via Riemann-Hurwitz applied to the base-change of X -> X(1) by Y -> X(1).  There are also functions to determine the ramification profile over infinity for any congruence subgroup of prime power level.  It will take less than 10 seconds to run and the file ram_data.txt in the data subfolder must be present.


Additionally, the file check_Theorem_3-7.py verifies the classification of lift-exceptional subgroups given in Theorem 3.7 in the special case that ell is 3 and n is at most 5.


The data files are in the subfolder data:

- The files gl2_N.txt contains data of all subgroups of GL_2(Z/NZ) in the form [label]: [generators] for N = ell^n with ell = 2 or 3 and n<=5.  The data file for ell = 2 and n=6 is too large to store on GitHub, and so must be generated to use.  This data is generated from Drew Sutherland's data on subgroups of GL_2(Z/ell^nZ) that can be found at http://math.mit.edu/~drew/gl_2_full.tar.  All labeling of groups is consistent with those labels. The file used to generate these data files is parse_drew_data.py from the auxiliary code folder.

- The files label_N_cosets.txt contain lists of coset representative for the group with that label inside of GL_2(Z/NZ).  The file used to generate these data files is coset_reps.py from the auxiliary code folder.

- The file ram_data.txt contains ramification data for all of the maximal exceptional subgroups in Appendix A with 2-power level, as well as the exceptional subgroups H_{ell^n, exc} for ell = 5, 7 from Proposition 7.4.  This data is in the form label: n: level: index: generators for Gammabar: nu2: nu3: genus.  The function WriteRam from the file WorkBads.m in the auxiliary code folder makes this data file.

The files in the auxiliary folder are used to generate data files as well as to do other manipulations with groups not directly needed in the verification of Proposition 5.1 or Theorem 1: 

- The file parse_drew_data.py creates the gl2_N.txt data file from http://math.mit.edu/~drew/gl_2_full.tar.

- The file parse_bads.py contains various functions needed to read output files and to create input for other functions.

- The file coset_reps.py contains a function which computes a set of coset representative for a subgroup H in G.  This is optimized for the case when H is large compared to G, and is applied to maximal subgroups.

- The file WorkBads.m contains the functions necessary to match exceptional subgroups with those in the work of Rouse-Zureick-Brown and Sutherland-Zywina.  The data of exceptional subgroups is in bad2.m and the data from the work of Rouse-Zureick-Brown is in rzbData.m.  The function WriteRam makes the file ram_data.txt which contains the ramification data necessary for the functions in ram.py.
