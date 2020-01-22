This is an outline of computations performed by various files in this folder.

Drew Sutherland's data on subgroups of GL_2(Z/2^nZ) for n <= 6 can be found at http://math.mit.edu/~drew/gl_2_full.tar.  All labeling of groups is consistent with those labels.

- The file check_subgroups_2.py contains sage functions written to check every subgroup of GL_2(Z/2^nZ) and determine if it is exceptional, as in Definition 5.3. The output is a list of exceptional groups, from which Prop 5.4 and Table 1 contain the maximal subgroups.

- The file parse_bads.py contains various functions needed to read output files and to create input for other functions.

- The file coset_reps.py contains a function which computes a set of coset representative for a subgroup H in G.  This is optimized for the case when H is large compared to G, and is applied to maximal subgroups.

- The file sieve_max_subgroups contains a function which sieves out the subgroups of some maximal subgroup from a list. The function which checks if group A is conjugate to a subgroup of group B uses additional input of coset representatives of B in the whole group G.  This function is applied to the list of exceptional subgroups, with the maximal ones being those of largest order.  In all cases studied (n <= 6) these suffice to contain all exceptional subgroups. 

- The file ram.py contains the functions necessary to compute the genera of fiber products of modular curves X \times_{X(1)} Y used in the proof of Theorem 1 part (2).  The method is via Riemann-Hurwitz applied to the base-change of X -> X(1) by Y -> X(1).  There are also functions to determine the ramification profile over infinity for any congruence subgroup of prime power level.

- The file WorkBads.m contains the functions necessary to match exceptional subgroups with those in the work of Rouse-Zureick-Brown and Sutherland-Zywina.  The data of exceptional subgroups is in bad2.m and the data from the work of Rouse-Zureick-Brown is in rzbData.m.  The function WriteRam makes the file ram_data.txt which contains the ramification data necessary for the functions in ram.py.
