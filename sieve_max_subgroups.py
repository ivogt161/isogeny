import parse_bads
from check_2 import *

from sage.all import *

import ast, sys



#################################################################
#This file contains functions to check if all subgroups listed in
#a given file are contained in a list of (hypothetically complete)
#maximal subgroups
#################################################################



ell = 2

R = []
for i in xrange(7):
	R.append(Integers(ell**i))

G = []
for i in xrange(7):
	G.append(GL(2, R[i]))

def reduce_element(h, n):
	R = Integers(ell**n)
	return [R(h[i]) for i in xrange(len(h))]

def reduce_list(L, n):
	return [reduce_element(h, n) for h in L]

#List is a list of coset reps of H2 in G.
#Tests to see if H1 (with generators H1list) is a subgroup of H2 (with
#elements H2list) by checking that for all g in G, and h in H1list,
#the conjugate of h by g is in H2.
def is_subgroup_from_list_cosets(label1, H1list, label2, H2list, List):
    H2list = set([tuple(i) for i in H2list])
    for g in List:
        is_sub = True
        for h in H1list:
            if not( tuple(mmmult(minv(g), mmmult(h, g))) in H2list):
                is_sub = False
                break
        if is_sub:
            #print str(g.matrix()) + "\n"
            return True

    return False

#Maxs is a list of hypothetically maximal subgroups (label, list of elements) 
#of all subgroups of GL_2(Z/ell^nZ) given in input_file, which are presented 
#in gensList with their generators.
#cosreps is a dictionary of coset reps of the maximal subgroups in GL_2(Z/ell^nZ)
#Sieves out subgroups of the groups in Maxs and writes remaining subgroups
#to file not_in_(labels of maxs)
def sieve_subs(gensList, Maxs, n, cosreps, input_file, out):
	N = 2**n
	string = out + str(N)
	groups = str()
	for label, group in Maxs:
		groups = groups + " " + str(label)

	fp = open("not_in_" + string + ".txt", "a")
	fp.write("Subgroups from file " + str(input_file) + " that are not contained in one of the groups:" + groups + ". \n")
	fp.close()
	
	bad = []

	for l, Hgens in gensList:
		print l
		Hgens = reduce_list(Hgens,n)
		subgroup = False
		for lm, Hm_elems in Maxs:
			if is_subgroup_from_list_cosets(l, Hgens, lm, Hm_elems, cosreps[lm]):
				subgroup = True
				fl = open("subgroups_" + string + ".txt", "a")
				fl.write(str(l) + ":" + str(lm) + "\n")
				fl.close()
				break
		if not subgroup:
			f = open("not_in_" + string + ".txt", "a")
			f.write(str(l)+ "\n")
			f.close()
			bad.append([l, Hgens])	

	h = open("not_in_" + string + ".txt", "a")
	h.write("All Done")
	h.close()

	return bad

#Prepares everything from the list of labels of maximal subgroups to be put into
#the sieve_max function.  Uses functions from parse_bads.py to prepare all groups
#from file listString from which subgroups of the groups in Max will be removed.
def do_sieve(Labels, n, listString, search_file, out):
	N = ell**n
	Maxs = []
	d_cosreps = {}
	for l in Labels:
		gens = make_gens_from_label(n, str(l), search_file)
		#print len(gens)
		Maxs.append([l, [g for g in enumerate_elements(gens, n)]])
		#print len([g for g in enumerate_elements(gens, n)])
	
		repList = make_list_elements(l + "_" + str(N) + "_cosets_new")
		#print repList
		d_cosreps[l] = reduce_list(repList, n)

	genList = make_gens(n, listString)

	return sieve_subs(genList, Maxs, n, d_cosreps, listString, out)


#Matrix multiplication functions as in check_subgroups_2.py:
def mvmult(M,v):
    return [M[0]*v[0] + M[1]*v[1], M[2]*v[0]+M[3]*v[1]]

def mmmult(M,N):
    return [M[0]*N[0] + M[1]*N[2], M[0]*N[1]+M[1]*N[3], M[2]*N[0]+M[3]*N[2], M[2]*N[1]+M[3]*N[3]]

def mdet(M):
    return M[0]*M[3]-M[1]*M[2]

def minv(M):
    return [mdet(M)**(-1)*M[3], -mdet(M)**(-1)*M[1], -mdet(M)**(-1)*M[2], mdet(M)**(-1)*M[0]]

