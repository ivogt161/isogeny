from check_2 import *
from sieve_max import *

from sage.all import *

import ast, sys


#################################################################
#This file contains a function to make coset representatives of
#a subgroup in a larger group.  It is optimized for large subgroups.
#################################################################

ell = 2



#Matrix manipulation functions
def mvmult(M,v):
	return [M[0]*v[0] + M[1]*v[1], M[2]*v[0]+M[3]*v[1]]

def mmmult(M,N):
	return [M[0]*N[0] + M[1]*N[2], M[0]*N[1]+M[1]*N[3], M[2]*N[0]+M[3]*N[2], M[2]*N[1]+M[3]*N[3]]

def mdet(M):
	return M[0]*M[3]-M[1]*M[2]

def minv(M):
	return [mdet(M)**(-1)*M[3], -mdet(M)**(-1)*M[1], -mdet(M)**(-1)*M[2], mdet(M)**(-1)*M[0]]


def better_tuple(h):
	return tuple([tuple(i) for i in h.list()])

#makes a list of coset representative of H in G.  
#algorithm is optimized for #H large compared to #G.
#output is written to file
def make_coset_reps(H, G, write=False, string=None):
	reps = []
	out = []
	index = G.order()/H.order()
	#print "---------->" + str(H)
	for g in G:
		#print g
		found = False
		for g1 in reps:
			#print H
			#print g, g1
			if g.inverse() * g1 in H:
				found = True
				break
		if not found:
			reps.append(g)
			out.append(nice_list(g.list()))
			if len(reps) == index:
				break
	if write:
		for elem in out:
			f = open(string +"_cosets"+ ".txt", "a")
			f.write(str(elem)+ '\n')
			f.close()
	return out

def nice_list(L):
	return [L[0][0], L[0][1], L[1][0], L[1][1]]

