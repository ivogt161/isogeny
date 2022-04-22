
from sage.all import *

import ast, sys


#################################################################
#This file contains a function to make coset representatives of
#a subgroup in a larger group.  It is optimized for large subgroups.
#################################################################

ell = 2

G = []
for i in range(7):
    RR = Integers(ell**i)
    G.append(GL(2, RR))


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


def make_reps(n, Hgens):
    if n == 1:
        return make_coset_reps(G[1].subgroup([G[1](v) for v in Hgens]), G[1])
    if n > 1:
        H = G[n].subgroup([G[n](v) for v in Hgens])
        reps_nminus1 = make_reps(n-1, Hgens)
        kernel_reps = []
        for x in range(ell):
            for y in range(ell):
                for z in range(ell):
                    for w in range(ell):
                        M = G[n]([1 + ell**(n-1)*x, ell**(n-1)*y, ell**(n-1)*z, 1+ ell**(n-1)*w])
                        found = False
                        for g1 in kernel_reps:
                            if M.inverse()*g1 in H:
                                found = True
                                break
                        if not found:
                            kernel_reps.append(M)
        out = []
        for g in reps_nminus1:
            g = G[n](g)
            for h in kernel_reps:
                out.append(nice_list((g*h).list()))
        assert(len(out) == G[n].order()/H.order())
        return out

def write_reps(n, Hgens, label):
    for elem in make_reps(n, Hgens):
        f = open(str(label) + "_" + str(ell**n) +"_cosets"+ ".txt", "a")
        f.write(str(elem)+ '\n')
        f.close()
    return
        

def nice_list(L):
	return [L[0][0], L[0][1], L[1][0], L[1][1]]

