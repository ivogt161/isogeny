from check_2 import *
from sieve_max import *

from sage.all import *

import ast, sys


#################################################################
#This file contains functions written to manipulate the output of
#check_exceptional_subgroups and similar functions that write 
#labels or labels and generators to a text file.
#################################################################

ell = 2

R = []
for i in xrange(7):
	R.append(Integers(ell**i))

G = []
for i in xrange(7):
	G.append(GL(2, R[i]))


def make_gens_from_label(n, search_string, search_file):
	N = ell**n
	for line in open(search_file + str(N) + ".txt", "r"):
		label = line.split(":")[0][1:-1]
		if str(label) == str(search_string):
			gens = line.split(":")[1]
			gens = ast.literal_eval(gens)
			return gens
	return False

def make_group_from_label(n, search_string, search_file):
	N = ell**n
	for line in open(search_file + str(N) + ".txt", "r"):
		label = line.split(":")[0][1:-1]
		print label
		if str(label) == str(search_string):
			gens = line.split(":")[1]
			gens = ast.literal_eval(gens)
			return G[n].subgroup([G[n](v) for v in gens])
	return False

def make_groups(n, string="Badgens_"):
	subgroups = []
	N = ell**n
	for line in open(string + str(N) + ".txt", "r"):
		label = line.split(":")[0]
		gens = line.split(":")[1]
		gens = ast.literal_eval(gens)
		L = []
		for v in gens:
			L.append(G[n](v))
		H = G[n].subgroup(L)
		subgroups.append([label, H])
	return subgroups


def make_gens(n, string="Badgens_"):
	subgroup_gens = []
	N = ell**n
	for line in open(string + str(N) + ".txt", "r"):
		label = line.split(":")[0]
		gens = line.split(":")[1]
		gens = ast.literal_eval(gens)
				#L = []
				#for v in gens:
				#		L.append(G[n](v))
		subgroup_gens.append([label, gens])
	return subgroup_gens


def make_gens_file(n, search_string="gl2", string="bad_gl2"):
	subgroup_gens = []
	N = ell**n
	for line in open(string + "_" + str(N) + ".txt","r"):
		label = line.split(":")[0][:-1]
		for gen_line in open(search_string + "_" + str(N) + ".txt","r"):
			if str(label) == str(gen_line.split(":")[0]):
				gens = gen_line.split(":")[1]
				f = open(string + "_gens_" + str(N) + ".txt", "a")
				f.write(str(label)+":"+str(gens))
				f.close()
	return


def make_list_elements(string):
	List = []
	for line in open(string + ".txt", "r"):
		List.append(ast.literal_eval(line))
	return List

				
def make_group_gen_file(n, search_string="gl2", string="bad"):

	subgroups = []

	N = ell**n
	
	for line in open(string + "_" + str(N) + ".txt","r"):
		label = line.split(":")[0][:-1]
		for gen_line in open(search_string + "_"+ str(N) + ".txt", "r"):
			#print gen_line.split(":")[0]
			#print str(label)
			if str(label) == str(gen_line.split(":")[0]):
				#print "yo"
				gens = gen_line.split(":")[1]
				f = open(string + "_gens_" + str(N) + ".txt", "a")
				f.write(str(label)+ ":" + str(gens))
				f.close()
				gens = ast.literal_eval(gens)
				#print gens
				L = []
				for v in gens:
					L.append(G[n](v))
				H = G[n].subgroup(L)
				subgroups.append([label,H])
				break
	
	return subgroups
