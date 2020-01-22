from sage.all import *

import ast, sys

#################################################################
#This file contains the function check_exceptional_subgroups
#which finds all exceptional subgroups from a file given generators
#of the subgroup in GL_2(Z/ell^nZ).

#All matrix/group computations are performed as tuples of entries
#to avoid initializing matrix groups in Sage.
#################################################################

ell = 2
XXM1 = {}
P1 = {}


#Precomputes the groups G[i] = GL(2, Integers(2**i)), 
#the set P1[i] = P1(Integers(2**i)), and the solution set 
#XXM1[i] used in has_solution function

G = []
for i in xrange(7):
	RR = Integers(ell**i)
	G.append(GL(2, RR))
	
	if ell == 2:
		S = set([])
		L = []
		for x in RR:
			S.add(x * (x - 1))
			L.append([1, x])
			if ZZ(x) % ell ==  0:
				L.append([x,1])
									  
		XXM1[i] = S
		P1[i] = L


#multiplies the 2x2 matrix M represented as a list by the vector v
def mvmult(M,v):
	return [M[0]*v[0] + M[1]*v[1], M[2]*v[0]+M[3]*v[1]]

#multiplies the 2x2 matrices M and N
def mmmult(M,N):
	return [M[0]*N[0] + M[1]*N[2], M[0]*N[1]+M[1]*N[3], M[2]*N[0]+M[3]*N[2], M[2]*N[1]+M[3]*N[3]]


#checks for ALL exceptional subgroups -- not just minimally exceptional ones.
#subgroup labels and generators are fed in from the file stringN.txt
#bad subgroups are output to file bad_all_stringN.txt
def check_exceptional_subgroups(n, string="gl2_"):
	N = ell**n
	for line in open(string + str(N) + ".txt", "r"):
		label = line.split(":")[0]
		gens = line.split(":")[1]
		Hgens = ast.literal_eval(gens)

		#does it satisfy the conclusions?
		if not check_conclusions(n, label, Hgens):

			#does it satisfy the hypotheses?
			if check_elements(Hgens, n):

				f = open("bad_all_"+ string + str(N) + ".txt", "a")
				f.write(str(label)+":"+str(Hgens) + "\n")
				f.close()
				print str(label) + " is bad for " + str(N) + "-torsion"

			else:
				print str(label) + " is good for " + str(N) + "-torsion"
		else:
			print str(label) + " is good for " + str(N) + "-torsion"

		g = open("done_all_" + string + str(N) + ".txt", "a")
		g.write(str(label) + "\n")
		g.close()
	
	return


#checks a subgroup of the ell^m torsion which is Borel mod ell^(n-1) to see if it
#satisfies the conclusions of the lg principle.  If it satisfies the hypotheses
#but fails this test, then the group is
#exceptional (not necessarily minimally exceptional) 
def check_conclusions(m, label, Hgens, n=1):
	N = ell**m
	
	previous = False
	previousB = None

	#need to check to see if there exist j >= k and j + k = m such that
	#the group is Borel mod ell^j and Cartan mod ell^k 
	#We already know that it is Borel mod ell^n, so we just need to make up m-n
	for i in xrange(n, m+1):
		#print i
		B = borels(modell(i, Hgens), i, previous, previousB)
		if len(B) != 0:
			if (i == m) or in_cartan(modell(m-i, Hgens), m-i):
				#print str(label) + " is good with "+ str((i, m-i))
				return True
			previousB = B
			previous = True
		else:
			#print str(label) + " is bad.  Its not Borel mod ell to the " +  str(i)
			return False

#enumerates elements of H as lists using a FIXME algorithm
def enumerate_elements(Hgens, n):

	R = Integers(ell**n)
	yield reduce_element([1,0,0,1],n)

	traced = set([tuple(reduce_element([1,0,0,1],n))])
	todo = [reduce_element([1,0,0,1],n)]

	while todo:
		g = todo.pop()
		for h in Hgens:
			hred = reduce_element(h,n)
			hg = mmmult(h,g)
			if not tuple(hg) in traced:
				yield hg
				traced.add(tuple(hg))
				todo.append(hg)

	return


#reduces mod ell^n a matrix represented as a list 
#with entries in Z/ell^mZ for some m>=n
def reduce_element(h, n):
	R = Integers(ell**n)
	return [R(h[i]) for i in xrange(len(h))]


#computes matrix trace
def mtrace(M):
	return M[0] + M[3]

#computes matrix determinant
def mdet(M):
	return M[0]*M[3]-M[1]*M[2]


#checks elements of H to see that their characteristic 
#polynomial has a rational root
def check_elements(Hgens, n):

	RR = Integers(2**n)
	for h in enumerate_elements(Hgens, n):
		if not has_solution(mtrace(h), mdet(h), RR, n):
			#print h
			return False

	return True


#determines if x^2 - ax + b = 0 has a solution in RR
def has_solution(a, b, RR, n):

	if ell != 2:
		return (a**2 - 4 * b).is_square()

	if a % 2 == 0:
		alpha = RR(ZZ(a) / 2) # a = 2 * alpha --> x^2 - 2 * alpha * x + alpha**2 = alpha**2 - b
		return (alpha**2 - b).is_square()
	else:
		alpha = RR(ZZ(a - 1) / 2) # a = 2 * alpha + 1 --> x^2 - 2 * alpha * x + alpha**2 - (x - alpha) = alpha**2 + alpha - b
		return (alpha**2 + alpha - b) in XXM1[n]


#determines if v is an eigenvector of all of H
#by checking that it is an eigenvector for each
#h in Hgens
def is_eigenvector(Hgens, v):
	
	assert(v[0] % ell != 0 or v[1] % ell != 0)

	for h in Hgens:
		hv = mvmult(h, v)
		if (v[0] * hv[1] != v[1] * hv[0]):
			return False
	return True


#v is a vector with entries in Z/ell^(n-1)Z,
#which is nonzero mod ell
#computes all lifts of v to Z/ell^nZ up to rescaling
def lifts(v, n):
	R = Integers(ell**n)
	v = [R(v[0]), R(v[1])]
	if v[0] % ell != 0:
		for i in xrange(ell):
			yield [v[0], v[1] + ell**(n - 1) * i]
	else:
		for i in xrange(ell):
			yield [v[0] + ell**(n - 1) * i, v[1]]
	return


#computes the set of borel subgroups of GL_2(Z/ell^nZ)
#containing H.  Borels are indexed by a representative
#of the stable line.
#previousB is the list of Borels mod n-1
def borels(Hgens, n, previous=False, previousB = None):
	output = []
	if n == 1:
		for v in P1[n]:
			if is_eigenvector(Hgens, v):
				output.append(v)
		return output
	if previous:
		B = previousB
	else:
		B = borels(modell(n - 1, Hgens), n - 1)
	for v in B:
		for w in lifts(v, n):
			if is_eigenvector(Hgens, w):
				output.append(w)
	return output


#checks to see if the group generated by Hgens is 
#contained in a Cartan subgroup of GL_2(Z/ell^nZ)
def in_cartan(Hgens, n):
	vecs = []
	for v in P1[n]:
		if len(vecs) == 1:
			w = vecs[0]
			if (v[0] * w[1] - w[0] * v[1]) % ell == 0:
				continue
		if is_eigenvector(Hgens, v):
			vecs.append(v)
		if len(vecs) == 2:
			return True
	return False

#reduces every element of Hgens (presented as a list)
#mod ell^n
def modell(n, Hgens):
	R = Integers(2**n)
	return [[R(h[0]), R(h[1]), R(h[2]), R(h[3])] for h in Hgens]
