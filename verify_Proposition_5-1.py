from sage.all import *

import ast, sys

#################################################################
#This file verifies Proposition 5.1 parts (1) and (2) of
#''A local-global principle for isogenies of composite degree''
#by Isabel Vogt

#Part (3) is verified by the computing the genus of each of the
#subgroups in Table A.1

#The following files must be in the subfolder data:
###'gl_2_N.txt' (for N = 4, 8, 16, 32, and 64)
###'label_N_cosets.txt' (for each of the groups H in appendix A.1)
#################################################################

ell = 2

#################################################################
#Functions to create groups from their label by searchin the file
#'gl2_N.txt'
#################################################################

def make_gens_from_label(n, search_string, search_file):
    N = ell**n
    for line in open(search_file + str(N) + ".txt", "r"):
        label = line.split(":")[0][1:-1]
        if str(label) == str(search_string):
            gens = line.split(":")[1]
            gens = ast.literal_eval(gens)
            return gens
    return False
    
def make_list_elements(string):
    List = []
    for line in open(string + ".txt", "r"):
        List.append(ast.literal_eval(line))
    return List

#################################################################
#All matrix/group computations are performed as tuples of entries
#to avoid initializing matrix groups in Sage.
#################################################################


XXM1 = {}
P1 = {}
V = {}

#Precomputes the groups G[i] = GL(2, Integers(2**i)),
#the set P1[i] = P1(Integers(2**i)),
#the set V[i] = (Integers(2**i))^2
#and the solution set XXM1[i] used in has_solution function

G = []
for i in range(7):
    RR = Integers(ell**i)
    G.append(GL(2, RR))
    
    if ell == 2:
        S = set([])
        L = []
        M = []
        for x in RR:
            S.add(x * (x - 1))
            L.append([RR(1), RR(x)])
            if ZZ(x) % ell ==  0:
                L.append([RR(x),RR(1)])
            for y in RR:
                M.append([RR(x),RR(y)])
        XXM1[i] = S
        P1[i] = L
        V[i] = M


#multiplies the 2x2 matrix M represented as a list by the vector v
def mvmult(M,v):
    return [M[0]*v[0] + M[1]*v[1], M[2]*v[0]+M[3]*v[1]]

#multiplies the 2x2 matrices M and N
def mmmult(M,N):
    return [M[0]*N[0] + M[1]*N[2], M[0]*N[1]+M[1]*N[3], M[2]*N[0]+M[3]*N[2], M[2]*N[1]+M[3]*N[3]]

#computes the inverse of a 2x2 matrix M
def minv(M):
    return [mdet(M)**(-1)*M[3], -mdet(M)**(-1)*M[1], -mdet(M)**(-1)*M[2], mdet(M)**(-1)*M[0]]

#reduces mod ell^n a matrix represented as a list
#with entries in Z/ell^mZ for some m>=n
def reduce_element(h, n):
    R = Integers(ell**n)
    return [R(h[i]) for i in range(len(h))]

def reduce_list(L, n):
    return [reduce_element(h, n) for h in L]

#computes matrix trace
def mtrace(M):
    return M[0] + M[3]

#computes matrix determinant
def mdet(M):
    return M[0]*M[3]-M[1]*M[2]
    
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

  

#################################################################
#check_exceptional_subgroups is the main function that determines
#if each subgroup in 'data/gl2_N.txt' is exceptional.  The function
#then verifies that any exceptional subgroup is contained in one
#of the maximal exceptional subgroups listed in appendix A.1.
#################################################################


#checks for ALL exceptional subgroups -- not just minimally exceptional ones.
#subgroup labels and generators are fed in from the file stringN.txt
#exceptional subgroups are returned as a list
def check_exceptional_subgroups(n, string="gl2_"):
    N = ell**n
    bads = []
    for line in open(string + str(N) + ".txt", "r"):
        label = line.split(":")[0]
        gens = line.split(":")[1]
        Hgens = ast.literal_eval(gens)
        #print('checking ', label)

        #does it satisfy the conclusions?
        if not check_conclusions(n, label, Hgens):

            #does it satisfy the hypotheses?
            if check_elements(Hgens, n):

                bads.append([label, Hgens])

                #f = open("bad_all_"+ string + str(N) + ".txt", "a")
                #f.write(str(label)+":"+str(Hgens) + "\n")
                #f.close()


        #g = open("done_all_" + string + str(N) + ".txt", "a")
        #g.write(str(label) + "\n")
        #g.close()
        
    return bads


#################################################################
#Functions to check the hypotheses of the LGP
#################################################################


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


#################################################################
#Functions to check the conclusions of the LGP
#################################################################

#checks to see if subgroup of ell^n torsion generated by Hgens
#satisfies the conclusions of the LGP, i.e., if its conjugate
#to a subgroup of A_{a,c,b}(ell^n) for a <= b <= c with b + c = n.
def check_conclusions(n, label, Hgens):
    return in_A(Hgens,n)


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
        for i in range(ell):
            yield [v[0], v[1] + ell**(n - 1) * i]
    else:
        for i in range(ell):
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

def independent_vectors(v, n):
    ind_vecs = []
    for w in V[n]:
        if v[0]*w[1] - v[1]*w[0] % ell != 0:
            ind_vecs.append(w)
    return ind_vecs

#computes the coordinates for Hgens in new_basis = [v, w]
def change_basis(Hgens, new_basis):
    COB = [new_basis[0][0], new_basis[1][0], new_basis[0][1], new_basis[1][1]]
    return [mmmult(minv(COB), mmmult(h, COB)) for h in Hgens]
    
#finds the ell-adic valuation of b in Z/ell^m Z
def val(b,m):
    for r in range(1, m+1):
        if b % ell**r != 0:
            return r-1
    return m


#determines if the group generated by Hgens is contained in a group
#of the form A_{a, c, b}(ell^n) for some choice of a, c, b with
#a <= c <= b, and b+c = n; this forces a+b >= ceiling(n/2)
def in_A(Hgens, n):
    
    #m = a+b, the power of ell mod which it is Borel
    previous = False
    previousB = None
    for m in range(1, n+1):
        Hgens_m = modell(m,Hgens)
        B = borels(Hgens_m, m, previous, previousB)
        
        #first check condition of being Borel mod ell^m and Cartan mod ell^{n-m}
        if len(B) != 0:
            if (m == n) or in_cartan(modell(n-m, Hgens), n-m):
                #print('borel mod ell to the ',m,', cartan mod ell to the ',n-m)
                return True
                
        #then check if contained in another A_{a,c,b} group that's not Borel + Cartan
        for v in B:
            iv = independent_vectors(v, m)
            for w in iv:
                CB_Hgens = change_basis(Hgens_m, [v, w])
                beta_min_val = m
                delta_minus_alpha_min_val = m
                for h in CB_Hgens:
                    assert h[2] % ell**m == 0
                    beta_min_val = min(beta_min_val, val(h[1],m))
                    delta_minus_alpha_min_val = min(delta_minus_alpha_min_val, val(h[3] - h[0],m))
                a = delta_minus_alpha_min_val - beta_min_val
                b = m - a
                c = n - b
                if (c > b) or (a > c) or (a < 0):
                    break
                contained = True
                for h in CB_Hgens:
                    if (h[3] - h[0] + ell**a * h[1]) % ell**c != 0:
                        contained = False
                        break
                if contained:
                    #print('a,c,b,[v,w] are ', a, c, b, [v,w])
                    return True
        previous = True
        previousB = B
    
    return False
    
#################################################################
#Function to check if exceptional group is contained in a maximal
#exceptional group
#################################################################


#cosreps_list is a list of coset reps of H2 in G.
#Tests to see if H1 (with generators H1list) is a subgroup of H2 (with
#elements H2_elem_list) by checking that for all g in G, and h in H1list,
#the conjugate of h by g is in H2.
def is_subgroup_from_list_cosets(n, label1, H1_gen_list, label2, H2_elem_list, cosreps_list):
    #H2_elem_list = [g for g in enumerate_elements(H2_gen_list, n)]
    H2_elem_list = set([tuple(i) for i in H2_elem_list])
    for g in cosreps_list:
        is_sub = True
        for h in H1_gen_list:
            if not( tuple(mmmult(minv(g), mmmult(h, g))) in H2_elem_list):
                is_sub = False
                break
        if is_sub:
            #print str(g.matrix()) + "\n"
            return True

    return False
    
#gens is a list of generators of the group with label 'label'
#Max_data is a tuple [max_label, elements, cosreps], where
#elements is a list of all elements of the group with label
#'max_label'
#checks to see if the group with label 'label' is a subgroup
#of the group specified by Max_data
def is_subgroup_of_max(n, label, gens, Max_data):
    return is_subgroup_from_list_cosets(n, label, gens, Max_data[0], Max_data[1], Max_data[2])


#ExceptionalList is a list of tuples [label, generators],
#where label is the label of an exceptional group, and generators
#are generators of the group
#MaxList is a list of tuples [max_label, generators, cosreps], where
#elements is a list of all elements of the group with label
#'max_label', and cosreps is a list of coset reps in GL_2(Z/ell^nZ)
def sieve_maxs(n, ExceptionalList, MaxList):
    ToDo = ExceptionalList
    for Max_data in MaxList:
        print('Sieving out subgroups of ', Max_data[0])
        Max_data = [Max_data[0], [g for g in enumerate_elements(Max_data[1], n)], Max_data[2]]
        NotFound = []
        for l, l_gens in ToDo:
            if not is_subgroup_of_max(n, l, l_gens, Max_data):
                NotFound.append([l, l_gens])
        ToDo = NotFound
    if len(ToDo) == 0:
        return True
    else:
        return ToDo

#################################################################
#The maximal exceptional subgroups of Appendix A.1
#################################################################

Max3_labels = [
 '2147',
 '2177',]
 
Max5_labels = [
 '189551',
 '189605',
 '189621',
 '189785',
 '189892',
 '189979',
 '189981',
 '189995',
 '190318',
 '190435',
 '190487',
 '190525',]
 
Max6_labels = [
 '876594',
 '878116',
 '881772',
 '885865',
 '890995',
 '891525',
 '891526',
 '891735',
 '891737',
 '893009',
 '893011',
 '893326',
 '894711']

#Precompute a list consisting of tuples [label, gens, cosreps], where
#label is one of the maximal exceptional labels, and gens is a list
#of generators of the corresponding subgroup, and cosreps is a list of
#coset representatives of the group in GL2(Z/NZ)

Max3 = []
Max5 = []
Max6 = []

print('Precomputing the data for maximal subgroups mod 8.')
for l in Max3_labels:
    gens = make_gens_from_label(3, l, 'data/gl2_')
    Max3.append([l, gens, reduce_list(make_list_elements("data/" + l + "_" + str(8) + "_cosets"),3)])

print('Precomputing the data for maximal subgroups mod 32.')
for l in Max5_labels:
    gens = make_gens_from_label(5, l, 'data/gl2_')
    Max5.append([l, gens, reduce_list(make_list_elements("data/" + l + "_" + str(32) + "_cosets"),5)])

print('Precomputing the data for maximal subgroups mod 64.')
for l in Max6_labels:
    gens = make_gens_from_label(6, l, 'data/gl2_')
    Max6.append([l, gens, reduce_list(make_list_elements("data/" + l + "_" + str(64) + "_cosets"),6)])

Maxs = {2: [], 3: Max3, 4: [], 5: Max5, 6:Max6}

#################################################################
#Verify Proposition 5.1 parts (1) and (2)
#################################################################

print('This file verifies the computational assertions in Proposition 5.1 of "A local-global principle for isogenies of composite degree".  It will take approximately 3 hours to run.')
for n in range(2,7):
    N = ell**n
    print('Checking which subgroups are exceptional for ', N, 'isogenies')
    Max_list = Maxs[n]
    Max_labels = [Max_list_item[0] for Max_list_item in Max_list]
    exceptionals = check_exceptional_subgroups(n, string="data/gl2_")
    print('All done checking for exceptional groups mod ', N, '.')
    New = sieve_maxs(n, exceptionals, Max_list)
    if New:
        print('Every exceptional group is contained in one of: ', Max_labels, '.')
    else:
        print('The exceptional groups: ', New, ' are not contained in one of ', Max_labels, '.')
    


