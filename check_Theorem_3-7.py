from sage.all import *

import ast, sys


ell = 3

#################################################################
#This file confirms Theorem 3.7 for ell = 3 and n <= 5 of
#''A local-global principle for isogenies of composite degree''
#by Isabel Vogt (references are to arXiv version 2)

#The following files must be in the subfolder data:
###'gl_2_N.txt' (for N = 9, 27, 81, and 243)
###'label_N_cosets.txt' (for the groups labeled '2788' and '81376')
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
    
def make_group_from_label(n, search_string, search_file):
    N = ell**n
    for line in open(search_file + str(N) + ".txt", "r"):
        label = line.split(":")[0][1:-1]
        print(label)
        if str(label) == str(search_string):
            gens = line.split(":")[1]
            gens = ast.literal_eval(gens)
            return G[n].subgroup([G[n](v) for v in gens])
    return False

#################################################################
#All matrix/group computations are performed as tuples of entries
#to avoid initializing matrix groups in Sage.
#################################################################


#XXM1 = {}
P1 = {}
V = {}

#Precomputes the groups G[i] = GL(2, Integers(ell**i)),
#the set P1[i] = P1(Integers(ell**i)),
#the set V[i] = (Integers(ell**i))^2

G = []
for i in range(6):
    RR = Integers(ell**i)
    G.append(GL(2, RR))
    
    
    L = []
    M = []
    for x in RR:
        L.append([RR(1), RR(x)])
        if ZZ(x) % ell ==  0:
            L.append([RR(x),RR(1)])
        for y in RR:
            M.append([RR(x),RR(y)])
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


#checks for lift-exceptional subgroups as in Definition 3.5
#subgroup labels and generators are fed in from the file stringN.txt
def check_lift_exceptional_subgroups(n, string="gl2_"):
    N = ell**n
    bads = []
    for line in open(string + str(N) + ".txt", "r"):
        label = line.split(":")[0]
        gens = line.split(":")[1]
        Hgens = ast.literal_eval(gens)
        #print('checking ', label)

        #is it contained in a borel mod ell^{n-1}?
        borels_nminus1 = borels(Hgens, n-1)
        if len(borels_nminus1) > 0:
        
            #does it satisfy the conclusions?
            if not check_conclusions(n, label, Hgens, borels_nminus1):

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

    RR = Integers(ell**n)
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

#checks to see if subgroup of ell^n torsion satisfies
#Definition 3.5 parts (i), (iii), (iv), (v)
def check_conclusions(n, label, Hgens, borels_nminus1):

    if in_cartan(Hgens, 1):
        return True

    borels_n = borels(Hgens, n, previous = True, previousB = borels_nminus1)
    if len(borels_n) > 0:
        return True
            
    if in_A_minimal(Hgens, n, borels_nminus1):
        return True
    
    return False


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
    R = Integers(ell**n)
    return [[R(h[0]), R(h[1]), R(h[2]), R(h[3])] for h in Hgens]

def independent_vectors(v, n):
    for w in V[n]:
        if v[0]*w[1] - v[1]*w[0] % ell != 0:
            break
    ind_vecs = []
    for i in range(ell):
        for j in range(1,ell):
            ind_vecs.append([i*v[0] + j*w[0], i*v[1] + j*w[1]])
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
#of the form A_{a, a+1, n-a-1}(ell^n) as in Definition 3.5 part (v)
#the input borels_nminus1 is a list of fixed lines mod ell^{n-1}
def in_A_minimal(Hgens, n, borels_nminus1):
    
    Hgens_nminus1 = modell(n-1,Hgens)
    for v in borels_nminus1:
        iv = independent_vectors(v, n-1)
        for w in iv:
            CB_Hgens = change_basis(Hgens_nminus1, [v, w])
            beta_min_val = n-1
            delta_minus_alpha_min_val = n-1
            for h in CB_Hgens:
                assert h[2] % ell**(n-1) == 0
                beta_min_val = min(beta_min_val, val(h[1],n-1))
                delta_minus_alpha_min_val = min(delta_minus_alpha_min_val, val(h[3] - h[0],n-1))
            a = delta_minus_alpha_min_val - beta_min_val
            b = n - 1 - a
            c = a + 1
            #print(a,b,c)
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

    return False


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
#The subgroups R(ell^3) and R(ell^5)
#################################################################


#Precompute a list consisting of tuples [label, gens, cosreps], where
#label is one of the maximal lift-exceptional labels, and gens is a list
#of generators of the corresponding subgroup, and cosreps is a list of
#coset representatives of the group in GL2(Z/NZ)

Max3 = [['2788', make_gens_from_label(3, '2788', 'data/gl2_'), reduce_list(make_list_elements("data/2788_27_cosets"),3)]]

Max5 = [['81376', make_gens_from_label(5, '81376', 'data/gl2_'), reduce_list(make_list_elements("data/81376_243_cosets"),5)]]

Maxs = {2: [], 3: Max3, 4: [], 5: Max5}

#Verify that 2788 is R(3^3) and 81376 is R(5^3)

R27_gens = [[2,0,0,25],[1,1,9,1], [2,3,0,2], [1,0,0,10]]
R27 = G[3].subgroup([G[3](v) for v in R27_gens])
assert(R27.order() == 2**2*3**6)
assert(is_subgroup_from_list_cosets(3, 'R27', R27_gens, '2788', [g for g in enumerate_elements(Max3[0][1], 3)], Max3[0][2]))


R243_gens = [[2,0,0,241],[1,1,81,1], [2,3,0,2], [1,0,0,28]]
R243 = G[5].subgroup([G[5](v) for v in R243_gens])
assert(R243.order() == 2**2*3**11)
assert(is_subgroup_from_list_cosets(5, 'R243', R243_gens, '81376', [g for g in enumerate_elements(Max5[0][1], 5)], Max5[0][2]))

#################################################################
#Verify Theorem 3.7 for ell = 3 and n <= 5
#################################################################

print('This file verifies Theorem 3.7 of "A local-global principle for isogenies of composite degree" in the cases where ell is 3 and n is at most 5.  It will take approximately 7 minutess to run.')
for n in range(2,6):
    N = ell**n
    print('Checking which subgroups are lift-exceptional mod ', N)
    Max_list = Maxs[n]
    Max_label = [Max_list_item[0] for Max_list_item in Max_list]
    lift_exceptionals = check_lift_exceptional_subgroups(n, string="data/gl2_")
    print('All done checking for lift-exceptional groups mod ', N, '.')
    New = sieve_maxs(n, lift_exceptionals, Max_list)
    if New:
        print('Every exceptional group is contained in: ', Max_label, '.')
    else:
        print('The exceptional groups: ', New, ' are not contained in ', Max_label, '.')
  
