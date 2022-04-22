from sage.all import *

import ast, sys

from sage.modular.arithgroup.congroup_generic import CongruenceSubgroup_constructor as CS

#################################################################
#This file verifies Theorem 1 part (1) of
#''A local-global principle for isogenies of composite degree''
#by Isabel Vogt

#This file contains the functions necessary to compute genera of
#fiber products of modular curves X_Gamma1 and X_Gamma2, given 
#nu2, nu3, and ramification profiles over infinity for the modular 
#maps X_Gamma1 —> X(1) and X_Gamma2 —> X(1).  This is contained in
#the function fiber_genus.

#The function cusp_ramification computes the ramification
#profile over infinity for Gamma from a set of generators of 
#Gammabar = subgroup of SL_2(Z/NZ) determining Gamma as a congruence
#subgroup.  This is done by considering the orbits of cusps of Gamma(N)
#under the action of Gamma.  This is contained in the function
#profile cusps.

#The main application of these functions is to determine for which M
#such that X_0(M) is genus 0 or 1, is the fiber product 
# X_Gamma \times_{X(1)} X_0(M)
#genus 0 or 1, where Gamma is the congruence subgroup corresponding
#to an exceptional subgroup.  This application is achieved by
#by the function test_bad_file.
#################################################################



#PRECOMPUTE SOME DATA FOR GAMMA0(N):
DataDictionary = {}
#Data4 = [Gamma0(4).index(), Gamma0(4).nu2(), Gamma0(4).nu3(), cusp_ram_profile(cusp_ramification([[1,1,0,1],[3,0,0,3]], 2, 2))]
DataDictionary[4] = [6, 0, 0, [1, 4, 1]]

#Data9 = [Gamma0(9).index(), Gamma0(9).nu2(), Gamma0(9).nu3(),cusp_ram_profile(cusp_ramification([[1,1,0,1],[2,0,0,5]], 3, 2))]
DataDictionary[9] = [12, 0, 0, [1, 9, 1, 1]]

#Data27 = [Gamma0(27).index(), Gamma0(27).nu2(), Gamma0(27).nu3(),cusp_ram_profile(cusp_ramification([[1,1,0,1],[2,0,0,14]], 3, 3))]
DataDictionary[27] = [36, 0, 0, [1, 27, 1, 3, 3, 1]]

#Data25 = [Gamma0(25).index(), Gamma0(25).nu2(), Gamma0(25).nu3(),cusp_ram_profile(cusp_ramification([[1,1,0,1],[2,0,0,13]], 5, 2))]
DataDictionary[25] = [30, 2, 0, [1, 25, 1, 1, 1, 1]]

#Data15 = [Gamma0(15).index(), Gamma0(15).nu2(), Gamma0(15).nu3(),fiber_ram_profile([1,3], [1,5])]
DataDictionary[15] = [24, 0, 0, [1, 5, 3, 15]]



#This is the main function that inputs a list of
#congruence subgroups corresponding to exceptions for
#ell^n-isogenies, and outputs the set of all integers
#M*ell^n such that the fiber product modular curve 
# X_Gamma \times_{X(1)} X_0(M)
#has genus 0 or 1.
def test_bad_file(file_name):
    L = set([])
    small_gps = []
    data_list = make_list_of_data(file_name)
    #print data_list
    for d in data_list:
        if d[-1] <= 1:
            small_gps.append(d)
            ell = d[2].prime_factors()[0]
            L.add(ell**d[1])
    for badgp in small_gps:
        #print("Working on group " +str(badgp[0]))
        still_low_genus = False
        small =[]
        ell = ZZ(badgp[2]).prime_factors()[0]
        #print(ell, badgp)
        gpData = makeData(badgp)
        
        #first check that no other badgp has fiber product with this one
        for badgp_prime in small_gps:
            ell_prime = badgp_prime[2].prime_factors()[0]
            if ell_prime != ell:
                gp_primeData = makeData(badgp_prime)
                if fiber_genus(gpData, gp_primeData) <= 1:
                    print('Warning: found small genus fiber product from', badgp[0], 'and', badgp_prime[0],'.')
                    L.append(ell**badgp[1] * ell_prime**badgp_prime[1])
        
        #now check fiber products with X_0(p)
        for p in prime_range(100):
            if (p != ell) and (Gamma0(p).genus() <= 1):
                #print(p)
                if fiber_genus(gpData,[p+1, Gamma0(p).nu2(), Gamma0(p).nu3(), [1,p]]) <= 1:
                    small.append(p)
                    still_low_genus = True
        new = small[:]
        while still_low_genus == True:
            
            #print(new)
            prod = []
            small_genus = []
            for x in new:
                for y in new:
                    if not x*y in prod:
                        prod.append(x*y)
                    while prod:
                        N = prod.pop()
                        if Gamma0(N).genus() <= 1:
                            DataN = DataDictionary[N]
                            if fiber_genus(gpData, DataN) <= 1:
                                small_genus.append(N)
                                small.append(N)
            new = small_genus[:]
            still_low_genus = bool(len(new) >= 1)

        print("For group " + str(badgp[0]) + " the following N give infinitely many counter-examples: " + str(small))
        for z in small:
            L.add(z*(ell**badgp[1]))

    return L




#Given a file with the data of label, n, level, index,
#generators for Gammabar, nu2, nu3, and genus for
#congruence groups Gamma, puts this into a list.
def make_list_of_data(file_name):
    data = []
    for line in open(file_name + ".txt", "r"):
        label = line.split(":")[0]
        nn = ZZ(line.split(":")[1])
        level = ZZ(line.split(":")[2])
        index = ZZ(line.split(":")[3])
        Gammabar_gens = ast.literal_eval(line.split(":")[4])
        nu2 = ZZ(line.split(":")[5])
        nu3 = ZZ(line.split(":")[6])
        genus = ZZ(line.split(":")[7])
        data.append([label, nn, level, index, Gammabar_gens, nu2, nu3, genus])
    return data



#makes the list called Data to be inputted to fiber_cusps function below
def makeData(L):
    ell = ZZ(L[2]).prime_factors()[0]
    n = ZZ(L[2]).valuation(ell)
    return [L[3], L[5], L[6], cusp_ram_profile(cusp_ramification(L[4], ell, n))]



#computes a set of representatives [a,b] with a and b relatively prime
#for the cusps of Gamma(N=ell^n)
def cusp_reps(ell, n):
    N = ell**n
    cusps =[]
    if (ell, n) == (2,1):
        return [(1,0), (0,1), (1,1)]
    if ell == 2:
        B = ell**(n-1)
    else:
        B = (N+1)/2
    for a in range(1,B):
        if a % ell != 0:
            for b in range(N):
                cusps.append(rel_prime_lift(a,b,N))
        else:
            for b in range(N):
                if b % ell != 0:
                    cusps.append(rel_prime_lift(a,b,N))
    for b in range(B):
        if b % ell !=0:
            cusps.append(rel_prime_lift(0,b,N))
    if ell == 2:
        for b in range(B):
            if b % ell != 0:
                cusps.append(rel_prime_lift(2**(n-1),b,N))
    return cusps

#given integers a,b,N relatively prime
#returns a lift of (a,b) mod N to the integers, which is relatively prime
def rel_prime_lift(a,b,N):
    if a !=0:
        t = 1
        for p in ZZ(a).prime_divisors():
            if b % p != 0:
                t *= p
        return (a, b + t*N)
    else:
        return (N, b)




#Gives a lift of the matrix M =[a,b,c,d] in SL_2(Z/NZ) to SL_2(Z)
def matrix_lift(M,NN):

    a,b,c,d = ZZ(M[0]), ZZ(M[1]), ZZ(M[2]), ZZ(M[3])

    #print(a,b,c,d, NN)

    cp, dp = rel_prime_lift(c, d, NN)
    g, x, y = xgcd(cp, dp)
    assert g == 1

    current_det = a*dp - b*cp

    diff = (current_det-1)/NN

    newM = [a-NN*y*diff, b+NN*x*diff, cp, dp]
    return newM

#Two rational numbers are equivalent under Gamma(N)
#iff they are congruent up to a sign.
def cusp_equiv(v1, v2, N):
    return ((v1[0]-v2[0])%N ==0 and (v1[1]-v2[1])%N ==0) or ((v1[0]+v2[0])%N ==0 and (v1[1]+v2[1])%N ==0)
    
#matches a rational number with a cusp representative
#for the action of Gamma(N) in a given list L
def match_cusp_index(v, L, N):
    for i in range(len(L)):
        if cusp_equiv(v, L[i], N):
            return i

#genList is a list of generators of Gammabar subset SL_2(Z/NZ)
#at the moment N must be a prime power
#returns the list of orbits of cusps of X(N) under the action
#of Gamma, with multiplicity equal to the number of cusps
#identified to that orbit.
def profile_orbits(genList, ell, n):

    N = ell**n
    genLifts = [matrix_lift(g, N) for g in genList]
    GammaN_orbits = cusp_reps(ell, n)
    #print "Orbits are", GammaN_orbits

    OrbitIndices = [i for i in range(len(GammaN_orbits))]

    #we'll keep track of cusp reps by their position in the list GammaN_orbits

    for h in genLifts:
        #print "h is", h
        for oi in range(len(OrbitIndices)):
            #print "Checking index", oi
            ho = mvmult(h, GammaN_orbits[oi])
            #print "The h translate is", ho
            ho_index = match_cusp_index(ho, GammaN_orbits, N)
            #print "The index is", ho_index
            if oi != ho_index:
                #print "Indices are", OrbitIndices
                #print (oi, ho_index)
                #k = max(oi,ho_index)
                new = [min(oi, ho_index)]
                to_do = [oi, ho_index]
                traced = []
                #same = False
                while to_do:
                    k = to_do.pop()
                    traced.append(k)
                    #print "traced is", traced
                    j =OrbitIndices[k]
                    if j != k:
                        new.append(j)
                        to_do.append(j)
                for ii in traced:
                    #print "hi"
                    OrbitIndices[ii] = min(new)
    for i in range(len(OrbitIndices)):
        r = OrbitIndices[i]
        if r != i:
            s = i
            value_same = False
            values = [r]
            while not value_same:
                s = OrbitIndices[r]
                values.append(s)
                if s != r:
                    r = s
                else:
                    value_same = True
            OrbitIndices[i] = min(values)
                
    return OrbitIndices


#multiplies the 2x2 matrix M represented as a list by the vector v
def mvmult(M,v):
    return [M[0]*v[0] + M[1]*v[1], M[2]*v[0]+M[3]*v[1]]
            
#genLst generates GammaBar is the subgroup of SL_2(Z/NZ) 
#determining the congruence group Gamma
#determines the size of Gamma/Gamma(N = ell^n)
def ambient_size(genList, ell, n):

    S = SL(2, Integers(ell**n))
    GammaBar = S.subgroup([S(m) for m in genList])
    if S([-1,0,0,-1]) in GammaBar:
        f = 2
    else:
        f = 1

    return GammaBar.order()/f


#manipulates the output of profile_orbits to return
#a list (preimage of infinity, ramification index) for
#the branched cover H/Gamma ——> H/Gamma(1)
#indexed by their position in the list of preimages of 
#infinity of H/Gamma(N) ——> H/Gamma(1)
#genList is a list of generators of Gammabar
def cusp_ramification(genList, ell, n):

    ram_data = []
    
    N = ell**n

    cusps = profile_orbits(genList, ell, n)
    for c in Set(cusps):
        ctr = 0
        for i in cusps:
            if i == c:
                ctr += 1
        ram_data.append([c,ZZ((N*ctr)/ambient_size(genList, ell, n))])

    return ram_data 


#Converts the list output L of cusp_ramification to a list
#of e_i giving the ramification index at every point over
#infinity 
def cusp_ram_profile(L):
    return [l[1] for l in L]


#L1 is a ramification profile of X_1 over a point b in B
#L2 is a ramification profile of X_2 over a point b in B
#returns the ramification profiles of X_1 \times_B X_2 over 
#the points above b in X_1
def fiber_cusps(L1, L2):
    profile = []
    for k in L1:
        k_profile = []
        for j in L2:
            j = ZZ(j)
            #want to add to k_profile gcd(j,k) points with e = j/gcd(j,k)
            d = gcd(j,k)  
            for i in range(d):
                k_profile.append(j/d)
        profile.append(k_profile)
    return profile


#L1 is a ramification profile of X_1 over a point b in B
#L2 is a ramification profile of X_2 over a point b in B
#return the cusp ramification profile of the points over infinity for the map
#from the fiber product to X(1)
def fiber_ram_profile(L1, L2):
    profile = []
    for k in L1:
        for j in L2: 
            d = gcd(j,k)
            for i in range(d):
                profile.append((j*k)/d)
    return profile


#Data is the list [degree, nu2, nu3, [ram profile]]
#We'll compute the genus using the Riemann-Hurwitz formula 
#g = 1 + (1/2)*(d2 * (2g1-2) + Sum ( e_i - 1))
def fiber_genus(Data1, Data2):
    #first compute g1 using the Hurwitz formula, which simplifies to 
    #g1 = 1 + d1/12 - (nu2_1/4) - (nu3_1/3) - (nuinfty_1/2)
    g1 = 1 + Data1[0]/12 - Data1[1]/4 - Data1[2]/3 - ZZ(len(Data1[3]))/2
    #print("g1 is", g1)
    g = 1 + Data2[0]*(2*g1-2)/2

    #now we need to add in the ramification contribution
    #start by the ramification contribution of points over 1728
    g += (Data1[1]*(Data2[0]-Data2[1])/4)
    #print "The number of ramified for X2 is",ZZ((Data2[0]-Data2[1])/2)
    #print "The contribution over 1728 is", (Data1[1]*(Data2[0]-Data2[1])/4)

    #now the ramification contribution of points over 0 
    g += (Data1[2]*(Data2[0]-Data2[2])/3)
    #print "The contribution over 0 is", (Data1[2]*(Data2[0]-Data2[2])/3)


    #finally we add in the contributions of the cusps
    cusp_ram_profiles = fiber_cusps(Data1[3], Data2[3])
    for l in cusp_ram_profiles:
        for e in l:
            g += (e-1)/2

    g = ZZ(g)
    
    return g


#################################################################
#Verify Proposition 5.1 parts (1) and (2)
#################################################################
print('This file verifies the computational assertions in Theorem 1 of "A local-global principle for isogenies of composite degree".  It will take approximately 5 seconds to run. For each maximal exceptional subgroup, the function returns the possible N for which the fiber product of that subgroup with X_0(N) has genus 0 or 1.  The function also checks that no fiber product of two exceptional modular curves can have genus 0 or 1.', '\n')
L = list(test_bad_file("data/ram_data"))
L.sort()
print('\n', 'For N not in the finite list', L, 'the list Sigma(K, N) is always finite.')
