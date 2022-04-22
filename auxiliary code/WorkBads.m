//This file contains several functions used to match subgroups 
//to those in the work of Rouse-Zureick-Brown and Sutherland-Zywina.  
//For that reason it is necessary to have the files GL2Invariants.m
//and g1groups.m from https://math.mit.edu/~drew/SZ16/ in the same
//directory.  The file rzbData.m is a table with relevant information
//from gl2data.gz available at http://users.wfu.edu/rouseja/2adic/

//bad2.m is a table of information on exceptional subgroups of GL_2(Z/2^nZ)
//for n <= 6.

load "bad2.m";
load "GL2Invariants.m";
load "rzbData.m";
load "g1groups.m";

//Constructs the subgroup of GL(2, ZZ/level) with given label
function GL2Group(tab, label)
	return sub<GL(2, Integers(tab[label]`level))|tab[label]`gens>;
end function;

//Constructs the transpose of subgroup with given label
//So as to match RZB setup.
function GL2TransposeGroup(tab, label)
	return sub<GL(2, Integers(tab[label]`level))|tab[label]`tgens>;
end function;

function szGroup(label)
	return sub<GL(2, Integers(gl2tab[label]`gl2level))|gl2tab[label]`gens>;
end function;

function rzbMatch(tab, label)
	return [k:k in Keys(rzbtab)|rzbtab[k]`level eq tab[label]`level and rzbtab[k]`index eq tab[label]`index and IsConjugate(GL(2, Integers(rzbtab[k]`level)), GL2Group(rzbtab, k), GL2TransposeGroup(tab, label))];
end function;

function szMatch(tab, label)
	return [k:k in Keys(gl2tab)|gl2tab[k]`gl2level eq tab[label]`level and gl2tab[k]`index eq tab[label]`index and IsConjugate(GL(2, Integers(tab[label]`level)), szGroup(k), GL2Group(tab, label))];
end function;

procedure analyzeBads(tab)
	for k in Keys(tab) do
		H := GL2Group(tab, k);
		indexDet := GL2DetIndex(H);
		genus := GL2Genus(H);
		printf "The group %o has determinant index %o, genus %o, rzb label %o, sz label %o, and ", k, indexDet, genus, rzbMatch(tab, k), szMatch(tab, k);
		if GL2ContainsCC(H) then
			print "it contains complex conjugation. \n";
		else
			print "it does not contain complex conjugation. \n";
		end if;
	end for;
end procedure;



function sm(L)
	return "\\" * "sm{" * Sprint(L[1]) * "}{" * Sprint(L[2]) * "}{" * Sprint(L[3]) * "}{" * Sprint(L[4]) * "}";
end function;

function minimizeGens(N, L)
	Gn := GL(2, Integers(N));
	H := sub< Gn | L >;
	Hord := Order(H);
	for l in L do
		if Order(sub< Gn | Exclude(L, l)>) eq Hord then Exclude(~L, l); end if;
	end for;
	return L;
end function; 

procedure TexFormat(tab)
    for k in Keys(tab) do
        H := GL2Group(tab, k);
        indexDet := GL2DetIndex(H);
        genus := GL2Genus(H);
		min_gens := minimizeGens(2^tab[k]`n, tab[k]`gens);
		sm_gens := "";
		for v in min_gens do
			sm_gens := sm_gens * "," * Sprint(sm(v));
		end for;

		R := rzbMatch(tab, k);
		if #R eq 1 then rl := R[1]; else rl := ""; end if;
		S := szMatch(tab, k);
		if #S eq 1 then sl := S[1]; else sl := ""; end if;

		Write("bad_tex.txt",  Sprint(k) * " & " * Sprint(tab[k]`n) * " & " * Sprint(tab[k]`level) * " & " * Sprint(genus) * " & $" *  Sprint(sm_gens) * "$ & " * Sprint(sl) *" & "* Sprint(rl) *"\\"*"\\"*"[5pt]");

	end for;
end procedure;

procedure InvestigateIndices(tab)
	for k in Keys(tab) do
		H := GL2Group(tab, k);
		S := SL(2, Integers(tab[k]`level));
		Gammabar := sub<S | H meet S>;
		Gam_gens := [g :  g in Generators(Gammabar)];
		sl2Index := Order(S)/Order(Gammabar);
		if not sl2Index eq tab[k]`index then printf "Indices for group %o are different: %o and %o \n", tab[k]`label, sl2Index, tab[k]`index; end if;
	end for;
end procedure;


procedure GroupRam(H)
	R := BaseRing(H);
	S := SL(2, R); 
	Gammabar := sub< S | H meet S>;
	Gam_gens := [g : g in Generators(Gammabar)];
	sl2Index := Order(S) / Order(Gammabar);
	gb_gens := "[";
	for v in Prune(Gam_gens) do
		gb_gens := gb_gens * "[" * Sprint(v[1][1])*","*Sprint(v[1][2])*","*Sprint(v[2][1])*","*Sprint(v[2][2])*"],";
	end for;
	w :=Gam_gens[#Gam_gens];
	gb_gens := gb_gens * "[" * Sprint(w[1][1])*","*Sprint(w[1][2])*","*Sprint(w[2][1])*","*Sprint(w[2][2])*"]]";
	print Sprint(#R)*":"*Sprint(sl2Index)*":"* gb_gens * ":" * Sprint(GL2Nu2(H)) * ":" * Sprint(GL2Nu3(H)) * ":" *  Sprint(GL2Genus(H));
end procedure;

procedure WriteRam(tab)
	for k in Keys(tab) do
		H := GL2Group(tab, k);
		S := SL(2, Integers(tab[k]`level));
		Gammabar := sub<S | H meet S>;
		Gam_gens := [g :  g in Generators(Gammabar)];
		sl2Index := Order(S)/Order(Gammabar);
		if not sl2Index eq tab[k]`index then printf "Indices are different: %o and %o", sl2Index, tab[k]`index; end if;
		gb_gens := "[";
		for v in Prune(Gam_gens) do
			gb_gens := gb_gens * "[" * Sprint(v[1][1])*","*Sprint(v[1][2])*","*Sprint(v[2][1])*","*Sprint(v[2][2])*"],";
		end for;
		w :=Gam_gens[#Gam_gens];
		gb_gens := gb_gens * "[" * Sprint(w[1][1])*","*Sprint(w[1][2])*","*Sprint(w[2][1])*","*Sprint(w[2][2])*"]]";

		Write("ram_data.txt", Sprint(k)*":"*Sprint(tab[k]`n)*":"*Sprint(tab[k]`level)*":"*Sprint(sl2Index)*":"* gb_gens * ":" * Sprint(GL2Nu2(H)) * ":" * Sprint(GL2Nu3(H)) * ":" *  Sprint(GL2Genus(H)));
	end for;
end procedure;

procedure LowerGenus(Gens,N,  CosetReps)
	GN := GL(2, Integers(N));
	genus := GL2Genus(sub< GN|Gens>);
	for g in CosetReps do
		if GL2Genus(sub<GN | Gens , g>) in [2..(genus-1)] then
			//for h in Exclude(CosetReps,g) do
			//	if GL2Genus(sub<GN | Gens , g, h>) eq 2 then
			//		print g, h;
			//	end if;
			//end for;
			print g;
		end if;
	end for;
end procedure;






