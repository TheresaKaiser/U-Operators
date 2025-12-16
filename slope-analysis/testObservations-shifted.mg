load "allslopes.txt";

listToMultiset := function(L)
	res := {* *};
	for tup in L do
		Include(~res, tup[1]^^tup[2]);
	end for;
	return res;
end function;

shiftValues := function(q, type, i, L)
	correction := Rationals() ! - (3-i)*((type-1) mod (q-1));
	for j -> tup in L do
		if tup[1] ne -1 then
			L[j][1] := tup[1]+correction;
		end if;
	end for;
	return L;
end function;



for current in AllSlopes do

q := current[1];
type := current[2];
if type eq 0 then
	n := q-1;
else
	n := type;
end if;
L1 := current[3];
L2 := current[4];

print "==========================================================================";
printf "Analyzing q=%o, type=%o. \nThe weight k goes up to %o.\n \n", q, type, L1[#L1][1];

fprintf "python1mult.txt", "\n],\n";
fprintf "python2mult.txt", "\n],\n";
fprintf "python1mult.txt", "(%o, %o, 1): [ \n    ", q, type;
fprintf "python2mult.txt", "(%o, %o, 2): [ \n    ", q, type;

assert #L1 eq #L2;


// I know the following, which helps the analysis: 
// * both lists contain the same weights in the same order
// * each list is sorted from smaller to bigger slopes
// * all finite slopes are >= 0, the slope infinity is given as -1
// * each slope appears at most once in each list
// I might use and assume this knowledge without further comment.
// Use lists when we need to make use of the ordering, multisets for easier comparison of multisets 

//compare nth smallest Ui slopes as k varies!
C := car< Rationals(), Integers() >;
smallest1 := [C | <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>];
smallest2 := [C | <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>, <0, -1>];

// find largest denominators
denominators1 := {* *};
denominators2 := {* *};

for counter := 1 to #L1 do
	tup1 := L1[counter];
	tup2 := L2[counter];
		
	k := tup1[1];
	assert k eq tup2[1];
	
	ListT1 := shiftValues(q, type, 1, tup1[2]);
	T1 := listToMultiset(ListT1);
	ListU1P0 := shiftValues(q, type, 1, tup1[3]);
	U1P0 := listToMultiset(ListU1P0);
	ListU1P2 := shiftValues(q, type, 1, tup1[4]);
	U1P2 := listToMultiset(ListU1P2);
	ListU1 := shiftValues(q, type, 1, tup1[5]);
	U1 := listToMultiset(ListU1);
	
	ListT2 := shiftValues(q, type, 2, tup2[2]);
	T2 := listToMultiset(ListT2);
	ListU2P0 := shiftValues(q, type, 2, tup2[3]);
	U2P0 := listToMultiset(ListU2P0);
	ListU2P2 := shiftValues(q, type, 2, tup2[4]);
	U2P2 := listToMultiset(ListU2P2);
	ListU2 := shiftValues(q, type, 2, tup2[5]);
	U2 := listToMultiset(ListU2);
	
	assert #T1 eq #T2;
	assert #U1P0 eq #U1P2 and #U2P0 eq #U1P2 and #U2P0 eq #U2P2;
	assert #U1 eq #U2;
	
	// Export slopes so we can plot them in python (ignoring multiplicities for now)
	for tup in ListU1 do
		fprintf "python1mult.txt", "(" cat Sprint(k) cat ", " cat Sprint(tup[1]) cat ", " cat Sprint(tup[2]) cat "), ";
	end for;
	for tup in ListU2 do
		fprintf "python2mult.txt", "(" cat Sprint(k) cat ", " cat Sprint(tup[1]) cat ", " cat Sprint(tup[2]) cat "), ";
	end for;	
	
	dimG := #T1;
	dimP := #U1P0;
	dimD := #U1;
	
	// check that zero (slope infty = -1 here) does not appear for Ti
	assert Multiplicity(T1, -1) eq 0;
	assert Multiplicity(T2, -1) eq 0;
	
	// check containment and multiplicities
	assert (T1 join {* (-1)^^(2*dimG) *}) subset U1P0;
	assert (T1 join T1 join {* (-1)^^(dimG) *}) subset U1P2;
	assert (T2 join T2 join {* (-1)^^(dimG) *}) subset U2P0;
	assert (T2 join {* (-1)^^(2*dimG) *}) subset U2P2;
	
	assert (T1 join T1 join {* (-1)^^(4*dimG) *}) subset U1;
	assert (T2 join T2 join {* (-1)^^(4*dimG) *}) subset U2;
	
	Old1 := (U1P0 join U1P2 join {* (-1)^^dimP *}) diff (T1 join {* (-1)^^(2*dimG) *});
	Old2 := (U2P0 join U2P2 join {* (-1)^^dimP *}) diff (T2 join {* (-1)^^(2*dimG) *});
	assert Old1 subset U1;
	assert Old2 subset U2;
	
	new1 := 2*(k/3 + 1-n);
	new2 := k/3 + 1-n;
	
	if new1 in Old1 then
		printf "Newform eigenvalue detected in oldforms for k=%o, i=1.\nMultiplicities in T1:%o, U1P0:%o, U1P2: %o.\n", k, Multiplicity(T1, new1), Multiplicity(U1P0, new1), Multiplicity(U1P2, new1);
	end if;
	if new2 in Old2 then
		printf "Newform eigenvalue detected in oldforms for k=%o, i=2.\nMultiplicities in T2:%o, U2P0:%o, U2P2: %o.\n", k, Multiplicity(T2, new2), Multiplicity(U2P0, new2), Multiplicity(U2P2, new2);
	end if;

	assert IsEmpty(U1 diff Old1) or MultisetToSet(U1 diff Old1) eq {new1};
	assert IsEmpty(U2 diff Old2) or MultisetToSet(U2 diff Old2) eq {new2};
	
	// count all denominators appearing in Ui-slopes (all slopes)
	for pair in ListU1 do
		Include(~denominators1, Denominator(pair[1]));
	end for;
	for pair in ListU2 do
		Include(~denominators2, Denominator(pair[1]));
	end for;
	
	// fails for q=3
	// compare largest finite Ui slope
	if #ListU1 gt 0 then
		max1 := ListU1[#ListU1][1];
		max2 := ListU2[#ListU2][1];
		if max1 lt 0 and #ListU1 gt 1 then
			max1 := ListU1[#ListU1-1][1];
		end if;
		if max2 lt 0 and #ListU2 gt 1 then
			max2 := ListU2[#ListU2-1][1];
		end if;
		if max1 ge 0 and max2 ge 0 then
			//fprintf "coords.txt", "(" cat Sprint(k) cat ", " cat Sprint(k-max2) cat "), ";
			if not max1 eq 2*max2 then
				print "not max1 eq 2*max2", k, max1, max2;
			end if;
			if k gt 5 then
				assert k - max2 ge 4;
			else
				assert max2 le k+1;
			end if;
			//k-max2;
		end if;
	end if;
	
	
	//compare nth smallest Ui slopes (as k varies!)
	for n := 1 to #smallest1 do
		if #ListU1 ge n then
			if (ListU1[n][1] gt smallest1[n][1]) or (ListU1[n][1] eq smallest1[n][1] and smallest1[n][2] eq -1) then // then it is automatically >= 0
				smallest1[n][1] := ListU1[n][1];
				smallest1[n][2] := k;
			end if;
		end if;
		if #ListU2 ge n then
			if (ListU2[n][1] gt smallest2[n][1]) or (ListU2[n][1] eq smallest2[n][1] and smallest2[n][2] eq -1) then // then it is automatically >= 0
				smallest2[n][1] := ListU2[n][1];
				smallest2[n][2] := k;
			end if;
		end if;
	end for;
end for;

printf "Assertions regarding the containment and multiplicities of slopes were successful for all weights. \n";

printf "\nUpper bounds for the n'th smallest slopes are the following. \n";
print "i=1: ";
for n := 1 to #smallest1 do
	printf "n=%o: %o, first appeared at k=%o.\n", n, smallest1[n][1], smallest1[n][2];
end for;
print "i=2: ";
for n := 1 to #smallest2 do
	printf "n=%o: %o, first appeared at k=%o.\n", n, smallest2[n][1], smallest2[n][2];
end for;

printf "\nMultiset containing all appearing denominators: \n";
denominators1;
denominators2;


end for;
quit;
