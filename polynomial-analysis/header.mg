// getSlopesList: A function calculating the slopes from the coefficients of 
// the given polynomial. The built-in Magma function for slopes does not
// return the slope infinity, so we calculate the slopes ourselves from
// the newton polygon that Magma generates.
// Since we know that all our slopes are non-negative, we use -1 to mark 
// the slope infinity.
// Input: 
	// f the polynomial
    // A the polynomial ring whose t-adic valuation we use
// Output: 
	// Slopes a list containing 2-tuples of slopes with multiplicities
getSlopesList := function(f, A)
	Slopes := [];
	V := [];
	for j := 0 to Degree(f) do
		c := Coefficient(f, j);
		if c ne 0 then
			Append(~V, <j, Valuation(A ! c)>);
		end if;
	end for;    
	// Magma gives us a list of the lower vertices of the Newton polygon
	p := NewtonPolygon(V : Faces := "Lower");
	lv := LowerVertices(p);
	// Compute the slope of each segment
	// Need to change sign due to the order we listed our points in
	for j := 0 to #lv-2 do
		current := lv[#lv - j];
		next := lv[#lv - j - 1];
		length := current[1] - next[1];
		slope := (next[2] - current[2])/length;
		Append(~Slopes, <slope, length>);
	end for;
	// The remaining segment belongs to the slope infinity. Modeled as -1
	if lv[1][1] ne 0 then
		Append(~Slopes, <-1, lv[1][1]>);
	end if;
	return Slopes;
end function;


getSlopesLaTeX := function(f, A, i, k : shift:=0)
	Slopes := "";
	if not IsZero(f) then
		V := [];
		for j := 0 to Degree(f) do
				c := Coefficient(f, j);
				if c ne 0 then
						Append(~V, <j, Valuation(A ! c)>);
				end if;
		end for;
		// Magma gives us a list of the lower vertices of the Newton polygon
		p := NewtonPolygon(V : Faces := "Lower");
		lv := LowerVertices(p);
		// Compute the slope of each segment
		// Need to change sign due to the order we listed our points in
		for j := 0 to #lv-2 do
				current := lv[#lv - j];
				next := lv[#lv - j - 1];
				length := current[1] - next[1];
				slope := (next[2] - current[2])/length;
				if (i eq 1 and (Rationals() ! slope) eq (Rationals() ! (2*k)/3)) or (i eq 2 and (Rationals() ! slope) eq (Rationals() ! k/3)) then
						Slopes := Slopes cat "\\textcolor{fullblue}{$" cat Sprint(slope+shift) cat "^{\\mathbf{" cat Sprint(length) cat "}}$}, ";
				else
						Slopes := Slopes cat "$" cat Sprint(slope+shift) cat "^{\\mathbf{" cat Sprint(length) cat "}}$, ";
				end if;
		end for;
		// The remaining segment belongs to the slope infinity
		if lv[1][1] ne 0 then
				Slopes := Slopes cat "$\\infty^{\\mathbf{" cat Sprint(lv[1][1]) cat "}}$";
		else // if there is no infinite slope, delete the last comma and space
				Slopes := Substring(Slopes, 1, #Slopes - 2);
		end if;
	end if;
	return Slopes;
end function;


matrixsize := function(k,q,i)
	return Floor((k-2-i)/(q-1))+1;
end function;

UMatrix := function(FF, k, q, j)
	i := j mod (q-1);
	S<t> := PolynomialRing(FF);
	sm := matrixsize(k,q,i);
	A := ZeroMatrix(S, sm, sm);
	
	for i1 := 0 to sm-1 do
		for i2 := 0 to sm-1 do
			A[i2+1,i1+1] := (-t)^(i+i1*(q-1))*Binomial(k-2-i-i2*(q-1),i+i1*(q-1));
			if i2 lt i1 then
				A[i2+1,i1+1] := A[i2+1,i1+1] +(-1)^(i+1)*(-t)^(i+i1*(q-1))*Binomial(k-2-i-i2*(q-1),(i1-i2)*(q-1));
			end if;
		end for;
	end for;
	
	return A;
end function;

charPolRank2 := function(q, k, type) 
	FF := GF(q);
	Pol := CharacteristicPolynomial(UMatrix(FF, k+2,q,type-1));
	r<X> := Parent(Pol);
	return Pol;
end function;

printLine := function(pair, A, i, q, type, filename)
	k := pair[1];
	polTup := pair[2]; // <Pol1D, Pol1P0, Pol1P2, Pol1G>
	
	correction := Rationals() ! - (3-i)*((type-1) mod (q-1));
	
	stringLaTeX := Sprint(k) cat " & ";

	slopeString := getSlopesLaTeX(polTup[4], A, i, k : shift := correction);
	stringLaTeX := stringLaTeX cat slopeString cat " & ";
	slopeString := getSlopesLaTeX(polTup[2], A, i, k : shift := correction);
	stringLaTeX := stringLaTeX cat slopeString cat " & ";
	slopeString := getSlopesLaTeX(polTup[3], A, i, k : shift := correction);
	stringLaTeX := stringLaTeX cat slopeString cat " & ";
	slopeString := getSlopesLaTeX(polTup[1], A, i, k : shift := correction);
	stringLaTeX := stringLaTeX cat slopeString cat " \\\\ \n ";

	fprintf filename, stringLaTeX cat "\\hline \n";
	return stringLaTeX;
end function;


testPair := function(pair, q, i, type, A)
	k := pair[1];
	polTup := pair[2]; // <Pol1D, Pol1P0, Pol1P2, Pol1G>
	
	r<X> := Parent(polTup[1]);
	
	assert Degree(polTup[2]) eq Degree(polTup[3]);
	
	if i eq 1 and (type-1) mod (q-1) eq 0 then
		polRank2 := r! charPolRank2(q, k, type);

		if not IsDivisibleBy(polTup[1], polRank2) then
			print "not IsDivisibleBy(polTup[1], polRank2) at k=", k;
		end if;
		if not IsZero(polTup[2]) then
			g := GCD(polTup[3], polRank2);
			if g ne polRank2 and not X*g eq polRank2 then
				print "g ne polRank2 and not X*g eq polRank2 at k=", k;
			end if;
			g := GCD(polTup[2], polRank2);
			
			commonSlopes := [s[1] : s in getSlopesList(g, A)];
			
			missingSlopes := [s[1] : s in getSlopesList(ExactQuotient(polRank2,g), A)];
			if not IsSubsequence(missingSlopes, [0]) then
				print "not IsSubsequence(missingSlopes, [0]) at k=", k;
				print missingSlopes;
			end if;
			if not IsZero(polTup[4]) then 
				g := GCD(polTup[4], polRank2);
				commonSlopes := [s[1] : s in getSlopesList(g, A)];
				// Need to use ExactQuotient-function to get correct data type for result
				missingSlopes := [s[1] : s in getSlopesList(ExactQuotient(polRank2,g), A)];
				if not IsSubsequence(missingSlopes, [0, k/2, -1]) then
					print "not IsSubsequence(missingSlopes, [0, k/2, -1]) at k=", k;
					print missingSlopes;
				end if;
			end if;
		end if;
		
	end if;
	
	if i eq 2 and q eq 2 and k mod 2 eq 0 then
		polRank2 := r! charPolRank2(q, ExactQuotient(k, 2), type);
		
		g := GCD(polTup[1], polRank2);
		print "missing Slopes Gamma0(t) at k=", k;
		print getSlopesList(ExactQuotient(polRank2,g), A);

		if not IsZero(polTup[2]) then
			g := GCD(polTup[3], polRank2);
			print "missing Slopes P2 at k=", k;
			print getSlopesList(ExactQuotient(polRank2,g), A);

			g := GCD(polTup[2], polRank2);
			print "missing Slopes P0 at k=", k;
			print getSlopesList(ExactQuotient(polRank2,g), A);

			if not IsZero(polTup[4]) then 
				g := GCD(polTup[4], polRank2);
				print "missing Slopes GL3(A) at k=", k;
				print getSlopesList(ExactQuotient(polRank2,g), A);
			end if;
		end if;
		
	end if;

	if not IsZero(polTup[2]) then
		assert IsDivisibleBy(polTup[1], polTup[2]);
		assert IsDivisibleBy(polTup[1], polTup[3]);
		assert IsDivisibleBy(polTup[1], X^Degree(polTup[3]));
		if not IsZero(polTup[4]) then 
			if Degree(polTup[4]) eq 1 then
				printf "q=%o, i=%o, k=%o, type=%o. T_i-Charpol (needs to be renormalized): %o. \n", q, i, k, type, polTup[4];
			end if;
			assert IsDivisibleBy(polTup[1+i], polTup[4]*X^(2*Degree(polTup[4])));
			assert IsDivisibleBy(polTup[4-i], polTup[4]*polTup[4]*X^Degree(polTup[4]));
			assert IsDivisibleBy(polTup[1], polTup[4]*polTup[4]*X^(4*Degree(polTup[4])));
		end if;
	end if;
	return "Test completed successfully for k = " cat Sprint(k);
end function;

