// Multinomialcoefficients with "bad" entries
myMultinomial := function(n, Q)
    for r in Q do
        if r lt 0 or r gt n then
            return 0;
        end if;
    end for;
    return Multinomial(n, Q);
end function;


getDims := function(F, q, k, type)  
	// Initialize cartesian product. Necessary to define tuples of indices lying in s
	C := car< Integers(), Integers() >;
	// List with tuples of indices describing basis of V
	s := [ C | <l, m> : l in [0..k], m in [0..k] | l+m le k ];
    dimV := #s;
	V := VectorSpace(F, dimV);
	sD := [];
	listBasisVectors := [];

	// Initialize other vector spaces to make sure all return values are defined
	VD := VectorSpace(F, 0);
	dimVD := #sD;
	VP0 := VectorSpace(F, 0);
	dimVP0 := 0;
	VP2 := VectorSpace(F, 0);
	dimVP2 := 0;
	VG := VectorSpace(F, 0);
	dimVG := 0;

	for i -> pair in s do
		l := pair[1];
		m := pair[2];
		if IsDivisibleBy(l+1-type, q-1) and IsDivisibleBy(m+1-type, q-1) and IsDivisibleBy(k-l-m+1-type, q-1) then
			Append(~sD, <l,m, i>);
			Append(~listBasisVectors, V.i);
		end if;
	end for;

	dimVD := #sD;

    if dimVD gt 0 then // should always be the case
		VD := VectorSpace(F, dimVD);
		inclVD := Matrix(F, dimVD, dimV, listBasisVectors);
		projVD := Transpose(inclVD);
		assert inclVD*projVD eq IdentityMatrix(F, dimVD);
		
		// Initialize matrices
		MId := IdentityMatrix(F, dimVD);
		M12 := ZeroMatrix(F, dimVD, dimVD);
		M23 := ZeroMatrix(F, dimVD, dimVD);

		MsumA := ZeroMatrix(F, dimVD, dimVD);
		MsumC := ZeroMatrix(F, dimVD, dimVD);

		// Dual iteration on s, allowing to get indices for Magma basis and our theoretical basis.
		for a -> tupA in sD do
			for b -> tupB in sD do
				// Careful: Magma works with row vectors.
				// la, mu: Index in domain -> a is index of row of resulting matrix
				// l, m: Index in codomain -> b is index of column of resulting matrix 
				la := tupA[1];
				mu := tupA[2];
				l := tupB[1];
				m := tupB[2];

				if l eq la then
					if k-(l+m) eq mu then
						M23[a, b] := (-1)^(1-type);
					end if;
					
					MsumC[a, b] := - Binomial(m, mu);
				end if;
				
				if l eq mu then
					if m eq la then
						M12[a, b] := (-1)^(1-type);
					end if;
				end if;
				
				if l+m-la eq mu then
					MsumA[a, b] := - Binomial(l, la);
				end if;		
			end for;    
		end for;

		Mstab12 := MId + MsumA;
		Mstab23 := MId + MsumC;

		// Initialize subspaces for VP0 and VP2
		VP0 := Nullspace(M12 + Mstab12);
		VP2 := Nullspace(M23 + Mstab23);
		
		// Compute subspaces V^{H_i} of V that contain VP0, VP2
		MEta0 := ZeroMatrix(F, dimV, dimV);
		MEta2 := ZeroMatrix(F, dimV, dimV);
		// Dual iteration on s, allowing to get indices for Magma basis and our theoretical basis.
		for a -> tupA in s do
			for b -> tupB in s do
				// Careful: Magma works with row vectors.
				// la, mu: Index in domain -> a is index of row of resulting matrix
				// l, m: Index in codomain -> b is index of column of resulting matrix 
				la := tupA[1];
				mu := tupA[2];
				l := tupB[1];
				m := tupB[2];
				
				i := la - l;
				if mu eq m-i then
					MEta0[a, b] := Binomial(m, i);
				end if;
				
				if l eq la then
					i := mu - m;
					MEta2[a, b] := Binomial(k-l-m, i);
				end if;
			end for;    
		end for;
		VP0 meet:= Kernel(inclVD*(MEta0 - IdentityMatrix(F, dimV)));
		VP2 meet:= Kernel(inclVD*(MEta2 - IdentityMatrix(F, dimV)));
		
		dimVP0 := Dimension(VP0);
		dimVP2 := Dimension(VP2);
		// We expect the two spaces to be isomorphic, because the groups
		// \Gamma^P_0 and \Gamma^P_2 are conjugated to each other
		assert dimVP0 eq dimVP2;

		// Test if P0- and P2- space are nonzero, only then need following computations
		if dimVP0 gt 0 then
			VG := VP0 meet VP2;
			
			M132 := ZeroMatrix(F, dimVD, dimVD);
			M123 := ZeroMatrix(F, dimVD, dimVD);
			M13 := ZeroMatrix(F, dimVD, dimVD);

			MsumB := ZeroMatrix(F, dimVD, dimVD);
			MsumAB := ZeroMatrix(F, dimVD, dimVD);
			MsumBC := ZeroMatrix(F, dimVD, dimVD);
			MsumAC := ZeroMatrix(F, dimVD, dimVD);
			MsumABC := ZeroMatrix(F, dimVD, dimVD);

			// Dual iteration on s, allowing to get indices for Magma basis and our theoretical basis.
			for a -> tupA in sD do
				for b -> tupB in sD do
					// Careful: Magma works with row vectors.
					// la, mu: Index in domain -> a is index of row of resulting matrix
					// l, m: Index in codomain -> b is index of column of resulting matrix 
					la := tupA[1];
					mu := tupA[2];
					l := tupB[1];
					m := tupB[2];
					
					if l eq mu then
						if k-(l+m) eq la then
							M132[a, b] := 1;
						end if;
						
					end if;
					
					if m eq la then
						if k-(l+m) eq mu then
							M123[a, b] := 1;
						end if;
					end if;
					
					if m eq mu then
						if k-(l+m) eq la then
							M13[a, b] := (-1)^(1-type);
						end if;
						
						MsumB[a, b] := - Binomial(l, la);
					end if;	

					MsumAB[a, b] := myMultinomial(l, [la, mu - m, l - la - mu + m]);
					MsumBC[a, b] := Binomial(l, la)*Binomial(m, mu);
					MsumAC[a, b] := Binomial(l, la)*Binomial(m, mu+la-l);
					
					sum := 0;
					for g := 0 to Minimum([m, mu]) do
						i := mu - g; // => g+i eq mu
						if 0 le i and i le l and IsDivisibleBy(i, q-1) then
							sum := sum + myMultinomial(l, [la, i, l-la-i]) * Binomial(m, g);
						end if;                                          
					end for;
					MsumABC[a, b] := -sum;
					
				end for;    
			end for;

			Mstab132 := MId + MsumB + MsumC + MsumBC;
			Mstab123 := MId + MsumA + MsumB + MsumAB;
			Mstab13 := MId + MsumA + MsumB + MsumC + MsumAB + MsumBC + MsumAC + MsumABC;

			VG meet:= Nullspace(M132 - Mstab132);
			VG meet:= Nullspace(M123 - Mstab123);
			VG meet:= Nullspace(M13 + Mstab13);
						
			dimVG := Dimension(VG); 
		end if;	
	end if;
	
	return <k, dimV, dimVD, dimVP0, dimVG>;
end function; 

// Values of q we want to consider (as pairs <p, exponent>)
Q := [<2,1>, <3,1>, <2,2>, <5,1>, <7,1>, <2,3>, <3,2>, <11,1>, <13,1>, <2,4>, <17,1>, <19,1>, <23,1>, <5,2>, <3,3>, <29,1>, <31,1>, <2,5>];

// fileName := "dimension-list.txt";
fileName := "dimension-list-type1.txt";
fprintf fileName, "allDimensions := [";

for j -> tupQ in Q do
	p := tupQ[1];
	F := GF(p);
	q := p^tupQ[2];
	
	//types := [1 .. (q-1)]; // Allow all possible types
	types := [1]; // Only comput for type 1 so we can get further

	print "=================";
	q;
	types;
	
	results := [];
	for type in types do
		Append(~results, <q, type, [car< Integers(), Integers(), Integers(), Integers(), Integers()> | ]>);
	end for;
	
	max_time := 0;
	k := 0;
	
	while max_time le 900 do // Allow 15 minutes for the calculation at each weight k. reduce this when considering more than one type!
		for n := 1 to #types do
			type := types[n];
			if IsDivisibleBy(k+3-3*type, q-1) then 
				secV := Cputime();
				dimTuple := getDims(F, q, k, type);
				Append(~results[n][3], dimTuple);
				sec := Cputime(secV);
				max_time := Maximum(sec, max_time);
			end if;
		end for;
	k := k+1;
	end while;
	
	k;
	
	resString := Sprint(results);
	if j ne #Q then
		fprintf fileName, Substring(resString, 2, #resString-3) cat ", \n";
	else
		fprintf fileName, Substring(resString, 2, #resString-3) cat "\n";
	end if;
	
end for;

fprintf fileName, "\n];"; 
quit;
