// Kronecker-Delta
delta := function(x, y)
    if x eq y then 
        return 1;
    else 
        return 0;
    end if;
end function;

// Multinomialcoefficients with "bad" entries
myMultinomial := function(n, Q)
    for r in Q do
        if r lt 0 or r gt n then
            return 0;
        end if;
    end for;
    return Multinomial(n, Q);
end function;



// getVectorSpaces: A function returning characteristic polynomials for the
// Hecke-action on the space of Gamma_0-invariant harmonic cocycles VD, 
// the spaces of \Gamma^P_0- and \Gamma^P_2-invariant harmonic cocycles VP0 and VP2,
// and the space of GL_3(A)-invariant harmonic cocycles VG. 
// Input:
    // F the current finite field F_p
    // q the current prime power q = p^...
    // k the current weight,
    // type the current type
// Output: 
	// A 3-tuple containing two 4-tuples with the characteristic polynomials 
	// of the U_1^\Gamma, resp. U_2^\Gamma operator acting on VD, VP0, VP2 and VG,
	// and a 4-tuple with the dimensions of V=V_{k, type}, VD, VP0 (isom. to VP2), VG 
	// If one of the spaces is zero, a zero is returned as the characteristic polynomial
getCharPols := function(F, q, k, type)  
	secV := Cputime();
	
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
	
	// Output the dimensions of the subspaces to the terminal,
	// so that we know this step is completed
	printf "Dim V = %o, Dim VD = %o, Dim VP0 = %o, Dim VP2 = %o, Dim VG = %o. \n", dimV, dimVD, dimVP0, dimVP2, dimVG;
	printf "Finding the vector spaces took %o seconds.\n", Cputime(secV);
	
	secABC := Cputime();

	A<t> := PolynomialRing(F); // Need to remind Magma about t

	// Initialize results //////////////////////////////////////////////
	Pol1G := 0;
	Pol1P0 := 0;
	Pol1P2 := 0;
	Pol1D := 0;
	Pol2G := 0;
	Pol2P0 := 0;
	Pol2P2 := 0;
	Pol2D := 0;
	
	if dimVD gt 0 then
		// The operators A, B, C ///////////////////////////////////////
		
		// A 
		A1 := ZeroMatrix(F, dimVD, dimVD);
		UA1 := ZeroMatrix(A, dimVD, dimVD);
		A2 := ZeroMatrix(F, dimVD, dimVD);
		UA2 := ZeroMatrix(A, dimVD, dimVD);
		
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

				// U1 
				// Initialize sums
				sum1 := 0;
				sum2 := 0;
				
				for g := 0 to Minimum([la, m]) do
					h1 := la - l - g;
					h2 := la - g;
					if h1 ge 0 then
						sum1 := sum1 + myMultinomial(k-l-m, [h1, mu, k-l-m - (h1+mu)]) * Binomial(m, g);
					end if;
					if h2 ge 0 then
						sum2 := sum2 + myMultinomial(k-l-m, [h2, mu, k-l-m - (h2+mu)]) * Binomial(m, g);
					end if;                                       
				end for;
				
				A1[a, b] := myMultinomial(k-l-m, [la - l, mu - m, k-l-m - ((la-l)+(mu-m))])
							- (-1)^l * myMultinomial(k-l-m, [la, mu - m, k-l-m - (la+(mu-m))])
							- (-1)^m * sum1
							+ (-1)^(l+m) * sum2
							+ delta(la, l) * (- Binomial(k-l-m, mu - m) + (-1)^m * Binomial(k-l-m, mu))
							+ delta(mu, m) * (- Binomial(k-l-m, la - l) + (-1)^l * Binomial(k-l-m, la) + delta(la, l));
				//UA1[a, b] := A1[a, b]*t^(la+mu); // already renormalized with T^(k+1-type) to get positive valuations
				
				// We changed the renormalization in the article since 
				// executing the calaculation. Now it should be as follows:
				UA1[a, b] := A1[a, b]*t^(la+mu - 2*((type-1) mod (q-1))); 

				//U2
				A2[a, b] := Binomial(m, m - mu) * Binomial(k-l-m, la + mu - (l+m))
						   - (-1)^l * Binomial(m, l + m - mu) * Binomial(k-l-m, la + mu - (l+m))
						   - (-1)^l * Binomial(m, mu) * Binomial(k-l-m, la)
						   + (-1)^(2*l + m - mu) * Binomial(m, mu - l) * Binomial(k-l-m, la)
						   + delta(la + mu, l + m) * (- Binomial(m, m - mu) + (-1)^l * Binomial(m, la))
						   + delta(mu, m) * (- Binomial(k-l-m, la - l) + (-1)^l * Binomial(k-l-m, la) + delta(la, l));
				//UA2[a, b] := A2[a, b]*t^(la);  // renormalized with T^(k+2*(1-type)) to get positive valuations
								
				// We changed the renormalization in the article since 
				// executing the calaculation. Now it should be as follows:
				UA2[a, b] := A2[a, b]*t^(la - ((type-1) mod (q-1)));
				
			end for;    
		end for;
		
		
		if dimVP0 gt 0 then
			B1 := ZeroMatrix(F, dimVD, dimVD);
			UAB1 := ZeroMatrix(A, dimVD, dimVD);
			B2 := ZeroMatrix(F, dimVD, dimVD);
			UAB2 := ZeroMatrix(A, dimVD, dimVD);
		
			// B
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

					// U1 
					if (k-l-m) eq mu then
						B1[a, b] := (-1)^(1-type) * (Binomial(m, la-l) - (-1)^l * Binomial(m, la));
					end if;
					B1[a, b] := B1[a, b] - M23[a, b]; // factor (-1)^(1-type) was already included in M23
					//UAB1[a, b] := UA1[a, b] + B1[a, b]*t^(la+mu);
					
					// We changed the renormalization in the article since 
					// executing the calaculation. Now it should be as follows:
					UAB1[a, b] := UA1[a, b] + B1[a, b]*t^(la+mu - 2*((type-1) mod (q-1))); 

					//U2
					if l eq mu then
						B2[a, b] := (-1)^(1-type) * (Binomial(k-l-m, la-m) - (-1)^m * Binomial(k-l-m, la));
					end if;
					B2[a, b] := B2[a, b] - M12[a, b]; // factor (-1)^(1-type) was already included in M12
					//UAB2[a, b] := UA2[a, b] + B2[a, b]*t^(la);

					// We changed the renormalization in the article since 
					// executing the calaculation. Now it should be as follows:
					UAB2[a, b] := UA2[a, b] + B2[a, b]*t^(la - ((type-1) mod (q-1)));
				end for;    
			end for;

			if dimVG gt 0 then
				C1 := ZeroMatrix(F, dimVD, dimVD);
				UABC1 := ZeroMatrix(A, dimVD, dimVD);
				C2 := ZeroMatrix(F, dimVD, dimVD);
				UABC2 := ZeroMatrix(A, dimVD, dimVD);
			
				// C
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
						// C1 = +M123
						//UABC1[a, b] := UAB1[a, b] + M123[a, b]*t^(la+mu);
						
						// We changed the renormalization in the article since 
						// executing the calaculation. Now it should be as follows:
						UABC1[a, b] := UAB1[a, b] + M123[a, b]*t^(la+mu - 2*((type-1) mod (q-1))); 

						// C2 = +M132
						//UABC2[a, b] := UAB2[a, b] + M132[a, b]*t^(la);
						
						// We changed the renormalization in the article since 
						// executing the calaculation. Now it should be as follows:
						UABC2[a, b] := UAB2[a, b] + M132[a, b]*t^(la - ((type-1) mod (q-1)));
					end for;    
				end for;
			end if;
		end if;
	end if;
	printf "Building the operator matrices took %o seconds.\n", Cputime(secABC);
	
	secR := Cputime();
	// The restriction to the correct subspaces ////////////////////////	
	if dimVP0 gt 0 then
		if dimVG gt 0 then
			BSeq := ExtendBasis(VG, VD);
			BM := Matrix(F, dimVD, dimVD, BSeq); 
			BMinv := BM^(-1);

			// Need to force coercion of matrix entries into A
			T1 := Matrix(A, RowSubmatrix(BM, 1, dimVG)) * UABC1 * Matrix(A, ColumnSubmatrix(BMinv, 1, dimVG)); 	
			T2 := Matrix(A, RowSubmatrix(BM, 1, dimVG)) * UABC2 * Matrix(A, ColumnSubmatrix(BMinv, 1, dimVG)); 		
		end if;
		
		BSeq0 := ExtendBasis(VP0, VD);
		BM0 := Matrix(F, dimVD, dimVD, BSeq0);
		BM0inv := BM0^(-1);
		BSeq2 := ExtendBasis(VP2, VD);
		BM2 := Matrix(F, dimVD, dimVD, BSeq2);
		BM2inv := BM2^(-1);
		
		U1P0 := Matrix(A, RowSubmatrix(BM0, 1, dimVP0)) * UA1 * Matrix(A, ColumnSubmatrix(BM0inv, 1, dimVP0)); 
		U1P2 := Matrix(A, RowSubmatrix(BM2, 1, dimVP2)) * UAB1 * Matrix(A, ColumnSubmatrix(BM2inv, 1, dimVP2)); 
		U2P0 := Matrix(A, RowSubmatrix(BM0, 1, dimVP0)) * UAB2 * Matrix(A, ColumnSubmatrix(BM0inv, 1, dimVP0)); 
		U2P2 := Matrix(A, RowSubmatrix(BM2, 1, dimVP2)) * UA2 * Matrix(A, ColumnSubmatrix(BM2inv, 1, dimVP2)); 
	end if;		
	printf "Restricting to subspaces took %o seconds.\n", Cputime(secR);
	
	secP := Cputime();
	// The characteristic polynomials //////////////////////////////////	
	if dimVD gt 0 then
		if dimVP0 gt 0 then
			if dimVG gt 0 then
				Pol1G := CharacteristicPolynomial(T1);
				Pol2G := CharacteristicPolynomial(T2);
			end if;
			Pol1P0 := CharacteristicPolynomial(U1P0);
			Pol1P2 := CharacteristicPolynomial(U1P2);
			Pol2P0 := CharacteristicPolynomial(U2P0);
			Pol2P2 := CharacteristicPolynomial(U2P2);
		end if;		
		Pol1D := CharacteristicPolynomial(UA1);
		Pol2D := CharacteristicPolynomial(UA2);
		
		// This changes the name of the variable from $.1 to X for all 
		// the characteristic polynomials (for printing)
		r<X> := Parent(Pol1D);
	end if;
	printf "Computing the charpols took %o seconds.\n", Cputime(secP);

	return <<Pol1D, Pol1P0, Pol1P2, Pol1G>, <Pol2D, Pol2P0, Pol2P2, Pol2G>, <dimV, dimVD, dimVP0, dimVG>>;
end function;


// getSlopes: A function calculating the slopes from the coefficients of 
// the given polynomial. The built-in Magma function for slopes does not
// return the slope infinity, so we calculate the slopes ourselves from
// the newton polygon that Magma generates.
// Input: 
	// f the polynomial
    // A the polynomial ring whose t-adic valuation we use
    // i either 1 or 2 depending on which U_i-operator the current polynomial belongs to
    // k current weight so that knowing i we can mark the expected newform slopes in blue
// Output: 
	// Slopes a string containing the slopes with multiplicities as bold exponents
	// formatted so that we copy the output directly into LaTeX
getSlopes := function(f, A, i, k)
	Slopes := "";
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
			Slopes := Slopes cat "\\textcolor{fullblue}{$" cat Sprint(slope) cat "^{\\mathbf{" cat Sprint(length) cat "}}$}, ";
		else
			Slopes := Slopes cat "$" cat Sprint(slope) cat "^{\\mathbf{" cat Sprint(length) cat "}}$, ";
		end if;
	end for;
	// The remaining segment belongs to the slope infinity
	if lv[1][1] ne 0 then
		Slopes := Slopes cat "$\\infty^{\\mathbf{" cat Sprint(lv[1][1]) cat "}}$";
	else // if there is no infinite slope, delete the last comma and space 
		Slopes := Substring(Slopes, 1, #Slopes - 2);
	end if;
	return Slopes;
end function;


////////////////////////////////////////////////////////////////////////
// Here is the code that will actually be executed,                   //
// calling the other functions:                                       //
////////////////////////////////////////////////////////////////////////

// Values of q we want to apply the program to. Each tuple should be of the form
// <p, exponent for p, type, lower limit for weight k, upper limit for k>
Q := [<2,1, 0, 0,50>, <3,1, 0, 0,80>, <3,1, 1, 0,80>, <2,2, 0, 0,100>, <2,2, 1, 0,100>, <2,2, 2, 0,100>, <5,1, 0, 0, 130>, <5,1, 1, 0, 130>, <5,1, 2, 0, 130>, <5,1, 3, 0,130>];


// Remember: Magma uses 1-based indexing
for n := 1 to #Q do
	// Initialize a timer
	secQ := Cputime();
    tupQ := Q[n];
    p := tupQ[1];
    m := tupQ[2];
    q := p^m;
    printf "\n ============= \n q: %o \n", q;
    type := tupQ[3];
    printf "\n Type : %o \n", type;
    
    lower := tupQ[4];
    upper := tupQ[5];
    
    // Will save results continously in the following files
    filenamesLaTeX := <"slopes1LaTeXq" cat Sprint(q) cat "k" cat Sprint(lower) cat "-" cat Sprint(upper) cat "type" cat Sprint(type) cat ".txt", "slopes2LaTeXq" cat Sprint(q) cat "k" cat Sprint(lower) cat "-" cat Sprint(upper) cat "type" cat Sprint(type) cat ".txt">;
    filenamesTuples := <"tuples1q" cat Sprint(q) cat "k" cat Sprint(lower) cat "-" cat Sprint(upper) cat "type" cat Sprint(type) cat ".txt", "tuples2q" cat Sprint(q) cat "k" cat Sprint(lower) cat "-" cat Sprint(upper) cat "type" cat Sprint(type) cat ".txt">;
    
    // Add Headers for files
    for i in [1,2] do
		fprintf filenamesLaTeX[i], "$k$ & $T_%o$-Slopes & $U_%o^{\G^P_0}$-Slopes & $U_%o^{\G^P_2}$-Slopes & $U_%o^{\G_0(t)}$-Slopes \\\\ \n \\hline \n", i, i, i, i;
    end for;

	 for k in [lower .. upper] do
        printf "\nk: %o - ", k;
        
        // Save computation time by ignoring cases where we already know that there are no Gamma_0-invariant harmonic cocycles
        if IsDivisibleBy(k+3-3*type, q-1) then 
			// Define structures we are working over
			F<g> := FiniteField(p);
			A<t> := PolynomialRing(F);

			try
				iTuple := getCharPols(F, q, k, type);
				
				dimTuple := iTuple[3];
				
				for i in [1,2] do
					polTuple := iTuple[i];
					
					fprintf filenamesTuples[i],  "Append(~listTuples, <" cat Sprint(k) cat ", " cat Sprint(polTuple) cat ">); \n \n ";
					stringLaTeX := Sprint(k) cat " & ";

					secA := Cputime();
					if Type(polTuple[4]) eq RngIntElt then
						stringLaTeX := stringLaTeX cat " & ";
					else
						slopeString := getSlopes(polTuple[4], A, i, k);
						stringLaTeX := stringLaTeX cat slopeString cat " & ";
					end if;
					if Type(polTuple[2]) eq RngIntElt then
						stringLaTeX := stringLaTeX cat " & ";
					else
						slopeString := getSlopes(polTuple[2], A, i, k);
						stringLaTeX := stringLaTeX cat slopeString cat " & ";
					end if;
					if Type(polTuple[3]) eq RngIntElt then
						stringLaTeX := stringLaTeX cat " & ";
					else
						slopeString := getSlopes(polTuple[3], A, i, k);
						stringLaTeX := stringLaTeX cat slopeString cat " & ";
					end if;
					if Type(polTuple[1]) eq RngIntElt then
						stringLaTeX := stringLaTeX cat "  \\\\ \n";
					else
						slopeString := getSlopes(polTuple[1], A, i, k);
						stringLaTeX := stringLaTeX cat slopeString cat " \\\\ \n ";
					end if;
					printf "Analyzing the charPols for i=%o took %o seconds.\n", i, Cputime(secA);
					
					fprintf filenamesLaTeX[i], stringLaTeX cat "\\hline \n";
				end for;
			catch e
				printf "\n An Error occured for q = %o, type = %o and k = %o \n", q, type, k;
				print "Error: ", e`Object;
				// Close Magma process
				quit;
			end try;
        else
            print "The space of Gamma_0 invariant cocycles is zero.";
        end if;
    end for; 
end for;
// Close Magma process
quit;
