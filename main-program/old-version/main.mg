// mToM: This stands for "matrix to matrix". The function returns a 
// transformation matrix MV describing the action of a matrix M on V.
// In our case, the entries will be in A=F_q[t] or K=F_q(t).
// Input: 
    // R the field or ring we are working over, 
    // s list with tuples of indices, allowing to switch between Magma basis for V and our basis with double indices
    // k the current weight
    // M the matrix acting on V
// Output: 
	// MV a #s x #s matrix with entries in R
mToM := function(R, s, k, M) 
    
    // Dimension of resulting matrix
    d := #s;
    
    // Initialize resulting matrix
    MV := ZeroMatrix(R, d, d);

    // Dual iteration on s, allowing to get indices for Magma basis and our theoretical basis.
    for a1 -> r in s do
        for a2 -> t in s do
            // Careful: Magma works with row vectors.
            // l, m: Index in codomain -> index of column of resulting matrix 
            // la, mu: Index in domain -> index of row of resulting matrix
            l := r[1];
            m := r[2];
            la := t[1];
            mu := t[2];
            
            // 6-fold-sum
            sum := 0;
            for b := 0 to Minimum(la, l) do
                for c := 0 to Minimum(mu, l-b) do
                    for e := 0 to Minimum(la, m) do
                        for f := 0 to Minimum(mu, m-e) do
                            for h := 0 to Minimum(la, k-l-m) do
                                for i := 0 to Minimum(mu, k-l-m-h) do
                                    if b+e+h eq la and c+f+i eq mu then
                                        d := l-b-c;
                                        g := m-e-f;
                                        j := k-l-m-h-i;
                                        
                                        koeff1 := Multinomial(l, [b,c,d]);
                                        koeff2 := Multinomial(m, [e,f,g]);
                                        koeff3 := Multinomial(k-l-m, [h,i,j]);
                                        
                                        sum := sum + koeff1*koeff2*koeff3 * M[1,1]^b * M[1,2]^c * M[1,3]^d * M[2,1]^e * M[2,2]^f * M[2,3]^g * M[3,1]^h * M[3,2]^i * M[3,3]^j;
                                    end if;                        
                                end for;                    
                            end for;                        
                        end for;                    
                    end for;                       
                end for;                    
            end for;
            
            // Here, another exponent for the determinant could be used, in order to apply the programm to a different type
            MV[a2, a1] := Determinant(M)*sum;
        end for;    
    end for;
    
    return MV;
end function;


// getVectorSpaces: A function returning the space of Gamma_0-invariant harmonic cocycles VD, 
// the spaces of \Gamma^P_0- and \Gamma^P_2-invariant harmonic cocycles VP0 and VP2,
// and the space of GL_3(A)-invariant harmonic cocycles VG. 
// Input:
    // K the field or ring we are working over, 
    // F the underlying finite field
    // k the current weight,
    // s list with tuples of indices, allowing to switch between Magma basis for V and our basis with double indices
    // sD list with tuples of indices describing basis of V^{D(F_q)}
// Output: 
	// a 5-tuple containing the vector spaces V, VD, VP0, VP2, VG
getVectorSpaces := function(K, F, k, s, sD)  
    // Dimension and initialization of V
    d := #s;
    V := VectorSpace(K, d);
    
    // Initialize other vector spaces to make sure all return values are defined
    VD := VectorSpace(K, 0);
    dimD := #sD;
    VP0 := VectorSpace(K, 0);
    dimP0 := 0;
	VP2 := VectorSpace(K, 0);
	dimP2 := 0;
	VG := VectorSpace(K, 0);
    dimG := 0;

    if dimD gt 0 then
		// List of base elements of V describing basis of V^{D(F_q)}
		sDbase := [];
		I := IdentityMatrix(K, d);
		q := #F;
		for i := 1 to d do
			l := s[i][1];
			m := s[i][2];
			if IsDivisibleBy(l+1, q-1) and IsDivisibleBy(m+1, q-1) then
				Append(~sDbase, I[i]);
			end if;
		end for; 
		
		// Define subspace V^{D(F_q)} of V
		VD := sub<V|sDbase>;
		
		// Initialize subspaces for VP0 and VP2
		VP0 := VD;
		VP2 := VD;
		
		// Compute subspaces V^{H_i} of V that contain VP0, VP2
		MEta0 := mToM(K, s, k, Matrix(3, [1,0,0, 1,1,0, 0,0,1]));
		MEta2 := mToM(K, s, k, Matrix(3, [1,0,0, 0,1,0, 0,1,1]));
		VP0 meet:= Eigenspace(MEta0, 1);
		VP2 meet:= Eigenspace(MEta2, 1);
		
		// Now apply condition for one permutation matrix each
		M0 := mToM(K,s,k, Matrix(3, [0,1,0, 1,0,0, 0,0,1]));
		for a in F do
			M0 := M0 + mToM(K,s,k, Matrix(3, [1,a,0, 0,1,0, 0,0,1]));
		end for;
		VP0 meet:= Nullspace(M0);  
		dimP0 := Dimension(VP0);
		
		M2 := mToM(K,s,k, Matrix(3, [1,0,0, 0,0,1, 0,1,0]));
		for a in F do
			M2 := M2 + mToM(K,s,k, Matrix(3, [1,0,0, 0,1,a, 0,0,1]));
		end for;
		VP2 meet:= Nullspace(M2);  
		dimP2 := Dimension(VP2);
		
		// We expect the two spaces to be isomorphic, because the groups
		// \Gamma^P_0 and \Gamma^P_2 are conjugated to each other
		assert dimP0 eq dimP2;
		
		if dimP0 gt 0 then
			VG := VP0 meet VP2;
			
			// Apply conditions for remaining three non trivial permutation matrices			
			M3 := -mToM(K,s,k, Matrix(3, [0,1,0, 0,0,1, 1,0,0]));
			for a in F do
				for b in F do
					M3 := M3 + mToM(K,s,k, Matrix(3, [1,0,b, 0,1,a, 0,0,1]));
				end for;
			end for;
			VG meet:= Nullspace(M3);
			
			M4 := -mToM(K,s,k, Matrix(3, [0,0,1, 1,0,0, 0,1,0]));
			for a in F do
				for b in F do
					M4 := M4 + mToM(K,s,k, Matrix(3, [1,a,b, 0,1,0, 0,0,1]));
				end for;
			end for;
			VG meet:= Nullspace(M4);
			
			M5 := -mToM(K,s,k, Matrix(3, [0,0,1, 0,1,0, 1,0,0]));
			for a in F do
				for b in F do
					for c in F do
						M5 := M5 - mToM(K,s,k, Matrix(3, [1,a,b, 0,1,c, 0,0,1]));
					end for;
				end for;
			end for;
			VG meet:= Nullspace(M5);
						
			dimG := Dimension(VG);
		end if;
	end if;
    
    // Output the dimensions of the subspaces to the terminal,
    // so that we know this step is completed
    printf "Dim V = %o, Dim VD = %o, Dim VP0 = %o, Dim VP2 = %o, Dim VG = %o. \n", d, dimD, dimP0, dimP2, dimG;
    return <V, VD, VP0, VP2, VG>;
end function;

// getCharPols: A function returning the characteristic polynomial of 
// the U_i, resp. T_i operators acting on VD, VP0, VP2, VG. 
// Since we know that the operators stay in VD, it should suffice to use  
// the entries of sD for the mToM function, but in that case, we have 
// difficulties with ExtendBasis.
// So, the code expects you to use s as the list when calling this.
// Input: 
	// i must be 1 or 2, indexing the U_i and T_i operators
	// K the field we are working over
	// F the current finite field
	// k the current weight
	// s the list of index-tuples for the basis of the vector space
	// VTuple 5-tuple containing the vector spaces V, VD, VP0, VP2, VG
// Output: 
	// a 4-tuple containing the characteristic polynomials of U_i resp. 
	//     T_i acting on VD, VP0, VP2, VG, or 0 if the corresponding 
	//     space is zero
getCharPols := function(i, K, F, k, s, VTuple)
	K<t> := K; // Need to remind Magma about t

	// Get vector spaces and dimensions ////////////////////////////////
	V := VTuple[1];
	d := #s;
	VD := VTuple[2];
	dimD := Dimension(VD);
	VP0 := VTuple[3];
	dimP0 := Dimension(VP0);
	VP2 := VTuple[4];
	dimP2 := Dimension(VP2);
	VG := VTuple[5];
	dimG := Dimension(VG);
	
	// Initialize results //////////////////////////////////////////////
	PolG := 0;
	PolP0 := 0;
	PolP2 := 0;
	PolD := 0;
	
	if dimD gt 0 then
		// The action of delta_i^{-1} //////////////////////////////////
		if i eq 1 then
			dInv := Matrix(3, [1,0,0, 0,1,0, 0,0,1/t]);
		else
			dInv := Matrix(3, [1,0,0, 0,1/t,0, 0,0,1/t]);
		end if;
		MdInv := mToM(K,s,k, dInv);
		
		// The operators A, B, C ///////////////////////////////////////
		
		// A 
		ListAab := [];
		ListAa := [];
		if i eq 1 then // i=1
			for a in F do
				if a ne 0 then			
					for b in F do
						if b ne 0 then
							Append(~ListAab, mToM(K,s,k, Matrix(3, [1,0,0, 0,1,0, a,b,1])) - mToM(K,s,k, Matrix(3, [0,0,-a, 0,1,0, 1/a,b/a,1])) 
										- mToM(K,s,k, Matrix(3, [1,0,0, -a*b,0,-a, b,1/a,1])) + mToM(K,s,k, Matrix(3, [0,0,-a, -b/a,0,-b, 1/a,1/b,1])));
						end if;
					end for;
					
					Append(~ListAa, mToM(K,s,k, Matrix(3, [1,0,0, 0,1,0, a,0,1])) + mToM(K,s,k, Matrix(3, [1,0,0, 0,1,0, 0,a,1])) 
							- mToM(K,s,k, Matrix(3, [0,0,-a, 0,1,0, 1/a,0,1])) - mToM(K,s,k, Matrix(3, [1,0,0, 0,0,-a, 0,1/a,1])));
				end if;
			end for;
		else // i = 2
			for a in F do
				if a ne 0 then			
					for b in F do
						if b ne 0 then
							Append(~ListAab, mToM(K,s,k, Matrix(3, [1,0,0, a,1,0, b,0,1])) - mToM(K,s,k, Matrix(3, [0,-a,0, 1/a,1,0, b,0,1])) 
										- mToM(K,s,k, Matrix(3, [0,0,-a, 0,1,b, 1/a,0,1])) + mToM(K,s,k, Matrix(3, [0,-a,0, 0,1,-b/a, 1/b,0,1])));
						end if;
					end for;
					
					Append(~ListAa, mToM(K,s,k, Matrix(3, [1,0,0, a,1,0, 0,0,1])) + mToM(K,s,k, Matrix(3, [1,0,0, 0,1,0, a,0,1])) 
							- mToM(K,s,k, Matrix(3, [0,-a,0, 1/a,1,0, 0,0,1])) - mToM(K,s,k, Matrix(3, [0,0,-a, 0,1,0, 1/a,0,1])));
				end if;
			end for;
		end if;
		
		A := IdentityMatrix(K, d);
		for m in ListAab do
			A := A + m;
		end for;
		for m in ListAa do
			A := A + m;
		end for;
		
		// If q > p: Check if entries are in F_p
		// Only print something if this fails
		if Degree(F) gt 1 then
			E := BaseField(F); 
			R := MatrixRing(E, d);
			 try
				MA := R ! A;
			catch e
				print "The entries of A are not in F_p";
			end try;
		end if;
		
		// Careful: We want to first apply d_i^{-1}, then A, B, C. 
		// Magma works with row vectors!!!
		U := MdInv * A * t^(k+i); // need to renormalize
		
		
		if dimP0 gt 0 then
		
			// B
			ListBa := [];
			if i eq 1 then // i = 1
				B := mToM(K,s,k, - Matrix(3, [1,0,0, 0,0,1, 0,1,0]));
				for a in F do
					if a ne 0 then			
						Append(~ListBa, mToM(K,s,k, Matrix(3, [0,0,-a, 1/a,0,1, 0,1,0])) -  mToM(K,s,k, Matrix(3, [1,0,0, a,0,1, 0,1,0])));
					end if;
				end for;
			else // i = 2
				B := mToM(K,s,k, - Matrix(3, [0,1,0, 1,0,0, 0,0,1]));
				for a in F do
					if a ne 0 then			
						Append(~ListBa, mToM(K,s,k, Matrix(3, [0,1,0, 0,0,-a, 1/a,0,1])) - mToM(K,s,k, Matrix(3, [0,1,0, 1,0,0, a,0,1])));
					end if;
				end for;
			end if;
			
			for m in ListBa do
				B := B + m;
			end for;
			// If q > p: Check if entries are in F_p
			// Only print something if this fails
			if Degree(F) gt 1 then
				E := BaseField(F); 
				R := MatrixRing(E, d);
				 try
					MB := R ! B;
				catch e
					print "The entries of B are not in F_p";
				end try;
			end if;
			
			UAB := MdInv * (A + B) * t^(k+i); // need to renormalize
			
			if dimG gt 0 then
				// C
				if i eq 1 then // i = 1
					C := - mToM(K,s,k, Matrix(3, [0,0,1, 1,0,0, 0,1,0]));
				else  // i = 2
					C := - mToM(K,s,k, Matrix(3, [0,1,0, 0,0,1, 1,0,0]));
				end if;
				// If q > p: Check if entries are in F_p
				// Only print something if this fails
				if Degree(F) gt 1 then
					E := BaseField(F); 
					R := MatrixRing(E, d);
					 try
						MC := R ! C;
					catch e
						print "The entries of C are not in F_p";
					end try;
				end if;
				T := MdInv * (A + B + C) * t^(k+i); // need to renormalize

	// The restriction to the correct subspaces ////////////////////////			
				BSeq := ExtendBasis(VG, V);
				BM := Matrix(K, d, d, BSeq);
				BMinv := BM^(-1);
				TG := RowSubmatrix(BM, 1, dimG) * T * ColumnSubmatrix(BMinv, 1, dimG); 	
				PolG := CharacteristicPolynomial(TG);
			end if;
			
			BSeq0 := ExtendBasis(VP0, V);
			BM0 := Matrix(K, d, d, BSeq0);
			BM0inv := BM0^(-1);

			BSeq2 := ExtendBasis(VP2, V);
			BM2 := Matrix(K, d, d, BSeq2);
			BM2inv := BM2^(-1);
			
			if i eq 1 then
				UP0 := RowSubmatrix(BM0, 1, dimP0) * U * ColumnSubmatrix(BM0inv, 1, dimP0); 
				UP2 := RowSubmatrix(BM2, 1, dimP2) * UAB * ColumnSubmatrix(BM2inv, 1, dimP2); 
			else
				UP0 := RowSubmatrix(BM0, 1, dimP0) * UAB * ColumnSubmatrix(BM0inv, 1, dimP0); 
				UP2 := RowSubmatrix(BM2, 1, dimP2) * U * ColumnSubmatrix(BM2inv, 1, dimP2); 
			end if;
			PolP0 := CharacteristicPolynomial(UP0);
			PolP2 := CharacteristicPolynomial(UP2);
		end if;		
		
		BSeq := ExtendBasis(VD, V);
		BM := Matrix(K, d, d, BSeq);
		BMinv := BM^(-1);
		
		UD := RowSubmatrix(BM, 1, dimD) * U * ColumnSubmatrix(BMinv, 1, dimD);
		
		PolD := CharacteristicPolynomial(UD);
		
		// This changes the name of the variable from $.1 to X for all 
		// the characteristic polynomials (at least for printing)
		r<X> := Parent(PolD);
	end if;

	return <PolD, PolP0, PolP2, PolG>;
end function;


// getSlopes: A function calculating the slopes from the coefficients of 
// the given polynomial. The built-in Magma function for slopes does not
// return the slope infinity, so we calculate the slopes ourselves from
// the newton polygon that Magma generates.
// (Also because we want to copy the output into LaTeX at the end)
// Input: 
	// f the polynomial
    // A the polynomial ring whose t-adic valuation we use
// Output: 
	// Slopes a string containing the slopes with multiplicities as bold exponents
getSlopes := function(f, A)
	Slopes := "";
	V := [];
	for i := 0 to Degree(f) do
		c := Coefficient(f, i);
		if c ne 0 then
			Append(~V, <i, Valuation(A ! c)>);
		end if;
	end for;    
	// Magma gives us a list of the lower vertices of the Newton polygon
	p := NewtonPolygon(V : Faces := "Lower");
	lv := LowerVertices(p);
	// Compute the slope of each segment
	// Need to change sign due to the order we listed our points in
	for i := 0 to #lv-2 do
		current := lv[#lv - i];
		next := lv[#lv - i - 1];
		length := current[1] - next[1];
		slope := (next[2] - current[2])/length;
		Slopes := Slopes cat "$" cat Sprint(slope) cat "^{\\mathbf{" cat Sprint(length) cat "}}$, ";
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
// calling all the other functions:                                   //
////////////////////////////////////////////////////////////////////////

// Initialize a timer
sec := Cputime();

// Values of q we want to apply the program to
Q := [2,3];
// For each q in Q: Lower and upper bounds for values of k we want to apply the program to
// need to have same length as Q
lower := [0,0];
upper := [16,16];

// Remember: Magma uses 1-based indexing
for n := 1 to #Q do
	// Initialize a timer
	secQ := Cputime();
    q := Q[n];
    printf "\n ============= \n q: %o \n", q;

    //Initialize lists to save results
    weights := ["$k$;\n"];
    listDimD := ["$\dim \Char(\G_0(t), V)$;\n"];
    listDimP := ["$\dim \Char(\G^P_j, V)$;\n"];
    listDimG :=["$\dim \Char(GL_3(A), V)$;\n"];
    slopesU := <["$U_1^{\G_0(t)}$-Slopes;\n"], ["$U_2^{\G_0(t)}$-Slopes;\n"]>;
    slopesUP0 := <["$U_1^{\G^P_0}$-Slopes;\n"], ["$U_2^{\G^P_0}$-Slopes;\n"]>;
    slopesUP2 := <["$U_1^{\G^P_2}$-Slopes;\n"], ["$U_2^{\G^P_2}$-Slopes;\n"]>;
    slopesT := <["$T_1$-Slopes;\n"], ["$T_2$-Slopes;\n"]>;
    factorsU := <["$U_1^{\G_0(t)}$-Factors;\n"], ["$U_2^{\G_0(t)}$-Factors;\n"]>;
    factorsUP0 := <["$U_1^{\G^P_0}$-Factors;\n"], ["$U_2^{\G^P_0}$-Factors;\n"]>;
    factorsUP2 := <["$U_1^{\G^P_2}$-Factors;\n"], ["$U_2^{\G^P_2}$-Factors;\n"]>;
    factorsT := <["$T_1$-Factors;\n"], ["$T_2$-Factors;\n"]>;
    
    for k in [lower[n] .. upper[n]] do
        printf " \n \n k: %o - ", k;
        
        // Save computation time by ignoring cases where we already know that there are no Gamma_0-invariant harmonic cocycles
        if IsDivisibleBy(k+3, q-1) then 
			// Define structures we are working over
			F<g> := FiniteField(q);
			A<t> := PolynomialRing(F);
			K<t> := FieldOfFractions(A);
			
			// Initialize cartesian product. Necessary to define tuples of indices lying in s
			C := car< Integers(), Integers() >;
			// List with tuples of indices describing basis of V
			s := [ C | <l, m> : l in [0..k], m in [0..k] | l+m le k ];
			// List with tuples of indices describing basis of V^{D(F_q)}
			sD := [ C | <l, m> : l in [0..k], m in [0..k] | l+m le k and IsDivisibleBy(l+1, q-1) and IsDivisibleBy(m+1, q-1) and IsDivisibleBy(k-l-m+1, q-1)];
			
            VTuple := getVectorSpaces(K, F, k, s, sD);
            printf "Finding the vector spaces took %o seconds.\n", Cputime(secQ);
            
            if Dimension(VTuple[2]) gt 0 then
				Append(~weights, Sprint(k) cat ";\n");
				Append(~listDimD, Sprint(Dimension(VTuple[2])) cat ";\n");
				Append(~listDimP, Sprint(Dimension(VTuple[3])) cat ";\n");
				Append(~listDimG, Sprint(Dimension(VTuple[5])) cat ";\n");
				for i in [1,2] do
					secP := Cputime();
					PolTuple := getCharPols(i, K, F, k, s, VTuple);
					printf "Finding the charPols for i=%o took %o seconds.\n", i, Cputime(secP);
					secA := Cputime();
					if Type(PolTuple[1]) eq RngIntElt then
						Append(~slopesU[i], " ;\n");
						Append(~factorsU[i], " ;\n");
					else
						Append(~slopesU[i], getSlopes(PolTuple[1], A) cat ";\n");
						Append(~factorsU[i], Sprint(Factorisation(PolTuple[1])) cat ";\n");
					end if;
					if Type(PolTuple[2]) eq RngIntElt then
						Append(~slopesUP0[i], " ;\n");
						Append(~factorsUP0[i], " ;\n");
					else
						Append(~slopesUP0[i], getSlopes(PolTuple[2], A) cat ";\n");
						Append(~factorsUP0[i], Sprint(Factorisation(PolTuple[2])) cat ";\n");
					end if;
					if Type(PolTuple[3]) eq RngIntElt then
						Append(~slopesUP2[i], " ;\n");
						Append(~factorsUP2[i], " ;\n");
					else
						Append(~slopesUP2[i], getSlopes(PolTuple[3], A) cat ";\n");
						Append(~factorsUP2[i], Sprint(Factorisation(PolTuple[3])) cat ";\n");
					end if;
					if Type(PolTuple[4]) eq RngIntElt then
						Append(~slopesT[i], " ;\n");
						Append(~factorsT[i], " ;\n");
					else
						Append(~slopesT[i], getSlopes(PolTuple[4], A) cat ";\n");
						Append(~factorsT[i], Sprint(Factorisation(PolTuple[4])) cat ";\n");
					end if;
					printf "Analyzing the charPols for i=%o took %o seconds.\n", i, Cputime(secA);
				end for;
			end if;
        else
            print "The space of Gamma_0 invariant cocycles is zero.";
        end if;

    end for; 
    // open filestream to save result
    filenameSlopes := "slopesq" cat Sprint(q) cat "k" cat Sprint(lower[n]) cat "-" cat Sprint(upper[n]) cat ".txt";
    filenameFactors := "factorsq" cat Sprint(q) cat "k" cat Sprint(lower[n]) cat "-" cat Sprint(upper[n]) cat ".txt";
	fileSlopes := Open(filenameSlopes, "w");
	fileFactors := Open(filenameFactors, "w");
	// Print results to files. Use "!" as limiter between lists 
	fprintf fileSlopes, "!%o!", &cat weights;
	fprintf fileFactors, "!%o!", &cat weights;
	fprintf fileSlopes, "!%o!", &cat listDimD;
    fprintf fileSlopes, "!%o!", &cat listDimP;
    fprintf fileSlopes, "!%o!", &cat listDimG;
    fprintf fileFactors, "!%o!", &cat listDimD;
    fprintf fileFactors, "!%o!", &cat listDimP;
    fprintf fileFactors, "!%o!", &cat listDimG;
	for i in [1, 2] do
		fprintf fileSlopes, "!%o!", &cat slopesU[i];
		fprintf fileSlopes, "!%o!", &cat slopesUP0[i];
		fprintf fileSlopes, "!%o!", &cat slopesUP2[i];
		fprintf fileSlopes, "!%o!", &cat slopesT[i];

		fprintf fileFactors, "!%o!", &cat factorsU[i];
		fprintf fileFactors, "!%o!", &cat factorsUP0[i];
		fprintf fileFactors, "!%o!", &cat factorsUP2[i];
		fprintf fileFactors, "!%o!", &cat factorsT[i];
	end for;
	fprintf fileSlopes, "\n \n The computation took %o seconds.", Cputime(secQ);
	fprintf fileFactors, "\n \n The computation took %o seconds.", Cputime(secQ);
    // Close filestreams (this does not delete the files)
    delete fileSlopes;
    delete fileFactors;
end for;

printf "\n \n The program took %o seconds to run.", Cputime(sec);

// Close Magma process
quit;


