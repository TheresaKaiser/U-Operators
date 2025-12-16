load "main.mg";

////////////////////////////////////////////////////////////////////////
// Here is the code that will actually be executed,                   //
// calling the other functions:                                       //
////////////////////////////////////////////////////////////////////////

// Values of q we want to apply the program to. Each tuple should be of the form
// <p, exponent for p>
Q := [<7,1>, <2,3>];



// Remember: Magma uses 1-based indexing
for n := 1 to #Q do
	// Initialize a timer
	secQ := Cputime();
    tupQ := Q[n];
    p := tupQ[1];
    m := tupQ[2];
    q := p^m;
    printf "\n ============= \n q: %o \n", q;
    type := 1;
	printf "\n Type : %o \n", type;
    
	// Will save results continously in the following files
    filenamesLaTeX := <"slopes1LaTeXq" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt", "slopes2LaTeXq" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt">;
    filenamesTuples := <"tuples1q" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt", "tuples2q" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt">;
	// Add Headers for file
    for i in [1,2] do
		fprintf filenamesLaTeX[i], "$k$ & $T_%o$-Slopes & $U_%o^{\G^P_0}$-Slopes & $U_%o^{\G^P_2}$-Slopes & $U_%o^{\G_0(t)}$-Slopes \\\\ \n \\hline \n", i, i, i, i;
    end for;
    
    maxTime := 0;
    k := 0;    
    while maxTime le 12000 do
		printf "\nk: %o - ", k;
        
        // Save computation time by ignoring cases where we already know that there are no Gamma_0-invariant harmonic cocycles
        if IsDivisibleBy(k+3-3*type, q-1) then 
			// Define structures we are working over
			F<g> := FiniteField(p);
			A<t> := PolynomialRing(F);

			try
				secK := Cputime();
				iTuple := getCharPols(F, q, k, type);
				maxTime := Maximum(maxTime, Cputime(secK));
				
				dimTuple := iTuple[3];
				
				for i in [1,2] do
					polTuple := iTuple[i];
					
					fprintf filenamesTuples[i],  "Append(~listTuples, <" cat Sprint(k) cat ", " cat Sprint(polTuple) cat ">); \n\n";
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
        
        k := k+1;
    end while;
    
    upper := k;

	for type in [2 .. (q-2)] cat [0] do
		printf "\n Type : %o \n", type;
		
		// Will save results continously in the following files
		filenamesLaTeX := <"slopes1LaTeXq" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt", "slopes2LaTeXq" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt">;
		filenamesTuples := <"tuples1q" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt", "tuples2q" cat Sprint(q) cat "type" cat Sprint(type) cat ".txt">;
		// Add Headers for file
		for i in [1,2] do
			fprintf filenamesLaTeX[i], "$k$ & $T_%o$-Slopes & $U_%o^{\G^P_0}$-Slopes & $U_%o^{\G^P_2}$-Slopes & $U_%o^{\G_0(t)}$-Slopes \\\\ \n \\hline \n", i, i, i, i;
		end for;
		
		for k in [0 .. upper] do
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
						
						fprintf filenamesTuples[i],  "Append(~listTuples, <" cat Sprint(k) cat ", " cat Sprint(polTuple) cat ">); \n\n";
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
end for;
// Close Magma process
quit;
