// The following three functions assume that the type 0 has been replaced by the type q-1.
checkDimVG := function(k, q, type)
	A := 1;
	B := q+1;
	C := q^2+q+1;

	S := {A, B, C};

    b, n := IsDivisibleBy(k+3-type*C, q-1);
    if b and n ge 0 then
        return #RestrictedPartitions(n, S);
    else
        return 0;
    end if;
end function;

checkDimVP := function(k, q, type)
    A := 1;
    B := q+1;
    
    K := k-(type-1)*(q+2);
    b, n := IsDivisibleBy(K, q-1);
    if b and n ge 0 then
		sum := 0;
		for i := 0 to n-1 do
            assert #RestrictedPartitions(i, {A, B}) eq (Floor(i/(q+1))+1);
			sum := sum + Floor(i/(q+1))+1;
		end for;
        return sum;
    else
        return 0;
    end if;
end function;

checkDimVD := function(k, q, type)
    b, n := IsDivisibleBy(k+3*(1-type), q-1);
    if b and n ge 0 then
        return Binomial(n+2,2);
    else
        return 0;
    end if;
end function;

analyzeTuple := function(current)
    q := current[1];
    type := current[2];
        if type eq 0 then
        type := q-1;
    end if;
    dims := current[3];

    kmax := dims[#dims][1];

    print "==========================================================================";
    printf "Analyzing q=%o, type=%o. \nThe weight k goes up to %o.\n \n", q, type, kmax;


    for counter := 1 to #dims do
        tup := dims[counter];
        
        k := tup[1];
        dimV := tup[2];
        dimVD := tup[3];
        dimVP := tup[4];
        dimVG := tup[5];

        assert dimV eq Binomial(k+2, 2);
        assert dimVD eq checkDimVD(k, q, type);
        assert dimVP eq checkDimVP(k, q, type);
        assert dimVG eq checkDimVG(k, q, type);
    end for;
    return "All tests successfull.";
end function;


// Values of q we want to consider (as pairs <p, exponent>)
Q := [<2,1>, <3,1>, <2,2>, <5,1>, <7,1>, <2,3>, <3,2>, <11,1>, <13,1>, <2,4>, <17,1>, <19,1>, <23,1>, <5,2>, <3,3>, <29,1>, <31,1>, <2,5>];

load "dimension-list.txt";
lower := 0;
upper := 0;
for j -> tupQ in Q do
	p := tupQ[1];
	q := p^tupQ[2];
    
    types := [1 .. (q-1)]; // Allow all possible types
    
    lower := upper + 1;
    upper := upper + #types;
    
    for n := lower to upper do
        analyzeTuple(allDimensions[n]);
    end for;
    
end for;


load "dimension-list-type1.txt"; 
for j -> tupQ in Q do
	p := tupQ[1];
	q := p^tupQ[2];
    analyzeTuple(allDimensions[j]);
end for;

// End Magma process
quit;