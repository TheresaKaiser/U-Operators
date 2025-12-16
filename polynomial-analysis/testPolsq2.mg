load "header.mg";

q := 2; 
q;
F := GF(q); A<t> := PolynomialRing(F); K<t> := FieldOfFractions(A); KK := AlgebraicClosure(K); r<X> := PolynomialRing(KK); 
C := car<Integers(), car<r, r, r, r> >;

LaTeX := false;
filenamesLaTeX := <"slopes1LaTeXq2k0-50.txt", "slopes2LaTeXq2k0-50.txt">;
if LaTeX then
	// Add Headers for files
	for i in [1,2] do
		fprintf filenamesLaTeX[i], "$k$ & $T_%o$-Slopes & $U_%o^{P_0}$-Slopes & $U_%o^{P_2}$-Slopes & $U_%o^{\G_0(t)}$-Slopes \\\\ \n \\hline \n", i, i, i, i;
	end for;
end if;

////////////////////////////////////////////////////////////////////////
i := 1;
print "i: ", i;
////////////////////////////////////////////////////////////////////////
type := 0;
print "type: ", type;
listTuples := [C|];  
load "pols-q2-i1.txt";
//load "pols-test-copy.txt";


for pair in listTuples do
	testPair(pair, q, i, type, A);
	
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;


////////////////////////////////////////////////////////////////////////
i := 2;
print "i: ", i;
////////////////////////////////////////////////////////////////////////
type := 0;
print "type: ", type;
listTuples := [C|]; 
load "pols-q2-i2.txt";

for pair in listTuples do
	testPair(pair, q, i, type, A);
	
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

quit;
