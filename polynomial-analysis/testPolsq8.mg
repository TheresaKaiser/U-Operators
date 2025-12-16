load "header.mg";

q := 8; 
q;
F := GF(q); A<t> := PolynomialRing(F); K<t> := FieldOfFractions(A); KK := AlgebraicClosure(K); r<X> := PolynomialRing(KK); 
C := car<Integers(), car<r, r, r, r> >;

LaTeX := false;
filenamesLaTeX := <"slopes1LaTeXq8k0-176.txt", "slopes2LaTeXq8k0-176.txt">;
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
load "tuples1q8type0.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 1;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type1.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;


type := 2;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type2.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 3;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type3.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 4;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type4.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 5;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type5.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 6;
print "type: ", type;
listTuples := [C|];  
load "tuples1q8type6.txt";
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
load "tuples2q8type0.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 1;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type1.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 2;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type2.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 3;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type3.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 4;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type4.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 5;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type5.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;

type := 6;
print "type: ", type;
listTuples := [C|]; 
load "tuples2q8type6.txt";
for pair in listTuples do
	testPair(pair, q, i, type, A);
	if LaTeX then
		temp := printLine(pair, A, i, q, type, filenamesLaTeX[i]);
	end if;
end for;


quit;
