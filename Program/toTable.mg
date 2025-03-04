filename := "factorsq2k0-16";
input := Read(filename cat ".txt");
inputLists := Split(input, "!");

// Choose numbers of columns in the order they should appear in the resulting table
// 1: weights k
// 2: Dimensions of VD
// 3: Dimensions of VP0 and VP2
// 4: Dimensions of VG
// 5: Results for U1 = U1^{\Gamma_0(t)}
// 6: Results for U1P0
// 7: Results for U1P2
// 8: Results for T1
// 9: Results for U2 = U2^{\Gamma_0(t)}
// 10: Results for U2P0
// 11: Results for U2P2
// 12: Results for T2
// (13: The sentence "\n \n The program took ... seconds to run.")

numbers := [1, 12, 10, 11, 9];
columns := [];

for i in numbers do
	Append(~columns, Split(inputLists[i], ";\n"));
end for;

l := #columns[1];
output := "";

for i := 1 to l do
	line := "";
	for j := 1 to #columns do
		line := line cat columns[j][i] cat " & ";
	end for;
	line := Substring(line, 1, #line - 3);
	output := output cat line cat "\n \\\\ \\hline \n";
end for;

// open filestream to save result
// Change to "table1_" when using numbers 5-8 and "table2_" when using numbers 9-12
// or choose another name
file := Open("table2_" cat filename cat ".txt", "w");
fprintf file, output;
delete file;

quit;
