// mToM: This stands for "matrix to matrix". The function returns a
// transformation matrix MV describing the action of a matrix M on V.
// In our case, the entries will be in A=F_q[t] or K=F_q(t).
// Input:
    // R the field or ring we are working over,
    // s list with tuples of indices, allowing to switch between Magma basis for V and our basis with double indices
    // k the current weight
    // type the current type
    // M the matrix acting on V
// Output:
        // MV a #s x #s matrix with entries in R
mToM := function(R, s, k, type, M)

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

            // This is where the type comes in
            MV[a2, a1] := Determinant(M)^(1-type)*sum;
        end for;
    end for;

    return MV;
end function;
