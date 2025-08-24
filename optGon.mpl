# -------------------------
# Input:
#     m — an integer m ≥ 3 with an odd factor
# Output:
#     K — an integer list K = [k_1,...,k_r] corresponding to a certain 
#         optimal convex m-gon, where r denotes the smallest odd 
#         prime factor of m, and k_1 = ... = k_r = m / r
giveAnOptIntlist := proc(m)
    local mm, r, k, i, K;
    if not type(m, integer) or m < 3 then
        error "The input m must be a positive integer greater than or equal to 3!";
    end if;

    mm := m;
    while mm mod 2 = 0 do
        mm := mm / 2;
    end do;

    if mm = 1 then
        error "The input m must contain an odd factor!";
    end if;

    r:=mm;
    for k from 3 by 2 to isqrt(mm) do
        if mm mod k = 0 then
            r:=k;
        end if;
    end do;

    K := [seq(m /r, i = 1..r)];
    return K;
end proc:

# -------------------------
# Input:
#     m — an integer m ≥ 3 with an odd factor
#     K — an integer list K = [k_1,...,k_r] corresponding to an 
#         optimal convex m-gon
# Output:
#     angInt — 1×m matrix, the counterclockwise angles of the 
#              optimal convex polygon associated with the list K,
#              from e_1 to e_i, i=1,...,m, normalized by (π/m)
giveAngles := proc(m :: posint, K :: list)
    local r, angInt, i, j, startIdx, endIdx, prefixSum, prefixIndex, m2, tmp, idx;
    r := nops(K);
    if add(K[j], j=1..r) <> m then error "giveAngles: sum(K) must equal m"; end if;
    angInt := [];
    for i from 1 to m do angInt := [op(angInt), 0]; end do;

    startIdx := 1; endIdx := K[1];
    for i from startIdx to endIdx do angInt[i] := i - 1; end do;

    if r >= 2 then
        startIdx := K[1] + 1; endIdx := K[1] + K[2];
        for i from startIdx to endIdx do angInt[i] := m + K[1] + (i - startIdx); end do;
    end if;

    if r >= 3 then
        for prefixIndex from 2 to (r-1)/2 do
            prefixSum := 0;
            for i from 1 to (2*prefixIndex-2) do prefixSum := prefixSum + K[i]; end do;
            startIdx := prefixSum + 1; endIdx := prefixSum + K[2*prefixIndex-1];
            for i from startIdx to endIdx do angInt[i] := prefixSum + (i - startIdx); end do;
            prefixSum := 0;
            for i from 1 to (2*prefixIndex-1) do prefixSum := prefixSum + K[i]; end do;
            startIdx := prefixSum + 1; endIdx := prefixSum + K[2*prefixIndex];
            for i from startIdx to endIdx do angInt[i] := m + prefixSum + (i - startIdx); end do;
        end do;
    end if;

    if r > 1 then
        prefixSum := 0;
        for i from 1 to (r-1) do prefixSum := prefixSum + K[i]; end do;
        startIdx := prefixSum + 1; endIdx := m;
        for i from startIdx to endIdx do angInt[i] := prefixSum + (i - startIdx); end do;
    end if;

    # modulize and sort
    m2 := 2*m;
    for idx from 1 to m do
        tmp := angInt[idx] mod m2;
        angInt[idx] := tmp;
    end do;
    angInt := sort(angInt);
    return angInt;
end proc;

# -------------------------
# Find all distinct optimal convex m-gons obtained by "reversing"
# all edge vectors in a zero-sum edge vector subset of 
# a given optimal convex m-gon (see Lemma 4.2)
#
# Note: Polygons that differ only by translations, rotations, flips, 
#   or their combination are considered equivalent.
#
# Input:
#     m — an integer m ≥ 3 with an odd factor
#     tol — Tolerance for detecting a zero-sum vector set (optional, default 1e-8)
#
# Output:
#   angleSets  — the counterclockwise angles of all optimal polygons, 
#            s×m matrix. Each row corresponds to the set of 
#            the counterclockwise angles of one optimal polygon, 
#            from e_1 to e_i, i=1,...,m.
#   vectorSets — all optimal tight frames, 
#             2×m×s matrix. Each slice (a 2×m matrix) corresponds to 
#             one optimal tight frame.
allOptPolygons := proc(m :: posint, tol :: numeric := 1.0e-12)
    local m2, K, theta0, theta0Int, i, j, mask, Ilist, sumComplex, angleSetsInt, thetaInt;
    local existsFlag, jRow, idx, nullflag, is_all_int, kround, row, theta_C, theta_Cflip, a;
    local tmpList, s_count, angleSets, vectorSetsSym, jj;
    m2 := 2*m;

    # Generate an initial optimal convex m-gon, normalize, add it to the result list
   K := giveAnOptIntlist(m):
   theta0 := giveAngles(m, K):

    is_all_int := true;
    for i from 1 to m do
        if not type(theta0[i], integer) then is_all_int := false; break; end if;
    end do;

    theta0Int := [];
    if is_all_int then
        for i from 1 to m do theta0Int := [op(theta0Int), theta0[i] mod m2]; end do;
    else
        for i from 1 to m do
            kround := round( evalf( (theta0[i] * m) / Pi, 20 ) );
            theta0Int := [op(theta0Int), kround mod m2];
        end do;
    end if;

    tmpList := [];
    for i from 1 to m do tmpList := [op(tmpList), (theta0Int[i] - theta0Int[1]) mod m2]; end do;
    theta0Int := sort(tmpList);

    angleSetsInt := [ theta0Int ];

    # Enumerate all subsets I (excluding the empty and the full set)
    for mask from 1 to 2^m - 2 do
        Ilist := [];
        for j from 1 to m do
            if floor(mask / 2^(j-1)) mod 2 = 1 then
                Ilist := [op(Ilist), j];
            end if;
        end do;

        if nops(Ilist) = 0 then next; end if;

        # sum_{j in I} exp(i * k * pi / m)
        sumComplex := 0;
        for j in Ilist do
            sumComplex := sumComplex + exp( I * Pi * theta0Int[j] / m );
        end do;

        # Verify if the sum of the subset equals the zero vector
        if simplify(Re(sumComplex)) = 0 and simplify(Im(sumComplex)) = 0 then
            nullflag := true;
        else
            nullflag := evalf( abs(evalf(sumComplex, 40)) ) < tol;
        end if;

        if nullflag then
	    printf("I=%a \n", Ilist);
            # reversing all edge vectors in I
            thetaInt := [];
            for i from 1 to m do thetaInt := [op(thetaInt), theta0Int[i]]; end do;
            for j in Ilist do thetaInt[j] := (thetaInt[j] + m) mod m2; end do;

            tmpList := [];
            for i from 1 to m do tmpList := [op(tmpList), (thetaInt[i] - thetaInt[1]) mod m2]; end do;
            thetaInt := sort(tmpList);

            # Remove polygons that are equivalent under translations, rotations, flips, or any combination thereof
            existsFlag := false;
            for jRow from 1 to nops(angleSetsInt) while not existsFlag do
                row := angleSetsInt[jRow];
                for idx from 1 to m while not existsFlag do
                    theta_C := [];
                    for a from 1 to m do theta_C := [op(theta_C), (thetaInt[a] - thetaInt[idx]) mod m2]; end do;
                    theta_C := sort(theta_C);

                    theta_Cflip := [];
                    for a from 1 to m do theta_Cflip := [op(theta_Cflip), (-theta_C[a]) mod m2]; end do;
                    theta_Cflip := sort(theta_Cflip);

                    if theta_C = row or theta_Cflip = row then
                        existsFlag := true;
                    end if;
                end do;
            end do;

            if not existsFlag then angleSetsInt := [op(angleSetsInt), thetaInt]; end if;
        end if;
    end do;

    angleSets := [];
    for jRow from 1 to nops(angleSetsInt) do
        tmpList := [];
        for i from 1 to m do tmpList := [op(tmpList), angleSetsInt[jRow][i] * Pi / m]; end do;
        angleSets := [op(angleSets), tmpList];
    end do;

    s_count := nops(angleSets);
    vectorSetsSym := Array(1..2, 1..m, 1..s_count, datatype = anything);
    for jj from 1 to s_count do
        for i from 1 to m do
            vectorSetsSym[1,i,jj] := 'cos'( angleSets[jj][i] / 2 );
            vectorSetsSym[2,i,jj] := 'sin'( angleSets[jj][i] / 2 );
        end do;
    end do;

    return angleSets, vectorSetsSym;
end proc;

# -------------------------
# Example: 
# Input: m --- an integer m ≥ 3 with an odd factor
m := 12:
angleSets, vectorSetsSym := allOptPolygons(m):
s := nops(angleSets):
printf("\n");
printf("Final results:\n");
printf("Found %d optimal sets\n", s);
printf("\n");

printf("All the optimal convex %d-gons represented by edge vector direction angles: \n", m);
for i from 1 to s do
   angleSets[i];
end do;
printf("\n");

for i from 1 to s do
   printf("The %d th optimal vector set (tight frame): \n", i);
   for j from 1 to m do
      printf("v_%d = (%a, %a) \n", j, vectorSetsSym[1, j, i], vectorSetsSym[2, j, i]);
   end do;
   printf("\n");
end do;
