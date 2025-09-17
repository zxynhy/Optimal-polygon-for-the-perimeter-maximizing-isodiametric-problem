# ----------------------- Function 1 -----------------------------------
# Input:
#     m — an integer m ≥ 3
#     fmtDecimals — Tolerance for detecting zero-sum (optional, default 1e-12)
# Output:
#     All vectors eps = (eps_0,...,eps_{m-1})^T with eps_j in {+1,-1}
#     such that g(eps) = sum_{j=0}^{m-1} eps_j * zeta_m^j = 0,
#     where zeta_m = exp(1i*pi/m)
#
# Notes:
# - The code supports two modes: for small m (default m<=22) it uses brute-force
#   enumeration (2^m). For larger m it uses a meet-in-the-middle (MITM)
#   strategy to reduce time complexity to roughly O(2^{m/2}).
# - Found solutions are merged according to the equivalence relations given in
#   the paper. The canonical representative of each equivalence class is chosen 
#   as the lexicographically smallest binary string obtained by encoding +1 -> '1', -1 -> '0'.
#
# Usage example:
# partition_sols = findEqPartition(m, fmtDecimals);
#
# Caution:
# - For large m (e.g. m>30) the search space can still be very large and may
#   require substantial time and memory. Adjust your resources accordingly.

# findEqPartition (Maple) — corrected with explicit local declarations
findEqPartition := proc(m :: posint, fmtDecimals :: posint := 5)
    local mm, zeta, idx, method;
    local sols, cnt;
    local n1, n2, idx1, idx2;
    local map, rightList, leftList, eps1, eps2, s1, s2, key, keyt;
    local i, k, listk, eps, eps_temp, keysList, S, solutions;
    local j, prev, res, v, s, n, cc, re, im;
    local gen_signs, partial_sum, complex_key, re_s, im_s;
    local vec_to_binstr, R_op, S_k, canonical_rep, crep, canonical_map;
    ####################################################################
    # input checks
    ####################################################################
    if m < 3 then
        error "The input parameter m must be >= 3";
    end if;

    zeta := exp(I * Pi / m);
    idx := [seq(j, j = 0 .. (m - 1))];

    if m <= 22 then
        method := "bruteforce";
    else
        method := "mitm";
    end if;

    ####################################################################
    # brute force (2^n)
    ####################################################################
    gen_signs := proc(n :: nonnegint)
        local prev, res, v;
        if n = 0 then
            return [[]];
        else
            prev := gen_signs(n - 1);
            res := [];
            for v in prev do
                res := [op(res), [op(v), -1], [op(v), 1]];
            end do;
            return res;
        end if;
    end proc;

    ####################################################################
    # sum_{j=1..n} v[j] * zeta^(j-1)
    ####################################################################
    partial_sum := proc(v :: list)
        local j, n, s;
        n := nops(v);
        s := 0;
        for j from 1 to n do
            s := s + v[j] * zeta^(j - 1);
        end do;
        return s;
    end proc;

    ####################################################################
    complex_key := proc(c, decimals :: posint)
        local cc, re, im;
        cc := evalf(c, decimals);
        re := evalf(Re(cc), decimals);
        im := evalf(Im(cc), decimals);
        return cat(convert(re, string), "_", convert(im, string));
    end proc;

    ####################################################################
    # method selection: brute force for small m, MITM otherwise
    ####################################################################
    sols := []; cnt := 0;

    if method = "bruteforce" then
        # printf("Using brute-force enumeration; m=%d\n", m);
        leftList := gen_signs(m);
        for eps in leftList do
        	s := partial_sum(eps);

    		re_s := evalf(Re(s), 2*fmtDecimals);
    		im_s := evalf(Im(s), 2*fmtDecimals);

    		if abs(re_s) < 10^(-fmtDecimals) and abs(im_s) < 10^(-fmtDecimals) then
			# print(eps);
        		sols := [op(sols), eps];
    		end if;
	end do;

    else
        # printf("Using meet-in-the-middle; m=%d\n", m);
        n1 := floor(m / 2);
        n2 := m - n1;
        idx1 := [seq(j, j = 0 .. (n1 - 1))];
        idx2 := [seq(j, j = n1 .. (m - 1))];

        rightList := gen_signs(n2);
        map := table();

        for eps2 in rightList do
            s2 := 0;
            for k from 1 to nops(eps2) do
                s2 := s2 + eps2[k] * zeta^(idx2[k]);
            end do;
            s2 := evalf(s2, fmtDecimals);
            key := complex_key(s2, fmtDecimals);
            if assigned(map[key]) then
                map[key] := [op(map[key]), eps2];
            else
                map[key] := [eps2];
            end if;
        end do;

        leftList := gen_signs(n1);
        for eps1 in leftList do
            s1 := 0;
            for k from 1 to nops(eps1) do
                s1 := s1 + eps1[k] * zeta^(idx1[k]);
            end do;
            s1 := evalf(s1, fmtDecimals);
            keyt := complex_key(-s1, fmtDecimals);
            if assigned(map[keyt]) then
                for eps2 in map[keyt] do
                    sols := [op(sols), [op(eps1), op(eps2)]];
                end do;
            end if;
        end do;
    end if;

    # printf("Found %d candidate sequences (before merging equivalence classes)\n", nops(sols));

    ####################################################################
    # Merge solutions according to the paper's equivalence relations 
    # and keep one canonical representative per equivalence class
    ####################################################################
    vec_to_binstr := proc(v :: list)
        local s, i;
        s := "";
        for i from 1 to nops(v) do
            if v[i] = 1 then
                s := cat(s, "1");
            else
                s := cat(s, "0");
            end if;
        end do;
        return s;
    end proc;

    R_op := proc(x :: list)
        local nloc, y, ii;
        nloc := nops(x);
        y := [seq(0, ii = 1 .. nloc)];
        for ii from 1 to nloc do
            y[ii] := x[nloc - ii + 1];
        end do;
        return y;
    end proc;

    S_k := proc(x :: list, k :: posint)
        local mloc, y, ii;
        mloc := nops(x);
        if k < 1 or k >= mloc then
            error "In S_k, k must satisfy 1 <= k <= m-1";
        end if;
        y := [seq(0, ii = 1 .. mloc)];
        for ii from 1 to (mloc - k) do
            y[ii] := x[ii + k];
        end do;
        for ii from (mloc - k + 1) to mloc do
            y[ii] := - x[ii - (mloc - k)];
        end do;
        return y;
    end proc;

    canonical_rep := proc(eps :: list)
        local mloc, bestStr, ks, k2, applyR, alpha, e, s2, C;
        mloc := nops(eps);
        bestStr := "";
        ks := [0, seq(i, i = 1 .. (mloc - 1))];
        for k2 in ks do
            for applyR from 0 to 1 do
                for alpha in [-1, 1] do
                    e := eps;
                    if k2 > 0 then
                        e := S_k(e, k2);
                    end if;
                    if applyR = 1 then
                        e := R_op(e);
                    end if;
                    e := [seq(alpha * e[i], i = 1 .. nops(e))];
                    s2 := vec_to_binstr(e);
                    if bestStr = "" then
                        bestStr := s2;
                    else
                        C := sort([bestStr, s2]);
                        if C[1] = s2 then
                            bestStr := s2;
                        end if;
                    end if;
                end do;
            end do;
        end do;
        return bestStr;
    end proc;

    canonical_map := table();
    for i from 1 to nops(sols) do
        crep := canonical_rep(sols[i]);
        if not assigned(canonical_map[crep]) then
            canonical_map[crep] := sols[i];
        end if;
    end do;

    keysList := [indices(canonical_map)];
    S := nops(keysList);
    solutions := [];
    for i from 1 to S do
	key := keysList[i];
	if type(key, list) and nops(key) = 1 then
		key := key[1];
	end if;
	if assigned(canonical_map[key]) then
        	eps_temp := canonical_map[key];
        	if eps_temp[1] = -1 then
            		eps_temp := [seq(-eps_temp[j], j = 1 .. nops(eps_temp))];
        	end if;
	end if;
        solutions := [op(solutions), eps_temp];
    end do;

    # printf("After merging, %d equivalence classes remain\n", nops(solutions));
    return solutions;
end proc;


# -------------------------- Function 2 -----------------------------------
# Input:
#     m — an integer m ≥ 3 with an odd factor
# Output:
#     (1) all optimal tight frames (m-vector sets) in R^2
#         that achieve the minimal condition number
#         for the stability of phase retrieval;
#     (2) all optimal convex m-polygons in R^2
#         for the perimeter-maximizing isodiametric problem

allOptPolygons := proc(m :: posint, tol :: numeric := 1.0e-12)
	local mm, sCount, th, partitionSols, eps, edgeVectorSets, frameVectorSets;
	local k, j;

	# input checks
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

	# Find all solutions of the partition of roots of unity problem
	partitionSols := findEqPartition(m, numeric);

	# Turn into opt gons and opt frames
	sCount := nops(partitionSols);
	printf("---------------------------- \n");
	printf("Found %d distinct optimal convex %d-gons.\n", sCount, m);
	printf("Found %d distinct optimal tight frames (%d-vector sets).\n", sCount, m);
	edgeVectorSets := Array(1..2, 1..m, 1..sCount, datatype = anything);
	frameVectorSets := Array(1..2, 1..m, 1..sCount, datatype = anything);
	for k from 1 to sCount do
		eps := partitionSols[k];
		th := [];
		for j from 1 to m do
			if eps[j] = 1 then
				th := [op(th), (j-1)*Pi/m];
			else
				th := [op(th), (j-1+m)*Pi/m];
			end if;
		end do;
		thSort := sort(th, (a, b) -> evalf(a,16) < evalf(b,16));
		for j from 1 to m do
			edgeVectorSets[1,j,k] := 'cos'( thSort[j] );
			edgeVectorSets[2,j,k] := 'sin'( thSort[j] );
			frameVectorSets[1,j,k] := 'cos'( thSort[j] / 2 );
			frameVectorSets[2,j,k] := 'sin'( thSort[j] / 2 );
		end do;
	end do;

	return edgeVectorSets, frameVectorSets;
end proc;

# ------------------------------- Example ------------------------------------
# Input: m --- an integer m ≥ 3 with an odd factor
m := 12:
edgeVectorSets,frameVectorSets := allOptPolygons(m):
S := upperbound(edgeVectorSets, 3):

printf("\n");
printf("---------------------------- \n");

for k from 1 to S do
	printf("The %d th edge vector set: \n", k);
	for j from 1 to m do
        	printf("v_%d = (%a, %a) \n", j, edgeVectorSets[1, j, k], edgeVectorSets[2, j, k]);
   	end do;
   	printf("\n");
end do;
printf("---------------------------- \n");

for k from 1 to S do
	printf("The %d th frame: \n", k);
        for j from 1 to m do
        	printf("v_%d = (%a, %a) \n", j, frameVectorSets[1, j, k], frameVectorSets[2, j, k]);
   	end do;
   	printf("\n");
end do;
# -----------------------------------------------------------


