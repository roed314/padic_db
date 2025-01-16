/********************************************
This file contains an intrinsic for counting p-adic extensions by Galois group,
as long as the inputs have 2-power order.
It uses the description of the maximal pro-2 quotient of Gal(Q2bar/Q2) as a free pro-2 group on three generators x,y,z modulo one relation: x^2 y^4 y^-1 z^-1 y z

It can be extended to other contexts by modifying the abelian case (Abelian2AdicSolutions) and the lifting check in LiftSolutions.  BuildQuotientGraph is generic.

EXAMPLES:

> T8 := Get8Ts("Roe/8T2pow.txt");
> X := Count2Adic2Pow(T8);
hashing groups
hashing complete
Finding quotients 1/21: 8T1 abelian
Finding quotients 2/21: 8T2 abelian
Finding quotients 3/21: 8T3 abelian
Finding quotients 4/21: 8T4 22 : 4.2 of size 4 with kernel of size 2
Finding quotients 5/22: 8T5 22 : 4.2 of size 4 with kernel of size 2
Finding quotients 6/22: 8T6 4 : 8T4 of size 8 with kernel of size 2
Finding quotients 7/22: 8T7 2 : 8T2 of size 8 with kernel of size 2
Finding quotients 8/22: 8T8 4 : 8T4 of size 8 with kernel of size 2
Finding quotients 9/22: 8T9 3 : 8T3 of size 8 with kernel of size 2
Finding quotients 10/22: 8T10 2 : 8T2 of size 8 with kernel of size 2
Finding quotients 11/22: 8T11 3 : 8T3 of size 8 with kernel of size 2
Finding quotients 12/22: 8T15 9 : 8T9 of size 16 with kernel of size 2
Finding quotients 13/22: 8T16 10 : 8T10 of size 16 with kernel of size 2
Finding quotients 14/22: 8T17 10 : 8T10 of size 16 with kernel of size 2
Finding quotients 15/22: 8T18 3 : 8T3 of size 8 with kernel of size 4
Finding quotients 16/22: 8T19 10 : 8T10 of size 16 with kernel of size 2
Finding quotients 17/22: 8T26 15 : 8T18 of size 32 with kernel of size 2
Finding quotients 18/22: 8T27 16 : 8T19 of size 32 with kernel of size 2
Finding quotients 19/22: 8T29 15 : 8T18 of size 32 with kernel of size 2
Finding quotients 20/22: 8T30 16 : 8T19 of size 32 with kernel of size 2
Finding quotients 21/22: 8T35 19 : 8T29 of size 64 with kernel of size 2
Finding quotients 22/22: 4.2 abelian
Abelian base case 4.2: 11 triples
Abelian base case 8T1: 24 triples
Abelian base case 8T2: 18 triples
Lifting with 11 triples below, 3 aut-costs, 2 kernel elements; to 8T4 (4) from 4.2 (22)...found 18 triples in 0.000s
Lifting with 11 triples below, 1 aut-costs, 2 kernel elements; to 8T5 (5) from 4.2 (22)...found 6 triples in 0.000s
Abelian base case 8T3: 1 triples
Lifting with 1 triples below, 21 aut-costs, 2 kernel elements; to 8T9 (9) from 8T3 (3)...found 9 triples in 0.000s
Lifting with 1 triples below, 28 aut-costs, 2 kernel elements; to 8T11 (11) from 8T3 (3)...found 16 triples in 0.010s
Lifting with 18 triples below, 1 aut-costs, 2 kernel elements; to 8T10 (10) from 8T2 (2)...found 12 triples in 0.000s
Lifting with 18 triples below, 2 aut-costs, 2 kernel elements; to 8T7 (7) from 8T2 (2)...found 36 triples in 0.000s
Lifting with 18 triples below, 1 aut-costs, 2 kernel elements; to 8T6 (6) from 8T4 (4)...found 16 triples in 0.000s
Lifting with 18 triples below, 2 aut-costs, 2 kernel elements; to 8T8 (8) from 8T4 (4)...found 36 triples in 0.010s
Lifting with 12 triples below, 4 aut-costs, 2 kernel elements; to 8T17 (14) from 8T10 (10)...found 48 triples in 0.000s
Lifting with 1 triples below, 28 aut-costs, 4 kernel elements; to 8T18 (15) from 8T3 (3)...found 4 triples in 0.010s
Lifting with 12 triples below, 2 aut-costs, 2 kernel elements; to 8T19 (16) from 8T10 (10)...found 24 triples in 0.000s
Lifting with 12 triples below, 1 aut-costs, 2 kernel elements; to 8T16 (13) from 8T10 (10)...found 12 triples in 0.010s
Lifting with 9 triples below, 8 aut-costs, 2 kernel elements; to 8T15 (12) from 8T9 (9)...found 38 triples in 0.010s
Lifting with 24 triples below, 1 aut-costs, 2 kernel elements; to 8T30 (20) from 8T19 (16)...found 24 triples in 0.010s
Lifting with 4 triples below, 12 aut-costs, 2 kernel elements; to 8T26 (17) from 8T18 (15)...found 24 triples in 0.020s
Lifting with 4 triples below, 8 aut-costs, 2 kernel elements; to 8T29 (19) from 8T18 (15)...found 16 triples in 0.020s
Lifting with 24 triples below, 2 aut-costs, 2 kernel elements; to 8T27 (18) from 8T19 (16)...found 48 triples in 0.020s
Lifting with 16 triples below, 6 aut-costs, 2 kernel elements; to 8T35 (21) from 8T29 (19)...found 48 triples in 0.080s
Success!  Quotient graph built in 0.770s, lifting step complete in 0.230s.

In the quotient stage, the total number of groups can increase (here from 21 to 22) in order to be able to recurse down to abelian groups.  Note that the counts are of Galois closed fields; counts of stem fields can be recovered by using the sibling data from the LMFDB.
********************************************/



/*
# Sage code for producing the files needed for Get16Ts and Get8Ts (just change 16 to 8 in query)

from lmfdb import db
X = list(db.gps_transitive.search({"n": 16}, ["label", "order", "siblings"]))
Y = []
seen = set()
for rec in X:
    if not rec["order"].is_power_of(2) or rec["label"] in seen:
        continue
    Y.append(rec["label"])
    for z in rec["siblings"]:
        seen.add(f"{z[0][0]}T{z[0][1]}")
with open("8T2pow.txt", "w") as F:
    for x in Y:
        _ = F.write(x+"\n")
*/

intrinsic Get16Ts(filename) -> SeqEnum[MonStgElt]
{}
    descs := Split(Read(filename), "\n");
    Gs := [StringToGroup(x) : x in descs];
    Gok := [i : i in [1..#Gs] | FrattiniQuotientRank(PCGroup(Gs[i])) lt 4];
    return [descs[i] : i in Gok];
end intrinsic;

intrinsic Get8Ts(filename) -> SeqEnum[MonStgElt]
{}
    descs := Split(Read(filename), "\n");
    Gs := [StringToGroup(x) : x in descs];
    Gok := [i : i in [1..#Gs] | FrattiniQuotientRank(PCGroup(Gs[i])) lt 4];
    return [descs[i] : i in Gok];
end intrinsic;

function CyclicGenerator(G)
    assert IsCyclic(G);
    gens := [x : x in Generators(G)];
    if #gens eq 1 then
        return gens[1];
    end if;
    // Pick maximal elements under divisibility ordering
    // Since G is cyclic, the sum/product of these will be a generator
    maxs := [];
    for g in gens do
        m := Order(g);
        if &or[x[2] mod m eq 0 : x in maxs] then
            //Already have a generator of order a multiple of this one
            continue;
        end if;
        maxs := [x : x in maxs | m mod x[2] ne 0];
        Append(~maxs, <g, m>);
    end for;
    if Type(G) eq GrpAb then
        return &+[x[1] : x in maxs];
    else
        return &*[x[1] : x in maxs];
    end if;
end function;

intrinsic Build2AdicOrdHsh(Gs::SeqEnum[Grp] : verbose:=true) -> Assoc
{}
    maxN := LCM([#G : G in Gs]);
    assert IsPrimePower(maxN) and Valuation(maxN, 2) gt 0;
    by_ordhsh := AssociativeArray();
    tp := 1;
    while tp le maxN do
        by_ordhsh[tp] := AssociativeArray();
        tp *:= 2;
    end while;
    if verbose then print "hashing groups"; end if;
    for i in [1..#Gs] do
        G := Gs[i];
        m := #G;
        hsh := hash(G);
        if not IsDefined(by_ordhsh[m], hsh) then
            by_ordhsh[m][hsh] := [];
        end if;
        Append(~by_ordhsh[m][hsh], i);
    end for;
    if verbose then print "hashing complete"; end if;
    return by_ordhsh;
end intrinsic;

intrinsic BuildQuotientGraph(descs::SeqEnum[MonStgElt], Gs::SeqEnum[Grp], by_ordhsh::Assoc : verbose:=true) -> Assoc, Assoc, SeqEnum, SeqEnum[Grp], Assoc
{-> quotients, auts, aut_imgs}
    quotients := AssociativeArray();
    auts := AssociativeArray();
    aut_imgs := []; // will only be defined for i when Gs[i] is nonabelian
    i := 1;
    while i le #Gs do
        G := Gs[i];
        if verbose then
            printf "Finding quotients %o/%o: %o ", i, #Gs, (i le #descs select descs[i] else GroupToString(G));
        end if;
        if i le #descs then
            A := AutomorphismGroup(G);
            Ambient, inj, proj := Holomorph(G, A);
            auts[i] := <A, Ambient, inj, proj>;
        else
            // We've stored this information when we found this as a quotient
            A, Ambient, inj, proj := Explode(auts[i]);
        end if;
        if IsAbelian(G) then
            if verbose then print "abelian"; end if;
        else
            //Ns := SetToSequence({ncl<Ambient|inj(N)> @@ inj : N in MinimalNormalSubgroups(G)});
            Ns := MinimalNormalSubgroups(G);
            Ns cat:= [ncl<Ambient|inj(N)> @@ inj : N in Ns];
            Ns := SetToSequence({N : N in Ns});
            // * Among minimal normal subgroups and their characteristic closures,
            // * we have two goals:
            //   * minimize (#N)^3 * aut-quo-index
            //   * use a quotient that we've already seen before or that's abelian.
            Qs := [];
            Qprojs := AssociativeArray();
            for j in [1..#Ns] do
                Q, Qproj := quo<G|Ns[j]>;
                Append(~Qs, Q);
                Qprojs[j] := Qproj;
            end for;
            AbQ := [j : j in [1..#Ns] | IsAbelian(Qs[j])];
            if #AbQ gt 0 then
                Ns := [Ns[j] : j in AbQ];
                Qs := [Qs[j] : j in AbQ];
                // Ugh.  I wish using a SeqEnum for Qprojs didn't give type errors
                rev := AssociativeArray();
                for k in [1..#AbQ] do
                    rev[k] := Qprojs[AbQ[k]];
                end for;
                Qprojs := rev;
            end if;
            AQs := [];
            AQimgs := [];
            quos := [];
            founds := [];
            isos := [];
            hshs := [];
            for j in [1..#Ns] do
                Q := Qs[j];
                Qproj := Qprojs[j];
                m := #Q;
                hsh := hash(Q);
                found := false;
                if IsDefined(by_ordhsh[m], hsh) then
                    for k in by_ordhsh[m][hsh] do
                        iso, phi := IsIsomorphic(Q, Gs[k]);
                        if iso then
                            founds[j] := k;
                            isos[j] := phi;
                            found := true;
                            break;
                        end if;
                    end for;
                end if;
                if not found then
                    founds[j] := 0;
                    hshs[j] := hsh;
                end if;
            end for;
            foundj := [j : j in [1..#Ns] | founds[j] ne 0];
            if #foundj gt 0 then
                Ns := [Ns[j] : j in foundj];
                Qs := [Gs[founds[j]] : j in foundj];
                // Ugh.  I wish using a SeqEnum for Qprojs didn't give type errors
                rev := AssociativeArray();
                for k in [1..#foundj] do
                    j := foundj[k];
                    rev[k] := Qprojs[j] * isos[j];
                end for;
                Qprojs := rev;
                founds := [founds[j] : j in foundj];
                AQs := [Explode(auts[k]) : k in founds];
            else
                AQs := [AutomorphismGroup(Q) : Q in Qs];
            end if;
            costs := []; // Cost heuristic: (#N)^3 #(AQ/AQimg)
            for j in [1..#Ns] do
                N := Ns[j];
                Q := Qs[j];
                AQ := AQs[j];
                Qproj := Qprojs[j];
                Qgens := [x : x in Generators(Q)];
                lifts := [x @@ Qproj : x in Qgens];
                Stab := Normalizer(Ambient, inj(Ns[j]));
                AQimgs[j] := sub<AQ| [hom<Q -> Q | [<Qgens[k], (inj(lifts[k])^n) @@ inj @ Qproj> : k in [1..#Qgens]]> : n in Generators(Stab)]>;
                costs[j] := (#N)^3 * #AQs[j] div #AQimgs[j];
            end for;
            minCost := Min(costs);
            cheapest := [j : j in [1..#Ns] | costs[j] eq minCost];
            costs := [];
            if #cheapest gt 1 and #AbQ gt 0 then
                // prefer larger rank and more balanced abelian invariants when the Q are abelian,
                minrank := Min([#AbelianInvariants(Qs[j]) : j in cheapest]);
                cheapest := [j : j in cheapest | #AbelianInvariants(Qs[j]) eq minrank];
                if #cheapest gt 1 then
                    minsum := Min([&+AbelianInvariants(Qs[j]) : j in cheapest]);
                    cheapest := [j : j in cheapest | &+AbelianInvariants(Qs[j]) eq minsum];
                end if;
            end if;
            if #cheapest gt 1 then
                // Now take the largest N, since this will make Q smaller
                minN := Min([#Ns[j] : j in cheapest]);
                cheapest := [j : j in cheapest | #Ns[j] eq minN];
            end if;
            // Finally, we just choose the first
            j := cheapest[1];
            Q := Qs[j];
            AmbientQ, injQ, projQ := Holomorph(Q, AQs[j]);
            if #foundj gt 0 then
                // Quotient group has already been added
                k := founds[j];
            else
                // Need to add a new group
                Append(~Gs, Q);
                k := #Gs;
                if not IsDefined(by_ordhsh[m], hshs[j]) then
                    by_ordhsh[m][hshs[j]] := [];
                end if;
                Append(~by_ordhsh[m][hshs[j]], k);
            end if;
            aut_imgs[i] := AQimgs[j];
            auts[k] := <AQs[j], AmbientQ, injQ, projQ>;
            quotients[i] := <k, Qprojs[j]>;
            if verbose then
                print k, ":", (k le #descs select descs[k] else GroupToString(Gs[k])), "of size", #Gs[k], "with kernel of size", #Ns[j];
            end if;
        end if;
        i +:= 1;
    end while;
    return quotients, auts, aut_imgs, Gs, by_ordhsh;
end intrinsic;

function CycSol(G)
    // We need at least one of x,y,z to generate on its own
    // It can't be x, by the relation.
    // In the case that z generates, we can set it to 1 using the aut group, removing the freedom
    // Now y can be anything, and x=-2y or N/2-2y
    // If z does not generate, y must, so we can set y to 1, x to -2 or N/2-2, and z to any non-unit.
    g := CyclicGenerator(G);
    N := #G;
    N2 := N div 2;
    ans := [];
    for y in G do
        Append(~ans, <g, y, y^-2>);
        Append(~ans, <g, y, y^(N2-2)>);
    end for;
    for k in [0..N2-1] do
        Append(~ans, <g^(2*k), g, g^-2>);
        Append(~ans, <g^(2*k), g, g^(N2-2)>);
    end for;
    //print "RANK 1", N, ans;
    return ans;
end function;

function AbRnk2Sol(G)
    // Z/2^a x Z/2^b with a <= b, A=2^a, B=2^b
    basis, AB := PrimaryAbelianBasis(G);
    A, B := Explode(AB);
    r, s := Explode(basis);
    twos := [Identity(G), r^(A div 2), s^(B div 2), r^(A div 2) * s^(B div 2)];
    ans := [];
    // Up to automorphism, there's a unique element of maximal order and it must be among our generators
    // With the notation above, this is s.
    // If z maps to s and y generates the rest, then we can act diagonally to send y to r or y to r*s.
    for t in twos do
        Append(~ans, <t*r^-2, r, s>);
        Append(~ans, <t*r^-2, r*s, s>);
    end for;
    // If y maps to s but z does not map to a generator of order B, we must have that a < b
    if A lt B then
        for t in twos do
            Append(~ans, <t*s^-2, s, r>);
        end for;
    end if;
    if A eq 2 then
        // In this case x can actually contribute to generating the group (e.g. if x is the first generator and y is the second, we're done)
        // Now we can send z to s, y to a power of s and x to r*y^-2
        // But some of these will be related by automorphisms, so we only need one value of x rather than the normal 2
        for k in [0..B-1] do
            y := s^k;
            Append(~ans, <r*y^-2, y, s>);
        // Not needed: Append(~ans, <r*y^-2*twos[3], y, s>);
        end for;
        // Or we can send y to s, z to a non-unit power of s, and x to r*s^-2
        for k in [0..(B div 2)-1] do
            z := s^(2*k);
            Append(~ans, <r*s^-2, s, z>);
        // Not needed: Append(~ans, <r*s^-2*twos[3], y, z>);
        end for;
    end if;
    //print "RANK 2", A, B, ans;
    return ans;
end function;

function AbRnk3Sol(G)
    basis, ABC := PrimaryAbelianBasis(G);
    A, B, C := Explode(ABC);
    // Now x must contribute to generating the group, so we must have A = 2
    ans := [];
    if A eq 2 then
        u, r, s := Explode(basis);
        // Have to handle the cases with more automorphisms first
        if B eq 2 then
            if C eq 2 then
                ans := [<u, r, s>];
            else
                // x has to help with u or r, and can take cases based on whether z maps to an element of maximal order or not
                ans := [<u, r, s>, <u, s, r>];
            end if;
        else
            twos := [u];
            for pair in [<r,B>, <s,C>] do
                twos cat:= [t * pair[1]^(pair[2] div 2) : t in twos];
            end for;
            // Still a unique element of maximal order
            // If z maps to such an element, we can take it to be s, and y can be conjugated to r
            for t in twos do
                Append(~ans, <t*r^-2, r, s>);
            end for;
            // If not, we must have B < C, and then we can map y to s and z to r instead
            for t in twos do
                Append(~ans, <t*s^-2, s, r>);
            end for;
        end if;
    end if;
    //print "RANK 3", A, B, C, ans;
    return ans;
end function;

intrinsic Abelian2AdicSolutions(G::GrpPC) -> SeqEnum
{}
    // Looking for triples x,y,z up to automorphism that generate G and with (using additive notation now)
    // 2x + 4y = 0
    // Break it up by cases based on the rank
    rank := FrattiniQuotientRank(G);
    if rank eq 1 then
        return CycSol(G);
    elif rank eq 2 then
        return AbRnk2Sol(G);
    elif rank eq 3 then
        return AbRnk3Sol(G);
    else
        return [];
    end if;
end intrinsic;

intrinsic LiftSolutions(i::RngIntElt, gens::Assoc, quotients::Assoc, auts::Assoc, aut_imgs::SeqEnum, descs::SeqEnum[MonStgElt] : verbose:=true) -> SeqEnum
{}
    t0 := Cputime();
    j, Qproj := Explode(quotients[i]);
    G := Domain(Qproj);
    N := Kernel(Qproj);
    Q := Codomain(Qproj);
    A, Hol, inj := Explode(auts[i]);
    phi := PermutationRepresentation(A);
    PA := Image(phi);
    Alist := [a @@ phi : a in PA];
    AQ, HolQ, injQ := Explode(auts[j]);
    AQimg := aut_imgs[i];
    phiQ := PermutationRepresentation(AQ);
    PAQ := Image(phiQ);
    // GrpAuto is missing some features which would make this nicer
    phiAQimg := sub<PAQ | [phiQ(g) : g in Generators(AQimg)]>;
    translators := [(f^-1) @@ phiQ : f in RightTransversal(PAQ, phiAQimg)]; // Using inverses gives a left transversal
    if verbose then
        printf "Lifting with %o triples below, %o aut-costs, %o kernel elements; to %o (%o) from %o (%o)...",
               #gens[j], #translators, #N, (i le #descs) select descs[i] else GroupToString(G), i, (j le #descs) select descs[j] else GroupToString(Q), j;
    end if;
    // Create the restriction map from a subgroup of A to AQ.
    // It is only defined on the stabilizer of N within A
    one := Identity(G);
    seen := {};
    ans := [];
    for trip in gens[j] do
        for f in translators do
            x0, y0, z0 := Explode(<f(t) @@ Qproj : t in trip>);
            for nx in N do
                x := x0 * nx;
                for ny in N do
                    y := y0 * ny;
                    for nz in N do
                        z := z0 * nz;
                        // Note that the lifts automatically generate,
                        // since they generate after projecting down to an abelian quotient
                        if x^2 * y^3 * z^-1 * y * z eq one then
                            // ?? Since N was characteristic, no need to remove duplicates
                            if not <x,y,z> in seen then
                                Append(~ans, <x,y,z>);
                                for a in Alist do
                                    Include(~seen, <a(x), a(y), a(z)>);
                                end for;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    if verbose then
        printf "found %o triples in %os\n", #ans, Cputime() - t0;
    end if;
    return ans;
end intrinsic;

intrinsic FindAllSolutions(quotients::Assoc, auts::Assoc, aut_imgs::SeqEnum, descs::SeqEnum[MonStgElt], Gs::SeqEnum[Grp], by_ordhsh::Assoc : verbose:=true) -> Assoc
{}
    gens := AssociativeArray();
    for m in Sort([x : x in Keys(by_ordhsh)]) do
        for hsh->gps in by_ordhsh[m] do
            for i in gps do
                G := Gs[i];
                // Either G is abelian or we have a map to a group where the images have been found
                gens[i] := [];
                if #G eq 1 then
                    one := Identity(G);
                    gens[i] := [<one, one, one>];
                elif IsAbelian(G) then
                    gens[i] := Abelian2AdicSolutions(G);
                    if verbose then
                        printf "Abelian base case %o: %o triples\n", (i le #descs) select descs[i] else GroupToString(G), #gens[i];
                    end if;
                else
                    gens[i] := LiftSolutions(i, gens, quotients, auts, aut_imgs, descs : verbose:=verbose);
                end if;
            end for;
        end for;
    end for;
    return gens;
end intrinsic;

intrinsic Count2Adic2Pow(descs::SeqEnum[MonStgElt] : verbose:=true) -> SeqEnum
{Descriptions should be of NON-ISOMORPHIC 2-groups}
    // Count number of solutions to x^2y^4[y,z] = 1, modulo the automorphism group
    // We start by iteratively finding quotients from which we can lift solutions.
    Gs := [PCGroup(StringToGroup(desc)) : desc in descs];

    by_ordhsh := Build2AdicOrdHsh(Gs : verbose:=verbose);

    t0 := Cputime();
    quotients, auts, aut_imgs, Gs, by_ordhsh := BuildQuotientGraph(descs, Gs, by_ordhsh : verbose:=verbose);

    // Preparations are done: now we iteratively build up the list of image generators
    t1 := Cputime();
    gens := FindAllSolutions(quotients, auts, aut_imgs, descs, Gs, by_ordhsh : verbose:=verbose);

    t2 := Cputime();
    if verbose then
        printf "Success!  Quotient graph built in %os, lifting step complete in %os.\n", t1 - t0, t2 - t1;
    end if;
    return [<descs[i], #gens[i]> : i in [1..#descs]];
end intrinsic;
