function GetByOrd(N)
    by_ord := AssociativeArray();
    for i in [1..#N] do
        H := N[i];
        if not IsDefined(by_ord, #H) then
            by_ord[#H] := {};
        end if;
        Include(~by_ord[#H], i);
    end for;
    return by_ord;
end function;

function GetSubs(N, by_ord)
    subs := AssociativeArray();
    for i in [1..#N] do
        H := N[i];
        subs[i] := {};
        for m in Divisors(#H) do
            if m eq #H then continue; end if;
            for j in by_ord[m] do
                K := N[j];
                if K subset H then
                    Include(~subs[i], j);
                end if;
            end for;
        end for;
    end for;
    return subs;
end function;

intrinsic WildFiltrations(P::Grp, p::RngIntElt : N:=[], by_ord:=[], subs:=[], subP:=[]) -> SeqEnum
{descending filtrations}
    // N is a list of normal subgroups of an ambient group
    // subs is an associative array giving indexes in N
    // subP is the list of indexes i so that N[i] is contained in P
    if #N eq 0 then
        N := [H`subgroup : H in NormalSubgroups(P)];
        by_ord := GetByOrd(N);
        subs := GetSubs(N, by_ord);
        subP := {1..#N};
    else
        for i in Keys(subs) do
            // Via a feature that usually annoys me about Magma, I don't think this will modify subs outside this intrinsic
            subs[i] meet:= subP;
        end for;
        for m in Keys(by_ord) do
            by_ord[m] meet:= subP;
        end for;
    end if;
    Pi := rep{x : x in by_ord[#P]};
    elem_subs := AssociativeArray();
    for i in subP do
        P := N[i];
        F := FrattiniSubgroup(P);
        Fi := rep{j : j in by_ord[#F] | N[j] eq F};
        elem_subs[i] := {j : j in subs[i] | j eq Fi or Fi in subs[j]};
    end for;
    from_here := function(i)
        if #N[i] eq 1 then return [[]]; end if;
        Filts := [];
        for j in elem_subs[i] do
            for sfilt in $$(j) do
                Append(~Filts, [N[j]] cat sfilt);
            end for;
        end for;
        return Filts;
    end function;
    return from_here(Pi);
end intrinsic;

intrinsic pAdicFiltrations(G::Grp, p::RngIntElt) -> SeqEnum
{descending filtrations}
    N := [H`subgroup : H in NormalSubgroups(G)];
    pprime := {q : q in PrimeDivisors(#G) | q ne p};
    by_ord := GetByOrd(N);
    subs := GetSubs(N, by_ord);
    Pmax := pCore(G, p);
    Pmaxi := 0;
    for i in by_ord[#Pmax] do
        if Pmax eq N[i] then
            Pmaxi := i;
            break;
        end if;
    end for;
    Pmin := DerivedSubgroup(Pmax);
    Pmini := 0;
    for i in by_ord[#Pmin] do
        if Pmin eq N[i] then
            Pmini := i;
            break;
        end if;
    end for;
    Filts := [];
    for i in {Pmaxi} join subs[Pmaxi] do
        P := N[i];
        if Pmini in {Pmini} join subs[i] and IsCyclic(quo<Pmax|P>) then
            Q, Qproj := quo<G | P>;
            if #pprime gt 0 then
                hall := Core(Q, HallSubgroup(Q, pprime));
            else
                hall := sub<Q|>;
            end if;
            // tame piece must be contained in hall, with cyclic quotient
            tames := [];
            for tame in [X`subgroup : X in CyclicSubgroups(hall)] do
                if not IsNormal(Q, tame) then continue; end if;
                U, Uproj := quo<Q | tame>;
                if not IsCyclic(U) then continue; end if;
                tau := Q!rep{x : x in tame | Order(x) eq #tame};
                sigma := rep{x : x in U | Order(x) eq #U};
                sigma := sigma@@Uproj;
                ts := tau^sigma;
                // I don't know how to do discrete log for PC groups, but e should be very small
                e := Order(tau);
                R := ResidueClassRing(e);
                if e eq 1 then
                    Append(~tames, tame@@Qproj);
                    continue;
                end if;
                m := rep{x : x in [1..e-1] | tau^x eq ts};
                x := 1;
                pe := p mod e;
                while true do
                    if x eq pe then
                        Append(~tames, tame@@Qproj);
                        break;
                    end if;
                    x := x*m mod e;
                    if x eq 1 then
                        break;
                    end if;
                end while;
            end for;
            if #tames gt 0 then
                wilds := WildFiltrations(P, p : N:=N, by_ord:=by_ord, subs:=subs, subP:={i} join subs[i]);
                for tame in tames do
                    for W in wilds do
                        Append(~Filts, [tame, P] cat W);
                    end for;
                end for;
            end if;
        end if;
    end for;
    return Filts;
end intrinsic;

intrinsic TransitiveSubs(G::Grp) -> SeqEnum
{The core-free subgroups of minimal index}
    S := [H`subgroup : H in Subgroups(G) | #Core(G, H`subgroup) eq 1];
    m := Min([Index(G, H) : H in S]);
    return [H : H in S | Index(G, H) eq m];
end intrinsic;

intrinsic pAdicFiltrationClasses(G::Grp, H::Grp, p::RngIntElt) -> SeqEnum
{}
    t0 := Cputime();
    Filts := pAdicFiltrations(G, p);
    print #Filts, "filtrations found in", Cputime() - t0;
    t0 := Cputime();
    rho := CosetAction(G, H);
    PG := Image(rho);
    Ambient, inj := Holomorph(G);
    print "Automorphism group of size", #Ambient div #G, "found in", Cputime() - t0;
    t0 := Cputime();
    T0 := t0;
    by_sizes := AssociativeArray();
    cnt := 1;
    for filt in Filts do
        if cnt mod 1000 eq 0 then
            print "filt", cnt, Cputime() - t0;
            t0 := Cputime();
        end if;
        cnt +:= 1;
        svec := [#Gi : Gi in filt];
        if not IsDefined(by_sizes, svec) then
            by_sizes[svec] := {};
        end if;
        Include(~by_sizes[svec], [inj(Gi) : Gi in filt]);
    end for;
    //F := {[inj(Gi) : Gi in filt] : filt in Filts};
    print #by_sizes, "Inj done", Cputime() - T0;
    t0 := Cputime();
    T0 := t0;
    ans := [];
    cnt := 1;
    for svec->F in by_sizes do
        if cnt mod 20 eq 0 then
            print cnt, Cputime() - t0;
            t0 := Cputime();
        end if;
        cnt +:= 1;
        FxA := CartesianProduct(F, Ambient);
        action := map<FxA -> F | x :-> [Gi^x[2] : Gi in x[1]]>;
        //print "Action done", Cputime() - t0;
        GS := GSet(Ambient, F, action);
        //print "GSet created in", Cputime() - t0;
        //t0 := Cputime();
        orbs := Orbits(Ambient, GS);
        //print #orbs, "orbits found in", Cputime() - t0;
        //t0 := Cputime();
        ans cat:= [[Gi@@inj : Gi in orb[1]] : orb in orbs];
    end for;
    print "Orbits finished in", Cputime() - T0;
    return ans;
end intrinsic;

intrinsic Get16Ts(filename) -> SeqEnum[MonStgElt]
{}
    // this file contains transitive groups of degree 16 with order a power of 2 and having siblings removed
    // /home/roed/16T2pow.txt
    descs := Split(Read(filename), "\n");
    Gs := [StringToGroup(x) : x in descs];
    Gok := [i : i in [1..#Gs] | FrattiniQuotientRank(PCGroup(Gs[i])) lt 4];
    return [descs[i] : i in Gok];
end intrinsic;

intrinsic Get8Ts(filename) -> SeqEnum[MonStgElt]
{}
    // this file contains transitive groups of degree 16 with order a power of 2 and having siblings removed
    // /home/roed/8T2pow.txt
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

intrinsic Count2Adic2Pow(descs::SeqEnum[MonStgElt]) -> SeqEnum
{Descriptions should be of NON-ISOMORPHIC 2-groups}
    // Count number of solutions to x^2y^4[y,z] = 1, modulo the automorphism group
    // We start by iteratively finding quotients from which we can lift solutions.
    Gs := [PCGroup(StringToGroup(desc)) : desc in descs];
    maxN := LCM([#G : G in Gs]);
    assert IsPrimePower(maxN) and Valuation(maxN, 2) gt 0;
    by_ordhsh := AssociativeArray();
    tp := 1;
    while tp le maxN do
        by_ordhsh[tp] := AssociativeArray();
        tp *:= 2;
    end while;
    print "hashing groups";
    for i in [1..#Gs] do
        G := Gs[i];
        m := #G;
        hsh := hash(G);
        if not IsDefined(by_ordhsh[m], hsh) then
            by_ordhsh[m][hsh] := [];
        end if;
        Append(~by_ordhsh[m][hsh], i);
    end for;
    print "hashing complete";
    quotients := AssociativeArray();
    auts := AssociativeArray();
    aut_imgs := []; // will only be defined for i when Gs[i] is nonabelian
    i := 1;
    while i le #Gs do
        G := Gs[i];
        //if i mod 20 eq 0 then
        printf "Finding quotients %o/%o: %o ", i, #Gs, (i le #descs select descs[i] else GroupToString(G));
        //end if;
        if i le #descs then
            A := AutomorphismGroup(G);
            Ambient, inj, proj := Holomorph(G, A);
            auts[i] := <A, Ambient, inj, proj>;
        else
            // We've stored this information when we found this as a quotient
            A, Ambient, inj, proj := Explode(auts[i]);
        end if;
        if IsAbelian(G) then
            print "abelian";
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
            print k, ":", (k le #descs select descs[k] else GroupToString(Gs[k])), "of size", #Gs[k], "with kernel of size", #Ns[j];
        end if;
        i +:= 1;
    end while;

    // Preparations are done: now we iteratively build up the list of image generators
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
                    // Looking for triples x,y,z up to automorphism that generate G and with (using additive notation now)
                    // 2x + 4y = 0
                    // Break it up by cases based on the rank
                    rank := FrattiniQuotientRank(G);
                    if rank eq 1 then
                        // We need at least one of x,y,z to generate on its own
                        // It can't be x, by the relation.
                        // In the case that z generates, we can set it to 1 using the aut group, removing the freedom
                        // Now y can be anything, and x=-2y or N/2-2y
                        // If z does not generate, y must, so we can set y to 1, x to -2 or N/2-2, and z to any non-unit.
                        g := CyclicGenerator(G);
                        N := #G;
                        N2 := N div 2;
                        for y in G do
                            Append(~gens[i], <g, y, y^-2>);
                            Append(~gens[i], <g, y, y^(N2-2)>);
                        end for;
                        for k in [0..N2-1] do
                            Append(~gens[i], <g^(2*k), g, g^-2>);
                            Append(~gens[i], <g^(2*k), g, g^(N2-2)>);
                        end for;
                        print "RANK 1", N, gens[i];
                    elif rank eq 2 then
                        // Z/2^a x Z/2^b with a <= b, A=2^a, B=2^b
                        basis, AB := PrimaryAbelianBasis(G);
                        A, B := Explode(AB);
                        r, s := Explode(basis);
                        twos := [Identity(G), r^(A div 2), s^(B div 2), r^(A div 2) * s^(B div 2)];
                        // Up to automorphism, there's a unique element of maximal order and it must be among our generators
                        // With the notation above, this is s.
                        // If z maps to s and y generates the rest, then we can act diagonally to send y to r or y to r*s.
                        for t in twos do
                            Append(~gens[i], <t*r^-2, r, s>);
                            Append(~gens[i], <t*r^-2, r*s, s>);
                        end for;
                        // If y maps to s but z does not map to a generator of order B, we must have that a < b
                        if A lt B then
                            for t in twos do
                                Append(~gens[i], <t*s^-2, s, r>);
                            end for;
                        end if;
                        if A eq 2 then
                            // In this case x can actually contribute to generating the group (e.g. if x is the first generator and y is the second, we're done)
                            // Now we can send z to s, y to a power of s and x to r*y^-2
                            // But some of these will be related by automorphisms, so we only need one value of x rather than the normal 2
                            for k in [0..B-1] do
                                y := s^k;
                                Append(~gens[i], <r*y^-2, y, s>);
                                // Not needed: Append(~gens[i], <r*y^-2*twos[3], y, s>);
                            end for;
                            // Or we can send y to s, z to a non-unit power of s, and x to r*s^-2
                            for k in [0..(B div 2)-1] do
                                z := s^(2*k);
                                Append(~gens[i], <r*s^-2, s, z>);
                                // Not needed: Append(~gens[i], <r*s^-2*twos[3], y, z>);
                            end for;
                        end if;
                        print "RANK 2", A, B, gens[i];
                    elif rank eq 3 then
                        basis, ABC := PrimaryAbelianBasis(G);
                        A, B, C := Explode(ABC);
                        // Now x must contribute to generating the group, so we must have A = 2
                        if A eq 2 then
                            u, r, s := Explode(basis);
                            // Have to handle the cases with more automorphisms first
                            if B eq 2 then
                                if C eq 2 then
                                    gens[i] := [<u, r, s>];
                                else
                                    // x has to help with u or r, and can take cases based on whether z maps to an element of maximal order or not
                                    gens[i] := [<u, r, s>, <u, s, r>];
                                end if;
                            else
                                twos := [u];
                                for pair in [<r,B>, <s,C>] do
                                    twos cat:= [t * pair[1]^(pair[2] div 2) : t in twos];
                                end for;
                                // Still a unique element of maximal order
                                // If z maps to such an element, we can take it to be s, and y can be conjugated to r
                                for t in twos do
                                    Append(~gens[i], <t*r^-2, r, s>);
                                end for;
                                // If not, we must have B < C, and then we can map y to s and z to r instead
                                for t in twos do
                                    Append(~gens[i], <t*s^-2, s, r>);
                                end for;
                            end if;
                        end if;
                        print "RANK 3", A, B, C, gens[i];
                    end if;
                else
                    // Not abelian
                    j, Qproj := Explode(quotients[i]);
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
                    print "Lifting", i, "from", j, ":", #translators, "aut-cosets,", #N, "kernel elements";
                    // Create the restriction map from a subgroup of A to AQ.
                    // It is only defined on the stabilizer of N within A
                    one := Identity(G);
                    seen := {};
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
                                                Append(~gens[i], <x,y,z>);
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
                end if;
            end for; // i in gps
        end for; // hsh->gps in by_ordhsh[m]
    end for; // m in Sort(Keys(by_ordhsh))
    return [<descs[i], #gens[i]> : i in [1..#descs]];
end intrinsic;
