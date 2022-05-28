/*********************
This file implements three intrinsics:

* pAdicFiltrations takes as input a group and a prime p, and outputs the possible ramification filtrations as chains of subgroups (without any information about the break locations).  The first two terms are the inertia subgroup and the wild inertia subgroup.
* WildFiltrations does the same, but ignoring the tame and unramified part.  As a consequence, the input must be a p-group.
* pAdicFiltrationClasses finds orbits of filtrations under the action of the automorphism group.  It is significantly slower than the other two.

EXAMPLES:

> G := SmallGroup(32, 7);
> W := WildFiltrations(G, 2);
> #W;
35
> W[25];
[
GrpPC of order 16 = 2^4
PC-Relations:
$.1^2 = $.3,
$.2^2 = $.4,
$.3^2 = $.4,
$.2^$.1 = $.2 * $.4,

GrpPC of order 8 = 2^3
PC-Relations:
$.2^2 = $.3,

GrpPC of order 2
PC-Relations:
$.1^2 = Id($),

GrpPC of order 1
PC-Relations:
]

> F := pAdicFiltrations(G, 2);
> #F;
73
> [#H : H in F[25]]; // this is wildly ramified since the first two terms are the whole group
[ 32, 32, 16, 8, 2, 1 ]
> C4 := [x : x in F | #x[1] eq 8];
> #C4; // There are 8 filtrations with an unramified start of degree 4
8

> Hs := TransitiveSubs(G);
> #Hs; // there are two subgroups realizing G as a transitive group of minimal degree
2
> TransitiveGroupIdentification(CosetImage(G, Hs[1]));
16 8
> TransitiveGroupIdentification(CosetImage(G, Hs[1])); // Two siblings 8T16
16 8
> H := Hs[1];
> F := pAdicFiltrationClasses(G, H, 2);
> #F; The 73 filtrations above fall into 44 classes under the action of the automorphism group
44

*********************/

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
