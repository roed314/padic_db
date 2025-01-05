// Usage: ls /scratch/lf/poly_2_8/ | parallel -j100 --timeout 600 "magma -b suff:=2_8 fname:={0} ComputePolRedPadicRelative.m > /scratch/lf/err_2_8/{0}"
// Usage: ls /scratch/lf/poly_t4/ | parallel -j120 "magma -b suff:=t4 fname:={0} ComputePolRedPadicRelative.m > /scratch/lf/err_t4/{0}"

Attach("../Pauli/polredabs.m");
SetColumns(0);

function remove_whitespace(X)
    return Join(Split(Join(Split(X," "),""),"\n"),"");
end function;
function sprint(X)
    return remove_whitespace(Sprint(X));
end function;

infile := Sprintf("/scratch/lf/poly_%o/" * fname, suff);
outfile := Sprintf("/scratch/lf/out_%o/" * fname, suff);
//infile := Sprintf("/Users/roed/Downloads/poly_%o_%o/" * fname, p, n);
//outfile := Sprintf("/Users/roed/Downloads/out_%o_%o/" * fname, p, n);
polys := Split(Read(infile), "\n");

pieces := Split(fname, ".");
p := StringToInteger(pieces[1]);
f := StringToInteger(pieces[2]);
e := StringToInteger(pieces[3]);
n := e * f;
/*if pieces[4][1] eq "0" then
    // There is a bug in PolRedPadic for unramified extensions
    PrintFile(outfile, Sprintf("%o|%o|%o", polys[1], polys[1], polys[1]));
    quit;
end if;*/

Zx<x> := PolynomialRing(Integers());
Zt<t> := PolynomialRing(Integers());
Ztx<x> := PolynomialRing(Zt);

prec := 4*n; // We use the same prec for all polynomials, since polynomials that need extra precision tend to arise together.
for poly in polys do
    while true do
        try
            k0 := pAdicRing(p, prec);
            R0<x> := PolynomialRing(k0);
            if f eq 1 then
                k := k0;
            else
                k<t> := UnramifiedExtension(k0, f);
            end if;
            R<x> := PolynomialRing(k);
            relpoly := eval poly;
            relbundle := PolRedPadic(relpoly);
            relbest := Distinguished(relbundle);
            relbest := Ztx![Zt!Coefficients(u) : u in Coefficients(relbest)];
            B<b> := TotallyRamifiedExtension(k, relpoly);
            Bx<x> := PolynomialRing(B);
            y := Roots(Bx!DefiningPolynomial(k) - b, B)[1][1];
            abspoly := Zx!MinimalPolynomial(y,k0);
            bundle := PolRedPadic(R0!abspoly);
            best := Distinguished(bundle);
            break;
        catch err
            print poly, prec;
            prec *:= 2;
        end try;
    end while;
    PrintFile(outfile, Sprintf("%o|%o|%o|%o", poly, sprint(Zx!best), sprint(relbest), Join([sprint(Zx!g) : g in bundle], ", ")));
end for;
quit;

