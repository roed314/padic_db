// Usage: parallel -j112 --timeout 86400 -a /scratch/lf/todo_rp.txt "magma -b fname:={0} ComputeAllExtensions.m > /scratch/lf/err_rp/{0}"

Attach("../Pauli/AllExtensions.m");
Attach("../Pauli/polredabs.m");
SetColumns(0);

function remove_whitespace(X)
    return Join(Split(Join(Split(X," "),""),"\n"),"");
end function;
function sprint(X)
    return remove_whitespace(Sprint(X));
end function;

infile := Sprintf("/scratch/lf/poly_rp/" * fname);
outfile := Sprintf("/scratch/lf/out_rp/" * fname);
pieces := Split(fname, ".");
p := StringToInteger(pieces[1]);
f := StringToInteger(pieces[2]);
e := StringToInteger(pieces[3]);
c := StringToInteger(pieces[4]);
RP := eval Read(infile);

Zx<x> := PolynomialRing(Integers());
prec := 64;
while true do
    try
        k0 := pAdicRing(p, prec);
        if f eq 1 then
            k := k0;
        else
            k := UnramifiedExtension(k0, f);
        end if;
        RK, res := ResidueClassField(k);
        RKz<z> := PolynomialRing(RK);
        phi0s := DistinctConstantCoefficients(k, e);
        for phi0 in phi0s do
            A := AllResidualPolynomials(k, RP, phi0);
            for a in A do
                for relpoly in AllTotallyRamifiedExtensions(k, RP, a, res(phi0/Prime(k))) do
                    if f eq 1 then
                        abspoly := Zx!relpoly;
                    else
                        B<b> := TotallyRamifiedExtension(k, relpoly);
                        Bx<x> := PolynomialRing(B);
                        y := Roots(Bx!DefiningPolynomial(k) - b, B)[1][1];
                        abspoly := Zx!Minpoly(y,k0);
                    end if;
                    R<x> := PolynomialRing(k0);
                    bundle := PolRedPadic(R!abspoly);
                    best := Distinguished(bundle);
                    PrintFile(outfile, Sprintf("%o|%o|%o", sprint(abspoly), sprint(Zx!best), Join([sprint(Zx!g) : g in bundle], ", ")));
                end for;
            end for;
        end for;
        break;
    catch err
        System("rm " * outfile);
        prec *:= 2;
    end try;
end while;
quit;

