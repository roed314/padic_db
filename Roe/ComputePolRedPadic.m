// Usage: ls /scratch/lf/poly_2_8/ | parallel -j100 --timeout 600 "magma -b p:=2 n:=8 fname:={0} ComputePolRedPadic> /scratch/lf/err_2_8/{0}"

Attach("../Pauli/polredabs.m");
SetColumns(0);

function remove_whitespace(X)
    return Join(Split(Join(Split(X," "),""),"\n"),"");
end function;
function sprint(X)
    return remove_whitespace(Sprint(X));
end function;

infile := Sprintf("/scratch/lf/poly_%o_%o/" * fname, p, n);
outfile := Sprintf("/scratch/lf/out_%o_%o/" * fname, p, n);
p := StringToInteger(p);
n := StringToInteger(n);
polys := Split(Read(infile), "\n");

pieces := Split(fname, ".");
if pieces[3][1] eq "0" then
    // There is a bug in PolRedPadic for unramified extensions
    PrintFile(outfile, Sprintf("%o|%o", poly, poly));
    quit;
end if;

Zx<x> := PolynomialRing(Integers());
prec := 4*n; // We use the same prec for all polynomials, since polynomials that need extra precision tend to arise together.
for poly in polys do
    while true do
        try
            R<x> := PolynomialRing(pAdicRing(p, prec));
            f := eval poly;
            bundle := PolRedPadic(f);
            break;
        catch err
            print poly, prec;
            prec *:= 2;
        end try;
    end while;
    PrintFile(outfile, Sprintf("%o|%o", poly, Join([sprint(Zx!g) : g in bundle], ", ")));
end for;
quit;
// x^8+x^4+x^3+x^2+1
