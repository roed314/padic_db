// Usage: ls /scratch/lf/poly_2_8/ | parallel -j100 --timeout 600 "magma -b suff:=2_8 fname:={0} ComputePolRedPadic.m > /scratch/lf/err_2_8/{0}"
// Usage: ls /scratch/lf/poly_db/ | parallel -j100 --timeout 600 "magma -b suff:=db fname:={0} ComputePolRedPadic.m > /scratch/lf/err_db/{0}"

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
n := StringToInteger(pieces[2]) * StringToInteger(pieces[3]);
if pieces[4][1] eq "0" then
    // There is a bug in PolRedPadic for unramified extensions
    PrintFile(outfile, Sprintf("%o|%o|%o", polys[1], polys[1], polys[1]));
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
            best := Distinguished(bundle);
            break;
        catch err
            print poly, prec;
            prec *:= 2;
        end try;
    end while;
    PrintFile(outfile, Sprintf("%o|%o|%o", poly, sprint(Zx!best), Join([sprint(Zx!g) : g in bundle], ", ")));
end for;
quit;

