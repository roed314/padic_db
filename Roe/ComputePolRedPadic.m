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
polys := Split(Read(infile), "\n");
R<x> := PolynomialRing(pAdicRing(p, 30));
Zx := PolynomialRing(Integers());
AssignNames(~Zx, ["x"]);
for poly in polys do
    f := eval poly;
    bundle := PolRedPadic(f);
    PrintFile(outfile, Sprintf("%o|%o", poly, Join([sprint(Zx!g) : g in bundle], ", ")));
end for;
quit;
