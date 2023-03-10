// Create a folder Roe/lf.todo/ with one file for each p-adic field (named for the label of the field, with contents the defining polynomial as a list of integers)
// Create folders Roe/lf.errors and Roe/lf.out for error output and results respectively
// Usage: ls lf.todo | parallel -j128 "magma -b label:={1} compute_columns.m > lf.errors/{1}"
// Output: Writes to lf.out/<label>, with a | separated line of computed columns

SetColumns(0);
AttachSpec("../spec");
print label;
p, n, c, num := Explode([StringToInteger(c) : c in Split(label, ".")]);
f, coeffs := Explode([* eval c : c in Split(Read("lf.todo/" * label), "|") *]);
//f, eispol := Explode(Split(Read("lf.todo/" * label), "|"));
//f := StringToInteger(f);
if n eq f then
    // Unramified, so we already know the output
    PrintFile("lf.out/" * label, Sprintf("%o|{}|{}|{}", label));
    exit;
end if;
// Taken from SuggestedPrecision
prec := Maximum([2, 2*c]);

while true do
    R := pAdicRing(p, prec);
    if f gt 1 then
        R<t> := UnramifiedExtension(R, f);
    end if;
    S<x> := PolynomialRing(R);
    pol := S!coeffs;
    try
        fac, prec_loss, ext_info := Factorization(pol : Extensions:=true);
        eispol := DefiningPolynomial(ext_info[1]`Extension);
    catch err
        prec := 2 * prec;
        continue;
    end try;
    break;
end while;
//y := x;
//eispol := eval eispol;

polygon := Vertices(RamificationPolygon(eispol));
// Remove the first vertex, since it has y-coordinate infinity
polygon := polygon[2..#polygon];
respols := ResidualPolynomials(eispol);
k := CoefficientRing(respols[1]);
AssignNames(~k, ["t"]);
associated_inertia := [LCM([Degree(he[1]) : he in Factorization(respol)]) : respol in respols];

// Now prepare for printing to file
polygon := "{" * Join([Sprintf("{%o,%o}", pair[1], pair[2]) : pair in polygon], ",") * "}";
respols := "{" * Join([Sprintf("\"%o\"", respol) : respol in respols], ",") * "}";
associated_inertia := "{" * Join([Sprint(c) : c in associated_inertia], ",") * "}";

PrintFile("lf.out/" * label, Sprintf("%o|%o|%o|%o", label, polygon, respols, associated_inertia));
exit;
