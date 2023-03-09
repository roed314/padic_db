// Create a folder Roe/lf.todo/ with one file for each p-adic field (named for the label of the field, with contents the defining polynomial as a list of integers)
// Create folders Roe/lf.errors and Roe/lf.out for error output and results respectively
// Usage: ls lf.todo | parallel -j128 "magma -b label:={1} compute_columns.m > lf.errors/{1}"
// Output: Writes to lf.out/<label>, with a | separated line of computed columns

SetColumns(0);
AttachSpec("../spec");
print label;
p, n, c, num := Explode([StringToInteger(c) : c in Split(label, ".")]);
f, coeffs := Explode([* eval c : c in Split(Read("lf.todo/" * label), "|") *]);
if n eq f then
    // Unramified, where the code below doesn't work
    PrintFile("lf.out/" * label, Sprintf("%o|{}|{}|{}", label));
    exit;
end if;
// Taken from SuggestedPrecision
prec := Maximum([2, 2*c]);

while true do
    R := pAdicRing(p, prec);
    S := PolynomialRing(R);
    pol := S!coeffs;

    try
        fac, prec_loss, t := Factorization(pol : Extensions:=true);
        eispol := DefiningPolynomial(t[1]`Extension);
    catch err
        prec := 2*prec;
        continue;
    end try;
    break;
end while;
AssignNames(~S, ["x"]);

polygon := Vertices(RamificationPolygon(eispol));
// Remove the first vertex, since it has y-coordinate infinity
polygon := polygon[2..#polygon];
respols := ResidualPolynomials(eispol);
k0X := Parent(respols[1]);
k0 := CoefficientRing(respols[1]);
for j in [2..#respols] do
    assert k0 eq CoefficientRing(respols[j]);
end for;
// Need to change the coefficient ring to one defined by a Conway polynomial
k<a> := ext<GF(p)|ConwayPolynomial(p, Degree(k0))>;
SetPowerPrinting(k, false);
rts := [pair[1] : pair in Roots(DefiningPolynomial(k0), k)];
for rt in rts do
    tok := hom<k0 -> k | rt>;
    kX<x>, tokX := ChangeRing(k0X, k, tok);
    PrintFile("lf.check/" * label, "{" * Join([Sprintf("\"%o\"", tokX(respol)) : respol in respols], ",") * "}");
end for;
exit;

associated_inertia := [LCM([Degree(he[1]) : he in Factorization(respol)]) : respol in respols];

// Now prepare for printing to file
polygon := "{" * Join([Sprintf("{%o,%o}", pair[1], pair[2]) : pair in polygon], ",") * "}";
respols := "{" * Join([Sprintf("\"%o\"", respol) : respol in respols], ",") * "}";
associated_inertia := "{" * Join([Sprint(c) : c in associated_inertia], ",") * "}";

PrintFile("lf.out/" * label, Sprintf("%o|%o|%o|%o", label, polygon, respols, associated_inertia));
exit;
