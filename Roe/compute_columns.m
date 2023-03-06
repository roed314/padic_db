// Create a folder Roe/lf.todo/ with one file for each p-adic field (named for the label of the field, with contents the defining polynomial as a list of integers)
// Create folders Roe/lf.errors and Roe/lf.out for error output and results respectively
// Usage: ls lf.todo | parallel -j128 "magma -b label:={1} compute_columns.m > lf.errors/{1}"
// Output: Writes to lf.out/<label>, with a | separated line of computed columns

AttachSpec("../spec");
p, n, c, num := Explode([StringToInteger(c) : c in Split(label, ".")]);
f, coeffs := Explode([* eval c : c in Split(Read("lf.todo/" * label), "|") *]);
// Taken from Suggested precision
prec := Maximum([4, 4*c]);
R := pAdicRing(p, prec);
S := PolynomialRing(R);
AssignNames(~S, ["x"]);
f := S!coeffs;

if f eq 1 then
    eispol := f;
else
    fac, prec_loss, t := Factorization(f : Extensions:=true);
    eispol := DefiningPolynomial(t[1]`Extension);
end if;

polygon := Vertices(RamificationPolygon(eispol));
// Remove the first vertex, since it has y-coordinate infinity
polygon := polygon[2..#polygon];
respols := ResidualPolynomials(eispol);
k := CoefficientRing(respols[1]);
AssignNames(~k, ["a"]);
associated_inertia := [LCM([Degree(he[1]) : he in Factorization(respol)]) : respol in respols];

// Now prepare for printing to file
polygon := "{" * Join([Sprintf("{%o,%o}", pair[1], pair[2]) : pair in polygon], ",") * "}";
respols := "{" * Join([Sprintf("\"%o\"", respol) : respol in respols], ",") * "}";
associated_inertia := "{" * Join([Sprint(c) : c in associated_inertia], ",") * "}";

PrintFile("lf.out/" * label, Sprintf("%o|%o|%o|%o", label, polygon, respols, associated_inertia));
exit;
