/////////////////////////////////
// our implemntations are based on different definitions
//
// AllExtensions.m is based on [PS] S. Pauli and B. Sinclair, Enumerating extensions of (pi)-adic fields with
given invariants, Int. J. Number Theory 13 (2017)
// PolRedPadic is based on [GJK+] J. Guardia, J.W. Jones, K. Keating, S. Pauli, D.P. Roberts, and
D. Roe, Distinguished defining polynomials for extensions of p-adic fields

AttachSpec("spec");

Z3 := pAdicRing(3,30);
K := UnramifiedExtension(Z3,2);
"K is unramified extension of Z3 of degree 2";
kappa<a> := ResidueClassField(K);
kappaz<z> := PolynomialRing(kappa);
J := PossibleDiscriminants(K,18);
R := AllRamificationPolygons(K,18,16);
"Choose ramification polygon";
R[2];
A := AllResidualPolynomials(K,R[3],3);
"Choose restidual polynomials";
A[2];
psis := AllTotallyRamifiedExtensions(k,R[3],A[2],1);

psi := psis[10];
"Considering the Eisenstein polynomial",psi;

"ramification polygons differ in order of segments";
"RamificationPolygon from AllExtensions: ramification polygon";
R1 := RamificationPolygon(psi);  
"Vertices", LowerVertices(R1);
"Slopes",LowerSlopes(R1);
"RamificationPoly from PolRedPadic/rampoly.m: ramification polygon and ramification polynomial";
R2, rho := RamificationPoly(psi); 
"Vertices", LowerVertices(R2);
"Slopes",LowerSlopes(R2);
"because the order of segments differs, the order of the residual polynomials differ";
"ResidualPolynomials from AllExtensions";
A1 := ResidualPolynomials(psi);        
A1;
"ResidualPolys from PolRedPadic/rampoly.m";
A2 := ResidualPolys(f);                 
A2;
"Residual from GenResPack.mi see GenResExamples";

"RamificationPoly and ResidualPolys also work for Oystein polynomials";
phi, nu, alpha := OysteinPoly(psi,Z3);
"Considering the Oystein polynomial",String(phi:wherenu);
"RamificationPoly from PolRedPadic/rampoly.m";
Rphi, rhophi := RamificationPoly(phi,nu,alpha);
"Vertices", LowerVertices(Rphi);
"Slopes",LowerSlopes(Rphi);
"ResidualPolys from PolRedPadic/rampoly.m";
Aphi := ResidualPolys(phi, nu, alpha);
Aphi;


