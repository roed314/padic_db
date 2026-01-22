///////////////////////////////

Z3 := pAdicRing(3,50);
Z3x<x> := PolynomialRing(Z3);

Lur<xi> := UnramifiedExtension(Z3,2);
Lury<y> := PolynomialRing(Lur);
lambda<xib> := ResidueClassField(Lur);
lambdaz<z> := PolynomialRing(lambda);

psi0 := y^9 + (6*xi + 12)*y^6 + 9*y + 24;
R0, rho0 := RamificationPoly(psi0);
"ramificaion polygon",LowerVertices(R0),"with slopes",LowerSlopes(R0),"and residual polynomials",ResidualPolys(psi0);

phi0 := PolRedPadic(psi0);
"Distinguished defining polynomial over Z3[x]/(x^2+2x+2)"; phi0;
phia := PolRedPadic(psi0,Z3);
"Distinguished defining polynomial over Z3"; String(phia:wherenu)," or "; phia;

psi1 := x^18+9*x^14+3*6*x^13+3*3*x^12+3;
R1, rho1 := RamificationPoly(psi0);
"ramificaion polygon",LowerVertices(R1),"with slopes",LowerSlopes(R1),"and residual polynomials",ResidualPolys(psi1);
phi1 := PolRedPadic(psi1);
"Distinguished defining polynomial over Z3"; phi1;


///////////////////////////////

Z3 := pAdicRing(3,30);
Z3x<x> := PolynomialRing(Z3);
Lur<xi> := UnramifiedExtension(Z3,3);
Lury<y> := PolynomialRing(Lur);
 
psi2 := y^6+3*y^5+3*y^4+9*6*y^3+3+9*3;
R2, rho2 := RamificationPoly(psi2);
"ramificaion polygon",LowerVertices(R2),"with slopes",LowerSlopes(R2),"and residual polynomials",ResidualPolys(psi2);

phi2 := PolRedPadic(psi);
"Distinguished defining polynomial over Z3[x]/(x^3+2x+1)", phi0;
phib := PolRedPadic(psi2,Z3);
"Distinguished defining polynomial over Z3"; String(phib:wherenu)," or ";phib;



