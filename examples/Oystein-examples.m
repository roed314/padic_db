///////////////////////////////////////////
// Examples for OysteinPoly and PolRedPadic

AttachSpec("spec");

////////////////////////////

Z3 := pAdicRing(3,50);
Z3x<x> := PolynomialRing(Z3);
nu := Z3x!ConwayPolynomial(3,4);

f := nu^9 + 9*nu^7+3*x*nu^4 + 3+99;
phi,nu,alpha := OysteinPoly(f); // A defining Oystein polynomial of K=Z3[x]/(f)
String(phi:wherenu);            // nu expansion of phi
Valuation(Evaluate(nu,alpha));  // nu(alpha) is uniformizer

psi := PolRedPadic(phi);        // The distinguished defining (Oystein) polynomial of K
String(psi:wherenu);

f1 := CharacteristicPoly(f,3*x^2+2*x); 
phi1, nu1, alpha1 := OysteinPoly(f1);  // A defining Oystein polynomial of Z3[x]/(f1)
String(phi1:wherenu);                  // nu expansion of phi   
Valuation(Evaluate(nu1,alpha1));       // nu1(alpha1) is uniformizer

psi1 := PolRedPadic(phi1);
String(psi1:wherenu);
psi eq psi1;             


/////////////////////////////////////////////////////////

Z5 := pAdicRing(5,20);
Z5x<x> := PolynomialRing(Z5);
nu := Z5x!ConwayPolynomial(5,3);

f := nu^25 + 50*nu^7+5*x*nu^4 + 55;
phi,nu,alpha := OysteinPoly(f);
String(phi:wherenu);
Valuation(Evaluate(nu,alpha));

f1 := CharacteristicPoly(f,25*x^2+5*x+2);
phi1, nu1, alpha1 := OysteinPoly(f1);
String(phi:wherenu);
Valuation(Evaluate(nu1,alpha1));

////////////////////////////////////////////

Z2 := pAdicRing(2,1000);
Z2x<x> := PolynomialRing(Z2);
nu := Z2x!ConwayPolynomial(2,5);
f := nu^8 + 8*nu^7+4*x*nu^4 + 90;

phi, nu, alpha := OysteinPoly(f);
String(phi:wherenu);
Valuation(Evaluate(nu,alpha));

psi := PolRedPadic(phi);
String(psi:wherenu);

f1 := CharacteristicPoly(f,x^2+4*x+1);
phi1, nu1, alpha1 := OysteinPoly(f1);
String(phi1:wherenu);
Valuation(Evaluate(nu1,alpha1));

psi1 := PolRedPadic(phi1);
psi1 eq psi;

////////////////////////////////

Z:=Integers();
QX<X>:=PolynomialRing(Rationals());
ZX<X> := PolynomialRing(Integers());
nu := X^2 + 2*X + 2;
f1 := nu^9 + 3*nu^7 + 3+99;
phi:=f1;
k<a>:=NumberField(phi);
ff := ZX!CharacteristicPolynomial(3*k.1+9);

Z3 := pAdicRing(3,20);
Z3x<x> := PolynomialRing(Z3);
f := Z3x!ff;

phi, nu, alpha := OysteinPoly(f);
String(phi:wherenu);
alpha in [r[1]: r in Roots(Polynomial(Parent(alpha),phi))];
Valuation(Evaluate(nu,alpha));

psi := PolRedPadic(phi);
String(psi:wherenu);


