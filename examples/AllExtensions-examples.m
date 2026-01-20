//////////////////////
// Examples for AllExtensions and IndicesOfInseperability

AttachSpec("spec");

//////////////////////

z3 := pAdicRing(3,30);
J := PossibleDiscriminants(z3,9);
J;
R := AllRamificationPolygons(z3,9,11); R;
R;
A := AllResidualPolynomials(z3,R[3],3);
A;
L := AllTotallyRamifiedExtensions(z3,R[3],A[4],1);
L;

/////////////////////

z3:=pAdicRing(3,30);
z3x<x>:=PolynomialRing(z3);
Rs := AllRamificationPolygons(z3,9,14); 
R := Rs[1];
R;
As := AllResidualPolynomials(z3,R,3); 
A := As[1];
A;
Ls := AllTotallyRamifiedExtensions(z3,R,A,1);
for L in Ls do
  A,ResidualPolynomials(L);
  //R,RamificationPolygon(L);
end for;

for L in Ls do
  IndicesOfInseperability(L);
  //R,RamificationPolygon(L);
end for;

AllIndicesOfInseperability(z3,9,14);
