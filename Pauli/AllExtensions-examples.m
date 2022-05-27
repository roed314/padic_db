

// The output of AllTotallyRamifiedExtensions may return several (reduced) generating polynomials for isomorphic extensions.

/*   2 */   Attach("AllExtensions.m");
/*   3 */   z3 := pAdicRing(3,30);
/*   6 */   J := PossibleDiscriminants(z3,9);
J;
/*  12 */   R := AllRamificationPolygons(z3,9,11); R;
R;
/*  22 */   A := AllResidualPolynomials(z3,R[3],3);
/*  23 */   A;
/*  24 */   L := AllTotallyRamifiedExtensions(z3,R[3],A[4],1);
/*  25 */   L;
/*  33 */   z3x<x> := PolynomialRing(z3);
/*  34 */   [<f[1],f[2],f[1] eq f[2],IsIsomorphic(f[1],f[2])> : f in CartesianProduct(L,L)];

// So in general we need filtering.  We try hard to make this unnecessary though.

/////////////////////

Attach("AllExtensions.m");
z3:=pAdicRing(3,30);
_<x>:=PolynomialRing(z3);
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
