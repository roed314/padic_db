intrinsic Indices(f::RngUPolElt) -> []
{This computes the indices of inseparability of the extension L/K,
where L is generated over K by a root of f, which is an Eisenstein
polynomial over O_K.}
   Rx:=Parent(f);
   R:=BaseRing(f);
   if Type(R) eq FldPad then
      R:=IntegerRing(R);
      Rx:=PolynomialRing(R);
      f:=Rx!f;
   end if;
   Fq:=ResidueClassField(R);
   p:=Characteristic(Fq);
   e:=Valuation(R!p);
   n:=Degree(f);
   k:=Valuation(n,p);
   coeffs:=Coefficients(f);
   vals:=[0^^(n+1)];
   for i in [1..n+1] do
      vals[i]:=Valuation(coeffs[i]);
   end for;
   indicestilde:=[0^^(k+1)];
   indices:=[0^^(k+1)];
   Z:=Integers();
   for j in [1..k] do
      indicestilde[k-j+1]:=Minimum({i+n*vals[i+1]-n:i in [1..n-1]|
         Valuation(i,p) le k-j});
      indices[k-j+1]:=Z!Minimum({indicestilde[k-j+1],indices[k-j+2]+n*e});
   end for;
   return indices;
end intrinsic;

intrinsic RamPolygon(indices::[],n::RngIntElt,p::RngIntElt) -> [], []
{This computes the ramification polygon of a totally ramified
extension L/K using the indices of inseparability of L/K.  One
must also specify the degree n of the extension and the residue
characteristic p.}
   k:=Valuation(n,p);
   error if k+1 ne #indices, "Inconsistent inputs.";
   corners:=[<-p^k,0>];
   edges:=[];
   a:=k;
   while a gt 0 do
      newslopes:={(indices[i+1]-indices[a+1])/(p^a-p^i):i in [0..a-1]};
      smallslope:=Minimum(newslopes);
      newedge:=[corners[#corners]];
      for i in [1..a] do
         if (indices[a-i+1]-indices[a+1])/(p^a-p^(a-i)) eq smallslope then
            Append(~newedge,<-p^(a-i),indices[a-i+1]>);
         end if;
      end for;
      Append(~edges,newedge);
      Append(~corners,newedge[#newedge]);
      a:=Valuation(corners[#corners][1],p);
   end while;
   if n ne p^k then
     corners:=Insert(corners,1,<-n,0>);
     edges:=Insert(edges,1,[<-n,0>,<-p^k,0>]);
   end if;
   return corners,edges;
end intrinsic;

intrinsic Residual(f::RngUPolElt,FqAz::RngUPol) -> []
{Computes the residual polynomials of the generic polynomial
corresponding to a family of p-adic fields.  The input "f" should
be the generic polynomial, given as a polynomial over Z_p.  The
input "FqAz" is a polynomial ring in one variable whose
coefficients are rational functions over F_p in the parameters
used to define the family.  The output will be an element of
FqAz.}
   FqA:=BaseRing(FqAz);
   p:=Characteristic(FqA);
   n:=Degree(f);
   k:=Valuation(n,p);
   OKax:=Parent(f);
   OKa:=BaseRing(OKax);
   OK:=BaseRing(OKa);
   params:=Rank(OKa);
   OKx:=PolynomialRing(OK);
   f1:=OKx!0;
   for i in [0..n] do
      f1:=f1+Evaluate(Coefficient(f,i),[(OKx!1)^^params])*OKx.1^i;
   end for;
   indices:=Indices(f1);
   _,edges:=RamPolygon(indices,n,p);
   if n ne p^k then
      Remove(~edges,1);
   end if;
   coeffs:=Coefficients(f);
   OK2:=ChangePrecision(OK,2);
   OK2a:=PolynomialRing(OK2,params);
   dpibar:=OK2a!(coeffs[1]);
   dpi:=OKa!dpibar;
   pi:=UniformizingElement(OK);
   d:=dpi div pi;
   dbar:=FqA!d;
   respolys:=[];
   for E in edges do
      rpoly:=FqA!0;
      r:=#E;
      for j in [1..r] do
         pj:=-E[j][1];
         aj:=(E[j][2]-1) div n;
         bj:=E[j][2]-n*aj;
         newcoeff:=(Binomial(bj,pj)*coeffs[bj+1]) div pi^(aj+1);
         rpoly:=rpoly+FqAz!(newcoeff)*(-dbar)^(-aj-1)*FqAz.1^(pj);
      end for;
   Append(~respolys,rpoly);
   end for;
   if n ne p^k then
      tamepoly:=(FqAz.1+1)^n-1;
      Insert(~respolys,1,tamepoly);
   end if;
   return respolys;
end intrinsic;
