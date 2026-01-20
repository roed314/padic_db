//////////////////////////////////////////////
// Functions for computing Oystein polynomials
// see Guardia, Jones, Keating, Pauli, and Roe, Distinguished Defining Polynomials for Extensions of p-Adic Fields, 2025
//
// By Sebastian Pauli, December 2025
//

"Loading oystein.m";

declare verbose Oystein, 6;
//////////////////////////////////////
// Miscillaneous

intrinsic ResidueFactor(phi) -> .
{The irreducible factor of phi over residue class field.
0 if phi is not a power of an irreducible polynomial over residue class field}
  K := CoefficientRing(phi);
  RK, KtoRK := ResidueClassField(K);
  Rphi := Polynomial(RK,phi);
  if Rphi eq 0 then return 0; end if;
  Rfact := Factorization(Rphi);
  if #Rfact ne 1 then return 0; end if;
  nu := Rfact[1][1];
  return Polynomial(K,nu);
end intrinsic;

intrinsic IsMonomial(f::RngUPolElt) -> .
{true if the polynomial f is a monomial}
  R := CoefficientRing(f);
  mono := 1 eq &+[1: a in Eltseq(f) | a ne 0];
  if not mono then 
    return false;
  elif R eq BaseRing(R) then
    return mono;
  else 
    coeff := &+[a : a in Eltseq(f)];
    ret := 1 eq &+[1: a in Eltseq(coeff) | a ne 0];
    return ret;
  end if;
end intrinsic;

///////////////////////////////////////
// Expansions

intrinsic Expansion(f::RngUPolElt,nu::RngUPolElt) -> .
{   
  The coefficients of the nu-expansion of f as a list.
}      
  expansion := [];
  l := 0; 
  repeat  
    l +:= 1;
    a := f mod nu;
    Append(~expansion,a);
    f := (f-a) div nu;
  until f eq 0;
  return expansion;
end intrinsic;
    
intrinsic Contraction(L::SeqEnum,nu::RngUPolElt) -> .
{ 
  Contraction(Expansion(f,nu),nu) = f
} 
  return &+[L[i]*nu^(i-1) : i in [1..#L]];
end intrinsic;

intrinsic Expansion2(f::RngUPolElt,nu::RngUPolElt:limit:=0) -> .
{
  The nu-expansion of f such that its coefficients are given as p expansions and the nu-expansion of f.
}
  if limit eq 0 then limit := Precision(CoefficientRing(f)); end if;

  Zx<x> := PolynomialRing(Integers()); 

  nuexp := Expansion(f,nu);
  p := Prime(CoefficientRing(f));

  if Degree(nu) gt 1 then
    expansion := [Zx!a : a in nuexp];
  else
    expansion := [Zx!Eltseq(ConstantCoefficient(a)): a in nuexp];
  end if;

  expexp := [];
// TODO which representatives are we using over an unramified extension ? 
  for g in expansion do
    h := g;
    gel := [];
    c := 0;
    while h ne 0 and c le limit do 
      Append(~gel,h mod p);
      h := h div p;
      c := c+1;
    end while;
    Append(~expexp,gel);
  end for;
  maxlen := Maximum([#gel:gel in expexp] cat [limit]); 
  for i in [1..#expexp] do 
    expexp[i] := expexp[i] cat [0:k in [1..maxlen-#expexp[i]+1]];
  end for;
  return expexp, nuexp;
end intrinsic;

intrinsic Contraction2(L::SeqEnum,nu::RngUPolElt) -> .
{
  Contraction2(Expansion2(f,nu),nu) = f
}
  Rx<x> := Parent(nu);
  R := CoefficientRing(nu);
  p := Prime(R);
  if R eq PrimeRing(R) then 
    return Rx!(&+[ &+[ p^(j-1)*L[i][j] : j in [1..#L[i]] ]*nu^(i-1) : i in [1..#L]]);
  elif Degree(nu) eq 1 then
    return Rx!([ &+[ p^(j-1)*Evaluate(L[i][j],R.1) : j in [1..#L[i]] ]: i in [1..#L]]);
  else
    error "not implemented yet";
  end if;  
end intrinsic;

/////////////////////////////////////////////
// Printing

intrinsic String(f::RngUPolElt:nu:=0,wherenu:=false,Latex:=false) -> .
{
  The nu expansion of f as a string.
}
  if (nu eq 0 and IsEisenstein(f)) or IsEisenstein(f) then
    nu := Parent(f).1;
    if wherenu then wherenu := false; end if;
  elif nu eq 0 then
    nu := ResidueFactor(f);
  end if;
  
  withnu := wherenu and not Latex; 
  if withnu then wherenu := false; end if;
 
  R := CoefficientRing(f);
  QR<a> := quo<R|UniformizingElement(R)^Precision(R)>;
  QRx<x> := PolynomialRing(QR);
 
  function tex(g)
    t := Sprintf("%O",QRx!g,"Latex");
    if IsMonomial(g) then 
      return t;
    else 
      return "(" cat t cat ")";
    end if;
  end function;
  nuexp := Expansion(f,nu);
 
  if wherenu then
    s := Sprintf("\\(\\nu^{%o} + ",#nuexp-1);
  elif withnu then
    s := Sprintf("nu^%o + ",#nuexp-1);
  elif Latex then
    s := Sprintf("\\(%o^{%o} + ",tex(nu),#nuexp-1);
  else 
    s := Sprintf("(%o)^%o + ",nu,#nuexp-1);
  end if;
  for j in [1..#nuexp-2] do
    i := #nuexp-j;
    if nuexp[i] ne 0 then
      if i eq 2 then
        if wherenu then
          s cat:= Sprintf("%o\\nu + ",tex(nuexp[2]));
        elif Latex then
          s cat:= Sprintf("%o %o + ",tex(nuexp[2]),tex(nu));
        elif withnu then
          s cat:= Sprintf("(%o)*nu + ",tex(nuexp[2]));
        else
          s cat:= Sprintf("(%o)*(%o) + ",QRx!nuexp[2],QRx!nu);
        end if;
      else
        if wherenu then
          s cat:= Sprintf("%o\\nu^{%o} + ",tex(nuexp[i]),i-1);
        elif withnu then  
          s cat:= Sprintf("(%o)*nu^%o + ",QRx!(nuexp[i]),i-1);
        elif Latex then
          s cat:= Sprintf("(%o) %o^{%o} + ",QRx!(nuexp[i]),QRx!nu,i-1);
        else
          s cat:= Sprintf("(%o)*(%o)^%o + ",QRx!nuexp[i],QRx!nu,i-1);
        end if;
      end if;
    end if;
  end for;
  if wherenu then
    s cat:= Sprintf("%o\\) where \\(\\nu = %O\\)",tex(nuexp[1]),QRx!nu,"Latex");
  elif withnu then
    s cat:= Sprintf("%o where nu = %o",QRx!(nuexp[1]), QRx!nu);
  elif Latex  then
    s cat:= Sprintf("%o\\)",tex(nuexp[1]));
  else
    s cat:= Sprintf("%o",QRx!nuexp[1]);
  end if;
  return s;
end intrinsic;

///////////////////////////////////////
// Finite fields and Conway Polynomials
//
// in PARI-GP: John Jones
// Magma : Sebastian Pauli 

/* is pol irreducible and primitive? */
function is_prim_pol(pol,p)
  m := Degree(pol);
  Fpz<z> := PolynomialRing(GF(p));
  pol := Fpz!pol;
  if not IsIrreducible(pol) then return false; end if;
  xoo := p^m-1;
  pps := [a[1]:a in Factorization(p^m-1)];
  for j in [1..#pps] do
    Fq<a> := ext<GF(p)|pol>;
    if (a^(xoo div pps[j])) eq 1 then
      return false;
    end if;
  end for;
  return true;
end function;


function unram_pol_jr(m,p)
  Zx<x> := PolynomialRing(Integers());
  pol := x^m;
  done := false;
  while not done do 
    j:=0; s:=1;
    while Coefficient(pol,j) eq (p-1)*s do
      pol -:= s*(p-1)*x^j;
      s := -s;
      j +:= 1;
    end while;
    pol +:= s*x^j;
    if is_prim_pol(pol,p) then
      done := true;
    end if;
  end while;
  return Zx![a mod p:a in Coefficients(pol)];
end function;

intrinsic ConwayOrJrPolynomial(K,n) -> .
{The Conway polynomial of degree n over K if it is known, a canonical irreducible polynomial of degree n over
K otherwise.}
  if K eq PrimeRing(K) then
    p := Prime(K);
    if ExistsConwayPolynomial(p,n) then
      return Polynomial(Integers(),ConwayPolynomial(p,n));
    else
      return unram_pol_jr(n,p);
    end if;
  else
    RK := ResidueClassField(K);
    return Polynomial(K!IrreduciblePolynomial(RK,n));
  end if;
end intrinsic;

intrinsic IsConwayOrJr(nu) -> .
{true if nu is a Conway polynomial or a Jones-Roberts irreducible polynomial}
   return ConwayOrJrPolynomial(Parent(nu),Degree(nu)) eq nu;
end intrinsic;

////////////////////////
// Oystein polynomials
//
// by Jordi Guardia and Sebastian Pauli
//

intrinsic IsOystein(phi::RngUPolElt[RngPad]:Conway:=false) -> .
{
True, if phi is in Oystein.  
}
  K := CoefficientRing(phi);
  if K ne PrimeRing(K) then return false; end if;
  nu := ResidueFactor(phi);
  if nu eq 0 then return false; end if;
  if not IsMonic(nu) then return false; end if;
  nuexp := Expansion(phi,Polynomial(K,nu));
  if Min([Valuation(a):a in Coefficients(nuexp[1])]) ne 1 then return false; end if;
  for i in [1..#nuexp-1] do
    if &or[Valuation(a) lt 1:a in Coefficients(nuexp[i])]  then return false; end if;
  end for;
  if not Conway then
    return true;
  else
    return IsConwayOrJr(nu);
  end if;
end intrinsic;

function oystein_poly_om(phi)
    Z := Integers();
    Zp := CoefficientRing(phi);
    Qp := FieldOfFractions(Zp);
    Qpx<x> := PolynomialRing(Qp);
    p := Prime(Zp);
    QX<X>:=PolynomialRing(Rationals());
    ZX<X> := PolynomialRing(Z);
    phiZ := ZX!phi;
    phiQ := QX!phi;

    // OM
    K<a>:=NumberField(phiZ);
    Montes(K,p);
    if #K`PrimeIdeals[p] ne 1 then error "OysteinPoly: polynomial must be irreducible over pAdicRing(",p,")"; end if;
    P := K`PrimeIdeals[p,1];
    D := Valuation(Discriminant(K),p) - 2*(K`LocalIndex[p]);
    piK := P`LocalGenerator;
    pi := QX!Eltseq(piK);    
    kp := ResidueField(P);
    e := P`e; // ramification index
    f := P`f; // inertia degree
    // D = f*(J+e-1)
    newprec := Max(Precision(Zp),2*((D div f)-e+1));

    Zpp := ChangePrecision(Zp,newprec);
    phip := ChangePrecision(Polynomial(Zpp,phi),newprec);  
    if f eq 1 then
      psi := CharacteristicPoly(phip,pi);
      L := TotallyRamifiedExtension(Zpp,psi); 
      return psi, X, L.1;
    end if;
    nu := ConwayOrJrPolynomial(Zpp,f);
    nubar := PolynomialRing(kp)!nu;
    ro := Roots(nubar);
    roK := [Lift(r[1],P): r in ro];
    polgammas:= [QX!Eltseq(gamma): gamma in roK];
    // Newton lift a root of nu to a root of nu(x)-pi
    dnu := Derivative(nu);
    gamma := polgammas[1];
    den := QX!Evaluate(dnu, gamma);
    _, invden := XGCD(den, phiZ);
    alpha := (gamma-(Evaluate(nu, gamma)-pi)*invden) mod phiQ;
    Phi := CharacteristicPoly(phip,alpha);    
    U := UnramifiedExtension(Zpp,nu);
    Uy<y> := PolynomialRing(U);
    RU := ResidueClassField(U);
    RUz<z> := PolynomialRing(RU);
    Phy := Uy!Phi;
    psi0 := z-RU.1;
    nupsi0 := nu div psi0;
    psifactors := HenselLift(Phy,[(Uy!psi0)^e,(Uy!nupsi0)^e]);
    psi := psifactors[1];
    L := TotallyRamifiedExtension(U,ChangePrecision(Evaluate(psi,y+U.1),newprec));
    alpha := Roots(Polynomial(L,psi))[1][1];
    nu := Polynomial(Zpp,nu);
    return Phi,nu,alpha;
end function;

intrinsic OysteinPoly(L::RngPad,K::RngPad) -> .
{
  A defining Oystein polynomial phi in K[x] of L along with 
  the polynomial nu generating the unramified subextensions of 
  L/K and gamma with phi(gamma) = 0.
}
vprintf Oystein,1: "OysteinPoly: Defining %o over %o\n",L,K;
  Lt<t> := PolynomialRing(L); 
  pi := UniformizingElement(L);

  if InertiaDegree(L,K) eq Degree(L,K) then
    nu := Polynomial(K,ConwayOrJrPolynomial(PrimeRing(K),InertiaDegree(L,K)));
    alpha := Roots(Lt!nu-pi:Max:=1)[1][1];
    phi := nu;
  elif RamificationIndex(L,K) eq Degree(L,K) then
    // L is totally ramified over K
    phi := DefiningPolynomial(L,K);
    nu := PolynomialRing(K).1;
    gamma := L.1;
  else
    nu := Polynomial(K,ConwayOrJrPolynomial(PrimeRing(K),InertiaDegree(L,K)));
    gamma := Roots(Lt!nu-pi:Max:=1)[1][1];
    phi := CharacteristicPolynomial(gamma,K);
  end if;
  vprintf Oystein,2: "OysteinPoly is %o\n",String(phi:wherenu);
  vprintf Oystein,4: "OysteinPoly check %o\n",HasRoot(Lt!phi);
  return phi, nu, gamma;
end intrinsic;

intrinsic OysteinPoly(f::RngUPolElt[RngPad],K::RngPad) -> .
{
Given an irreducible f, return a nu-Oystein polynomial g in such that K[x]/(g) is isomorphic to 
K[x]/(f) along with a defining polynomial nu of the unreamified part of the extension 
and a root of g in K[x]/(g)
}
  L := CoefficientRing(f);
  RL := ResidueClassField(L);
  f := ChangePrecision(f,Precision(L));
  vprintf Oystein,1: "OysteinPoly: %o over %o\n",f,K;

  if not IsMonic(f) then error "OysteinPoly: polynomial",f,"must be monic"; end if;
  if IsEisenstein(f) then
    L := TotallyRamifiedExtension(L,f);
    return OysteinPoly(L,K);
  elif IsIrreducible(Polynomial(RL,f)) then
    U := UnramifiedExtension(L,f);
    return OysteinPoly(U,K);
  elif L eq PrimeRing(L) then
    return oystein_poly_om(f);
  else
    factors, _ ,Cs := Factorization(f:Certificates:=true); 
    if #factors gt 1 then error "OysteinPoly: polynomial must be irreducible"; end if;
    U := UnramifiedExtension(L,ConwayOrJrPolynomial(PrimeRing(L),Cs[1]`F));
    factors, _ ,Ls := Factorization(Polynomial(U,f):Extensions:=true); 
    return OysteinPoly(Ls[1]`Extension,K);
  end if;
end intrinsic;

intrinsic OysteinPoly(f::RngUPolElt[RngPad]) -> .
{
  Given f in K[x] irreducible, return a defining Oystein polynomial phi of L=K[x]/(f) 
  along with  the polynomial nu generating the unramified subextensions of 
  L/K and gamma with phi(gamma) = 0.
}
  K := CoefficientRing(f);
  return OysteinPoly(f,K);
end intrinsic;



