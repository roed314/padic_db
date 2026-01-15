
//////////////////////////////////////////////////////////////////////////
// compute the characteristic polynomials using newton relations
// see Henri Cohen, A course in Computational Number Theory, page 161
//
// By Sebastian Pauli, 2001 and 2025
//

"Loading charpoly.m";

traces_from_poly := function(Phi)
 
  n := Degree(Phi);
  R := CoefficientRing(Parent(Phi));
  
  t := function(i)
    if i lt 0 then
      return 0;
    else
      return Coefficient(Phi,i);
    end if;
  end function;
  
  S := [];

  for k := 1 to n do
    Sk := -k*t(n-k);
    for i := 1 to Minimum(n,k-1) do
      Sk := Sk - t(n-i)*S[k-i];
    end for;
    Append(~S,Sk);
  end for;
  return S;
end function;


poly_from_traces := function(S)

  n := #S;
  Tn := [];
  
  for k := 1 to n do
    Tnk := -S[k];
    for i := 1 to k-1 do
      Tnk := Tnk - Tn[i]*S[k-i];
    end for;
    Tnk := Tnk div k;
    Append(~Tn,Tnk);
  end for;
  Tn := Reverse(Tn);
  Append(~Tn,Parent(S[1])!1);
  return PolynomialRing(Parent(S[1]))!Tn;
  
end function;

chi := function(Phi, traces_Phi, theta)

  n := Degree(Phi);
  assert n gt 0;

  tmp_PR := PolynomialRing(Universe(traces_Phi));
  tmp_theta := tmp_PR!theta;
  tmp_Phi := tmp_PR!Phi;

  S := [];

  pow := 1;
  for i := 1 to n do
    pow := pow*tmp_theta mod tmp_Phi;
    Si := 0;
    Si := n*Coefficient(pow,0);
    for j := 1 to Degree(Phi)-1 do
      if Coefficient(pow,j) ne 0 then
        Si +:= traces_Phi[j]*Coefficient(pow,j);
      end if;
    end for;
    Append(~S,Si);
  end for;
  return Parent(Phi)!poly_from_traces(S);
end function;

chi_with_den := function(Phi, traces_Phi, theta,denval)
    vprint RoundFour,5:" -- char poly";
    vprint RoundFour,6:" -- of",theta;
    vprint RoundFour,6:" -- with denominator pi ^",denval;
    
    Zp := CoefficientRing(Phi);
    Zpy<y> := PolynomialRing(Zp);
    res := chi(Phi, traces_Phi, theta);
    pi := UniformizingElement(Zp);
    res := Evaluate(res, pi^denval * y);
    res div:= LeadingCoefficient(res);
    return res;
end function;

intrinsic CharacteristicPoly(Phi,theta) -> .
{Characteristic polynomial of the element represented by theta in the extension defined by Phi}

  Zpin := CoefficientRing(Phi);
  Qpin := FieldOfFractions(Zpin);
  thetain := Polynomial(Qpin,theta); 
  p := Prime(Zpin);
  denval := Valuation(Denominator(thetain));
  prec := Precision(Zpin)+denval*Degree(Phi)+(Ceiling(Degree(Phi)/(p-1)));
  Zp := ChangePrecision(Zpin,prec);
  Psi := ChangePrecision(Polynomial(Zp,Phi),prec);
  pi := UniformizingElement(Zp);
  Qp := FieldOfFractions(Zp);
  thetap := ChangePrecision(Polynomial(Qp,theta),prec);
  vals := [Valuation(t):t in Eltseq(thetap)];
  minval := Minimum(vals);
  traces_Psi := traces_from_poly(Psi);
  if minval ge 0 then
    chitheta := chi(Psi,traces_Psi,thetap);
  else
    denval := -minval;
    tau := Polynomial(Zp,(thetap*pi^denval));
   chitheta := chi_with_den(Psi,traces_Psi,tau,denval);
  end if;
  return ChangePrecision(Polynomial(Zpin,chitheta),Precision(Zpin));
end intrinsic;

