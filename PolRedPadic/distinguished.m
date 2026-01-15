//////////////////////////////////////////////////////////////////////
// Distinguished Defining Polynomials for Extensions of p-Adic Fields
//
// Main function is PolRedPadic
//
// see Guardia, Jones, Keating, Pauli, and Roe, Distinguished Defining Polynomials for Extensions of p-Adic Fields, 2025
//
// By Sebastian Pauli, December 2025
//


"Loading distinguished.m";

declare verbose Monge, 6;

function discrete_log(a)
  return Log(a);
end function;

//////////////////////////////////////
// Distinguished Residual Polynomials

intrinsic PolynomialCompareLog(f,g) -> .
{Compare two polynomials over a local field of the same degree by comparing the discrete logarithms of the coefficients and using lexicographic ordering with leading coefficient first.}
  for i := Max(Degree(f),Degree(g)) to 0 by -1 do
    a := Coefficient(f,i); b := Coefficient(g,i);
    if a eq 0 and b ne 0 then 
      return -1;
    elif b eq 0 and a ne 0 then 
      return 1;
    elif a ne 0 and b ne 0 and a ne b then
      c := discrete_log(a); d := discrete_log(b);
      return c-d;
    end if;
  end for;
  return 0;
end intrinsic;

intrinsic ResidialPolyCompare(A,B) -> .
{
Return 1 if A>B, -1 if A<B, 0 otherwise
}
  if #A ne #B then error "ResidialPolyCompare: Lists of residual polynomials must be of the same length."; end if;
  for i in [1..#A] do
    c := PolynomialCompareLog(A[i],B[i]);
    if c ne 0 then return c; end if;    
  end for;  
  return 0;
end intrinsic;

intrinsic ResidualPolyDistinguished(phi::RngUPolElt:conjugates:=false,constant_first:=true) -> .
{
The distinguished (minimal) representative of the residual polynomial class of an Eisenstein polynomial phi
along with Eisenstein polynomials psi that yield the distinguished representative.
If conjugates is true, the conjugates under AutomorphismGroup(CoefficientRing(phi)) are also considered.
If constant_first is true, it is first assured that the dicrete logarithm of the constant coefficiant of
psi is minimal and only afterwards the residual polynomials are minimized.
}
      if not IsEisenstein(phi) then
        conjugates := true;
        _ , nu, alpha := OysteinPoly(phi); 
        pi := Evaluate(nu,alpha);
        phi := CharacteristicPolynomial(pi);
      end if;

      K := CoefficientRing(phi);
      p := Prime(K);
      piK := UniformizingElement(K);
      Kx<x> := PolynomialRing(K);
      e := Degree(phi);
      phi := ChangePrecision(phi,Precision(K));
      L<alpha> := TotallyRamifiedExtension(K,phi);
      LX<X> := PolynomialRing(L);

      // TODO
      rho:=Evaluate(LX!phi,X + alpha);
      rp := NewtonPolygon(rho);
      //rp, rho := ramification_poly_raw(phi,alpha);
      
      slopes := Reverse([-m:m in LowerSlopes(rp)]);
      vertices := Reverse(LowerVertices(rp));
      Fq, KtoFq := ResidueClassField(K);
      Fqz<z> := PolynomialRing(Fq);
      q := #Fq;
      xi := PrimitiveElement(Fq);

      // tamely ramified is easy, since the residual polynomial only depends on the degree.
      if Valuation(Degree(phi),p) eq 0 then
        return ResidualPolys(phi), [phi];
      end if;

      vprintf Monge, 2: "ResidualPolyDistinguished: conjugates %o\n",conjugates;  
      vprintf Monge, 4: "ResidualPolyDistinguished: %o with slopes %o\n",phi,slopes;  

      function residual_polynomial_distinguished_sub(phi:constant_first:=true);
      // 
          if constant_first then
            // nice constant term takes precendence over nice residual polynomial class representative
            phi0 := ConstantCoefficient(phi);
            phi01 := KtoFq(phi0 div piK);
            a := Log(phi01);
            d, s0, _ := Xgcd(e,q-1);
            k, b  := Quotrem(a,d);
            t0 := (-s0 * k) mod (q-1);
            Delta := (q-1) div d;
            x_base := [t0];
          else
            Delta := 1;
            x_base := [0];
          end if;

          //L<alpha> := TotallyRamifiedExtension(K,phi);
          //LX<X> := PolynomialRing(L);
          //slopes := Reverse([-m:m in LowerSlopes(rp)]);
          //vertices := Reverse(LowerVertices(rp));
          A := ResidualPolys(phi);
          g := 0;
          for i in [1..#slopes] do
            n := Degree(A[i]);
            m := slopes[i];
            t := Numerator(m);
            d := Denominator(m);
            g := g+(d-t)*n;
            for j := n to 0 by -1 do
              Aij := Coefficient(A[i],j);
              vprintf Monge, 6: "ResidualPolyDistinguished: m=%o %o-th monomial %o\n",slopes[i],j,Aij;
              if Aij ne 0 then
                a := discrete_log(Aij) mod (q-1);
                D := (Delta*((t-d)*j+g)) mod (q-1);
                if D ne 0 then 
                  b,s,_ := Xgcd(D,q-1);
                  new_Delta := Lcm(Delta,(q-1) div b);
                  minexp := q;
                  for xij in x_base do
                      J := a+xij*((t-d)*j+g);                    
                      r := J mod b;
                      k := J div b;
                      x := (xij-k*s*Delta) mod (q-1);
                      vprintf Monge, 5: "ResidualPolyDistinguished: m=%o solutions %o+n*%o\n",slopes[i],x,Delta;  
                      if r lt minexp then
                        minexp := r;
                        new_x_base := [x];
                      elif r eq minexp then
                        Append(~new_x_base,x);
                      end if;
                  end for;
                  Delta := new_Delta;
                  x_base := new_x_base;
                end if;
              end if; 
            end for;
          end for;
          return x_base, Delta;
       end function;

       function residual_polynomial_phis(thisphi,s_base,s_diff)
          vprintf Monge, 5: "ResidualPolyDistinguished:     final difference: %o\n",s_diff;  
          minphis := [];
          for sb in s_base do
            s := sb;
            repeat
              s := (s+s_diff) mod (q-1);
              deltaK := K!(xi^s);
              phidelta := Kx!([Coefficient(thisphi,i)*deltaK^(Degree(thisphi)-i) : i in [0..Degree(thisphi)]]);
              Include(~minphis,<ResidualPolys(phidelta),phidelta>);
            until s eq sb; 
          end for;
          vprintf Monge, 5: "ResidualPolyDistinguished:     #phis=%o\n",#minphis;  
          return minphis;
       end function;

       if not conjugates then
         base, delta := residual_polynomial_distinguished_sub(phi:constant_first:=constant_first);
         As := residual_polynomial_phis(phi,base,delta);
         phis := [a[2]: a in As | a[1] eq As[1][1]];
         //return As[1][1],phis;
       else
         As := [];
         gaut, maut := AutomorphismGroup(K);
         aut := [ maut(tau) : tau in gaut ];
         for tau in aut do
           vprintf Monge,3: "ResidualPolyDistinguished: %o->%o\n",KtoFq(K.1),KtoFq(tau(K.1));
           tauphi := Kx![tau(c): c in Coefficients(phi)];
           base, delta := residual_polynomial_distinguished_sub(tauphi);
           As cat:= residual_polynomial_phis(tauphi,base,delta);
         end for;           
         Sort(~As,func<x,y|ResidialPolyCompare(x[1],y[1])>);
       end if;
       // needed when constant first ?
       philogs := [<discrete_log(KtoFq(ConstantCoefficient(a[2])/piK)),a[2]>: a in As | a[1] eq As[1][1]];
       minlog := Min([a[1] : a in philogs]);
       phis := [a[2]:a in philogs|a[1] eq minlog];
       return As[1][1],phis;
end intrinsic;


//////////////////////////
//

intrinsic Distinguished(M::SetEnum[RngUPolElt[RngPad]]:nu:=0) -> .
{
  Given a set of reduced nu-Oystein polynomials, return the distinguished polynomial.
}
  vprintf Monge,1: "Distinguished: among %o polynomials\n",#M; 
  vprintf Monge,2: "Distinguished: among \n" cat (&cat[String(phi:wherenu) cat "\n":phi in M]);
  if #M eq 1 then return Rep(M); end if;

  Kx<x> := Parent(Rep(M));
  K := CoefficientRing(Rep(M));
  eisen := IsEisenstein(Rep(M));
  p := Prime(K);
  Fp := ResidueClassField(K);
  Fpz<z> := PolynomialRing(Fp);
  if eisen then
    nu := x;
  elif nu eq 0 then
    nu := ResidueFactor(Rep(M));
  end if;
  L := SetToSequence(M);
  function dcompare(f,g)
    fg := f-g;
    fgseq := Eltseq(fg);
    fgval := [Valuation(a): a in fgseq];
    fgmin := Min(fgval);
    i := Min([j: j in [1..#fgval] | fgval[j] eq fgmin]);
    fi := Coefficient(f,i-1);
    gi := Coefficient(g,i-1);
    fexp := Expansion2(f,nu);
    gexp := Expansion2(g,nu);
    for k in [1..#fexp[1]] do
      for i in [1..#fexp] do
          fik := fexp[i][k];
          gik := gexp[i][k];
          if fik ne gik then
             if fik eq 0 then
               return 1;
             elif gik eq 0 then 
               return -1;
             elif eisen then
               fik1 := Evaluate(fik,K.1);
               gik1 := Evaluate(gik,K.1);
               return Log(Fp!fik1)-Log(Fp!gik1);
             else
	       ret := PolynomialCompareLog(Fpz!fik,Fpz!gik);
               return ret;
             end if;
          end if;
       end for;
    end for;
    error "Distinguished: compare trouble";
    return 0;
  end function;

  Sort(~L, dcompare);
  return L[1];
end intrinsic;



///////////////////////////////////////////////////////////////
function pol_red_padic_sub(Phi,nu,alpha,psi01)
// Phi in K[x]
// nu in K[x] generates unramified subextension of L = K[x]/(Phi) = K(alpha)  
// Phi(alpha) = 0
// psi01 desired constant coefficient mod pi^2
        n := Degree(Phi);
        f := Degree(nu);
        e := n div f;
        
        Kx<x> := Parent(Phi);
        K := CoefficientRing(Kx);

        p := Prime(K);
        Zp := PrimeRing(K);
        RK, KtoRK := ResidueClassField(K);
        RKs, RKstoRK := UnitGroup(RK);
        RKz<z> := PolynomialRing(RK);

        L := Parent(alpha);

        Lt<t> := PolynomialRing(L);
        RL, LtoRL := ResidueClassField(L);
        Fp := BaseField(RL);
        RLs, RLstoRL := UnitGroup(RL);
        psi01R := LtoRL(psi01);
        xi := RL.1; // primitive element
        Pi := Evaluate(Polynomial(L,nu),alpha);        
        A_phi := ResidualPolys(Phi);
       
        function is_iso(S)
          FB := Basis(RL,RK);  
          FL := [Eltseq(Evaluate(S,b)):b in FB];
          FM := Matrix(FL);
          return Rank(FM) eq #FB; 
        end function;

        // for m high enough we can set coefficients to 0       
       
        rp, rho := RamificationPoly(Phi,nu,alpha);

        slopes := LowerSlopes(rp); // slopes := slopes[2..#slopes]; // remove infinite slope
        slopes := [s : s in slopes | Abs(s) lt Precision(PrimeRing(K))];
        vprintf Monge,1: "PolRedPadic: Ramification polygon %o with slopes %o\n", LowerVertices(rp), slopes;
	maxslope := Max(slopes); 
	easystart := Floor(maxslope)+2;
        Smax, PHImax := AdditiveResidualPoly(Phi,nu,alpha,easystart-1);
	easylimit := PHImax div e + 1;
	
// TODO: actually use precision
        vprintf Monge,1:"PolRedPadic: easy reduction starts with m=%o, using a preision of %o\n",easystart,easylimit;
 
        function easyreduce(phi)
          m := easystart;
          nuexp := Expansion2(phi,nu : limit := easylimit);
          repeat
            //def, wm := IsDefined(Wm,m);
            //if not def then wm := Wm[max_m]+m-max_m; end if;
	    wm := PHImax+m-easystart; 
	    i := wm mod e;
            k := wm div e;
            nuexp[i+1][k+1] := 0;
            vprintf Monge,6:"PolRedPadic: easy m=%o setting phi*_(%o,%o)=%o to 0\n",m,i,k,nuexp[i+1][k+1];
            // printf "PolRedPadic:   m = %o : setting phi*_(%o,%o) = %o to 0\n",m,i,k,nuexp[i+1][k+1];
            vprintf Monge,6:"PolRedPadic: easy m=%o still isomorphic %o\n",m,HasRoot(Lt!Contraction2(nuexp,nu));
            m := m+1;
          until k gt easylimit or k ge Precision(Zp);
          newphi := Contraction2(nuexp,nu);
          return newphi;
        end function;
       
        // reduction of constant coefficient mod pi^2
        vprint Monge,4:"PolRedPadic: m=0 reducing",String(easyreduce(Phi):wherenu);
        vprint Monge,1:"PolRedPadic: m=0 reduction with alpha->alpha+theta*nu(alpha)";
        nuexp2, nuexp := Expansion2(Phi,nu : limit := easylimit);
        phi0 := nuexp[1];
        phi0alpha := Evaluate(Lt!phi0,alpha); 
        nualpha := Evaluate(Lt!nu,alpha);
        eta := LtoRL((nualpha^e) div p);
        S1, r1 := AdditiveResidualPoly(Phi,nu,alpha,0);
        S1eta := eta*S1;
        vprintf Monge,1: "PolRedPadic: m=0 Phi(0)=%o eta*A1=\n",r1,S1eta;
        if Valuation(alpha) eq 0 then
          gamma := alpha;
        else 
          gamma := L!xi;
        end if; 
        phi01 := RL!Evaluate(nuexp2[1][2],gamma);
        new_phis := {};
        Thetas := [r[1]:r in Roots(S1eta-(phi01-psi01R))];
        vprintf Monge,2:"PolRedPadic: m=0 transforming phi*_(0,1) from %o to %o\n",phi01,psi01R;
        vprintf Monge,2:"PolRedPadic: m=0 Thetas=%o\n",Thetas;
        if Thetas eq [] then
          error "PolRedPadic: reduction step slope 1 failed";
        end if;
        for theta in Thetas do
          vprintf Monge,2:"PolRedPadic: m=0 transformation alpha->alpha+(%o)*nu(alpha)\n",theta;
          new_beta := alpha+(L!theta)*nualpha; 
          new_phi := CharacteristicPolynomial(new_beta,K);
          if IsEisenstein(Phi) then // TODO needed ?
            if Parent(A_phi)!ResidualPolys(new_phi) eq A_phi then
              Include(~new_phis,<new_phi,new_beta>);
            end if;
          else
            Include(~new_phis,<new_phi,new_beta>);
          end if;
          vprintf Monge,4:"PolRedPadic: m=0 theta=%o now phi*_(0,1)=%o\n",theta,RL!Evaluate(Expansion2(new_phi,nu)[1][2],gamma);
          vprintf Monge,3:"PolRedPadic: m=0 theta=%o now %o\n",theta,String(easyreduce(new_phi):wherenu);
          if not RL!Evaluate(Expansion2(new_phi,nu)[1][2],gamma) eq psi01R then
             error "PolRedPadic: reduction step m=1 failed";
          end if;
        end for;
        M := new_phis;
        vprintf Monge,1:"PolRedPadic: m=%o reduced are isomorphic %o\n",1,[HasRoot(Polynomial(L,psi[1])):psi in M];

        // other levels
        for m in [1..easystart-1] do
           vprintf Monge,2:"PolRedPadic: m=%o reduction with alpha -> alpha+theta*nu(alpha)^%o\n",m,m+1;
           new_M := {};
           for phiandbeta in M do
              phi := phiandbeta[1]; beta := phiandbeta[2];
              nuexp2, nuexp := Expansion2(phi,nu : limit := easylimit);
              
              phi0 := -nuexp[1];
              phi0beta := Evaluate(phi0,beta); 
              nubeta := Evaluate(nu,beta);
              eta := LtoRL(nubeta^e/p);
              vprintf Monge,2:"PolRedPadic: m=%o reducing %o\n",m,String(easyreduce(phi):wherenu);
              Am, PHIm := AdditiveResidualPoly(phi,nu,beta,m);
              i := PHIm mod e;
              k := PHIm div e;
              vprintf Monge,1:"PolRedPadic: m=%o Phi(m)=%o eta=%o A_m=%o\n",m,PHIm,eta,Am;

              phisik := nuexp2[i+1][k+1];
              vprintf Monge,2:"PolRedPadic: m=%o improving phi*_(%o,%o)=%o\n",m,i,k,phisik;
              //G phisikbeta := LtoRL(Evaluate(phisik,beta));
              if Valuation(beta) eq 0 then
                gamma := LtoRL(beta);
              else
                gamma := RL.1;
              end if;
              phisikbeta := Evaluate(phisik,LtoRL(gamma));
              
              FB := Basis(RL,Fp); 
              FL := [Eltseq(Evaluate(eta^k*Am,b)):b in FB]; 
              FM := Matrix(FL);
              // reduce phisikbeta by the image of eta^j*Am 
              Mecho := EchelonForm(Matrix([Reverse(r): r in RowSequence(FM)]));
              vdelta := Vector(Reverse(Eltseq(phisikbeta)));
              jB := 1;
              iB := 1;
              done := false;
              while iB le #FB and not done do // row counter
                while Mecho[iB][jB] eq 0 and not done do
                  if jB lt #FB then jB := jB + 1; else done := true; end if;
                end while;
                if not done then
                  vb := Vector(Mecho[iB]);
                  ab := vdelta[jB]/vb[jB];
                  vdelta := vdelta - ab*vb;
                  iB := iB + 1;
                end if;
              end while;
              delta := RL!Reverse(Eltseq(vdelta));
              vprintf Monge,2:"PolRedPadic: m=%o want phi*_(%o,%o)=%o\n",m,i,k,delta;
              // find coefficient for reduction
              sol, kernel  := Solution(Matrix(FM),Vector(Eltseq(phisikbeta-delta)));
              theta := RL!Eltseq(sol);
              vprintf Monge,3:"PolRedPadic: m=%o theta=%o\n",m,theta;
              Thetas := [theta+RL!Eltseq(a):a in Set(Kernel(FM))];
              vprintf Monge,2:"PolRedPadic: m=%o Thetas=%o\n",m,Thetas;
              new_phis := {};
              for theta in Thetas do
                vprintf Monge,2:"PolRedPadic: m=%o theta=%o transformation alpha->alpha+(%o)*nu(alpha)^%o\n",m,theta,theta,m+1;
                new_beta := beta+(L!theta)*nubeta^(m+1);
                new_phi := CharacteristicPolynomial(new_beta,K);
                Include(~new_phis,<new_phi,new_beta>);
                vprintf Monge,4:"PolRedPadic: m=%o theta=%o now phi*_(%o,%o)=%o\n",m,theta,i,k,Expansion2(new_phi,nu)[i+1][k+1];
                vprintf Monge,3:"PolRedPadic: m=%o theta=%o now phi=%o\n",m,theta,String(easyreduce(new_phi):wherenu);
              end for;
              new_M join:= new_phis;
            end for;
            M := new_M;
            vprintf Monge,5:"PolRedPadic: candidates \n %o", &cat[String(easyreduce(psi[1]):wherenu) cat "\n":psi in M];
            vprintf Monge,4:"PolRedPadic: m=%o reduced are isomorphic %o\n",m,[HasRoot(Polynomial(L,psi[1])):psi in M];
          end for;
        M := {easyreduce(phibeta[1]): phibeta in M};
        vprintf Monge,4:"PolRedPadic: before final distinction\n %o", &cat[String(psi:wherenu) cat "\n":psi in M];
        vprintf Monge,5:"PolRedPadic: easy reduced are isomorphic %o\n",[HasRoot(Polynomial(L,psi)):psi in M];
        return M;

end function;


intrinsic PolRedPadicTame(phi::RngUPolElt) -> .
{Reduction of a tamely ramified polynomial }
  K := CoefficientRing(phi);
  p := Prime(K);
  e0 := Degree(phi);
  if (e0 mod p) eq 0 then
    error "PolRedPadicTame works for tamely ramified extensions only";
  elif not IsEisenstein(phi) then
    error "PolRedPadicTame works for Eisenstein polynomials only";
  end if;
  Kx<x> := PolynomialRing(K);
  pi := UniformizingElement(K);
  U, KtoU := ResidueClassField(K);
  xi := PrimitiveElement(U);
  phi0 := ConstantCoefficient(phi);
  phi01 := KtoU(phi0 div pi);
  l := discrete_log(phi01);
  b := Gcd(e0,#U-1);
  r := l mod b;
  psi := x^e0+pi*K!(xi^r);
  vprintf Monge,3:"PolRedPadicTame: reduced to %o\n",psi;
  vprintf Monge,5:"PolRedPadicTame: reduced is isomorphic %o\n",IsIsomorphic(phi,psi);
  return psi;
end intrinsic;

intrinsic PolRedPadicTame(Phi::RngUPolElt,nu::RngUPolElt,alpha:distinguished:=true,conjugates:="auto") -> .
{}
  K := CoefficientRing(Phi);
  Kx<x> := PolynomialRing(K);
  L := Parent(alpha);
  Ly<y> := PolynomialRing(L);

  if conjugates cmpeq "auto" then 
    if Degree(nu) eq 1 then 
      conjugates := false;
    else
      conjugates := true;
    end if;
  end if;

  pi := UniformizingElement(K);
  p := Prime(L);
  phi := DefiningPolynomial(L);
  Lur<a> := CoefficientRing(phi);
  Lurx<x> := PolynomialRing(Lur);
  U, LurtoU := ResidueClassField(Lur);
  e0 := Degree(phi);
  if conjugates and Degree(nu) ne 1 then
    gaut, maut := AutomorphismGroup(Lur);
    aut := [ maut(tau) : tau in gaut ];
    phis := {};
    for tau in aut do
      vprintf Monge,1: "PolRedPadicTame: %o->%o\n",Lur.1,tau(Lur.1);
      tauphi := Lurx![tau(c): c in Coefficients(phi)];
      Include(~phis,tauphi);  
    end for;
  else
    phis := {phi};
  end if;
  M := {};
  for tauphi in phis do
    psi := PolRedPadicTame(tauphi);
    if Degree(nu) eq 1 then
      Include(~M,psi);
    else
      psi0 := ConstantCoefficient(psi);
      psi01 := LurtoU(psi0 div pi);
      psi01coeffs := Eltseq(psi01);
      Psi01 := Kx!psi01coeffs;
      Psi := Kx!nu^e0+Psi01*p;
      vprintf Monge,3:"PolRedPadicTame: reduced to %o\n",String(Psi);
      vprintf Monge,5:"PolRedPadicTame: reduced is isomorphic %o\n",HasRoot(Polynomial(L,Psi));
      Include(~M,Psi);
    end if;
  end for;
  if distinguished then
    PSI := Distinguished(M:nu:=nu);
    return PSI;
  else
    return M;
  end if;
end intrinsic;


intrinsic PolRedPadic(Phi::RngUPolElt,nu::RngUPolElt,alpha:distinguished:=true,conjugates:="auto") -> .
        {Phi in Zp[x] in Eisenstein Form, Phi(alpha)=0, nu(alpha) uniformizer of Qp(alpha), 
return the Krasner- Monge reduction of Phi}
  Kx<x> := Parent(Phi); 
  K := CoefficientRing(Phi);

  L := Parent(alpha);
  Lt<t> := PolynomialRing(L); 

  if conjugates cmpeq "auto" then 
    if Degree(nu) eq 1 then 
      conjugates := false;
    else
      conjugates := true;
    end if;
  end if;

  RL, LtoRL := ResidueClassField(L);
  p := Prime(L);
  U := BaseRing(L);
  Uy<y> := PolynomialRing(U);  
  pi := UniformizingElement(L);
  psi := DefiningPolynomial(L);
  gamma := Roots(Lt!nu-pi)[1][1];
  phi := CharacteristicPolynomial(gamma,BaseRing(U));
  vprintf Monge,4:"PolRedPadic(Phi,%o,alpha): char poly isomorphic is %o\n",nu,HasRoot(Polynomial(L,phi));
  n := Degree(phi);
  A, psis := ResidualPolyDistinguished(psi:conjugates := conjugates,constant_first);
  //vprint Monge,2:"ResidualPolyDistinguished: ",A, [String(psi:wherenu):psi in psis];
  M := {};
  for psi in psis do
    //vprintf Monge,4: "ResidualPolyDistinguished: tau(phi) = %o\n" ,String(psi:wherenu);
    vprintf Monge,5:"PolRedPadic(Phi,%o,alpha): distinguished is isomorphic is %o\n",nu,HasRoot(Polynomial(L,psi));
    thisphi, nu, thisalpha := OysteinPoly(psi,K);
    //vprintf Monge,5:"PolRedPadic(Phi,%o,alpha): Oystein polynomial is isomorphic is %o\n",nu,HasRoot(Polynomial(L,phi));
    psi01 := Coefficient(psi,0) div p;
    newphis := pol_red_padic_sub(thisphi,Kx!nu,thisalpha,psi01);
    vprintf Monge,5:"PolRedPadic(Phi,%o,alpha): reduced are isomorphic %o\n",nu,[HasRoot(Polynomial(L,psi)):psi in newphis];
    M join:= newphis;
  end for;
  if distinguished then
    PSI := Distinguished(M);
    vprintf Monge,4:"PolRedPadic(Phi,%o,alpha): distinguished is isomorphic is %o\n",nu,HasRoot(Polynomial(L,PSI));
    return PSI;
  else
    return M; 
  end if;
end intrinsic;


intrinsic PolRedPadic(Phi::RngUPolElt,K::RngPad:distinguished:=true,conjugates:="auto") -> .
{For Phi in O_L irreducible return a Krasner-Monge reduced polynomial Psi such that L[x]/(Phi)=K[x]/(Psi).}
   p := Prime(K);
   vprintf Monge,2:"PolRedPadic: converting to Oystein polynomial\n";
   phi, nu, alpha := OysteinPoly(Phi,K);
   vprintf Monge,2:"PolRedPadic: ramification index is %o and inertia degree is %o\n",Degree(phi)/Degree(nu),Degree(nu);
   L := Parent(alpha);
   // psi := DefiningPolynomial(L);
   if RamificationIndex(L,K) mod p ne 0 then
     M := PolRedPadicTame(phi,nu,alpha:distinguished:=distinguished,conjugates:=conjugates);
   else
     M := PolRedPadic(phi,nu,alpha:distinguished:=distinguished,conjugates:=conjugates);
   end if;
   return M;
end intrinsic;

intrinsic PolRedPadic(Phi::RngUPolElt:distinguished:=true,conjugates:="auto") -> .
{For Phi in OK irreducible return a Krasner-Monge reduced polynomial Psi such that K[x]/(Phi)=K[x]/(Psi).}
  return PolRedPadic(Phi,CoefficientRing(Phi):conjugates:=conjugates);
end intrinsic;

intrinsic PolRedPadic(f::RngUPolElt[RngInt],p::RngIntElt:prec:=300,distinguished:=true) -> .
{
The distinguished reduced generating polynomial of the extension generated by f over Zp. 
}
  Zp := pAdicRing(p,prec);
  ZpX<X> := PolynomialRing(Zp);
  Phi := ZpX!f;
  prec := Max(SuggestedPrecision(Phi),prec);
  Zp := pAdicRing(p,prec);
  ZpX<X> := PolynomialRing(Zp);
  Phi := ZpX!f;
  Psi := PolRedPadic(Phi,Zp: distinguished:=distinguished, conjugates :=true);
  F := Polynomial(Integers(),Psi);
  return Psi;
end intrinsic;


