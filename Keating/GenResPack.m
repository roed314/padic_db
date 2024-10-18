intrinsic pRed(g::RngMPolElt,p::RngIntElt) -> RngMPolElt
{This is a utility for use with Residual.  It eliminates
any terms in g which are not of minimum valuation.  Here
x is regarded as a uniformizer for the extension field.}
   R:=Parent(g);
   gcdc:=GCD(Coefficients(g));
   pval:=Valuation(gcdc,p);
   I:=Ideal(R!(p^(pval+1)));
   RI:=quo<R|I>;
   gRI:=RI!g;
   g2:=R!gRI;
   n:=Degree(g,1);
   d:=Degree(GCD(g2,R.1^n));
   J:=Ideal(R.1^(d+1));
   RJ:=quo<R|J>;
   g2J:=RJ!g2;
   g3:=R!g2J;
   return g3;
end intrinsic;;

intrinsic Residual(f::RngMPolElt,vertices::[.],Fppolz::RngUPol) -> []
{Computes the residual polynomials of the generic polynomial
corresponding to a family of p-adic fields.  The input "f"
should be the generic polynomial, but given as a polynomial
over Z rather than Z_p.  The input "vertices" is the sequence
of vertices of the ramification polygon of the family.  The
input "Fppolz" is a polynomial ring in one variable whose
coefficients are polynomials over F_p in the parameters used
to define the family.  The output will be an element of Fppolz.}
   Fppol:=BaseRing(Fppolz);
   Fp:=BaseRing(Fppol);
   p:=Characteristic(Fp);
   n:=Degree(f,1);
   R:=Parent(f);
   I:=Ideal(f);
   RI:=quo<R|I>;
   rank:=Rank(R);
   n:=Degree(f,1);
   fcoeffs:=Coefficients(f,1);
   ramcoeffs:=[(R!0)^^(n+1)];
   for i in [0..n] do
      for j in [i..n] do
         ramcoeffs[i+1]:=ramcoeffs[i+1]+Binomial(j,i)*fcoeffs[j+1]*R.1^j;
      end for;
      ramcoeffbar:=RI!ramcoeffs[i+1];
      ramcoeffs[i+1]:=R!ramcoeffbar;
   end for;
   Fpvars:=[];
   for i in [1..rank-1] do
      Append(~Fpvars,Fppol.i);
   end for;
   v:=#vertices;
   respolys:=[];
   Z:=IntegerRing();
   for i in [1..v-1] do
      k:=vertices[i][1];
      l:=vertices[i+1][1]-k;
      m:=vertices[i+1][2]-vertices[i][2];
      slope:=m/l;
      h:=-Numerator(slope);
      e:=Denominator(slope);
      lequo:=l div e;
      kppow:=Valuation(Coefficients(pRed(ramcoeffs[k+1],p))[1],p);
      kxpow:=Degree(GCD(pRed(ramcoeffs[k+1],p),R.1^n));
      valk:=n*kppow+kxpow;
      newrespoly:=Fppolz!0;
      for j in [0..lequo] do
         ramcoeff:=pRed(ramcoeffs[j*e+k+1],p);
         xpowerbar:=RI!(R.1^(valk-j*h));
         xpower:=pRed(R!xpowerbar,p);
         rppow:=Valuation(Coefficients(ramcoeff)[1],p);
         xppow:=Valuation(Coefficients(xpower)[1],p);
         rxpow:=Degree(GCD(ramcoeff,R.1^n));
         xxpow:=Degree(GCD(xpower,R.1^n));
         vr:=n*rppow+rxpow;
         vx:=n*xppow+xxpow;
         error if vx gt vr, "Nonintegral quotient.",ramcoeff,xpower;
         if vx eq vr then
            dr:=Degree(ramcoeff,R.1);
            xcoeff:=Z!(Coefficients(xpower)[1]);
            rcoeff:=Z!(Coefficients(ramcoeff)[1]);
            frac:=rcoeff/xcoeff;
            Fpfrac:=Fp!frac;
            Zfrac:=Z!Fpfrac;
            quot:=Zfrac*(ramcoeff div(rcoeff*R.1^dr));
            zFpvars:=Insert(Fpvars,1,Fppol!0);
            rescoeff:=Evaluate(quot,zFpvars);
            newrespoly:=newrespoly+rescoeff*Fppolz.1^j;
         end if;
      end for;
      Append(~respolys,newrespoly);
   end for;
   return respolys;
end intrinsic;
