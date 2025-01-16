intrinsic PReduce(u::RngPadElt) -> .
{Returns x,y such that u=x*y^p, with v(x-1) maximized.}
   A:=Parent(u);
   p:=Prime(A);
   e:=Valuation(A!p);
   if e gt 1 then
      pi:=A.1;
   else
      pi:=A!p;
   end if;
   prec:=Precision(u);
   bprec:=Maximum({Ceiling(prec/p),prec-e});
   F:=ResidueClassField(A);
   deg:=Degree(F);
   v:=A!1;
   flag:=0;
   while flag eq 0 do
      if Valuation(u-1) eq Precision(u) then
         return u, v^(-1)+O(pi^bprec);
      end if;
      a:=Valuation(u-1);
      if a lt p*e/(p-1) and not IsDivisibleBy(a,p) then
         return u, v^(-1)+O(pi^bprec);
      elif a lt p*e/(p-1) and IsDivisibleBy(a,p) then
         r:=a div p;
         coeff:=(u-1) div pi^a;
         w:=1-coeff^(p^(deg-1))*pi^r;
         u:=u*w^p;
         v:=v*w;
      elif a gt p*e/(p-1) then
         s:=(u-1) div p;
         u:=u*(1-s)^p;
         v:=v*(1-s);
      elif a eq p*e/(p-1) then
         e0:=e div (p-1);
         lifts:=TeichmuellerSystem(A);
         flag:=1; i:=1;
         while i le p^deg and flag eq 1 do
            w:=1+lifts[i]*pi^e0;
            y:=u*w^p;
            if Valuation(y-1) gt p*e0 then
               u:=y; v:=v*w;
               flag:=0;
            end if;
            i:=i+1;
         end while;
      end if;
   end while;
   return u, v^(-1)+O(pi^bprec);
end intrinsic;

intrinsic pPowReduce(u::RngPadElt) -> .
{Returns x,y,r such that u=x*y^(p^r), with r as large
as possible allowing v(x-1)>v(u-1), and v(x-1) maximized
for r.}
   A:=Parent(u);
   p:=Prime(A);
   r:=0;
   c:=u;
   while 0 ne 1 do
      a,b:=PReduce(c);
      if Valuation(b-1) eq Precision(b) then
         return u*a^(-p^r), a, r;
      end if;
      c:=b;
      r:=r+1;
   end while;
end intrinsic;

intrinsic Zetap(A::RngPad) -> RngPadElt
{Returns a primitive pth root of unity in A.
If A does not have pth roots of unity, returns 0.}
   p:=Prime(A);
   e:=Valuation(A!p);
   if not IsDivisibleBy(e,p-1) then
      return 0;
   end if;
   if e gt 1 then
      pi:=A.1;
   else
      pi:=p;
   end if;
   e0:=e div (p-1);
   lifts:=TeichmuellerSystem(A);
   for t in lifts do
      z:=1+t*pi^e0;
      if z ne 1 and Valuation(z^p-1) gt p*e0 then
         a,b:=PReduce(z^p);
         return z*b^(-1);
      end if;
   end for;
   return 0;
end intrinsic;

intrinsic pPowRoots(jset::[],e::RngIntElt,p::RngIntElt) -> RngIntElt
{The input is a jump set of a p-adic field K whose
absolute ramification index is e.  The output is the
maximum nonnegative integer d such that K contains a
primitive p^d root of unity.}
   k:=#jset;
   if k eq 0 then
      return 0;
   end if;
   d:=1;
   while d lt k do
      if jset[k-d+1] eq Min({p*jset[k-d],jset[k-d]+e}) then
         d:=d+1;
      else
         return d;
         break;
      end if;
   end while;
   return d;
end intrinsic;

intrinsic JumpSet(A::RngPad) -> .
{Outputs the jump set associated to the 1-units of A.
If A does not contain a primitive pth root of unity,
returns the empty sequence.}
   p:=Prime(A);
   e:=Valuation(A!p);
   e0:=e div (p-1);
   n:=Valuation(e,p)+1;
   A:=ChangePrecision(A,n*e+e0+1);
   u:=Zetap(A);
   if u eq 0 then
      return [],0;
   end if;
   jumps:=[0^^n];
   while Valuation(u-1) lt Precision(u) do
      a,b,r:=pPowReduce(u);
      jumps[n-r]:=Valuation(b-1);
      u:=u*b^(-p^r);
   end while;
   d:=1;
   while d lt n do
      if jumps[n-d+1] eq 0 then
         d:=d+1;
      else
         break;
      end if;
   end while;
   for i:=2 to n do
      if jumps[i] eq 0 then
         jumps[i]:=Minimum({p*jumps[i-1],e+jumps[i-1]});
      end if;
   end for;
   return jumps,d;
end intrinsic;

intrinsic JumpSetJJ(data) -> .
{Outputs a list of jump sets for 1-units in p-adic fields
which are characterized by data in the format given in
https://hobbes.la.asu.edu/LocalFields/file-src/File-format }
   jumplist:=[];
   rootslist:=[];
   for field in data do
      p:=field[1];
      e:=field[3];
      if p ne 2 then
         r:=Valuation(e,p)+2;
      else
         r:=Valuation(e,p)+3;
      end if;
      Zp:=pAdicRing(p,r);
      R<t>:=PolynomialRing(Zp);
      UnramPoly:=R!(field[12][1]);
      Zq<t>:=ext<Zp|UnramPoly>;
      f:=hom<PolynomialRing(IntegerRing())->Zq|t>;
      S<y>:=PolynomialRing(Zq);
      g:=hom<PolynomialRing(PolynomialRing(IntegerRing()))->S|f,y>;
      EisPoly:=g(field[12][2]);
      A<pi>:=ext<Zq|EisPoly>;
      jumps,roots:=JumpSet(A);
      jumplist:=Append(jumplist,jumps);
      rootslist:=Append(rootslist,roots);
   end for;
   return jumplist,rootslist;
end intrinsic;
