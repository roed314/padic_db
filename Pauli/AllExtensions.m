///////////////////////////////
// Indices of inseparability

intrinsic AllIndicesOfInseperability(K,n,j) -> .
{Input:  Output: }
  js := PossibleDiscriminants(K,n);
  if not j in js then error "j+n-1 =",j,"+",n,"- 1 is not a valid discriminant exponent.  Valid values for j are",js; end if;
  p := Prime(K);
  nu := Valuation(n,p);
  ek := RamificationIndex(K);
  Is := [[0]];
  k := nu;
  while k gt 1 do
    Isk := [];
    for this_i in Is do
      //"this_i",this_i;
      if Valuation(this_i[#this_i],p) lt k then
        Append(~Isk, this_i cat [this_i[#this_i]]);
      else
        for ij1 in [this_i[#this_i]..Minimum(this_i[#this_i]+n*ek,j)] do
          if ij1 ge j-Valuation(ij1,p)*n*ek then 
            if ij1 eq this_i[#this_i] + n*ek then
              Append(~Isk, this_i cat [ij1]);  
            elif 1 le Valuation(ij1,p) and Valuation(ij1,p) le k-1 then
              Append(~Isk, this_i cat [ij1]);
            elif ij1 eq j then
              Append(~Isk, this_i cat [ij1]);
            end if;
          end if;
        end for;
      end if;
      Is := Isk;
    end for;
    //k,Is;
    k := k-1;
  end while;
  Isr := [];
  for this_i in Is do
    new_i := this_i cat [j];
    if (Valuation(new_i[#new_i],p) ge 1 and new_i[#new_i] eq new_i[#new_i-1]+n*ek) or Valuation(new_i[#new_i],p) eq 0 then
      Append(~Isr,Reverse(new_i));
    end if;
  end for;
  return Isr;
end intrinsic;


intrinsic IndicesOfInseperability(f) -> .
{Input: f - an Eisenstein polynomial (correctness of input is not checked).  Output: List of inices of Inseperability}

  n := Degree(f);
  K := CoefficientRing(f);
  nu := Valuation(K!n);
  p := Prime(K);
  i_0 := Valuation(Discriminant(f))-n+1;
  itilde := [];
  for j in [0..nu] do
    Append(~itilde, Minimum([n*Valuation(Coefficient(f,i))+i-n : i in [1..n] | Valuation(K!i) le j]));
  end for;
  i := [-1: k in [0..nu]];
  i[nu+1] := 0;
  for j in [nu-1..0 by -1] do
    i[j+1] := Minimum(itilde[j+1],i[j+2]+n*Valuation(K!p));
  end for;
  return i;
end intrinsic;

////////////////////////////////////////////////////////////////
///



///////////////////////////////////////////////////////////////


function RamificationPolygon_sub(phi);

        K:=CoefficientRing(phi);

        n:=Degree(phi);
        L:=ext<K|phi>;
        alpha:=L.1;

        Lx<x>:=PolynomialRing(L);
        rho:=Evaluate(phi,alpha*x + alpha) div (alpha^n); //Ramification Polynomial
        rho:=Lx!rho;
        ChangePrecision(~rho,Precision(L));

        ramification_polygon := NewtonPolygon(rho);

        return ramification_polygon,rho;
end function;


intrinsic RamificationPolygon(f::RngUPolElt[RngPad]) -> .
{
Returns the ramification polygon of f.

The ramification polygon is the Newton polygon of f(a*x+a)/a^n where a is a root of f and n is the degree of f.

        EXAMPLES:

        The vertices of a ramification polygon and the slopes of its segments:

            R := pAdicRing(3,20); 
            Rx<x> := PolynomialRing(R);
            f := x^27+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6;
            rp := RamificationPolygon(f);
            Vertices(rp);
            // [(1, 51), (3, 24), (9, 9), (27, 0)]
            Slopes(rp);
            // [-27/2, -5/2, -1/2]

        A ramification polygon with a horizontal segment:

            R := pAdicRing(3,20); 
            Rx<x> := PolynomialRing(R);
            f := x^108+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6;
            rp := RamificationPolygon(f);
            Vertices(rp);
            // [(1, 132), (3, 24), (9, 9), (27, 0), (108, 0)]
            Slopes(rp);
            // [-54, -5/2, -1/2, 0]

        AUTHORS:

}

        if not IsEisenstein(f) then
            Error("Ramification polynomials are only defined for Eisenstein polynomials");
        end if;


        return LowerVertices(RamificationPolygon_sub(f));

end intrinsic;


function vertices_slopes(ramification_polygon)

        vertices := LowerVertices(ramification_polygon);
        slopes := Slopes(ramification_polygon)[1..#vertices-1];

        if vertices[1][1] eq 0 then
                vertices:=vertices[2..#vertices];               //get rid of point with x-coordinate=0
                slopes:=slopes[2..#slopes];                     //get rid of segment with infinite slope
        end if;

        return vertices, slopes;

end function;

/*

for k in [10..3000] do
"k =",k;
  p := Random([2,3,5,7,11,13,17]);
"p =",p;
  n := Random([a : a in [p^2,2*p^2,6*p^2,3*p^2,p^4,p^3,p^4,p^5] | a le k] cat [12*p]);
"n =",n;
  zp := pAdicRing(p,30);
  J := (PossibleDiscriminants(zp,n));
  j := Random(J[1..Minimum(#J,20)]);
"j =",j;
  R := Random(AllRamificationPolygons(zp,n,j));
"R =",R;
  A := Random(AllResidualPolynomials(zp,R,p));
"A =",A;
  Ls := (AllTotallyRamifiedExtensions(zp,R,A,1:want_filter:=false));
  for kk in [1..Minimum(#Ls div 2,30)] do
    L := Random(Ls);
"R_L",RamificationPolygonWithColinearPoints(L);
"A_L",ResidualPolynomials(L);
    if R ne RamificationPolygonWithColinearPoints(L) then
       error p,R,L;
    elif A ne ResidualPolynomials(L) then
       error p,A,L;
    end if;
  end for;
end for;
*/


intrinsic RamificationPolygonWithColinearPoints(f::RngUPolElt[RngPad]) -> .
{
        Returns the ramification polygon of f as a list of points including all points on segments of the lower convex hull.
        The ramification polygon is the Newton polygon of f(a*x+a)/a^n where a is a root of self and n is the degree of f.
        We compare the output of this method with the output of `ramification_polygon`.

             R := pAdicRing(5,20); 
             Rx<x> := PolynomialRing(R);
             f := x^25 + 20*x^6 + 20*x^5 + 5;
             RamificationPolygonWithColinearPoints(f);
             // [(1, 6), (5, 5), (25, 0)]

             // The colinear point (5, 5) is missing on the ramification polygon returned by::

             RamificationPolygon(f)
             // Finite Newton polygon with 2 vertices: (1, 6), (25, 0)

        When the generated extensions has a tamely ramified subextension::

             f := x^100 + 20*x^6 + 20*x^5 + 5;
             RamificationPolygonWithColinearPoints(f);
             // [(1, 6), (5, 5), (25, 0), (50, 0), (75, 0), (100, 0)]
             RamificationPolygon(f);
             //Finite Newton polygon with 3 vertices: (1, 6), (25, 0), (100, 0)

        AUTHORS:

        - Brian Sinclair (2017-07-19): initial version
        - Sebastian Pauli (2018-03-13): Magma version
}
        
        if not IsEisenstein(f) then
            Error("the polynomial self must be Eisenstein");
        end if;

        //  First we find the ordinates of points above p^k
        verts := [];
        vv := [Valuation(cc) : cc in Coefficients(f)];
        k := BaseRing(f);
        p := Prime(k);
        n := Degree(f);
        su := Valuation(n,p);
        for i in [0..su-1] do
            abscissa := p^i;
            ordinate := Minimum([n * (Valuation(k!(Binomial(kk,abscissa))) + vv[kk+1] - 1) + kk :  kk in [abscissa..n]]);
            Append(~verts,<abscissa,ordinate>);
        end for;

        //  Now we add the tame segment
        for i in [p^su..n-1] do
            if Valuation(Binomial(n,i),p)  eq  0 then
                Append(~verts,<i,0>);
            end if;
        end for;

        //  Finally the point for the monic leading term
        Append(~verts,<n,0>);

        //  Next we need to take the lower convex hull of these points
        function cross(o, a, b)
            //  2D cross product of the vectors oa and ob.
            return (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1]);
        end function;
//"verts",verts;
        lower := [verts[1]];
        segments := [];
        for i in [2..#verts] do
//"i",i,#lower;
            //  We check cross < 0 since we want to retain points on the boundary.
//if #lower ge 2 then cross(lower[#lower-1], lower[#lower], verts[i]); end if;
            while #(lower) ge 2 and cross(lower[#lower-1], lower[#lower], verts[i]) lt 0 do
                lower := lower[1..#lower-1];
            end while;
            Append(~lower,verts[i]);
        end for;
        if #(lower) le 1 then
            Error("Not enough vertices");
        end if;
        return lower;
end intrinsic;

intrinsic HasseHerbrand(f::RngUPolElt[RngPad],m) -> .
    {
        Returns n times the (generalized) Hasse-Herbrand function of f evaluated at m.
        
        EXAMPLES:
             R := pAdicRing(2,200); 
             Rx<x> := PolynomialRing(R);
             f := x^16+2;
             RamificationPolygon(f);

             // Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)

        // We evaluate the Hasse-Herbrand function at various integers.

            [HasseHerbrand(f,m) : m in [0..18]];
            //[0, 16, 32, 40, 48, 52, 56, 60, 64, 66, 68, 70, 72, 74, 76, 78, 80, 81]

        Now a different degree:

             f := x^80+2;
             RamificationPolygon(f);
             //Finite Newton polygon with 6 vertices: (1, 320), (2, 240), (4, 160), (8, 80), (16, 0), (80, 0)
             [HasseHerbrand(f,m) : m in [0..16]];
             //[0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 168, 176, 184, 192, 200]

        AUTHORS:

        - Brian Sinclair (2017-07-20): initial version
        - Sebastian Pauli (2018-03-13): back to Magma
}
        rp := RamificationPolygon(f);
        return Minimum([pt[2]+m*pt[1] : pt in Vertices(rp)]);
end intrinsic;





function ResidualPolynomials_sub(R,rho)


        //rho is ramification polynomial.

        L:= CoefficientRing(rho);

        vertices, slopes:=vertices_slopes(R);
//"slopes",slopes;
        if vertices[1][1] eq 0 then
                vertices:=vertices[2..#vertices];               //get rid of point with x-coordinate=0
                slopes:=slopes[2..#slopes];                     //get rid of segment with infinite slope
        end if;

        a:=[Integers()!vertices[i][1]: i in [1..#vertices]];    //list of x-coordinates of vertices
        b:=[Integers()!vertices[i][2]: i in [1..#vertices]];    //list of y-coordinates of vertices

        l:=#slopes;

        pi_L:=UniformizingElement(L);

        e:=[];
        h:=[];

        for i in [1..l] do
                e[i]:=Denominator(-slopes[i]);          //list of (negative) slope denominators
                h[i]:=Numerator(-slopes[i]);                    //list of (negative) slope numerators
        end for;

        A:=[];

        RL,omega:=ResidueClassField(L);
        Ry<y>:=PolynomialRing(RL);

        for i in [1..l] do
                di:=a[i+1]-a[i];
                fin:=Integers()!(di/e[i]);
                Ai:=Polynomial(Ry,[&+[omega(Coefficient(rho,j*e[i]+a[i]) div pi_L^(-j*h[i]+b[i]))*y^j: j in [0..fin]]]);
                Ai:=Ry!Ai;
                A[i]:=Ai;
        end for;


        return A;

end function;





intrinsic ResidualPolynomials(f::RngUPolElt[RngPad]) -> .
        {
        Returns a list of the residual polynomials of the ramification polynomial of an Eisenstein polynomial self.
        
        EXAMPLES::

        The residual polynomials of the segments of the ramification polygon of an Eisenstein polynomial::

             R := pAdicRing(3,30); 
             Rx<x> := PolynomialRing(R);
             f := x^9+6*x^3+3;
             IsEisenstein(f);
             Slopes(RamificationPolygon(f));
             //Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
             ResidualPolynomials(f);
             // [z + 2, z^3 + 1]
             g := x^9 + 6*x^5 + 3; 
             Slopes(RamificationPolygon(g));
             ResidualPolynomials(g);



        //A ramfication polygon with a horizontal segment::

             f := x^90+6*x^3+3;
             RamificationPolygon(f);
             // Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
             Slopes(RamificationPolygon(f));
             // [-45, -1/2, 0]
             ResidualPolynomials(f);
             //[z^2 + 2, z^3 + 1, z^81 + z^72 + 1]

        //A ramification polygon with more segments::

             R := pAdicRing(2,300); 
             Rx<x> := PolynomialRing(R);
             f := x^16+2;
             NewtonPolygon(RamificationPolynomial(f));
             RamificationPolygon(f);
             // Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)
             ResidualPolynomials(f);
             // [z + 1, z^2 + 1, z^4 + 1, z^8 + 1]




        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version
        - Sebastian Pauli (2018-03-13): Magma
        }
        if not IsEisenstein(f) then
            Error("the polynomial self must be Eisenstein");
        end if;

        
        R,rho:=RamificationPolygon_sub(f);
        A:=ResidualPolynomials_sub(R,rho);

        return <a:a in A>;


end intrinsic;



intrinsic ResidualPolynomialOfComponent(f::RngUPolElt[RngPad],m::RngIntElt) -> .
        {
        Return the residual polynomial S_m of the (-m)-component of the ramifation polygon of polynomials self, which must be Eisenstein.

        Let N be the ramification polygon then [(k,w) in N | (-m)k + w := min[(-m)l+u|(l,u) in N]] is the (-m)-component of N.

        INPUT:

            A natural number m

        OUTPUT:

            The residual polynomial of the (`-m`)-component of the ramificaton polygon of self.

        EXAMPLES:

        In our first example, we have a polynomial over the degree 2 unramified extension of Q_2
        which has a ramification polygon with two segments of integral slope::

             R<g> := UnramifiedExtension(pAdicRing(2,30),2); 
             Rx<x> := PolynomialRing(R);
             f := x^8 + 2*g*x^6 + 4*g*x + 2;
             RamificationPolygon(f);
             //  Finite Newton polygon with 3 vertices: (1, 9), (2, 6), (8, 0)
             Slopes(RamificationPolygon(f)); 
             //  [-3, -1]
             ResidualPolynomials(f);
             //  [g0*z + g0, z^6 + g0]
             [ResidualPolynomialOfComponent(f,m) : m in [1..9]];
             //  [z^8 + g0*z^2, z^2, g0*z^2 + g0*z, z, z, z, z, z, z]

        Here we have a nonic polynomial over Q_3 whose ramification polygon has no segments
        of integral slope::

             R := pAdicRing(3,30); 
             Rx<x> := PolynomialRing(R);
             f := x^9 + 6*x^3 + 3;
             RamificationPolygon(f);
             //  Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
             Slopes(RamificationPolygon(f));
             //  [-9/2, -1/2]
             ResidualPolynomials(f);
             //  [z + 2, z^3 + 1]
             [ResidualPolynomialOfComponent(f,m) : m in [0..10]];
             //  [z^9, z^3, z^3, z^3, z^3, z, z, z, z, z]

        Here the ramification polygon has a horizontal segment::

             f := x^90 + 6*x^3 + 3;
             RamificationPolygon(f);
             //  Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
             Slopes(RamificationPolygon(f));
             //  [-45, -1/2, 0]
             ResidualPolynomials(f);
             //  [z^2 + 2, z^3 + 1, z^81 + z^72 + 1]
             [ResidualPolynomialOfComponent(f,m) : m in [0..10]];
             //  [z^90 + z^81 + z^9, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3]

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        REFERENCES:

        [PS17] S. Pauli and B. Sinclair, "Enumerating Extensions of (pi)-adic Fields with Given Invariants", International Journal of Number Theory
          (2017)
        }

        if not IsEisenstein(f) then
            Error("residual polynomials are only defined for Eisenstein polynomials");
        end if;

        Rx<x> := Parent(f);
        R := BaseRing(Rx);
        F := ResidueClassField(R);
        Fz<z> := PolynomialRing(F);

        rp := RamificationPolygon(f);

        if -m in Slopes(rp) then
            i := Position(Slopes(rp),-m);
            respol := ResidualPolynomials(f);
//"verti1",Vertices(rp)[i][1];            
            return (respol[i])*z^Integers()!(Vertices(rp)[i][1]);
        else
            L := [v[2]+v[1]*m : v in Vertices(rp)];
            mini := Minimum(L);
            mindex := Position(L,mini);
            return z^Integers()!(Vertices(rp)[mindex][1]);
        end if;
end intrinsic;



intrinsic HasResidualPolynomialClass(f,A) -> .
       { 
        Checks whether a list `A` of polynomials is in the same residual polynomial class as the residual polynomials of `f`.

        INPUT::

            A list `A` of polynomials over the residual class field of the coefficient ring of the polynomial `f`.

        OUTPUT::

            true if the polynomials in `A` are in thee same residual polynomial class as the residual polynomials of f.

        EXAMPLES::

        Clearly the residual polynomials of `f` are in the residual polynomial class of f.

             R := pAdicRing(3,30); Rx<x> := PolynomialRing(R);
             F := ResidueClassField(R);  
             f := x^9+6*x^3+3; 
             A := ResidualPolynomials(f);
             A;
             //  [z + 2, z^3 + 1]
             HasResidualPolynomialClass(f,A);
             //  true

        In the following the polynomial g generates an extension isomorphic to the extension generated by f.
        So the residual polynomials of g are in the residual polynomial class of g.

             g := f(2*x)/512
             ResidualPolynomials(g);
             //  [2*z + 2, z^3 + 2]
             HasResidualPolynomialClass(g,A);
             //  true

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-21): initial version

        REFERENCES:

        [PS17] S. Pauli and B. Sinclair, "Enumerating Extensions of (\pi)-adic Fields with Given Invariants",
        International Journal of Number Theory (2017)
}

        if not IsEisenstein(f) then
            Error("residual polynomials are only definesd for Eisenstein polynomials");
        end if;

        B := ResidualPolynomials(f);

        // quickly check whether the two residual polynomial classes have the same length
        if #A ne #B then
            return false;
        end if;

        // quickly check whether corresponding residual polynomials have same degree
        for i in [1..#B] do
            if Degree(B[i]) ne  Degree(A[i]) then
                return false;
            end if;
        end for;

        slopes := Slopes(RamificationPolygon(f));

        Rx<x> := Parent(f);
        R := BaseRing(Rx);
        F := ResidueClassField(R);
        Fz<z> := PolynomialRing(F);

        h := [Numerator(-s) : s in slopes];

        for delta in F do
            if delta ne 0 then
                g:=[delta^(-h[#h]*Degree(B[#B]))];
                for i in Reverse([1..#B-1]) do
                    g := [g[1]*delta^(-h[i]*Degree(B[i]))] cat g;
                end for;
                Bdelta :=  [g[i]*Evaluate(B[i],delta^h[i]*z) : i in [1..#B]];
                if A  eq  Bdelta then
                    return true;
                end if;
            end if;
        end for;
        return false;
end intrinsic;




/////////////////////////////////////////////////////////////////
//
// ENUM.M - Enumerate all extensions given additional invariants
//
//////////////////////////////////////////////////////////////////
// load "rampol-enum.m";
//////////////////////////////////////////////////////////////////


//
// RAMPOL-ENUM.M - Enumerate possible ramificaiton polygons
//

//////////////////////////////////////////////////////////////////
//load "rampol.m";

// Computes the lower convex hull of a set of two-dimensional points.
//
// Input: an enuemrated sequence of <x, y> pairs representing the points.
// Output: a list of vertices of the lower convex hull.
// Implements Andrew's monotone chain algorithm. O(n log n) complexity.
lower_convex_hull := function(points)
    // Sort the points lexicographically (tuples are compared lexicographically).
    // Remove duplicates to detect the case we have just one unique point.
    points := Sort(points);
 
    // Boring case: no points or a single point, possibly repeated multiple times.
    if #points le 1 then
        return points;
    end if;
 
    // 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    // Returns a positive value, if OAB makes a counter-clockwise turn,
    // negative for clockwise turn, and zero if the points are collinear.
    cross := func<o,a,b | (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])>;
 
    // Build lower hull 
    lower := [];
    for p in points do
        //while #lower ge 2 and cross(lower[#lower-1], lower[#lower], p) le 0 do  // removes points on lines
        while #lower ge 2 and cross(lower[#lower-1], lower[#lower], p) lt 0 do    // does not remove points on lines
            Prune(~lower);
        end while;
        Append(~lower,p);
    end for;

    return lower;
end function;

// Create the Ramification Polygon of polynomial f
//
// Input:
//   f - an Eisenstein polynomial (correctness of input is not checked)
//   method (optional): "min" the minimum-based definition (default)
//                      "def" the formal definition
//   output (optional): "newton" a NewtonPolyon object (default)
//                      "list" a list of points on the Newton polygon
// Output: 
//   The Newton polygon of the ramificaiton polygon of f (as NewtonPolygon or List)
intrinsic RamificationPolygon(f:method:="min",output:="newton") -> .
{Input: f - an Eisenstein polynomial (correctness of input is not checked) method (optional): "min" the minimum-based definition (default) "def" the formal definition output (optional): "newton" a NewtonPolyon object (default) "list" a list of points on the Newton polygon. Output: The Newton polygon of the ramificaiton polygon of f (as NewtonPolygon or List)
}    
    
    if method eq "def" then
      n := Degree(f);
      k := CoefficientRing(f);
      K<pi> := ext<k|f>;
      KX<X> := PolynomialRing(K);
      F := Evaluate(KX!f,pi*X+pi) div pi^n;
      rp := NewtonPolygon(F);
      c := Vertices(rp)[2..#Vertices(rp)]; // Removes the point at (0,+inf)
      return NewtonPolygon(lower_convex_hull(c));
    end if;
    if method eq "min" then
      n := Degree(f);
      p := Prime(BaseRing(f));
      e := RamificationIndex(BaseRing(f));
      a := Eltseq(f); a := a[2..#a];
      va := [Valuation(ai) : ai in a];
      c := [<1,Min([n*(e*Valuation(k,p)+va[k]-1)+k : k in [1..#va]])>];
      for i in [1..Valuation(n,p)] do
        Append( ~c,<p^i,Min([n*(e*Valuation(Binomial(k,p^i),p)+va[k]-1)+k : k in [p^i..#va]])> );
      end for;
	  
	  // append the horizontal segment if there is tame ramification
      if n gt c[#c][1] then Append(~c,<n,0>); end if;
	  if output eq "list" then
          return lower_convex_hull(c);
      end if;
      return NewtonPolygon(lower_convex_hull(c));
    end if;
	error "Invalid method: use 'min' or 'def'";
end intrinsic;

// Create Ramification Polygon from sequence of valuations of f_{p^i}
//
// Input:
//  va - sequence of valuations of f_{p^i}
//   p - prime of our p-adic field
//   j - from discriminant = p^{j+n-1}   (optional)
//   e - ramification index of base ring (optional)
// Output: NewtonPolygon
RamificationPolygonSeq := function(va,p:j:=0,e:=1)
      // If va is the output of most functions, then
      // it lacks the 0 for the monic leading term
      if va[#va] ne 0 then va cat:= [0]; end if;
      n := #va;
      if j eq 0 then
        c := [<1,Min([n*(e*Valuation(k,p)+va[k]-1)+k : k in [1..#va]])>];
      else
        c := [<1,j>];
      end if;
      for i in [1..Valuation(n,p)] do
        Append( ~c,<p^i,Min([n*(e*Valuation(Binomial(k,p^i),p)+va[k]-1)+k : k in [p^i..#va]])> );
      end for;

	  // append the horizontal segment if there is tame ramification
      if n gt c[#c][1] then Append(~c,<n,0>); end if;
	  
      return NewtonPolygon(lower_convex_hull(c));
end function;



// Compute the Ramification Polynomial of f
intrinsic RamificationPolynomial(f:RngPolElt) -> .
{The ramification polynomial of an Eisienstein  polynomial over a local ring.}
  n := Degree(f);
  k := CoefficientRing(f);
  K<pi> := ext<k|f>;
  KX<X> := PolynomialRing(K);
  F := Evaluate(KX!f,pi*X+pi) div pi^n;
  return F;
end intrinsic;


//////////////////////////////////////////////////////////////////
// Current Ramification Polygon
// used to hold data for building ramification polygons
CRP := recformat<
  K : RngPad,
  p : Integers(),   // prime of our p-adic field
  rpolyg : SeqEnum, // Current Vertices of our ramificaiton polygon
  lowerv : SeqEnum, // lower bound for valuations of f_{i}
  fixedv : SeqEnum, // list of fixed valuations of f_{p^i} (as points <i,v(f_p^i)>
  points : BoolElt  // true if we list all points, false for vertices only
                    // (lowerv will be higher to avoid unlisted points on lines if true)
>;

// Compute lower bound for valuations of all f_{p^i} due to point (p^s,y)
// Input:
//  s - abcissa p^s to consider
//  y - valuation of ramification polygon above p^s : v(c_{p^s})
//  p - prime of p-adic field
//  n - degree of polynomial
// Output: l(i,s) as SeqEnum
lfunc := function(s,y,K,n)
  as,bs := Quotrem(Integers()!y,Integers()!n);
  p := Prime(K);
  L := [];
  for i in [1..n-1] do
    if i mod p^s ne 0 then
      Append(~L,1);
    elif i lt bs then
      Append(~L, Max(2+as-Valuation(K!Binomial(i,p^s)),1) );
    else
      Append(~L, Max(1+as-Valuation(K!Binomial(i,p^s)),1) );
    end if;
  end for;
  return L;
end function;

// Compute lower bound for valuations due to no point above p^k
// Input:
//  k - abcissa p^k to consider
//  lvert,rvert - left and right endpoints of segment over p^k
//  p - prime of p-adic field
//  n - degree of polynomial
lnvfunc := function(k,lvert,rvert,K,n:colinear:=false)
  // Find the lower bound for valuations of f due to no point above p^k
  // above a segment with given left and right vertices
  slope := (lvert[2]-rvert[2])/(lvert[1]-rvert[1]);
  p := Prime(K);
  L := [];
  for i in [1..n-1] do
    if i mod p^k ne 0 then
      Append(~L,1);
    else
      mv := (1/n)*(slope*(p^k-lvert[1])+lvert[2]-i)+1-Valuation(K!Binomial(i,p^k));
      // So that we separate cases to include points on lines
      if colinear and mv in Integers() then
        Append(~L,mv+1);
      else
        Append(~L,Ceiling(mv));
      end if;
    end if;
  end for;
  return L;
end function;

//
// AllRamificationPolygonsSub
//

// Find all possible values for ramification polygon point above p^s given current polygon
// Input:
//  crp - current ramificaiton polygon
//   s  - abcissa p^s to consider
AllRamificationPolygonsSub := function(crp,s)
    // Get the degree of the polynomial from the lower bound length
    n := #crp`lowerv + 1;
    // Get j from the point <1,j>
    j := crp`rpolyg[1][2];

    // We know the point above p^0 to be <1,j> so we are almost done.
    if s eq 0 then
      // Check that missing vertices are acceptable
      missing := Valuation(crp`K!crp`rpolyg[#crp`rpolyg][1])-1;
      if missing gt 0 then      
        for k in [1..missing] do
          // Check that we can have a point above the line with current fixed valuations
          fbound := Min([n*(Valuation(crp`K!Binomial(crp`fixedv[iter][1],crp`p^k))+crp`fixedv[iter][2]-1)+crp`fixedv[iter][1] : iter in [1..#crp`fixedv] | crp`fixedv[iter][1] ge crp`p^k]);
          crossing := (j-crp`rpolyg[#crp`rpolyg][2])/(1-crp`rpolyg[#crp`rpolyg][1])*(crp`p^k-1)+j;
          if fbound lt crossing then
            // Upper bound is below the current segment
            return [];
          end if;
          
          // Find the lower bound for valuations of f due to no point above p^k
          Lnv := lnvfunc(k,<1,j>,crp`rpolyg[#crp`rpolyg],crp`K,n:colinear:=crp`points);
          for i in [1..#Lnv] do
            crp`lowerv[i] := Max([crp`lowerv[i],Lnv[i]]);
          end for;
        end for;
      end if;
      return [crp];
    end if;

    // Handle the case of s gt 0
    
    lastvert := crp`rpolyg[#crp`rpolyg];
    
    // Lower bound based on current minimum valutions
//"crp",#crp`lowerv,Parent(#crp`lowerv),crp`p,Parent(crp`p);
//[crp`p^s..#crp`lowerv];
    minbound := Min([n*(Valuation(crp`K!Binomial(k,crp`p^s))+crp`lowerv[k]-1)+k : k in [crp`p^s..#crp`lowerv]] cat [n*(Valuation(crp`K!Binomial(n,crp`p^s)))]);

    // Lower bound based on last added segment
    slope := func<a,b| (a[2]-b[2])/(a[1]-b[1])>;
    if crp`rpolyg[#crp`rpolyg][2] gt 0 then
      segbound := Floor( (lastvert[2]-crp`rpolyg[#crp`rpolyg-1][2])/(lastvert[1]-crp`rpolyg[#crp`rpolyg-1][1]) * (crp`p^s-lastvert[1])+lastvert[2] );
      segbound +:= (crp`p^s - segbound mod crp`p^s);
    else segbound := crp`p^s;
    end if;

    // Upper bound based on current fixed valuations (This assumes that <n,0> is present in fixedv)
    fbound := Min([n*(Valuation(crp`K!Binomial(crp`fixedv[k][1],crp`p^s))+crp`fixedv[k][2]-1)+crp`fixedv[k][1] : k in [1..#crp`fixedv] | crp`fixedv[k][1] ge crp`p^s]);

    // Upper bound based on segment <1,j> to lastvert <p^u,y_u> (only applies to vertices of polygon)
    ubound := j - ( (j-crp`rpolyg[#crp`rpolyg][2])/(crp`rpolyg[#crp`rpolyg][1]-1) )* (crp`p^s-1);
    // The following insures that we do not choose a y value on the segment
    if not crp`points and ubound in Integers() then ubound -:= 1;
    else ubound := Floor(ubound);
    end if;
    
    lowerbound := Max(minbound,segbound);
    upperbound := Min(ubound,fbound);
    possibley := [lowerbound..upperbound by crp`p^s];

    // Handle the trivial no vertex case
    
    if lowerbound gt ubound then
      // There is no possible vertex with the current valuation lower bounds.
      return $$(crp,s-1);
    end if;

    // Handle the non-trivial cases

    // Always pass the current polygon forward.
    crps := [crp];
    
    for y in possibley do
      as,bs := Quotrem(Integers()!y,Integers()!n);
      L := lfunc(s,y,crp`K,n);
      
      // CHECK THE VALIDITY OF THIS Y-VALUE
      
      polygon_is_okay := true;

      // Check that any previous "no vertex" choices work.
      missing := Valuation(crp`K!crp`rpolyg[#crp`rpolyg][1])-s-1;
      if missing gt 0 then
        slope := (y-crp`rpolyg[#crp`rpolyg][2])/(crp`rpolyg[#crp`rpolyg][1]-crp`p^s);
        for k in [s+1..s+missing] do
          // Upper bound based on current fixed valuations (This assumes that <n,0> is present in fixedv)
          fbound := Min([n*(Valuation(crp`K!Binomial(crp`fixedv[iter][1],crp`p^k))+crp`fixedv[iter][2]-1)+crp`fixedv[iter][1] : iter in [1..#crp`fixedv] | crp`fixedv[iter][1] ge crp`p^k]);
          crossing := (y-crp`rpolyg[#crp`rpolyg][2])/(crp`p^s-crp`rpolyg[#crp`rpolyg][1])*(crp`p^k-crp`p^s)+y;
          if fbound lt crossing then
            polygon_is_okay := false;
            break;
          end if;
          
          if polygon_is_okay then
            // Find the lower bound for valuations of f due to no point above p^k
            Lnv := lnvfunc(k,<crp`p^s,y>,crp`rpolyg[#crp`rpolyg],crp`K,n:colinear:=crp`points);

            // Check that lnv(i,s) doesn't conflict with fixed valuations
            for i in [1..#Lnv] do
              if Lnv[i] gt crp`lowerv[i] and Index([fx[1]:fx in crp`fixedv],i) gt 0 then
                // Increasing minimum valuation at fixed valuation - discard this polygon
                polygon_is_okay := false;
                break;
              end if;
            end for;
          end if;

          // If this missing point is okay, then increase lower bounds
          if polygon_is_okay then
            for i in [1..#Lnv] do
              L[i] := Max([crp`lowerv[i],Lnv[i],L[i]]);
            end for;
          end if;
        end for;
      end if;
      
      if bs gt 0 and crp`lowerv[bs] gt L[bs] then
        // L[bs] is less than current lower bound
        polygon_is_okay := false;
      end if;

      // Final check - Do the valuations prescribed by this point actually generate it?
      if polygon_is_okay then
        for i in [1..#L] do
          L[i] := Max([crp`lowerv[i],L[i]]);
        end for;
        ycheck := Min([n*(Valuation(crp`K!Binomial(k,crp`p^s))+L[k]-1)+k : k in [crp`p^s..#L]] cat [n*(Valuation(crp`K!Binomial(n,crp`p^s)))]);
        if y ne ycheck then
          // Prescribed valuations do not generate this y value - discard this polygon
          polygon_is_okay := false;
        end if;
      end if;
      
      // If y passes all checks, add the point <p^s,y>
      if polygon_is_okay then      
        // Add a fixed val (if bs ne 0), update our lower bounds, and point to RP
        newcrp := crp;
        if bs ne 0 then
          // If bs eq 0, then we do not add a new fixed valuation, as <n,0> defined it.
          Append(~newcrp`fixedv,<bs,L[bs]>);
        end if;
        newcrp`lowerv := L;
        Append(~newcrp`rpolyg,<crp`p^s,y>);
        Append(~crps,newcrp);
      end if;
    end for;

    // We need to return a flat list of all CRPs
    retrcrps := [];
    for i in [1..#crps] do
      retrcrps cat:= $$(crps[i],s-1);
    end for;
    return retrcrps;
end function;

//
// AllRamificationPolygons
//

// Find all possible values for ramification polygon point above p^i given current polygon
// Input:
//  n - degree of polynomial
//  p - prime of our p-adic field
//  j - discriminant p^{n+j-1}
//  output - (optional) "polygons" or "crps"
//  points - (optional) true to find all points, false for just vertices
intrinsic AllRamificationPolygons(K::RngPad,n,j:output:="polygons",points:=true) -> .
{Find all possible values for ramification polygon point above p^i given current polygon.  Input: n - degree of polynomial, p - prime of our p-adic field,j - discriminant p^(n+j-1). Output - (optional) "polygons" or "crps" points - (optional) true to find all points, false for just vertices}

    // Check that choice of j satisfies Ore's Conditions
    if not OreConditions(K,n,j) then
      error "Choice of discriminant fails Ore's Conditions.";
    end if;
    p := Prime(K);
    
    // Set up valuation lower bounds based on <1,j>
    a,b := Quotrem(Integers()!j,Integers()!n);
    L := lfunc(0,j,K,n);

    // Initialize fixed valuations based on <1,j>,<n,0>
    fixedv := [<n,0>];
    if b gt 0 then Append(~fixedv,<b,L[b]>); end if;

    // Handle the trivial case of n = p
    if n eq p then
      if output eq "crps" then
        crp := rec< CRP | K := K, p := p, rpolyg := [<1,j>,<n,0>], lowerv := L, fixedv := fixedv, points := true>;
        return [[crp]];
      else
        return [[<1,j>,<n,0>]];
      end if;
    end if;

    // Call the recursive step starting at p^{r-1}
    r := Valuation(n,p);
    if n eq p^r then
      crp := rec< CRP | K := K, p := p, rpolyg := [<1,j>,<n,0>], lowerv := L, fixedv := fixedv, points:=points>;
    else
      if points then
        rpolyg := [<1,j>,<n,0>] cat [<i,0> : i in [n-1..p^r+1 by -1] | Valuation(Binomial(n,i),p) eq 0] cat [<p^r,0>];
        crp := rec< CRP | K := K, p := p, rpolyg := rpolyg, lowerv := L, fixedv := fixedv, points:=true>;
      else  
        crp := rec< CRP | K := K, p := p, rpolyg := [<1,j>,<n,0>,<p^r,0>], lowerv := L, fixedv := fixedv, points:=false>;
      end if;
    end if;
    crps := AllRamificationPolygonsSub(crp,Valuation(n,p)-1);

    // Sort the polygons for presentation
    if output eq "polygons" then
      polygons := [];
      for i in [1..#crps] do
        Sort(~crps[i]`rpolyg);
        Append(~polygons,crps[i]`rpolyg);
      end for;
      return polygons;//crps;
    elif output eq "crps" then
      for i in [1..#crps] do
        Sort(~crps[i]`rpolyg);
        Sort(~crps[i]`fixedv);
      end for;
      return crps;
    end if;
end intrinsic;


//
//  DIAGNOSTIC FUNCTIONS
//


CheckCrp := function(crp)
  rp := RamificationPolygonSeq(crp`lowerv,crp`p);
  return Vertices(rp) eq crp`rpolyg;
end function;

// Create Example Polynomial f with a given ramification polygon (input as a CRP)
// Input: crp - a crp to create an example of
// Output: list of coefficients of f (ready to be coerced into a PolynomialRing
ExamplePolyFromCRP := function(crp)
    L := [1] cat crp`lowerv cat [0];
	return [crp`p ^ i : i in L];
end function;



//////////////////////////////////////////////////////////////////
// load "resseg-enum.m";
//////////////////////////////////////////////////////////////////

//
// RESSEG-ENUM.M - Enumerate possible residual polynomials of segments
//

// AllResidualPolynomials
// INPUT:  k   - a p-adic field
//         R   - a ramification polygon
//       phi_0 - constant term of generating polynomial ( can this be made optional ? only needed for linked coeffs )
// OUTPUT: A   - List of possible [_A_] given inputs

intrinsic AllResidualPolynomials(k,R,phi_0) -> .
{INPUT:  k is a p-adic field, R is a ramification polygon, phi_0 is the constant term of generating polynomial.
OUTPUT: A list of possible of representatives of residual polynomial classes [_A_] given inputs.}
    //print Valuation(k!phi_0);
    if Valuation(k!phi_0) ne 1 then
        error "Valuation of phi_0 must be 1";
    end if;
    
    // PHASE ZERO - INITIALIZE DATA
    
    K := FieldOfFractions(k);  // in case we pass a pAdicRing to the function
    ok := RingOfIntegers(K);
    rk := ResidueClassField(ok);
    rkz := PolynomialRing(rk);
    n := R[#R][1];
    
    phi_0 := K!(phi_0);
    
    // Initialize a and b lists
    a := []; b := [];
    for pt in R do
        q,r := Quotrem(pt[2],n);
        Append(~a,q);
        Append(~b,r);
    end for;
    
    // PHASE ONE - ASSIGN RESIDUES TO EACH POINT OF THE POLYGON
    
    // Initialize Aijs using the first point
    if b[1] eq 0 then
        // First point is linked to monic leading term
        binomt := Binomial(n,R[1][1]);
        deltap := ( K!(-phi_0) )^(-a[1]);
        Aijs := [[rk!(K!( binomt * deltap ))]];
    else
        // First point is free
        Aijs := [[i] : i in rk | i ne 0];
    end if;
    
    while Min([#A : A in Aijs]) lt #R do
        // Pop A of the front of Aijs
        A := Aijs[1];
        Aijs := Aijs[2..#Aijs];
        
        // find the index for the next point
        s := #A+1;

        // If this is a linked coeff, assign it's value,
        // Else it is free.
        if b[s] eq 0 then
            if R[s][1] eq n then
                // This is the monic leading term
                Append(~Aijs,Append(A,1));
            else
                // This term is linked to the monic leading term
                //print "linked to leading term";
                //binomt := Binomial(n,R[s][1])^(-1);  // WRONG
                //deltap := ( (-phi_0) )^(a[s]);
				binomt := Binomial(n,R[s][1]);         // Corrected
                deltap := ( K!(-phi_0) )^(-a[s]);
                Append(~Aijs,Append(A,rk!(k!( binomt * deltap ))));
            end if;
        elif b[s] in b[1..s-1] then
            // This term is linked to a previous term
            // Bin(b,pst)^(-1) * Bin(b,psq) * (d0)^(at-aq) * Aq  // THIS IS WRONG
            // Bin(b,pst) * Bin(b,psq)^(-1) * (d0)^(aq-at) * Aq  // Corrected		
            //print "linked coeff";
            at := a[s];
            aq := a[Index(b,b[s])];
            // binomt := Binomial(b[s],R[s][1])^(-1);  // WRONG
            // binomq := Binomial(b[s],R[Index(b,b[s])][1]);
            // deltap := ( (-phi_0) )^(at-aq);
            binomt := Binomial(b[s],R[s][1]);          // Corrected
            binomq := Binomial(b[s],R[Index(b,b[s])][1])^(-1);
            deltap := ( (-phi_0) )^(aq-at);
            LinkedA := rk!(k!( binomt * binomq * deltap )) * A[Index(b,b[s])];
            //print binomt,binomq,deltap,A[Index(b,b[s])];
            Append(~Aijs,Append(A,LinkedA));
        else
            // This term is free
            //print "free coeff";
            for i in [i : i in rk | i ne 0] do
                Append(~Aijs,Append(A,i));
            end for;
        end if;
    end while;

    //print "Aijs",Aijs;
    
    // PHASE TWO - CONSTRUCT POLYNOMIALS USING THOSE RESIDUES
    
    //S := [1] cat [R[i][1] : i in [2..#R-1] | allslopes[i] ne allslopes[i-1]] cat [R[#R][1]];
    allslopes := [ -( (R[i+1][2]-R[i][2]) / (R[i+1][1]-R[i][1]) ) : i in [1..#R-1]];
    Sindex := [1] cat [i : i in [2..#R-1] | allslopes[i] ne allslopes[i-1]] cat [#R];
    S := [R[i][1] : i in Sindex];
    slopes := [allslopes[i] : i in Sindex[1..#Sindex-1]];
    h := [-Numerator(s) : s in slopes];
    e := [Denominator(s) : s in slopes];

    //print "Sindex",Sindex;
    //print "S",S;
    //print slopes;

    Aflat := [];
    for aij in Aijs do
        ThisA := <>;
        for i in [1..#S-1] do
            // RKz![res(rhol[Integers()!(j*e[r]+rpv[r][1])] div (K.1)^(-Integers()!(j*h[r]-rpv[r][2]))) : j in [0..Integers()!(d[r]/e[r])]]);
            deg := Integers()!((S[i+1]-S[i])/e[i]);
            respoly := [Zero(rk) : cnt in [0..deg]];
            //for j in segs[i] do
            for j in [Sindex[i]..Sindex[i+1]] do
                term := Integers()!( (R[j][1]-S[i])/e[i] );
                //print "term",<i,j>,aij[j],"* z ^",term;
                respoly[term+1] := rk!aij[j];
            end for;
            Append(~ThisA,Polynomial(rk,respoly));
        end for;
        Append(~Aflat,ThisA);
    end for;
    
    return Aflat;
    
    // PHASE THREE - PARTITION POLYNOMIAL SETS INTO CLASSES  --- NOT IMPLEMENTED
    
end intrinsic;




//////////////////////////////////////////////////////////////////
// used to hold polynomials and invariants
//////////////////////////////////////////////////////////////////
FRP := recformat<
K : RngPad,       // base field
deg,              // degree pf extensions
j,                // discriminant exponent j
rampol,           // ramification polygon
phi0,             // constant coefficient modulo pi^2
A,                // residual polynomials
generators        // gernerating polynomials of extensions
>;

intrinsic DistinctConstantCoefficients(K::RngPad,n::RngIntElt) -> .
{The possible constant coefficients of extensions of K of degree n, that
yield distinct extensions}
    //K is base field, n is degree of polynomial.  n can be recovered from polygon.

    p:=Prime(PrimeRing(K));

    m:=Valuation(n,p);

    if (n eq p^m) then
        return {p};
    else
        _K_,vk:=ResidueClassField(K);    //vk is Mapping from: RngPad: K to FldFin: _K_

        Kx,vkx:=UnitGroup(_K_);        //vkx is Mapping from: GrpAb: Kx to FldFin: _K_

        if IsCyclic(Kx) then

            a:=Setseq(Generators(Kx))[1];
            den:=sub<Kx|n*a>;

            S,vS:=quo<Kx |den>;      //vS is Mapping from: GrpAb: Kx to GrpAb: S.  

            //Correct to here.

            Delta:={(vkx(s@@vS))@@vk  :s in S};   //Map each element of S to K (through _K_ and Kx)

            return {delta*p: delta in Delta};
        else
            
            error "Error: Something has gone terribly wrong.";
        end if;
        
    end if;        
end intrinsic;


Minval := function(k,R)
    // Find the minimum valuations prescribed by the ramification polygon R

    n := R[#R][1];
    K  := RingOfIntegers(k);
    p := Prime(K);
    xes := [pt[1] : pt in R];

    Lpt := [lfunc(Valuation(K!pt[1]),pt[2],K,n) : pt in R];
    for s in [1..Valuation(n,p)] do
        if p^s notin xes then
            left := #[x : x in xes | x lt p^s];
            Append(~Lpt,lnvfunc(s,R[left],R[left+1],K,n:colinear:=true));
        end if;
    end for;
    minval := [];
    for i in [1..n-1] do
        Append(~minval,Max([1] cat [l[i] : l in Lpt]));
    end for;
    return minval;
end function;

intrinsic IsomorphicInFields(f,Ls) -> .
{True if the extension generated by f is ismorphic to an extension in L}
  for L in Ls do
    Lz<z> := PolynomialRing(L);
    if HasRoot(Lz!f) then 
      return true;
    end if;
  end for;
  return false;    
end intrinsic;


intrinsic AllTotallyRamifiedExtensions(k::RngPad,R::SeqEnum,A,delta_0:CountOnly:=false,TemplateOnly:=false,want_filter:=true) -> .
{All extension of pi-adic ring or field k with give invariants:
R - points of a ramification polygon, as [<x,y>], 
A - list of residual polynomials for segments of R,
delta_0 - a representative of a class in RK*/(RK*)^n
}
    n  := R[#R][1];
    K  := RingOfIntegers(k);
    Kx := PolynomialRing(K);
    BK := BaseRing(K);
    
    // integral prime
    p := Prime(K);
    pi := UniformizingElement(K);

    // Define Resdidue Class Field of K
    RK := ResidueClassField(K);
    RKz:= PolynomialRing(RK);
    
    // Set up vector space to find quotient with images of Sm
    k_is_base := K eq BK;
    if not k_is_base then
        RBK := ResidueClassField(BK);
        V,Vm := VectorSpace(RK,RBK);
    end if;
    
    // Take R, and find xes[], a[], b[] s.t. <x,y> eq <x,an+b>
    a := []; b := [];
    xes := [pt[1] : pt in R];
    for pt in R do
        q,r := Quotrem(pt[2],n);
        Append(~a,q);
        Append(~b,r);
    end for;
    
    // Find the denominators of the slopes of segments
    allslopes := [ -( (R[i+1][2]-R[i][2]) / (R[i+1][1]-R[i][1]) ) : i in [1..#R-1]];
    Sindex := [1] cat [i : i in [2..#R-1] | allslopes[i] ne allslopes[i-1]] cat [#R];
    S := [R[i][1] : i in Sindex];
    slopes := [allslopes[i] : i in Sindex[1..#Sindex-1]];
    //slopes := [-((R[Index(xes,S[i+1])][2]-R[Index(xes,S[i])][2]) / (S[i+1]-S[i])) : i in [1..#S-1]];
    h := [Numerator(s) : s in slopes];
    e := [Denominator(s) : s in slopes];
    
    // Find the minimum valuations prescribed by the RP
    minval := Minval(K,R);

    // 1. Set Krasner Bound
    c := Ceiling(1+2*a[1]+(2*b[1])/n)-1;   // J_0 is a[1]n+b[1] (magma indexes from 1)

    // Precompute n*Hasse-Herbrand of RP
    m := 1;
    nhhr := [Min([pt[2] : pt in [<pt[1],pt[2]+pt[1]> : pt in R]])];
    while nhhr[#nhhr] lt n*c-1 do
        m +:= 1;
        Append(~nhhr,Min([pt[2] : pt in [<pt[1],pt[2]+m*pt[1]> : pt in R]]));
    end while;
    
    // 2. Initialize Template
    tau := [];
    for i in [1..n] do
        Append(~tau,[{Zero(RK)} : j in [1..c]]);
     end for;
    
    // 3. Set Free Choices (respecting ell)
    free_choices := [];
    //reducedzeros := ([<hhm mod n,(n-(hhm mod n)+hhm)/n> : hhm in nhhr | (n-(hhm mod n)+hhm) mod n eq 0]);
    // correcting for indexing from 1:
    reducedzeros := ([<(hhm mod n)+1,(n-(hhm mod n)+hhm)/n> : hhm in nhhr | (n-(hhm mod n)+hhm) mod n eq 0]);
    for i in [1..n] do
        for j in [1..c] do
            if <i,j> notin reducedzeros and (i eq 1 or j ge minval[i-1]) then
                tau[i][j] := Set(RK);
                if i gt 1 and j gt 1 then
                    Append(~free_choices,<i,j>);
                end if;
            end if;
        end for;
    end for;
    
    // 4. Set coefficients by Sm
    non_surjective_m := [];
    non_surjective_cell := [];
    surjective_s_places := [];
    for m in [1..Floor((R[1][2]-R[2][2])/(R[2][1]-1))] do
        i := nhhr[m] mod n;
        if (n-i+nhhr[m]) mod n eq 0 then
            j := (n-i+nhhr[m]) div n;
            if m notin slopes then
                tau[i+1][j] := {Zero(RK)};
                Append(~surjective_s_places,<i+1,j>);
            else
                image := {Evaluate(Parent(A[1]).1^S[Index(slopes,m)]*A[Index(slopes,m)],gamma) : gamma in RK};
                if image eq Set(RK) then
                    //print m,"is in slopes, Sm is surjective, so",<i+1,j>,"gets 0.",image;
                    tau[i+1][j] := {Zero(RK)};
                    Append(~surjective_s_places,<i+1,j>);
                else
                    //m,Index(slopes,m),slopes,S[Index(slopes,m)],A[Index(slopes,m)];
                    //print m,"is in slopes, Sm not surjective!",image;
                    if k_is_base then
                        tau[i+1][j] := Set(RK);
                    else
                        q,qm := quo<V|[Vm(i):i in image]>;
                        tau[i+1][j] := {Inverse(Vm)(Inverse(qm)(i)):i in q};
                    end if;
                    Append(~non_surjective_m,m);
                    Append(~non_surjective_cell,<i+1,j>);
					// FIX ME: Record the cokernal for possible reduction
                end if;
            end if;
        end if;
    end for;
    
    // 5. Set coefficients by A
    for i in [1..#R] do
        if b[i] ne 0 then
            seg := #[cnt : cnt in S | cnt le R[i][1]];
            j := a[i] + 1 - Valuation(Binomial(b[i],R[i][1]),p);
            tau[b[i]+1][j] := { Coefficient(A[seg],(R[i][1]-S[seg]) div e[seg]) * (-delta_0)^(a[i]+1) * RK!(Binomial(b[i],xes[i])^(-1) * pi^(Valuation(K!(Binomial(b[i],xes[i]))))) };
            // Remove from list of free choices
			Exclude(~free_choices,<b[i]+1,j>);
        end if;
    end for;
    
    // 6. Set tau_{0,1} to delta_0
    tau[1][1] := {delta_0};
    if TemplateOnly then
       return tau;
    end if;

    // 7. Construct list of polynomials
    if CountOnly then
        L := &*[ &*[#cell : cell in row] : row in tau];
    else
        cp := CartesianProduct([CartesianProduct(i):i in tau]);
        L := [Kx.1^n + &+[ Kx.1^(i-1)*&+[K!tup[i][j]*pi^j : j in [1..#tup[i]]] : i in [1..#tup]] : tup in cp];
    end if;

    // 7.5 Postprocess if reduction may be needed
    need_filter := true;
    if #non_surjective_cell eq 0 then
        need_filter := false;
    elif #non_surjective_cell eq 1 then
        if &and[f lt non_surjective_cell[1] : f in free_choices] then
            need_filter := false;
        end if;
    end if;
    //if #free_choices gt 0 and #non_surjective_m gt 0 then
        // Postprocess with reduction
		//print "free_choices", free_choices;
		//print "non_surjective_m",non_surjective_m;
        //print "surjective_s_places",surjective_s_places;
    //end if;

    // 8. Return list of polynomials
    if CountOnly then
        return L, need_filter;
    end if;
    if need_filter and want_filter then
      Ks := [TotallyRamifiedExtension(K,L[1])];
      fs := [L[1]];
      "Filtering";
      for f in L do
        if not IsomorphicInFields(f,Ks) then
          Append(~Ks,TotallyRamifiedExtension(K,f));
          Append(~fs,f);      
        end if;
      end for;
      return fs;
    else  
      return L;
    end if;
end intrinsic;

function print_template_old(tau)
    rtau := Reverse(tau);
    c := #rtau[1];
    //print "   " cat &cat[" x^" cat IntegerToString(i) : i in [#tau..0 by -1]];
    str := "   " cat &cat[" x^" cat IntegerToString(i) : i in [#tau..0 by -1]];
    for i in [c..1 by -1] do
        //print "p^" cat IntegerToString(i) cat "  0 " cat &cat[" " cat t[i]:t in rtau];
        str cat:= "\n" cat "p^" cat IntegerToString(i) cat "  0 " cat &cat[" " cat t[i]:t in rtau];
    end for;
    return str;
end function;

function print_template(tau:output:="text")
    rtau := Reverse(tau);
    c := #rtau[1];
    if output eq "text" then
        //print "   " cat &cat[" x^" cat IntegerToString(i) : i in [#tau..0 by -1]];
        str := "   " cat &cat[" x^" cat IntegerToString(i) : i in [#tau..0 by -1]];
        for i in [c..1 by -1] do
            //print "p^" cat IntegerToString(i) cat "  0 " cat &cat[" " cat t[i]:t in rtau];
            str cat:= "\n" cat "p^" cat IntegerToString(i) cat "  0 " cat &cat[" " cat t[i]:t in rtau];
        end for;
    elif output eq "latex" then
        str := "\\begin{center}\n\\small\n\\begin{tabular}{l|cc";
        str cat:= &cat["c":i in [1..#tau]];
        str cat:= "}\n";
        str cat:= "    " cat &cat["&$x^{" cat IntegerToString(i) cat "}$ " : i in [#tau..0 by -1]];
        str cat:= "\\\\\\hline";
        for i in [c..1 by -1] do
            str cat:= "\n" cat "$p^" cat IntegerToString(i) cat "$&$\\{0\\}$ " cat &cat["&$\\{" cat t[i] cat "\\}$ ":t in rtau] cat "\\\\";
        end for;
        str cat:= "\n\\end{tabular}\n\\end{center}";
        return str;
    else
        str := "Invalid output parameter";
    end if;
    return str;
end function;

intrinsic AllTotallyRamifiedExtensionsDemo(K,R:verbose:=false,output:="text") -> .
{}
    // INPUT:
    //     K - a pi-adic field
    //     R - points of a ramification polygon, as [<x,y>]

    // degree of polynomial
    n := R[#R][1];
    
    // integral prime
    p := Prime(K);
    pi := UniformizingElement(K);

    // Precompute n*Hasse-Herbrand of RP
    //nhhr := [Min([pt[2] : pt in rpolyplusm]) : rpolyplusm in [[<pt[1],pt[2]+m*pt[1]> : pt in R] : m in [1..2*n] ]];

    // Take R, and find xes[], a[], b[] s.t. <x,y> eq <x,an+b>
    a := []; b := [];
    xes := [pt[1] : pt in R];
    for pt in R do
        q,r := Quotrem(pt[2],n);
        Append(~a,q);
        Append(~b,r);
    end for;
    
    // Find the denominators of the slopes of segments
    allslopes := [ -( (R[i+1][2]-R[i][2]) / (R[i+1][1]-R[i][1]) ) : i in [1..#R-1]];
    Sindex := [1] cat [i : i in [2..#R-1] | allslopes[i] ne allslopes[i-1]] cat [#R];
    S := [R[i][1] : i in Sindex];
    slopes := [allslopes[i] : i in Sindex[1..#Sindex-1]];
    //slopes := [-((R[Index(xes,S[i+1])][2]-R[Index(xes,S[i])][2]) / (S[i+1]-S[i])) : i in [1..#S-1]];
    h := [Numerator(s) : s in slopes];
    e := [Denominator(s) : s in slopes];
    
    // Find the minimum valuations prescribed by the RP
    minval := Minval(K,R);

    // 1. Set Krasner Bound
    c := Ceiling(1+2*a[1]+(2*b[1])/n)-1;   // J_0 is a[1]n+b[1] (magma indexes from 1)

    // Precompute n*Hasse-Herbrand of RP
    m := 1;
    nhhr := [Min([pt[2] : pt in [<pt[1],pt[2]+pt[1]> : pt in R]])];
    while nhhr[#nhhr] lt n*c-1 do
        m +:= 1;
        Append(~nhhr,Min([pt[2] : pt in [<pt[1],pt[2]+m*pt[1]> : pt in R]]));
    end while;
    
    // 2. Initialize Template
    if verbose then print "2. Init Template"; end if;
    taustr := [];
    for i in [1..n] do
        Append(~taustr,[" 0 " : j in [1..c]]);
     end for;
     if verbose then print print_template(taustr:output:=output) cat "\n"; end if;
    
    // 3. Set Free Choices (respecting ell)
    if verbose then print "3. Set free choices (respecting ell)"; end if;
    reducedzeros := ([<(hhm mod n)+1,(n-(hhm mod n)+hhm)/n> : hhm in nhhr | (n-(hhm mod n)+hhm) mod n eq 0]);
    //print "nhhr",nhhr;
    //print "reducedzeros",reducedzeros;
    for i in [1..n] do
        for j in [1..c] do
            //print <i,j>,<i,j> notin reducedzeros,(i eq 1 or j ge minval[i-1]);
            if <i,j> notin reducedzeros and (i eq 1 or j ge minval[i-1]) then
                taustr[i][j] := "R_K";
            end if;
        end for;
    end for;
    if verbose then print print_template(taustr:output:=output) cat "\n"; end if;
    
    // 4. Set coefficients by Sm
    if verbose then print "4. Set coeffs by Sm"; end if;
    for m in [1..Floor((R[1][2]-R[2][2])/(R[2][1]-1))] do
        i := nhhr[m] mod n;
        if (n-i+nhhr[m]) mod n eq 0 then
            j := (n-i+nhhr[m]) div n;
            //m,<i,j>;
            taustr[i+1][j] := "S_" cat IntegerToString(m);
        end if;
    end for;
    if verbose then print print_template(taustr:output:=output) cat "\n"; end if;
    
    // 5. Set coefficients by A
    if verbose then print "5. Set coeffs by A"; end if;
    for i in [1..#R] do
        if b[i] ne 0 then
            seg := #[cnt : cnt in S | cnt le R[i][1]];
            //seg := 1;
            j := a[i] + 1 - Valuation(Binomial(b[i],R[i][1]),p);
            //print R[i],"(",a[i],"n +",b[i],") sets f_{",b[i],",",j,"} to A_{",seg,",",R[i][1]-S[seg],"/",e[seg],"}";
            taustr[b[i]+1][j] := "A" cat IntegerToString(seg) cat IntegerToString((R[i][1]-S[seg]) div e[seg]);
        end if;
    end for;
    if verbose then print print_template(taustr:output:=output) cat "\n"; end if;
    
    // 6. Set tau_{0,1} to delta_0
    if verbose then print "6. Set f_{0,1} to delta_0"; end if;
    taustr[1][1] := "d_0";
    return print_template(taustr:output:=output);
end intrinsic;

NumberOfExtensionsFromPoly := function(f)
    return Degree(f) / #Roots(Polynomial(ext<CoefficientRing(f)|f>,f));
end function;

IsExtensionIsomorphic := function(f,g)
    return HasRoot(g,ext<CoefficientRing(f)|f>);
end function;

intrinsic AllTotallyRamifiedExtensions(k::RngPad,n::RngIntElt:j:=0,phi0:=0,invariants:=false) -> .
{All extension of pi-adic ring or field k with given invariants:
R - points of a ramification polygon, as [<x,y>], 
A - list of residual polynomials for segments of R,
phi0 - the constant coefficient mod pi^2
}
    K     := RingOfIntegers(k);
    Kx<x> := PolynomialRing(K);
   
    // integral prime
    p := Prime(K);
    pi := UniformizingElement(K);

    // Define Resdidue Class Field of K
    RK,res := ResidueClassField(K);
    RKz<z> := PolynomialRing(RK);
    
    "AllRamificationPolygons"; 
    if j ne 0 then
        RR := AllRamificationPolygons(K,n,j:points:=true);
    else
        RR := &cat[AllRamificationPolygons(K,n,j:points:=true) : j in PossibleDiscriminants(k,n)];
    end if;
    " -- done";
    "AllResidualPolynomials";
    if phi0 ne 0 then
       AA := <AllResidualPolynomials(k,R,phi0) : R in RR>;
    //E := [[&cat+[AllTotallyRamifiedExtensions(k,RR[i],a,res(phi_0/Prime(k))) : a in A] : A in AA[i]] : i in [1..#RR]];
      " -- done";
       if invariants then
         E := Flat([[rec<FRP|K:=K,deg:=n,j:=RR[i][1][2],rampol:=RR[i],A:=a,generators:=AllTotallyRamifiedExtensions(k,RR[i],a,res(phi0/Prime(k)))> : a in AA[i]] : i in [1..#RR]]);
       else
         E := [&cat[AllTotallyRamifiedExtensions(k,RR[i],a,res(phi0/Prime(k))) : a in AA[i]] : i in [1..#RR]];
      end if;
    else
      phi0s := DistinctConstantCoefficients(K,n);      
    //E := [[&cat+[AllTotallyRamifiedExtensions(k,RR[i],a,res(phi_0/Prime(k))) : a in A] : A in AA[i]] : i in [1..#RR]];
      " -- done";
       if invariants then
         E := Flat([[[rec<FRP|K:=K,deg:=n,j:=RR[i][1][2],rampol:=RR[i],phi0:=phi0,A:=a,generators:=AllTotallyRamifiedExtensions(k,RR[i],a,res(phi0/Prime(k)))> : a in AllResidualPolynomials(k,RR[i],phi0)] : phi0 in phi0s]:i in [1..#RR]]);
       else
         E := [[&cat[AllTotallyRamifiedExtensions(k,RR[i],a,res(phi0/Prime(k))) : a in AllResidualPolynomials(k,RR[i],phi0)] :phi0 in phi0s]: i in [1..#RR]];
      end if;
    end if;
    return E;
end intrinsic;

