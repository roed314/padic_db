// Assumes the following are added as attributes of Q
// slopecache = associative array from file
// galdat = the triple from computing the Galois group
// dat16 = the degree 16 subfields of the Galois closure

intrinsic pdiscf(f,p) -> Any
{ Valuation of the field discriminant at a prime p}
  //v:= Valuation(Discriminant(pMaximalOrder(EquationOrder(f),p)), p);
  RR<x> := PolynomialRing(Integers());
  f := RR ! f;
  u,v:=pDiscriminant(f,p);
  return v;
end intrinsic;

intrinsic slope(p1,p2) -> Any
{ Slope between two points }
  return (p2[2]-p1[2])/(p2[1]-p1[1]);
end intrinsic;

intrinsic addone2slope(sl1, totdisc) -> Any
{ If we know all but one slope and disc of the Galois closure, add the last one }
  kk:=1;
  going := true;
  sl := Sort(MultisetToSequence(sl1));
  while going do
    s := &+[sl[j]*2^((j ge kk) select j else j-1) : j in [1..#sl]];
    newval := (totdisc-s)/2^(kk-1);
    if newval le sl[kk] and (kk eq 1 or newval ge sl[kk-1]) then
      going := false;
    else
      kk := kk+1;
    end if;
  end while;
  return sl1 join {* newval *};
end intrinsic;

intrinsic DelSpaces(s::MonStgElt) ->MonStgElt
  {Delete spaces from a string s}
  return &cat([x: x in Eltseq(Sprint(s)) | (x ne " ") and (x ne "\n")]);
end intrinsic;

intrinsic polredabs(f::Any) -> Any
  {Have gp compute polredabs}
  if Degree(f) gt 17 then
    return f;
  end if;
  R<x>:=PolynomialRing(Rationals());
  out := Sprintf("/tmp/polredabs%o.out", Random(10^30));
  txt := Sprintf("/tmp/polredabs%o.txt", Random(10^30));
  //f:=R!f * Denominator(VectorContent(Coefficients(f)));
  // Avoid hardwiring gp path
  write(txt,Sprintf("polredabs(%o)",f): rewrite:=true);
  System("which sage>"*out);
  gppath:= DelSpaces(Read(out));
  System("rm "* out);
  System(gppath*" -gp -f -q --default parisizemax=1G <"*txt*">"*out);
  //try
  f:=eval DelSpaces(Read(out));
  //catch e;
  //end try;
  System("rm "* out);
  System("rm "* txt);
  return f;
end intrinsic;


intrinsic dovisible(points, p) -> Any
{ Get the list of visible slopes of k for the prime p }

  Sort(~points);
  hull := [points[1]];
  cur := hull[1][1];
  for j:=2 to #points do
    if points[j][1] gt cur then
      Append(~hull, points[j]);
      cur:=points[j][1];
    end if;
  end for;
  cur := 1;
  sl := [];
  while cur le #hull -1 do
    pos := cur+1;
    msl := slope(hull[cur], hull[pos]);
    k := cur + 2;
    while k lt #hull do
      if slope(hull[cur], hull[k]) le msl then
        pos := k;
        msl := slope(hull[cur], hull[pos]);
      end if;
      k := k+1;
    end while;
    if msl gt -1 then
      degdiff:=Valuation(hull[pos][1]/hull[cur][1], p);
        for jj:=1 to degdiff do
          Append(~sl, slope(hull[cur],hull[pos]));
        end for;
    end if;
    cur := pos;
  end while;

  return sl;
end intrinsic;

intrinsic visible(k, p) -> Any
{ Get the list of visible slopes of k for the prime p }
  kk := NumberField(k);
  subs := Subfields(kk); // Includes k, but not Q
  subs := [DefiningPolynomial(z[1]) : z in subs];
  points := [[Degree(z), pdiscf(z,p)] : z in subs];
  points := [[1,0]] cat points;
  return dovisible(points, p);
end intrinsic

intrinsic minires(f,n,G,S) -> Any
{ my res without recomputing }
  order := Order(G) / n;
  R:=PolynomialRing(Integers());
  if order in Integers() then
    order := Integers() ! order;
    ss:= Subgroups(G : OrderEqual:=order);
    return <R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ss>;
  else
    return <>;
  end if;
end intrinsic;

intrinsic try16(f,p, looking) -> Any
{ Try to get slopes using just degree looking subfields.
  Assumes Gal data is in Q
  Puts degree 16 resolvents in Q if looking == 16}
  nf := NumberField(f);
  Q:= Rationals();
  galdat := Q`galdat;
  G,r,S := Explode(galdat);
  OG := Order(G);
  OGv := Valuation(OG,2);
  R<x>:=PolynomialRing(Rationals());
  vis := visible(f,p);
  slopes := Multiset(vis);
  neworder:= Integers() ! (OG/looking);
  ss:= Subgroups(G : OrderEqual:=neworder);
  "Got ", #ss, " subgroups";
  pols := <R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ss>;
  "Got pols";
  if looking eq 16 then
    polspra := <polredabs(z) : z in pols>;
    "Did polredabs";
    Q`dat16 := polspra;
    pols := polspra;
  end if;
  for po in pols do
    mysl := mini_gs(po,p);
    //"Adding ", mysl;
    slopes := merge_slopes(slopes, mysl);
    if #slopes eq OGv then 
      return slopes;
    end if;
  end for;
  return slopes;
end intrinsic;

intrinsic allslopes(f, p) -> Any
{ Get the list of slopes of the field defined by f and the prime p.
  This assumes that the Galois group is a 2-group}
  R:=PolynomialRing(Integers());
  G,r,S := GaloisGroup(f);
  Q :=Rationals();
  Q`galdat := <G, r, S>;
  slopes16 := try16(f,p,16);
  OG := Order(G);
  "Group order ", OG, " group ", TransitiveGroupIdentification(G), " with ", #slopes16, " slopes: ", slopes16;
  if #slopes16 eq Valuation(OG,2) then
    "Done at 16";
    return slopes16;
  elif OG lt 1050 and #slopes16 eq Valuation(OG,2)-1 then
    ns1:=NormalSubgroups(G : OrderEqual:=1); // Trivial subgroup
    "Got ns";
    po:=<R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ns1>;
    "Got pol";
    totdisc:=pdiscf(po[1],p);
    "Got disc";
    return addone2slope(slopes16, totdisc);
  else
    "Try 32";
    slopes32 := try16(f,p,32);
    if #slopes32 eq Valuation(OG,2) then
      return slopes32;
    end if;

    "Need ns";
    ns:=NormalSubgroups(G);
    po:=<R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ns>;
    points:=[[Degree(z), pdiscf(z,p)] : z in po];
    return dovisible(points, p);
  end if;
end intrinsic;

intrinsic merge_slopes(sl1, sl2) -> Any
{ Merge 2 multisets } 
  sl1 := {* Rationals()! z : z in sl1 *};
  sl2 := {* Rationals()! z : z in sl2 *};
  all:=MultisetToSet(sl1) join MultisetToSet(sl2);
  return {* z^^Max(Multiplicity(sl1, z), Multiplicity(sl2,z)) : z in all *};
end intrinsic;

intrinsic write(filename::MonStgElt,str::MonStgElt: console:=false, rewrite:=false)
  {Write str to file filename as a line
   rewrite:= true means we overwrite the file, default is to append to it
   console:= true means we echo the string as well.
   If the filename is the empty string, don't write it.}
  if console then str; end if;
  if filename ne "" then
    F:=Open(filename,rewrite select "w" else "a");
    WriteBytes(F,[StringToCode(c): c in Eltseq(Sprint(str)*"\n")]);
    Flush(F);
  end if;
end intrinsic;

intrinsic res(f::Any,n::Any) -> Any
  { resolvents of degree n}
  G,r,S := GaloisGroup(f);
  order := Order(G)/n;
  R<x> := PolynomialRing(Integers());
  if order in Integers() then
    order := Integers() ! order;
    ss:= Subgroups(G : OrderEqual:=order);
    return(<R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ss>);
  else
    return <>;
  end if;
end intrinsic;

intrinsic readslopecache() -> Any
  {Read the stored slope information}
  Q:=Rationals();
  a,b := HasAttribute(Q, "slopecache");
  if a then return b; end if;

  slopec:= AssociativeArray();
  try
    slopecstr:=Split(Read("Slopedata"));
  catch e;
    return slopec;
  end try;
  R<x>:=PolynomialRing(Rationals());
  for pdat in slopecstr do
    sllist:=Split(pdat, " ");
    slopec[eval(sllist[1])] := eval(sllist[2]);
  end for;
  Q`slopecache := slopec;
  return slopec;
end intrinsic;

intrinsic savefield(f, sl)
  {Save slope data to a file}
  slopec := readslopecache();
  if not IsDefined(slopec, f) then
    ff:= DelSpaces(Sprint(f));
    sll:= Sort([z : z in sl]);

    // Update internal slope cache
    slopec[f] := sll;
    Q := Rationals();
    Q`slopecache := slopec;

    sll:= DelSpaces(Sprint(sll));
    write("Slopedata", Sprintf("%o %o", ff, sll) : rewrite:= false);
  else
    "Already had ", f;
  end if;
end intrinsic;


intrinsic mini_gs(f,p) -> Any
{ Get slopes.  Lookup in dictionary, otherwise compute visible slopes.}
  slopec := readslopecache();
  if IsDefined(slopec, f) then
    return Multiset(slopec[f]);
  end if;
  return Multiset(visible(f,p));
end intrinsic;

/* Set for degree 16 fields */
intrinsic get_slopesold(f, p) -> Any
{ Naive function to get slopes for the field defined by f }
  slopec := readslopecache();
  if IsDefined(slopec, f) then
    return Multiset(slopec[f]);
  end if;
  nf := NumberField(f);
  G,r,S := GaloisGroup(f);
  OG := Order(G);
  R<x>:=PolynomialRing(Rationals());
  vis := visible(nf,p);
  slopes := Multiset(vis);
  looking:=8;
  while p^#slopes lt OG do
    looking:=looking*2;
    "Looking ", looking;
    neworder:= Integers() ! (OG/looking);
    ss:= Subgroups(G : OrderEqual:=neworder);
    pols := <R ! DefiningPolynomial(NumberField( GaloisSubgroup(S, z`subgroup))) : z in ss>;
    polspra := <polredabs(z) : z in pols>;
    allsl := <mini_gs(z,p) : z in pols>;
    for sl in allsl do
      slopes := merge_slopes(slopes, sl);
    end for;
  end while;

  savefield(f, slopes);
  return slopes;
end intrinsic;

intrinsic get_slopes(f, p) -> Any
{ Naive function to get slopes for the field defined by f }
  slopec := readslopecache();
  if IsDefined(slopec, f) then
    return Multiset(slopec[f]);
  end if;
  slopes:=Multiset(allslopes(f,p));
  Q := Rationals();
  secodecs := Q`dat16;
  OG:=Order(GaloisGroup(f));
  secodecs := [z : z in secodecs | Order(GaloisGroup(z)) eq OG];
  "Got ", #secodecs, " siblings";

  savefield(f, slopes);
  for z in secodecs do
    savefield(z, slopes);
  end for;
  return slopes;
end intrinsic;


/*

intrinsic topslope(f,p) -> Any
 { }
 sl:=visible(f,p);
 return sl[#sl];
end intrinsic;

intrinsic hasf2(f) -> Any
 { }
 sf := Subfields(NumberField(f), 2);
 sf := <z[1] : z in sf>;
 sf := <z : z in sf | pdiscf(z,2) eq 0 >;
 return #sf ne 0;
end intrinsic;

*/
