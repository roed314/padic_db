
intrinsic pdiscf(f,p) -> Any
{ Valuation of the field discriminant at a prime p}

  v:= Valuation(Discriminant(pMaximalOrder(EquationOrder(f),p)), p);
  return v;
end intrinsic;

intrinsic slope(p1,p2) -> Any
{ Slope between two points }

  return (p2[2]-p1[2])/(p2[1]-p1[1]);
end intrinsic;

intrinsic visible(k, p) -> Any
{ Get the list of visible slopes of k for the prime p }

  kk := NumberField(k);
  subs := Subfields(kk); // Includes k, but not Q
  subs := [z[1] : z in subs];
  points := [[Degree(z), pdiscf(z,p)] : z in subs];
  points := [[1,0]] cat points;
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

intrinsic merge_slopes(sl1, sl2) -> Any
{ Merge 2 multisets } 
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

intrinsic DelSpaces(s::MonStgElt) ->MonStgElt
  {Delete spaces from a string s}
  return &cat([x: x in Eltseq(Sprint(s)) | (x ne " ") and (x ne "\n")]);
end intrinsic;

intrinsic polredabs(f::Any) -> Any
  {Have gp compute polredabs}
  if Degree(f) gt 1000 then
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

intrinsic readslopecache() -> Any
  {Read the stored slope information}
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
  return slopec;
end intrinsic;

intrinsic savefield(f, sl)
  {Save slope data to a file}
  ff:= DelSpaces(Sprint(f));
  sll:= Sort([z : z in sl]);
  sll:= DelSpaces(Sprint(sll));
  write("Slopedata", Sprintf("%o %o", ff, sll) : rewrite:= false);
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
intrinsic get_slopes(f, p) -> Any
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
