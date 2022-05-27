/* Input comes as
   [label, p, eisen, c, f]
   */

valsofcoefs(pol,p)=
{
  my(v,v1,m);
  v1=Vecrev(pol);
  v=vector(#v1);
  for(j=1,#v1,
    m = 10^10;
    if(v1[j]!=0, 
      g=Vec(v1[j]);
      g=gcd(g);
      g=valuation(g,p);
      if(g<m,m=g);
    );
    v[j] = m;
  );
  return(v);
}

/* discrim is the valuation of the discriminant of the eisenstein
   polynomial.
   This equals c/f
   */


indofsep(eis,p,discrim)=
{
  my(vals,n,nu,i0,itilde=List(),ii,tmp,tmp2);
  n = poldegree(eis);
  vals = valsofcoefs(eis, p);  /* valuations of the coefficients */
  nu = valuation(n,p);
  i0 = discrim-n+1;
  tmp=vector(n,i,n*vals[i+1]+i-n);
  for(j=0, nu,
    tmp2=List();
    for(i=1,#tmp, if(valuation(i,p) <= j, listput(tmp2,tmp[i])));
    mm=tmp2[1];
    for(k=2,#tmp2,mm=min(tmp2[k],mm));
    listput(itilde,mm);
  );
  ii = vector(nu+1,h,-1);
  ii[nu+1]=0;
  forstep(j=nu-1, 0, -1,
    ii[j+1] = min(itilde[j+1], ii[j+2]+n);
  );
  return(ii);
}
