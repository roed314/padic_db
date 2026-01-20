177,0
V,montestalk,3
A,FldNum,8,Discriminant,FactorizedDiscriminant,FactorizedPrimes,IndexPrimeFactors,IntegralBasis,LocalIndex,PrimeIdeals,TreesIntervals
S,AdaptPrecision,,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,Cancel,Cancell the powers of p in the numerator and denominator of the fraction poly/p^vden,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,312,148,-38,-38,-38,-38
S,CompensateValue,tree is an interval [i..j] inside [1..#K`PrimeIdeals[p]] and exponents is a sequence of integers of length #tree. The output is an element beta such that v_P(beta) >= exponents[P] for all P in the tree,0,4,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,Construct,"This routine constructs a polynomial Ppol with integer coefficients such that: deg Ppol<m_i+1 and y^nu*R_i(Ppol)(y)=respol(y), where nu=ord_y(respol). For non-negative integers s,u, the variable point=[s,u] is the left endpoint of a segment of slope -type[i]`slope supporting N_i(Ppol)",0,6,0,0,1,0,0,1,1,-38,,0,0,-1,,0,0,-1,,0,0,-1,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,ConvertLogs,The class mod P of the product p^log[1]Phi_1^log[2] ... Phi_i^log[i+1],0,3,0,0,1,0,0,1,1,-38,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,CrossValues,"Compute a matrix of cross values Mat[P,Q]=w_Q(phi_P) for all prime ideals P,Q in the input interval tree of K`PrimeIdeals[p]",0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,177,-38,-38,-38,-38,-38
S,CRT,Compute x congruent to elements[j] mod ideals[j] for every j. Integrality of the given elements is not checked!,2,0,1,82,0,28,1,1,82,0,270,2,0,0,0,0,0,0,0,82,,0,0,82,,28,-38,-38,-38,-38,-38
S,Different,Valuation of the different ideal of the local extension of Qp given by the p-adically irreducible polynomial represented by the given prime ideal P,0,2,0,0,1,0,0,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,EqualizeLogs,Add zeros to the shorter first list to achieve the same length as the second list,0,2,0,0,1,0,0,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,FaithfulpAdicConversion,,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,312,82,-38,-38,-38,-38
S,Generators,Compute the generators of the prime ideals in K above the rational prime p,0,2,0,0,1,0,0,0,0,-1,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,IndexOfCoincidence,"The index is 0 if t1,t2 belong to different trees. Otherwise, it is the least index such that the triplets (t1[i]`Phi,t1[i]`slope,t1[i]`psi) and (t2[i]`Phi,t2[i]`slope,t2[i]`psi) are different. The types must correspond to different prime ideals",0,4,0,0,1,0,0,0,0,-1,,1,1,-38,,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,IndexOfCoincidence,The index of coincidence of two different types t1 and t2,0,2,0,0,0,0,0,0,0,270,,0,0,270,,148,-38,-38,-38,-38,-38
S,InitialPrimeIdeal,Initialize a prime IdealRecord with the given data,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,270,-38,-38,-38,-38,-38
S,Inversionloop,Apply one iteration of the classical p-adic Newton method to find and approximation xnum/xden to the inverse of A,0,6,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,-1,,1,1,-38,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,IsIntegralM,True iff the algebraic number alpha is integral,0,1,0,0,0,0,0,0,0,28,,36,-38,-38,-38,-38,-38
S,LastLevel,This routine is called when an irreducible p-adic factor is detected in the Montes algorithm,0,5,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,-1,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,Lift,,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,28,-38,-38,-38,-38,-38
S,Lift,The output is a lift of class to an integral element in K,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,28,-38,-38,-38,-38,-38
S,LocalCRT,,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,Localize,"output=den,exp,Pol such that alpha = (1/den)*Pol(theta)*p^exp, and den is coprime to p",0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,148,148,312,-38,-38,-38
S,LocalLift,"class should belong to the residue class field P`Type[r]`Fq. The output is a pair g,e such that g(theta)/p^e is a lift to a P-integral element in K and deg g(x)<n_P",0,4,0,0,1,0,0,1,1,-38,,1,1,-38,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,LocalLift,class should belong to the residue class field Z_K/P. The output is a lift to a P-integral element in K,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,28,-38,-38,-38,-38,-38
S,LocalLift,,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,28,-38,-38,-38,-38,-38
S,Montes,Triple output: OMreps= list of OM representations of the p-adic irreducible factors of poly; TreesIntervals: list of intervals with the position of the prime ideals in each disconnected tree of OM representations; totalindex= p-index of the input polynomial. The option NumberField:=true forces the computation of the psi polynomial at the last level of each type of each OM representation,0,2,0,0,0,0,0,0,0,148,,0,0,312,,82,82,148,-38,-38,-38
S,Montes,Apply the Montes algorithm to the number field K and the rational prime p,0,2,0,0,1,0,0,0,0,148,,0,0,27,,-38,-38,-38,-38,-38,-38
S,Montesloop,Main loop of the Montes algorithm. The iteration stops as soon as totalindex is greater than the given mahler bound,0,4,0,0,1,0,0,0,0,-1,,1,1,-38,,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,MultiplierLift,Compute an element mult in the number field P`Parent which is congruent to 1 modulo P^a_P and it has Q-value >= exponents[Q],0,3,0,0,1,0,0,1,1,-38,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,Multipliers,"Computes multipliers c_P, one for each prime ideal P|p, satisfying v_P(c_P)=0, v_Q(c_P)ge values[P,Q]",0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,MultiplyByInverse,"Divide alpha by a pseudo-inverse, so that after the routine, it is congruent to 1 modulo P^m",0,3,0,0,1,0,0,0,0,-1,,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,MultPiece,"Compute bp in Parent(P) which has P-value zero and v_Q(bp) >= expsTree[Q], for all Qne P in the tree",0,5,0,0,1,0,0,1,1,-38,,1,1,-38,,0,0,-1,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,Newton,"Given a type of order at least i, and the phiadic expansion of a polynomial, compute: - sides=list of sides of the r-th order Newton polygon w.r.t. the type; - devsEachSide[k]=list of multiadic phi_1,...,phi_i-1 expansion of the coefficients of the polynomial whose attached point lies on sides[k]",0,5,0,0,1,0,0,1,1,-38,,1,1,-38,,1,1,-38,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,pAdicFactors,"Computes a list of Okutsu approximations to the irreducible p-adic factors of the given polynomial, all of them correct modulo p^precision. The routine detects a non-squarefree input polynomial and outputs an empty list in this case",0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,312,,82,-38,-38,-38,-38,-38
S,PathOfPrecisions,,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,pDiscriminant,"Compute: -pdiscf: is the p-adic valuation of the discriminant of polynomial. -pdiscK: sum of the p-adic valuations of the discriminants of all local extensions of Q_p, given by the irreducible p-adic factors of the given polynomial",0,2,0,0,0,0,0,0,0,148,,0,0,312,,148,148,-38,-38,-38,-38
S,PolToFieldElt,"Equivalent to g(K.1), but more efficient",0,2,0,0,0,0,0,0,0,312,,0,0,27,,28,-38,-38,-38,-38,-38
S,PrescribedValue,"Compute a polynomial Phi=phi_1^a_1...phi_r^a_r such that v_P(p^a_0 Phi(theta))=value, where r is the Okutsu depth of P. The exponents a_0,...,a_r are stored in the list logphi",0,4,0,0,1,0,0,1,1,-38,,1,1,-38,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,pResultant,Compute the p-adic valuation of the resultant of the two given univariate polynomials,0,3,0,0,0,0,0,0,0,148,,0,0,312,,0,0,312,,148,-38,-38,-38,-38,-38
S,PValuation,"Compute the P-valuation v of alpha at the prime ideal P. If the option RED is set to true, compute also the class of alpha in P^v/P^(v+1)",0,2,0,0,0,0,0,0,0,270,,0,0,267,,148,85,-38,-38,-38,-38
S,PValuation,"Compute the P-valuation v of alpha at the prime ideal P. If the option RED is set to true, compute also the class of alpha in P^v/P^(v+1)",0,2,0,0,0,0,0,0,0,270,,0,0,148,,148,85,-38,-38,-38,-38
S,PValuation,"Compute the P-valuation v of alpha at the prime ideal P. If the option RED is set to true, compute also the class of alpha in P^v/P^(v+1)",0,2,0,0,0,0,0,0,0,270,,0,0,28,,148,85,-38,-38,-38,-38
S,mod,Compute the reduction map ZK--> ZK/P,0,2,0,0,0,0,0,0,0,270,,0,0,28,,85,-38,-38,-38,-38,-38
S,Reduction,Compute the reduction map ZK--> ZK/P,0,2,0,0,0,0,0,0,0,270,,0,0,28,,85,-38,-38,-38,-38,-38
S,Reduction,Compute the reduction map ZK--> ZK/P^m,0,3,0,0,0,0,0,0,0,148,,0,0,270,,0,0,28,,82,-38,-38,-38,-38,-38
S,Representative,Construction of a representative phi of type. We add phi and V at a new level of this type,0,2,0,0,1,0,0,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,ResidualPolynomial,Internal procedure to compute the i-th residual polynomial psi with respect to a side S of slope -P`Type[i]`slope of the Newton polygon of a certain polynomial Pol. The coefficients of Pol whose attached points in N_i(Pol) lie on S have multiadic expansions contained in the list devsSide,0,4,0,0,1,0,0,1,1,-38,,1,1,-38,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,SFL,The aim is P`Type[#P`Type]`slope>=slope,0,2,0,0,1,0,0,0,0,148,,1,1,270,,-38,-38,-38,-38,-38,-38
S,SFLInitialize,"Initialize certain values of the given type, which are necessary for the SFL routine",0,1,0,0,1,0,0,1,1,-38,,-38,-38,-38,-38,-38,-38
S,SFLprecision,The final polynomial P`Type[#P`Type]`Phi is congruent to the true p-adic irreducible factor determined by P modulo p^precision,0,2,0,0,1,0,0,0,0,148,,1,1,270,,-38,-38,-38,-38,-38,-38
S,Taylor,Compute the first omega+1 coefficients of the phi-expansion of pol,0,3,0,0,0,0,0,0,0,148,,0,0,312,,0,0,312,,82,-38,-38,-38,-38,-38
S,TreeInterval,Returns the interval of positions in K`PrimeIdeals[p] of the tree to which O belongs,0,2,0,0,1,0,0,1,1,-38,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,TrueDiscriminant,"Computes the discriminant of K/Q, the factorization of the discriminant of the defining polynomial of K and the list of prime numbers dividing the index of Z[theta] in the maximal order",0,1,0,0,1,0,0,0,0,27,,-38,-38,-38,-38,-38,-38
S,UpdateLastLevel,"Updates the values of slope, h, psi and logGamma in the last level of P`Type. Also, it updates P`sflPols[1]=a0 and P`sflPols[2]=a1",0,1,0,0,1,0,0,1,1,-38,,-38,-38,-38,-38,-38,-38
S,Value,"Given a level i, a type and a polynomial pol, compute: - devs=multiexpansion with respect to phi_1,...,phi_i-1 of the points in S_lambda_i-1(pol); - val=(i-1)-th valuation of pol with respect to P`Type",0,5,0,0,1,0,0,1,1,-38,,1,1,-38,,1,1,-38,,1,1,-38,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,Z,The primitive element z_i of the i-th residual finite field F_(i+1) of the type is stored in the variable z,0,3,0,0,1,0,0,1,1,-38,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,eq,"True iff the fractional ideals I,J are equal. They are both factored if their factorization is not kwown",0,2,0,0,0,0,0,0,0,270,,0,0,270,,36,-38,-38,-38,-38,-38
S,ideal,Define the fractional ideal generated by the elements of listgen,1,1,1,82,0,28,2,0,0,0,0,0,0,0,82,,0,0,27,,270,-38,-38,-38,-38,-38
S,ideal,Define the principal fractional ideal generated by a,0,2,0,0,0,0,0,0,0,28,,0,0,27,,270,-38,-38,-38,-38,-38
S,ideal,Define the principal fractional ideal generated by the integer a,0,2,0,0,0,0,0,0,0,148,,0,0,27,,270,-38,-38,-38,-38,-38
S,in,,0,2,0,0,0,0,0,0,0,270,,0,0,-1,,36,-38,-38,-38,-38,-38
S,IsIdealRecord,,0,1,0,0,0,0,0,0,0,270,,36,-38,-38,-38,-38,-38
S,IsIntegral,,0,1,0,0,0,0,0,0,0,270,,36,-38,-38,-38,-38,-38
S,IsOne,True iff I is the total ideal,0,1,0,0,0,0,0,0,0,270,,36,-38,-38,-38,-38,-38
S,IsPrimeIdeal,,0,1,0,0,0,0,0,0,0,270,,36,-38,-38,-38,-38,-38
S,IsZero,,0,1,0,0,0,0,0,0,0,270,,36,-38,-38,-38,-38,-38
S,*,"Compute the product of the fractional ideals I,J. They are both factored if their factorization is not yet known",0,2,0,0,0,0,0,0,0,270,,0,0,270,,270,-38,-38,-38,-38,-38
S,^,Compute the n-th power of the fractional ideal I. The ideal I is factored if its factorization is not known,0,2,0,0,0,0,0,0,0,148,,0,0,270,,270,-38,-38,-38,-38,-38
S,/,"Compute the quotient of the fractional ideals I,J. They are both factored if their factorization is not knwon",0,2,0,0,0,0,0,0,0,270,,0,0,270,,270,-38,-38,-38,-38,-38
S,+,"Compute the greatest common divisor of the fractional ideals I,J",0,2,0,0,0,0,0,0,0,270,,0,0,270,,270,-38,-38,-38,-38,-38
S,Factorization,Compute the decomposition of the fractional ideal I into prime ideals,0,1,0,0,0,0,0,0,0,270,,82,-38,-38,-38,-38,-38
S,Factorization,Compute the decomposition of the fractional ideal I into prime ideals,0,1,0,0,1,0,0,1,1,270,,-38,-38,-38,-38,-38,-38
S,FactorListToString,Write down a factorization in pretty form,0,1,0,0,0,0,0,0,0,-1,,298,-38,-38,-38,-38,-38
S,Norm,Compute the norm of the ideal I,0,1,0,0,0,0,0,0,0,270,,148,-38,-38,-38,-38,-38
S,PValuation,Compute the v_P-valuation of the ideal I,0,2,0,0,0,0,0,0,0,270,,0,0,270,,148,-38,-38,-38,-38,-38
S,RationalDenominator,The least positive integer a such that aI is included in the maximal order O,0,1,0,0,0,0,0,0,0,-1,,148,-38,-38,-38,-38,-38
S,RationalRadical,Compute the rational prime numbers dividing the norm of the ideal I,0,1,0,0,0,0,0,0,0,270,,82,-38,-38,-38,-38,-38
S,ResidueField,"Given a prime ideal P, return the residue field Z_K/P",0,1,0,0,0,0,0,0,0,270,,84,-38,-38,-38,-38,-38
S,subset,"True iff the fractional ideal J divides I, i.e., iff I/J is integral. Both ideals are factored if their factorization is not yet known",0,2,0,0,0,0,0,0,0,270,,0,0,270,,36,-38,-38,-38,-38,-38
S,TwoElement,Compute a pair of elements generating the ideal I,0,1,0,0,0,0,0,0,0,270,,82,-38,-38,-38,-38,-38
S,TwoElement,Compute a pair of elements generating the ideal,0,1,0,0,1,0,0,1,1,270,,-38,-38,-38,-38,-38,-38
S,InitializeType,"Initializa a typelevel record with the given data, and set two lists z, Y to store the primitive elements of the residual fields and the variables of the polynomial rings over these fields",0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,82,82,82,-38,-38,-38
S,EnlargeType,"Enlarge the given type t (and the lists Y, z) with the slope -h/e and residual polynomial psi",0,6,0,0,1,0,0,1,1,-38,,1,1,-38,,0,0,-1,,1,1,-38,,0,0,-1,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,CreateType,"Create a random type t whose invariants [h1,e1,f1,...,hr,er,fr] are specified in the list ll",0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,CreateRandomType,Create a random type of order r,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,RandomMultiplicityType,"Creates a random type of depth r and randomly combines s of its phi-polynomials, adding some random refinements. The full type is always included",0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,CreateRandomMultipleTypePolynomial,Compute a random irreducible polynomial with k types of order AT MOST r. The value of s is used to add a power of p to ensure irreducibility,0,4,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,CombineTypes,Compute and irreducible polynomial whose attached types are the given ones in the specified list,0,2,0,0,0,0,0,0,0,-1,,0,0,82,,312,-38,-38,-38,-38,-38
S,CombinePolynomialsWithDifferentPrimes,"Compute a polynomial which is congruent to the given polynomials f1, f2 modulo the specified powers of the primes p1, p2",0,5,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,GlobalExpansion,Compute the coefficients of the multi-phi-adic expansion of pol. They are stored in a recursive list,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,168,-38,-38,-38,-38,-38
S,Expand,This function is only useful to check the validity of expansions given by GlobalExpand,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,312,-38,-38,-38,-38,-38
S,ExpandTeX,Write in TeX form the multi-phi-adic expansion of pol,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandTeXAux,,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandPhiTeX,Write the phi-adic expansion of phi_k in TeX format,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandAllPhiTeX,Write in TeX format the phi-adic expansion of all the phi in the type,0,1,0,0,0,0,0,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandMagma,Write in Magma form the multi-phi-adic expansion of pol,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandMagmaAux,,0,3,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandPhiMagma,Write the phi-adic expansion of phi_k in Magma format,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,-54,-38,-38,-38,-38,-38
S,ExpandAllPhiMagma,Write in Magma format the phi-adic expansion of all the phi in the type,0,1,0,0,0,0,0,0,0,-1,,-54,-38,-38,-38,-38,-38
S,pIntegralBasis,Computes a triangular (Hermite) p-integral basis of the fractional ideal I,0,2,0,0,0,0,0,0,0,148,,0,0,270,,82,-38,-38,-38,-38,-38
S,pIntegralBasis,Returns a triangular (Hermite) p-basis of the maximal order,0,2,0,0,0,0,0,0,0,148,,0,0,27,,82,-38,-38,-38,-38,-38
S,reduceIdeal,Returns a new ideal J and exponent a such that the p-part of I is p^a J,0,2,0,0,0,0,0,0,0,148,,0,0,270,,270,148,-38,-38,-38,-38
S,HermiteFormBasis,,0,4,0,0,0,0,0,0,0,82,,0,0,82,,0,0,148,,0,0,270,,82,-38,-38,-38,-38,-38
S,IdealBasis,Compute a (Hermitian) basis of the ideal I,0,1,0,0,0,0,0,0,0,270,,82,-38,-38,-38,-38,-38
S,SIdealBasis,Compute an S-integral basis of I for the given set of primes S=primelist,0,2,0,0,0,0,0,0,0,82,,0,0,270,,82,82,267,-38,-38,-38
S,HermiteFormBasis,The input parameterizes a global triangular basis,0,2,0,0,0,0,0,0,0,82,,0,0,82,,82,-38,-38,-38,-38,-38
S,IntegralBasis,Compute a triangular basis of the maximal order ZK of K,0,1,0,0,0,0,0,0,0,27,,82,-38,-38,-38,-38,-38
S,MaxMinCore,"The core of the MaxMin algorithm. 	Input: 	 - The $Q$-value of every element of each Okutsu $P$-basis for all 	 $P$, $Q$ in the input set. 	 - 	Output: 	 - Indices of final basis elements as product of input bases elements 	 - Denominator exponents of each basis element 	 - The $P$-value of each basis element for all input $P$ 	 - The required $P$-value of each Montes approximation to $F_P$",0,2,0,0,0,0,0,0,0,82,,0,0,82,,82,82,168,168,-38,-38
S,MaxMin,,0,3,0,0,0,0,0,0,0,-1,,0,0,148,,0,0,27,,82,82,82,-38,-38,-38
S,liftMontesApproximations,Increase each $phi_P$ so that it's $P$-value is at least that of 	the corresponding entry in `req_phip_vals`,1,0,1,82,0,270,2,0,0,1,0,0,0,0,168,,1,1,82,,-38,-38,-38,-38,-38,-38
S,calculateOkutsuFramesValues,Calculate the primary and secondary values for the phi-polynomials of 	the Okutsu frames for all types,0,1,0,0,0,0,0,0,0,-1,,168,-38,-38,-38,-38,-38
S,updateOkutsuFrames,,0,2,0,0,1,0,0,0,0,82,,1,1,82,,-38,-38,-38,-38,-38,-38
S,computeOkutsuBasis,Efficiently computes the Okutsu basis for the given Okutsu frame. This is 	produced by the canonical product of the $phi$-polynomials from,0,1,0,0,0,0,0,0,0,-1,,82,-38,-38,-38,-38,-38
S,computeBasisValues,Efficiently computes the values of a basis,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,82,-38,-38,-38,-38,-38
S,computeBasesValues,Efficiently computes the values of all Oktusu bases,0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,168,-38,-38,-38,-38,-38
S,calculateBasisIndices,"Calculate the indices that represent basis elements. The 0th indiex is 	the exponent of x (only used if f_0 > 1) then the i-th index is the 	exponent of phi_i,P for the P associated with this type",0,1,0,0,0,0,0,0,0,270,,168,-38,-38,-38,-38,-38
S,computeLocalBasis,"Efficiently compute a local basis given the indices of which element 	from each $P$-basis is used to make up an element of the final basis. 	Note: We don't need to compute the Okutsu basis for each $P$ to do this, 	 just the Okutsu frame",0,4,0,0,0,0,0,0,0,-1,,0,0,-1,,0,0,-1,,0,0,-1,,168,-38,-38,-38,-38,-38
S,computeOkutsuBasisElement,,0,4,0,0,1,0,0,0,0,-1,,0,0,-1,,0,0,-1,,1,1,-38,,-38,-38,-38,-38,-38,-38
S,reductionExponents,Calculate the exponents for reduction modulo p^nu,0,2,0,0,0,0,0,0,0,82,,0,0,82,,82,-38,-38,-38,-38,-38
S,reducepBasis,Reduce all basis numerators mod their valuation,0,4,0,0,1,0,0,0,0,148,,0,0,82,,0,0,82,,1,1,82,,-38,-38,-38,-38,-38,-38
S,itertoolsProduct,The ugly implementation of the product function from python's itertools,0,1,0,0,0,0,0,0,0,-1,,168,-38,-38,-38,-38,-38
