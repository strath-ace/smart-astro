"""
    ellP(k,phi,m)
    computation of incomplete elliptic integral of third kind
    int_0^phi 1/(1+k*sin^2(theta))/sqrt(1-m sin^2(theta))
"""
function ellP(k,phi,m)
    if m<0.0
        out = ellPC(k,phi,m)
        return out;
    end


    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end
    nnn = 0.0;
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        PP = Elliptic.Pi(k,pi/2.0,m);
    end

    if  phi<pi/2.0
        out = Elliptic.Pi(k,phi,m);
    else 
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            theta=(nn+1.0)*pi/2.0-phi; 
            out = (nn+1.0)*PP-Elliptic.Pi(k,theta,m);
        else
            theta = phi-nn*pi/2.0;
            out = nn*PP+Elliptic.Pi(k,theta,m);
        end
    end

    if nnn>0.0
        out = out + nnn*4.0*PP;
    end

    out = out * ss;

    return out;
    # a = 0; 
    # b = phi; 
    # fa = 1/sqrt(1-m*sin(a)^2)/(1-k*sin(a)^2); 
    # fb = 1/sqrt(1-m*sin(b)^2)/(1-k*sin(b)^2);
    # out = (fa+fb)/2;
    # n = 2000;
    # for jj=1:n-1
    #  arg = a+jj*(b-a)/n; 
    #  out = out + 1/sqrt(1-m*sin(arg)^2)/(1-k*sin(arg)^2);
    # end
    # out = out*(b-a)/n;
    # return out;
end

function ellPC(k,phi,m)
    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end

    mmod = -m;
    kk   = mmod/(1+mmod);
    kkp  = 1/sqrt(1+mmod);
    k1   = (k+mmod)/(1+mmod);
    
    nnn = 0.0
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        PP = (kkp/k1)*(kk*Elliptic.K(kk) + kkp^2*k*Elliptic.Pi(k1,pi/2.0,kk));
    end

    if phi<pi/2.0
        theta  = asin(sqrt(1+mmod)*sin(phi)/sqrt(1+mmod*sin(phi)^2.0));
        out = (kkp/k1)*(kk*Elliptic.F(theta,kk) + kkp^2*k*Elliptic.Pi(k1,theta,kk));
    else
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            phimod=(nn+1.0)*pi/2.0-phi; 
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = (nn+1.0)*PP - (kkp/k1)*(kk*Elliptic.F(theta,kk) + kkp^2*k*Elliptic.Pi(k1,theta,kk));
        else
            phimod = phi-nn*pi/2.0;
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = nn*PP + (kkp/k1)*(kk*Elliptic.F(theta,kk) + kkp^2*k*Elliptic.Pi(k1,theta,kk));
        end       
    end

    if nnn>0.0
        out = out + nnn*4.0*PP;
    end

    out = ss*out;
    
    return out;
end

"""
    ellF(phi,m)
    computation of incomplete elliptic integral of first kind
    int_0^phi 1/sqrt(1-m sin^2(theta))
"""
function ellF(phi,m)
    if m<0.0
        out = ellFC(phi,m)
        return out;
    end

    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end
    
    nnn = 0.0
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        KK = Elliptic.K(m);
    end

    if  phi<pi/2.0 
        out = Elliptic.F(phi,m);
    else 
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            theta=(nn+1.0)*pi/2.0-phi; 
            out = (nn+1.0)*KK-Elliptic.F(theta,m);
        else
            theta = phi-nn*pi/2.0;
            out = nn*KK+Elliptic.F(theta,m);
        end        
    end

    if nnn>0.0
        out = out + nnn*4.0*KK;
    end

    out = out * ss;

    return out;

end

function ellFC(phi,m)
    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end

    mmod = -m;
    kk   = mmod/(1+mmod);
    kkp  = 1/sqrt(1+mmod);

    nnn = 0.0
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        KK = kkp*Elliptic.K(kk);
    end

    if phi<pi/2.0
        theta  = asin(sqrt(1+mmod)*sin(phi)/sqrt(1+mmod*sin(phi)^2.0));
        out = kkp*Elliptic.F(theta,kk);
    else
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            phimod=(nn+1.0)*pi/2.0-phi; 
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = (nn+1.0)*KK - kkp*Elliptic.F(theta,kk);
        else
            phimod = phi-nn*pi/2.0;
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = nn*KK + kkp*Elliptic.F(theta,kk);
        end       
    end

    if nnn>0.0
        out = out + nnn*4.0*KK;
    end

    out = ss*out;
    return out;


end

"""
    ellE(phi,m)
    computation of incomplete elliptic integral of second kind
    int_0^phi sqrt(1-m sin^2(theta))
"""
function ellE(phi,m)
    if m<0.0
        out = ellEC(phi,m)
        return out;
    end

    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end
    
    nnn = 0.0
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        EE = Elliptic.E(m);
    end

    if  phi<pi/2.0
        out = Elliptic.E(phi,m);
    else 
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            theta=(nn+1.0)*pi/2.0-phi; 
            out = (nn+1.0)*EE-Elliptic.E(theta,m);
        else
            theta = phi-nn*pi/2.0;
            out = nn*EE+Elliptic.E(theta,m);
        end
    end

    if nnn>0.0
        out = out + nnn*4.0*EE;
    end

    out = out * ss;

    # a = 0; 
    # b = phi; 
    # fa = sqrt(1-m*sin(a)^2); 
    # fb = sqrt(1-m*sin(b)^2);
    # out = (fa+fb)/2;
    # n = 2000;
    # for jj=1:n-1
    #  arg = a+jj*(b-a)/n; 
    #  out = out + sqrt(1-m*sin(arg)^2);
    # end
    # out = out*(b-a)/n;
    return out;
end

function ellEC(phi,m)
    ss = 1.0;
    if phi<0.0
        phi = -phi;
        ss = -1.0;
    end

    mmod = -m;
    kk   = mmod/(1+mmod);
    kkp  = 1/sqrt(1+mmod);

    nnn = 0.0
    if phi>2.0*pi
        nnn  = myfix(phi/(2.0*pi));
        phi = mod(phi,2.0*pi); 
    end

    if nnn>0.0 || phi>=pi/2.0
        EE = (1/kkp)*Elliptic.E(kk);
    end

    if phi<pi/2.0
        theta  = asin(sqrt(1+mmod)*sin(phi)/sqrt(1+mmod*sin(phi)^2.0));
        out = (1/kkp)*(Elliptic.E(theta,kk) - kk*sin(theta)*cos(theta)/sqrt(1-kk*sin(theta)^2.0))
    else
        nn = floor(phi/(pi/2.0)); 
        if isodd(nn)
            phimod=(nn+1.0)*pi/2.0-phi; 
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = (nn+1.0)*EE - (1/kkp)*(Elliptic.E(theta,kk) - kk*sin(theta)*cos(theta)/sqrt(1-kk*sin(theta)^2.0));
        else
            phimod = phi-nn*pi/2.0;
            theta = asin(sqrt(1+mmod)*sin(phimod)/sqrt(1+mmod*sin(phimod)^2.0))
            out = nn*EE + (1/kkp)*(Elliptic.E(theta,kk) - kk*sin(theta)*cos(theta)/sqrt(1-kk*sin(theta)^2.0));
        end       
    end

    if nnn>0.0
        out = out + nnn*4.0*EE;
    end

    out = ss*out;

    return out;
end




################################################
   
# """
#     lellipe3_mod(phi,m,errtol)
    
#     Exension of lellipe by Thomas Hoffend to calculate the Legendre elliptic
#     integral of the second kind for every value of phi

#     Compute Legendre's (incomplete) elliptic integral E(phi,k).
#     Uses a vectorized implementation of Carlson's Duplication Algorithms 
#     for symmetric elliptic integrals as found in "Computing Elliptic 
#     Integrals by Duplication," by B. C. Carlson, Numer. Math. 33, 1-16 (1979)
#     and also found in ACM TOMS Algorithm 577.  Section 4 in the paper cited
#     here describes how to convert between the symmetric elliptic integrals
#     and Legendre's elliptic integrals.

#     Returns NaN's for any argument values outside input range.

#     INPUT:
#         phi = input angle vector size 1xN. There are no constraints on phi
#         m = input parameter vector size 1 or 1xN (m = k^2)
#         errtol = error tolerance for Carlson's algorithms

#     OUTPUT:
#         e = value of the elliptic integral of the second kind

#     functions called: lellipe by Thomas Hoffend

#     evolution of: lellipe

#     issue - check for which value of m is valid

#     - Camilla Colombo - 20/02/2007
#     Modified by Marilena Di Carlo - 16/04/2015   

# """
# function lellipe3_mod(phi,m,errtol)
    
#     phi0 = phi;
#     npi  = myfix(phi0/(2*pi));
#     phi  = phi0 - npi*(2*pi); 
            
#     eval=0.0;
#     logi_a = phi>=-pi/2 && phi<=pi/2;
#     logi_b = ((phi>pi/2) && (phi<pi)) || ((phi<-pi/2) && (phi>-pi));
#     logi_c = ((phi>=pi) && (phi<=3/2*pi)) || ((phi<=-pi) && (phi>=-3/2*pi));
#     logi_d = !(logi_a || logi_b || logi_c);
    
#     if logi_a
#         eval = lellipe(phi,m,errtol);
#     end
#     if logi_b || logi_c ||logi_d
#         De = lellipe(pi/2,k,errtol);
#         if logi_b
#             phisup = -sign(phi)*2*pi + phi;
#             n      = myfix(phisup/(pi));
#             esup   = 2*n*De + lellipe(phisup-n*pi,m,errtol);
#             eval   = esup+sign(phi)*4*lellipe(pi/2,m,errtol);
#         end
#         if logi_c
#             n      = myfix(phi/(pi));
#             eval   = 2*n*De+lellipe(phi-n*pi,m,errtol);
#         end
#         if logi_d
#             n      = myfix(phi/(pi));
#             phi_pi = phi-n*pi;
#             phisup = -sign(phi_pi)*2*pi+phi_pi;
#             nsup   = myfix(phisup/(pi));
#             esup   = 2*nsup*De+lellipe(phisup-nsup*pi,k,errtol);
#             eval   = 2*n*De+esup+sign(phi_pi)*4*De;
#         end
#     elseif abs(npi)>0.0
#         De = lellipe(pi/2,m,errtol);
#     else
#         De = 0;
#     end
    
#     eval = eval + npi*4.0*De;
    
#     return eval;
# end

# function lellipe(phi, m, errtol)

#     snphi = sin(phi);
#     csphi = cos(phi);
#     snphi2 = snphi.^2;
#     csphi2 = csphi.^2;

#     y = 1.0 - m*snphi2;
#     f = snphi * ellrf(csphi2,  y, 1, errtol) - m*snphi*snphi2*ellrd(csphi2, y, 1.0, errtol)/3.0;
    
#     return f;
# end


# """
#     lellipf3_mod(phi,m,errtol)

#     Exension of lellipf by Thomas Hoffend to calculate the Legendre elliptic
#     integral of the first kind for every value of phi

#     f = lellipf2(phi,m,errtol)

#     Compute Legendre's (incomplete) elliptic integral F(phi, k).
#     Uses a vectorized implementation of Carlson's Duplication Algorithms 
#     for symmetric elliptic integrals as found in "Computing Elliptic 
#     Integrals by Duplication," by B. C. Carlson, Numer. Math. 33, 1-16 (1979)
#     and also found in ACM TOMS Algorithm 577.  Section 4 in the paper cited
#     here describes how to convert between the symmetric elliptic integrals
#     and Legendre's elliptic integrals.

#     Returns NaN's for any argument values outside input range.

#     INPUT:
#         phi = input angle vector size 1xN. There are no constraints on phi
#         m = input parameter vector size 1 or 1xN (m = k^2)
#         errtol = error tolerance for Carlson's algorithms

#     OUTPUT:
#         f = value of the elliptic integral of the first kind

#     functions called: lellipf by Thomas Hoffend

#     evolution of: lellipf

#     issue - check for which value of m is valid

#     - Camilla Colombo - 20/02/2007
#     Modified by Marilena Di Carlo - 16/04/2015
#     Scrivere dove
# """
# function lellipf3_mod(phi,m,errtol)
        
#     phi0 = phi;
#     npi  = myfix(phi0/(2*pi));
#     phi  = phi0 - npi*(2*pi); 
            
#     f = 0.0;
#     logi_a = (phi>=-pi/2) && (phi<=pi/2);
#     logi_b = ((phi>pi/2) && (phi<pi)) || ((phi<-pi/2) && (phi>-pi));
#     logi_c = ((phi>=pi)  && (phi<=3/2*pi)) || ((phi<=-pi) & (phi>=-3/2*pi));
#     logi_d = !(logi_a || logi_b || logi_c);
    
#     if logi_a
#         f = lellipf(phi,m,errtol);
#     end
#     if logi_b || logi_c || logi_d
#         Df = lellipf(pi/2,m,errtol);
#         if logi_b
#             phisup = -sign(phi)*2*pi+phi;
#             n      = myfix(phisup/(pi));
#             esup   = 2*n*Df+lellipf(phisup-n*pi,m,errtol);
#             f      = esup + sign(phi) * 4.0 * Df;
#         end
#         if logi_c
#             n = myfix(phi/(pi));
#             f = 2*n*Df + lellipf(phi-n*pi,m,errtol);
#         end
#         if logi_d
#             n      = myfix(phi/(pi));
#             phi_pi = phi-n*pi;
#             phisup = -sign(phi_pi)*2*pi+phi_pi;
#             nsup   = myfix(phisup/(pi));
#             esup   = 2*nsup*Df + lellipf(phisup-nsup*pi,m,errtol);
#             f      = 2*n*Df + esup + sign(phi_pi)*4*Df;
#         end
#     elseif npi>0.0
#         Df=lellipf(pi/2,m,errtol);
#     else 
#         Df=0;
#     end
    
#     f = f + npi*4*Df;
    
#     return   f;
# end

# function lellipf(phi, m, errtol)
    
#     snphi = sin(phi);
#     csphi = cos(phi);
#     csphi2 = csphi * csphi;
#     y = 1.0 - m * snphi * snphi;
#     f = snphi * ellrf(csphi2,  y, 1.0, errtol);

#     return f;
# end



# """
#     ellrf(x, y, z, errtol)

#     Inputs:

#     x       scalar
#     y       scalar
#     z       scalar
#     errtol  Error tolerance.

#     Julia function to compute Carlson's symmetric elliptic integral Rf.
#     Implementation of Carlson's Duplication Algorithm 1 in "Computing
#     Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
#     33, 1-16 (1979).

#     Returns NaN's for any argument values outside input range.

#     Algorithm is also from Carlson's ACM TOMS Algorithm 577.

#     This code is a complete rewrite of the algorithm in vectorized form.
#     It was not produced by running a FORTRAN to Matlab converter.

#     The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:

#     X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
#     %   ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
#     RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
#     16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).

#     SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
#                                 ERROR LESS THAN
#                     1.D-3      3.D-19
#                     3.D-3      2.D-16
#                     1.D-2      3.D-13
#                     3.D-2      2.D-10
#                     1.D-1      3.D-7

#     Note by TRH:

#     Absolute truncation error when the integrals are order 1 quantities
#     is closer to errtol, so be careful if you want high absolute precision.

#     Thomas R. Hoffend Jr., Ph.D.
#     3M Company
#     3M Center Bldg. 236-GC-26
#     St. Paul, MN 55144
#     trhoffendjr@mmm.com

# """
# function ellrf(x, y, z, errtol)

#     # Argument limits as set by Carlson:
#     LoLim = 5.0 * floatmin();
#     UpLim = 5.0 * floatmax();
    
#     # Check input arguments for acceptability:
#     mask = minimum([x, y, z]) >= 0 && minimum([(x + y), (x + z), (y + z)]) >= LoLim && maximum([x, y, z]) < UpLim;

#     if !mask 
#         return NaN;
#     end
    
#     # Define internally acceptable variable ranges for iterations:
#     Xi = x;
#     Yi = y;
#     Zi = z;
    
#     # Carlson's duplication algorithm for Rf:
#     Xn = deepcopy(Xi);
#     Yn = deepcopy(Yi);
#     Zn = deepcopy(Zi);
#     Mu = (Xn + Yn + Zn) / 3.0;
#     Xndev = 2.0 - (Mu + Xn) / Mu;
#     Yndev = 2.0 - (Mu + Yn) / Mu;
#     Zndev = 2.0 - (Mu + Zn) / Mu;
#     epslon = maximum( abs.([Xndev, Yndev, Zndev]) );
#     while (epslon >= errtol)
#         Xnroot = sqrt(Xn);
#         Ynroot = sqrt(Yn);
#         Znroot = sqrt(Zn);
#         lambda = Xnroot * (Ynroot + Znroot) + Ynroot * Znroot;
#         Xn = 0.25 * (Xn + lambda);
#         Yn = 0.25 * (Yn + lambda);
#         Zn = 0.25 * (Zn + lambda);
#         Mu = (Xn + Yn + Zn) / 3.0;
#         Xndev = 2.0 - (Mu + Xn) / Mu;
#         Yndev = 2.0 - (Mu + Yn) / Mu;
#         Zndev = 2.0 - (Mu + Zn) / Mu;
#         epslon = maximum( abs.([Xndev, Yndev, Zndev]) );
#     end
#     C1 = 1.0 / 24.0;
#     C2 = 3.0 / 44.0;
#     C3 = 1.0 / 14.0;
#     E2 = Xndev * Yndev - Zndev * Zndev;
#     E3 = Xndev * Yndev * Zndev;
#     S = 1.0 + (C1 * E2 - 0.1 - C2 * E3) * E2 + C3 * E3;
#     f = S ./ sqrt(Mu);
   
#     return f;
# end

# """
#     ellrd(x, y, z, errtol)

#     Inputs:

#     x       Input vector size 1xN.
#     y       Input vector size 1xN.
#     z       Input vector size 1xN.
#     errtol  Error tolerance.

#     Julia function to compute Carlson's symmetric elliptic integral Rd.
#     Implementation of Carlson's Duplication Algorithm 4 in "Computing
#     Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
#     33, 1-16 (1979).

#     Returns NaN's for any argument values outside input range.

#     Algorithm is also from Carlson's ACM TOMS Algorithm 577.

#     This code is a complete rewrite of the algorithm in vectorized form.
#     It was not produced by running a FORTRAN to Matlab converter.

#     The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:

#     X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).

#     ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
#     RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
#     16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).

#     SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
#                                 ERROR LESS THAN
#                     1.D-3      3.D-19
#                     3.D-3      2.D-16
#                     1.D-2      3.D-13
#                     3.D-2      2.D-10
#                     1.D-1      3.D-7

#     Note by TRH:

#     Absolute truncation error when the integrals are order 1 quantities
#     is closer to errtol, so be careful if you want high absolute precision.
    
#     Thomas R. Hoffend Jr., Ph.D.
#     3M Company
#     3M Center Bldg. 236-GC-26
#     St. Paul, MN 55144
#     trhoffendjr@mmm.com
    
    
# """
# function ellrd(x, y, z, errtol)

#     # Argument limits as set by Carlson:
#     LoLim = 5.0 * floatmin();
#     UpLim = 5.0 * floatmax();
    
#     # Check input arguments for acceptability:
#     mask = minimum([x, y]) >= 0 && minimum([(x + y), z]) >= LoLim && maximum([x, y, z]) < UpLim;

#     if !mask 
#         return NaN;
#     end
    
#     # Define internally acceptable variable ranges for iterations:
#     Xi = x;
#     Yi = y;
#     Zi = z;
    
#     # Carlson's duplication algorithm for Rf:
#     Xn = deepcopy(Xi);
#     Yn = deepcopy(Yi);
#     Zn = deepcopy(Zi);
#     sigma = 0.0;
#     power4 = 1.0;
    
#     Mu = (Xn + Yn + 3.0 * Zn) * 0.2;
#     Xndev = (Mu - Xn) / Mu;
#     Yndev = (Mu - Yn) / Mu;
#     Zndev = (Mu - Zn) / Mu;
#     epslon = maximum( abs.([Xndev, Yndev, Zndev]) );
#     while (epslon >= errtol)
#         Xnroot = sqrt(Xn);
#         Ynroot = sqrt(Yn);
#         Znroot = sqrt(Zn);
#         lambda = Xnroot * (Ynroot + Znroot) + Ynroot * Znroot;
#         sigma = sigma + power4 ./ (Znroot * (Zn + lambda));
#         power4 = 0.25 * power4;
#         Xn = 0.25 * (Xn + lambda);
#         Yn = 0.25 * (Yn + lambda);
#         Zn = 0.25 * (Zn + lambda);
#         Mu = (Xn + Yn + 3.0 * Zn) * 0.2;
#         Xndev = (Mu - Xn) / Mu;
#         Yndev = (Mu - Yn) / Mu;
#         Zndev = (Mu - Zn) / Mu;
#         epslon = maximum( abs.([Xndev, Yndev, Zndev]) );
#     end
#     C1 = 3.0 / 14.0;
#     C2 = 1.0 / 6.0;
#     C3 = 9.0 / 22.0;
#     C4 = 3.0 / 26.0;
#     EA = Xndev * Yndev;
#     EB = Zndev * Zndev;
#     EC = EA - EB;
#     ED = EA - 6.0 * EB;
#     EF = ED + EC + EC;
#     S1 = ED * (-C1 + 0.25 * C3 * ED - 1.50 * C4 * Zndev * EF);
#     S2 = Zndev * (C2 * EF + Zndev * (-C3 * EC + Zndev * C4 * EA));
#     f  = 3.0 * sigma + power4 * (1.0 + S1 + S2) / (Mu * sqrt(Mu));
    
#     return f;
# end