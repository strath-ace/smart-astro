
"""
    eq2deg(coeffs)
    finds the real roots of an equation of degree 2
    of the form 
            a0x^2 + a1x + a2 = 0
            
    INPUT 
    coeffs=[a0,a1,a2]

    OUTPUT
    nreal = number of real roots
    xsol  = vector of size (nreal x 1) containing the real roots

    AUTHOR:
    Irene Cavallari

"""
function eq2deg(coeffs)
    # a0x^2 + a1x + a2 = 0
    # coeffs=[a0,a1,a2,a3]

    if length(coeffs)!=3
        Base.error("wrong inputs");
    end

    a0 = coeffs[1];
    a1 = coeffs[2];
    a2 = coeffs[3];

    nreal = 0;
    xsol = [];
    if a0 == 0.0
        if a1 != 0.0
            nreal = 1;
            xsol = [-a2/a1];
        end
    else
        discr = a1^2.0 - 4.0*a0*a2;
        if discr<0.0 && abs(discr)<1e-15
            discr = 0.0;
        end
        if discr>=0.0            
            x1 = (- a1 + sqrt(discr))/2.0/a0;
            x2 = (- a1 - sqrt(discr))/2.0/a0;
            nreal = 2;
            xsol = [x1,x2];
        end
    end
    return nreal,xsol;
end

"""
    eq3deg(coeffs)
    finds the real roots of an equation of degree 2
    of the form 
                a0 x^3 + a1 x^2 + a2 x + a3 = 0
            
    INPUT 
    coeffs=[a0,a1,a2,a3]

    OUTPUT
    nreal = number of real roots
    xsol  = vector of size (nreal x 1) containing the real roots

    AUTHOR:
    Irene Cavallari

"""
function eq3deg(coeffs)

    if length(coeffs)!=4
        Base.error("wrong inputs");
    end

    a0 = coeffs[1];
    a1 = coeffs[2];
    a2 = coeffs[3];
    a3 = coeffs[4];

    nreal = 0;
    xsol = [];
    if a0==0.0
        nreal,xsol = eq2deg([a1,a2,a3]);
    else
        pp = a2/a0 - (a1/a0)^2.0/3.0 
        qq = a3/a0 - a1*a2/(a0^2.0)/3.0 + 2.0/27.0*(a1/a0)^3.0;
        discr = pp^3.0/27.0 + qq^2.0/4.0;

        if discr>0.0
            uu3 = -qq/2.0 - sqrt(discr);
            if uu3>0.0
                uu = uu3^(1.0/3.0)
            else
                uu = -(-uu3)^(1.0/3.0);
            end
            vv3 = -qq/2.0 + sqrt(discr);
            if vv3>0.0
                vv = vv3^(1.0/3.0)
            else
                vv = -(-vv3)^(1.0/3.0);
            end
            nreal = 1;
            xsol = [uu+vv-a1/a0/3.0];
        elseif discr == 0.0
            nreal = 2;
            if qq>0.0
                x1 = (qq/2.0)^(1.0/3.0);
            else
                x1 = -(-qq/2.0)^(1.0/3.0);
            end
            xsol = [x1-a2/a0/3.0, -2.0*x1-a1/a0/3.0];
        else
            pR = -qq/2.0
            pI = sqrt(-discr);
            theta = mod(atan(pI,pR),2.0*pi);
            nreal = 3;
            xsol = 2.0*sqrt(-pp/3.0)*[cos(theta/3.0), cos((theta+2.0*pi)/3.0), cos((theta+4.0*pi)/3.0)]-a1/a0/3.0*[1.0, 1.0, 1.0];
        end
    end
    return nreal,xsol;
end

"""
    eq4deg(coeffs)
    finds the real roots of an equation of degree 2
    of the form 
        a0 x^4 + a1 x^3 + a2 x^2 + a3 x + a4 = 0
            
    INPUT 
    coeffs=[a0,a1,a2,a3,a4]

    OUTPUT
    nreal = number of real roots
    xsol  = vector of size (nreal x 1) containing the real roots

    AUTHOR:
    Irene Cavallari

"""
function eq4deg(coeffs)

    if length(coeffs)!=5
        Base.error("wrong inputs");
    end

    a0 = coeffs[1];
    a1 = coeffs[2];
    a2 = coeffs[3];
    a3 = coeffs[4];
    a4 = coeffs[5];

    nreal = 0;
    xsol = [];
    if a0==0.0
        nreal,xsol = eq3deg([a1,a2,a3,a4]);
    else
        cc = ((256.0*a0^3.0*a2 - 96.0*a0^2.0*a1^2.0))/(256.0*a0^4.0)
        dd = ((256.0*a0^3.0*a3 - 128.0*a0^2.0*a1*a2 + 32.0*a0*a1^3.0))/(256.0*a0^4.0);
        ee = (256.0*a0^3.0*a4  - 64.0*a0^2.0*a1*a3  + 16.0*a0*a1^2.0*a2 - 3.0*a1^4.0)/(256.0*a0^4.0)
        if dd==0.0
            pp = 0.0;
            rr = 0.0;
            if cc==0.0 && ee==0.0
                nreal = 4
                xsol = - a1/a0/4.0*[1.0,1.0,1.0,1.0];
            else
                nqsol,qqsol = eq2deg([1.0,-cc,ee]);
                if nqsol>0
                    qq = qqsol[1];
                    ss = ee/qq;
                    nzz1sol,zz1sol = eq2deg([1.0,pp,qq]);
                    nzz2sol,zz2sol = eq2deg([1.0,rr,ss]);
                    nreal = nzz1sol+nzz2sol;
                    zzsol = append!(zz1sol,zz2sol);
                    xsol  = zzsol - a1/a0/4.0*ones(size(zzsol));
                end
            end
        else
            npreal,psol = eq3deg([1.0,2.0*cc,cc^2.0-4.0*ee,-dd^2.0]);
            if npreal>0
                kk = 0;
                for jj=1:npreal
                    if psol[jj]>0.0
                        kk = jj
                        break;
                    end
                end
                if kk>0
                    pp = sqrt(psol[kk]);
                    rr = -pp;
                    ss = (cc+pp^2.0+dd/pp)/2.0;
                    qq = (cc+pp^2.0-dd/pp)/2.0;
                    nzz1sol,zz1sol = eq2deg([1.0,pp,qq]);
                    nzz2sol,zz2sol = eq2deg([1.0,rr,ss]);
                    nreal = nzz1sol+nzz2sol;
                    zzsol = append!(zz1sol,zz2sol);
                    xsol  = zzsol - a1/a0/4.0*ones(size(zzsol));
                end
            end
        end
    end
    return nreal,xsol;
end


