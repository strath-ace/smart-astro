"""
    myfix(numb)
    rounds the input towards zero
"""
function myfix(numb)
    if numb>0.0
        return floor(numb);
    elseif numb <0.0
        return ceil(numb); 
    else
        return numb; 
    end
end

"""
    my_isnan
    determine whether the input is NaN
"""
function my_isnan(numb)
    if ((numb/numb==1) || numb==0)
       my_isnan = false;
    else
       my_isnan = true;
    end
    return my_isnan;
end 

"""
    simps1(fun::Function,a,b,mx)
    Simpson rule to evaluate integrals of the type 
    int_a^b f(x) dx.

    INPUTS
    fun = integrand function in two variables
    a = lower bound of the first integration interval
    b = upper bound of the first integration interval
    mx = natural number: 2mx-1 is the number of subintervals considered in [a,b] to perform the integration
    
    OUTPUTS
    out = value of the integral
"""
function simps1(fun::Function,a,b,mx)
    
    if b<a 
        Base.error("wrong inputs");
    end

    x  = collect(range(a,b,2*mx+1));
    hx = x[2]-x[1]
    
    # hx = (b-a)/2.0/mx;
    # x = collect(a:hx:b);
     
    I = fun(x[1])+fun(x[2*mx+1]);
            
    for k = 1:mx-1
        I = I + 2.0*fun(x[2*k+1]);
    end

    for k = 1:mx
        I = I + 4.0*fun(x[2*k]);
    end
        
    out = I*hx/3.0; 
    return out;    
end

"""
    simps2(fun::Function,a,b,c,d,mx,my)
    Simpson rule to evaluate integrals of the type 
    int_a^bint_c^d f(x,y) dx dy.

    INPUTS
    fun = integrand function in two variables
    a = lower bound of the first integration interval
    b = upper bound of the first integration interval
    c = lower bound of the second integration interval
    d = upper bound of the second integration interval
    mx = natural number: 2mx-1 is the number of subintervals considered in [a,b] to perform the integration
    my = natural number: 2my-1 is the number of subintervals considered in [c,d] to perform the integration

    OUTPUTS
    out = value of the integral
"""
function simps2(fun::Function,a,b,c,d,mx,my)
    
    if b<a || d<c
        Base.error("wrong inputs")
    end
    
    hx = (b-a)/2.0/mx;
    hy = (d-c)/2.0/my;

    x = collect(a:hx:b); 
    y = collect(c:hy:d); 
    
    I =   fun(x[1],y[1]) + fun(x[2*mx+1],y[1])+
          fun(x[1],y[2*my+1]) + fun(x[2*mx+1],y[2*my+1])+
          4.0*(fun(x[2*mx],y[1])+fun(x[2*mx],y[2*my+1])+
          fun(x[1],y[2*my])+fun(x[2*mx+1],y[2*my]))+
          16.0*fun(x[2*mx],y[2*my]);
        
    for k = 1:mx-1
        I = I + 2.0*(fun(x[2*k+1],y[1])+fun(x[2*k+1],y[2*my+1]))+
            4.0*(fun(x[2*k],y[1])+fun(x[2*k],y[2*my+1]))+
            8.0*fun(x[2*k+1],y[2*my]) + 16*fun(x[2*k],y[2*my]);
    end
    
    for i = 1:my-1
        I = I + 2.0*(fun(x[1],y[2*i+1])+fun(x[2*mx+1],y[2*i+1]))+
            4.0*(fun(x[1],y[2*i])+fun(x[2*mx+1],y[2*i]))+
            8.0*fun(x[2*mx],y[2*i+1])+
            16.0*fun(x[2*mx],y[2*i]);
    
        for k = 1:mx-1
            I = I + 4.0*fun(x[2*k+1],y[2*i+1]) + 8.0*(fun(x[2*k],y[2*i+1])+
                fun(x[2*k+1],y[2*i])) + 16.0*fun(x[2*k],y[2*i]);
        end
    end
    
    out = I*hx*hy/9.0; 
    return out;    
end

"""
    getintegral(fun::Function,A,B,npoints)
    interface for solver of int_A^B fun(x) dx
        
    Currently: trapezoidal rule implemented

    INPUTS:
    fun = integrand function in one variable
    A = lower bound of the integration interval
    B = upper bound of the integration interval
    npoints = number of subintervals considered in [A,B] to perform the integration
    OUTPUTS:
    out = value of the integral     

"""
function getintegral(fun::Function,A,B,npoints)
    # intval = simps1(fun,A,B,Int64(floor(npoints/2)));
    intval = trapezoidalrule(fun,A,B,npoints);
    return intval
end

"""
    trapezoidalrule(fun::Function,a,b,nn)
    Trapezoidal rule to evaluate integrals of the type 
    int_a^b f(x) dx.

    INPUTS:
        fun = integrand function in one variable
        a = lower bound of the integration interval
        b = upper bound of the integration interval
        nn = number of subintervals considered in [a,b] to perform the integration
    OUTPUTS:
        out = value of the integral        
"""
function trapezoidalrule(fun::Function,a,b,nn)
    fA = fun(a);
    fB = fun(b);
    intval = (fA+fB)/2
    cc = (b-a)/nn;
    for jj = 1:nn-1
        intval = intval + fun(a+jj*cc)
    end
    intval = intval*cc;
    return intval
end

"""
    trapezoidalrule_V2(yVect,xVect)
    Trapezoidal rule to evaluate integrals of the type 
    int_a^b f(x) dx.

    INPUTS:
    xVect = vector of points in the integration interval [a,b]
    yVect = vector of values of the integrand evaluated at xVect
    
    OUTPUTS:
    out = value of the integral

"""
function trapezoidalrule_V2(yVect,xVect)
    ll = size(yVect)[1]; 
    intval = sum((yVect[2:ll]+yVect[1:ll-1]).*(xVect[2:ll]-xVect[1:ll-1])/2.0)
    return intval
end



