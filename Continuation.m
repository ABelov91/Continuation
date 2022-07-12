TIME = 15;

problem = 3; % 1 -> одиночный полюс,
             % 2 -> полюсы 1-го пор, авт.,
             % 3 -> полюсы 3-го пор, авт, без несинг ОТ,
problem_list;

K0 = 1;

A = 5;
tol2 = 1e-1;
tol3 = 1e-4;

P = 10;
N0 = 100;
Ns = zeros(1,P);

for p = 1:P
    N = N0*2^(p-1);
    Ns(p) = N;

    tau = TIME/N;
    t = 0:tau:TIME;

    u  = zeros(1,N+1);
    v  = zeros(1,N+1);
    w = zeros(1,N+1);
    rh = zeros(1,N+1);
    
    flag = zeros(1,N+1);
    order_u  = zeros(1,N+1);
    order_u2 = zeros(1,N+1);
    tmp      = zeros(1,N+1);
    check = ones(1,N+1);
    
    u(1) = u0;
    v(1) = 1/u0;
    w(1) = 1;
    rh(1) = f( 0, u0 );
    u_ex = zeros(1,length(t));
    for n = 1:length(t)
        u_ex(n) = u_exact(t(n));
    end
    
    Num = double.empty;
    t_pole = double.empty;
    t_pole_error = double.empty;

    curr = u0;
    
    for n = 1:N
        if( flag(n) == 0 )
            u(n+1) = erk4(  f , t(n) , curr , tau );
            
            curr = u(n+1);
            flag(n+1) = 0;
            v(n+1) = 1/u(n+1);
            w(n+1) = sign(u(n+1))/abs(u(n+1))^(1/K0);
            rh(n+1) = f1( t(n+1), v(n+1) );
            if( abs(u(n+1)) > A )
                flag(n+1) = 1;
                curr = v(n+1);
                signum = sign( u(n+1) );
                f2 = @(t,w)( -1/K0*signum*w.^(1+K0) * f( t,signum/w.^K0 ) );
            end
            
            cond_diagn1 =  v(n) *  v(n+1);
            cond_diagn2 = rh(n) * rh(n+1);
            cond_diagn3 = v(n+1)*rh(n+1);
            cond_diagn4 = abs( v(n) ) - abs( v(n+1) );
            
            if( cond_diagn1 > 0 && cond_diagn2 > 0 && cond_diagn3 < 0 && cond_diagn4 > 0 )
                num = t(n+1) - t(n);
                denom = u(n)/f( t(n),u(n) ) - u(n+1)/f( t(n+1),u(n+1) );
                order_u(n+1) = num/denom;
                
                rel1 = rh(n)/rh(n+1);
                rel2 =  v(n)/ v(n+1);
                order_u2(n+1) = ( 1 - log(rel1) / log(rel2) )^(-1);
            end
        end

        if( flag(n) == 1 )
            
            v(n+1) = erk4(  f1 , t(n) , curr , tau );
            curr = v(n+1);
            flag(n+1) = 1;
            u(n+1) = 1/v(n+1);
            w(n+1) = sign(v(n+1))*abs(v(n+1))^(1/K0);
            rh(n+1) = f1( t(n+1), v(n+1) );
            
            check(n) = abs( v(n+1) );
            
            cond_diagn1 =  v(n) *  v(n+1);
            cond_diagn2 = rh(n) * rh(n+1);
            cond_diagn3 = v(n+1)*rh(n+1);
            cond_diagn4 = abs( v(n) ) - abs( v(n+1) );
            
            if( cond_diagn1 > 0 && cond_diagn2 > 0 && cond_diagn3 < 0 && cond_diagn4 > 0 )
                num = t(n+1) - t(n);
                denom = u(n)/f( t(n),u(n) ) - u(n+1)/f( t(n+1),u(n+1) );
                order_u(n+1) = num/denom;
                
                rel1 = rh(n)/rh(n+1);
                rel2 =  v(n)/ v(n+1);
                order_u2(n+1) = ( 1 - log(rel1) / log(rel2) )^(-1);
            end
            
            test = round( order_u(n+1) );
            cond1 = abs(test - order_u(n+1));
            cond2 = abs(order_u(n+1) - order_u( n ));
            cond3 = abs(order_u( n ) - order_u(n-1));
            if( test > 1 )
                if( cond1 < tol2 && cond2 < tol3 && cond3 < tol3 && test ~= 1 )
                    flag(n+1) = 2;
                    K0 = test;
                    w(n+1) = sign(v(n+1))*abs(v(n+1))^(1/K0);
                    curr = w(n+1);
                    signum = sign( u(n+1) );
                    f2 = @(t,w)( -1/K0*signum*w.^(1+K0) * f( t,signum/w.^K0 ) );
                end
            end
                
            if( check(n) > 1/A )
                flag(n+1) = 0;
                curr = u(n+1);
            end
        end
        
        if( flag(n) == 2 )
            w(n+1) = erk4(  f2 , t(n) , curr , tau );
            curr = w(n+1);
            flag(n+1) = 2;
            u(n+1) = 1/(w(n+1))^K0;
            v(n+1) = w(n+1)^K0;
            rh(n+1) = f1( t(n+1), v(n+1) );
            
            cond_diagn1 =  v(n) *  v(n+1);
            cond_diagn2 = rh(n) * rh(n+1);
            cond_diagn3 = v(n+1)*rh(n+1);
            cond_diagn4 = abs( v(n) ) - abs( v(n+1) );
            
            if( cond_diagn1 > 0 && cond_diagn2 > 0 && cond_diagn3 < 0 && cond_diagn4 > 0 )
                num = t(n+1) - t(n);
                denom = u(n)/f( t(n),u(n) ) - u(n+1)/f( t(n+1),u(n+1) );
                order_u(n+1) = num/denom;
                
                rel1 = rh(n)/rh(n+1);
                rel2 =  v(n)/ v(n+1);
                order_u2(n+1) = ( 1 - log(rel1) / log(rel2) )^(-1);
            end

            if( abs(v(n+1)) > 1/A )
                flag(n+1) = 0;
                curr = u(n+1);
            end
        end
        
        if( flag(n) == 3 )
            w(n+1) = erk4(  f2 , t(n) , curr , tau );
            curr = w(n+1);
            flag(n+1) = 3;
            u(n+1) = signum/(w(n+1))^K0;
            v(n+1) = signum*w(n+1)^K0;
            
            tmp1 = v(n+1);
            tmp2 = u(n+1);
            
            rh(n+1) = f1( t(n+1), v(n+1) );
            
            cond_diagn1 =  v(n) *  v(n+1);
            cond_diagn2 = rh(n) * rh(n+1);
            cond_diagn3 = v(n+1)*rh(n+1);
            cond_diagn4 = abs( v(n) ) - abs( v(n+1) );
            
            if( cond_diagn1 > 0 && cond_diagn2 > 0 && cond_diagn3 < 0 && cond_diagn4 > 0 )
                num = t(n+1) - t(n);
                denom = u(n)/f( t(n),u(n) ) - u(n+1)/f( t(n+1),u(n+1) );
                order_u(n+1) = num/denom;
                
                rel1 = rh(n)/rh(n+1);
                rel2 =  v(n)/ v(n+1);
                order_u2(n+1) = ( 1 - log(rel1) / log(rel2) )^(-1);
            end

            if( abs(v(n+1)) > 1/A )
                flag(n+1) = 0;
                curr = u(n+1);
            end
        end

    end
    
    counter = 1;
    for n = 1:N
        if( flag(n) == 2 || flag(n) == 3 )
            test = w(n)*w(n+1);
            if( test < 0 )
                Num(counter) = n;
                counter = counter + 1;
            end
        elseif( flag(n) == 1 )
            test = v(n)*v(n+1);
            if( test < 0 )
                Num(counter) = n;
                counter = counter + 1;
            end
        end
    end
    counter = counter - 1;
    
    for count = 1:counter
        
        T = zeros(1,4);
        V = zeros(1,4);
        for m = 1:4
            T(m) = t( Num(count) + m-2 );
            V(m) = w( Num(count) + m-2 );
        end
        d1 = zeros(1,3);
        d2 = zeros(1,2);
        for m = 1:3
            d1(m) = (T(m) - T(m+1))/(V(m) - V(m+1));
        end
        for m = 1:2
            d2(m) = (d1(m) - d1(m+1))/(V(m) - V(m+2));
        end
        d3 = (d2(1) - d2(2))/(V(1) - V(4));

        t_pole(count) = T(1) - V(1)*d1(1) + V(1)*V(2)*d2(1) - V(1)*V(2)*V(3)*d3(1);
    end
       
    
end

illustrations;