if( problem == 1 )
    K = 1;
    u0 = 1;
    u_exact = @(t)( u0./( 1 - t*u0 ) );

    f = @(t,u)( u^2 );
end
if( problem == 2 )
    K = 1;
    u0 = pi/4;
    a  = pi/4;
    b = -atan( u0 - a );
    u_exact = @(t)(a + tan(t-b));

    f  = @(t,u)( 1 + (u - a)^2 );
end
if( problem == 3 )
    K = 3;
    t_star = -pi/12;
    u_star = pi/4;
    u_exact = @(t)( u_star + tan(t - t_star) + tan(t - t_star).^3 );
    phi = @(u)( (1/3)*asinh( abs( 0.5*3^1.5*u ) ) );
    x = @(u)( -2*3^(-0.5)*sign( u )*sinh( phi(u) ) );
    u0 = u_exact(0);
    
    f = @(t,u)( (1 + x(u-u_star).^2).*( 1 + 3*x(u-u_star).^2 ) );
end

f1 = @(t,v)( -v^2 * f( t,1/v ) );