function muO = solve_for_muO(Ts,Press)
    % ochem is the model calculated chemical potential of O in eV

    syms P Ho So T
    conv = 96.487;  % 1 eV = 96.487 kJ/mol
    R_gas = 8.314462/1000/conv; %kJ/K*mol -> eV/K
    eqn = Ho/conv - (T).*So/conv + R_gas.*(T).*log(P);
    sP = matlabFunction(eqn);


    muO = zeros(size(Ts));
    for i = 1:length(Ts)
        T0 = Ts(i);
        P0 = Press(i);

        if ( T0< 700)
            a=31.32234;b=-20.23531;c=57.86644;d=-36.50624;e=-0.007374;f=-8.903471;g=246.7945;h=0;
        elseif (T0 >= 700 && T0 < 2000)
            % constats from NIST for O2 @ 700-2000 K
            a=30.03235;b=8.772972;c=-3.988133;d=0.788313;e=-0.741599;f=-11.32468;g=236.1663;h=0;
        elseif ( T0 >= 2000)
            a=20.91111;b=10.72071;c=-2.020498;d=0.146449;e=9.245722;f=5.337651;g=237.6185;h=0; 
        end
        t = T0/1000;
        H0 = a*t + b*t^2/2 + c*t^3/3 + d*t^4/4 - e/t + f - h; %kJ/mol
        S0 = (a*log(t) + b*t + c*t^2/2 + d*t^3/3 - e/(2*t^2) + g)/1000;% kJ/mol*K
        muO(i) = double(sP(H0,P0,S0,T0))/2;   
    end
end