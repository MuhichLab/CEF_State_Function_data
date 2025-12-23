function [Ho So] = get_O2_thermo(T)
    % These constats are based off 0.1MPa == 1 bar
    
    % constats from NIST for O2 @ 700-2000 K
    if ( T< 700)
        a=31.32234;b=-20.23531;c=57.86644;d=-36.50624;e=-0.007374;f=-8.903471;g=246.7945;h=0;
    elseif (T >= 700)
        % constats from NIST for O2 @ 700-2000 K
        a=30.03235;b=8.772972;c=-3.988133;d=0.788313;e=-0.741599;f=-11.32468;g=236.1663;h=0;
    end
    t = T/1000;
    Ho = a*t + b*t^2/2 + c*t^3/3 + d*t^4/4 - e/t + f - h; %kJ/mol
    So = (a*log(t) + b*t + c*t^2/2 + d*t^3/3 - e/(2*t^2) + g)/1000;% kJ/mol*K
end
