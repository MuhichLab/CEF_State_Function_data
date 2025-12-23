function [G,H,S,coeffs,S0,diffH] = fit_einstein(theta,H0,N)

%{ 
    Pass ins:
    Ho    == formation energy [J/mol] or Enthalpy at 0 K
    theta == Einstien Temperature to model Cp after
    plots == 'yes' produces various plots 

    Notes:
    einstein cp and fit cp are calcualted in J/mol
    cp parameters are fit to T/1000

    Outputs:
    Thermodynamic functions of G [J/mol], H [J mol] and S [J/(mol K)
%}

    syms T So Ho
    
    % Einstein Solid -- as T -> inf Cp = 3*N*kb

    %ABO3 => N = 5 atoms/mol
    kb = 8.314462; % J/(mol K)
    conv = 96.487; %  1 eV/K == 96.487 kJ/(mol K)
    
    % symbolic cp_ein
    cp_ein = 3*N*kb*(theta/T)^2*(exp(theta/T)/(exp(theta/T) -1)^2);
    fun_ein = matlabFunction(cp_ein);

    Temp = 0:5:2273;
    Cvtrue = fun_ein(Temp'); % J/(mol K)
    
    idx = Temp>298;
    
    ft = fittype('a + b*x/1000 + c*(x/1000)^2 + d*(x/1000)^3 + e/(x/1000)^2');
    % can get goodness of fit by replacing ~ with variable and printing
    [f,~] = fit(Temp(idx)',Cvtrue(idx),ft); 
    
    coeffs = coeffvalues(f);
    
    cp = coeffs(1) + coeffs(2)*T/1000 + coeffs(3)*(T/1000)^2 +...
    coeffs(4)*(T/1000)^3 + coeffs(5)/(T/1000)^2;

    cp_a = coeffs(1);

    % We need einstein model thermodynamics to correct integrated Cp
    % constant terms Ho and So
    dH_ein = int(cp_ein,T) + Ho; % J/mol
    dS_ein = int(cp_ein/T,T);    % J/(mol K)
    dG_ein = dH_ein - T*dS_ein;  % J/mol
    
    dH = int(cp,T) + Ho;    % J/mol
    dS = int(cp/T,T) + So;  % J/(mol K)

    dHa = int(cp_a,T) + Ho;    % J/mol
    dSa = int(cp_a/T,T) + So;  % J/(mol K)
    
    % Get corrected Ho 
    diffH = matlabFunction(subs(dH_ein,Ho,H0) - subs(dH,Ho,H0));
    diffH = mean(diffH(Temp(idx)));
    
    % Get corrected So
    S0 = matlabFunction(solve(dS - dS_ein,So));
    S0 = mean(S0(Temp(idx)));
    
    % Final Equations
    H = subs(dH,Ho,(Ho + diffH));
    S = subs(dS,So,S0);
    G = H - T*S;

end