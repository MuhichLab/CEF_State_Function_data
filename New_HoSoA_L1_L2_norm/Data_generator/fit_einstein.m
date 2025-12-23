function [G,H,S,coeffs,S0,diffH] = fit_einstein(theta,H0,T0,N,plots)

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

    syms T To So Ho
    
    % Einstein Solid -- as T -> inf Cp = 3*N*kb

    %ABO3 => N = 5 atoms/mol
    kb = 8.314462; % J/(mol K)
    conv = 96.487; %  1 eV/K == 96.487 kJ/(mol K)
    
    % symbolic cp_ein
    cp_ein = 3*N*kb*(theta/T)^2*(exp(theta/T)/(exp(theta/T) -1)^2);
    fun_ein = matlabFunction(cp_ein);

    Temp = 100:5:2273;
    Cvtrue = fun_ein(Temp'); % J/(mol K)
    
    idx = Temp>293;
    
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
    
    dH = int(cp,T) - subs(int(cp,T),T,To) + Ho;    % J/mol
    dS = int(cp/T,T) - subs(int(cp/T,T),T,To) + So;  % J/(mol K)

    dHa = int(cp_a,T) - subs(int(cp_a,T),T,To) + Ho;    % J/mol
    dSa = int(cp_a/T,T) - subs(int(cp_a/T,T),T,To) + So;  % J/(mol K)
    
    % Get corrected Ho 
    diffH = matlabFunction(subs(dH_ein,Ho,H0) - subs(dH,[To Ho],[T0 H0]));
    diffH = mean(diffH(Temp(idx)));
    
    % Get corrected So
    S0 = matlabFunction(solve(subs(dS,To,T0) - dS_ein,So));
    S0 = mean(S0(Temp(idx)));
    
    % Final Equations
    H = subs(dH,[To Ho],[T0 (Ho + diffH)]);
    S = subs(dS,[To So],[T0 S0]);
    G = H - T*S;

%%
    if plots == 1
        
        figure
        hold on
        fplot(subs(dH_ein,Ho,H0),[300 2000],'-k','linewidth',2.0)
        fplot(subs(H,Ho,H0),[300 2000],'--b','linewidth',2.0)
        fplot(subs(dHa,[To Ho],[T0 (H0 + diffH)]),[300 2000],':r','linewidth',2.0)
        fplot(subs(dHa,[To Ho],[T0 (H0)]),[300 2000],':g','linewidth',2.0)
        xlabel('Temperature [K]')
        ylabel('\DeltaH [J/mol]')
        box on
        grid on
        legend('Cp = Einstein','Full Cp','C_{p}=A','C_{p}=A no Ho correction')
        
        figure
        hold on
        fplot(dS_ein,[300 2000],'-k','linewidth',2.0)
        fplot(S,[300 2000],'--b','linewidth',2.0)
        fplot(subs(dSa,[To So],[T0 S0]),[300 2000],':r','linewidth',2.0)
        fplot(subs(dSa,[To So],[T0 0]),[300 2000],':g','linewidth',2.0)
        xlabel('Temperature [K]')
        ylabel('\DeltaS [J/(mol K)]')
        box on
        grid on
        legend('Cp = Einstein','Full Cp','C_{p}=A','C_{p}=A no Ho correction')

        figure 
        plot(Temp,Cvtrue,'-k','Linewidth',2.0)
        ylabel('Heat Capacity [J/(mol K)]')
        ylim([0 130])
        yyaxis right
        plot(Temp,Cvtrue/conv/1000,'-k','Linewidth',2.0)
        ylabel('Heat Capacity [eV/K]')
        ylim([0 130/conv/1000])
        xlabel('Temperature [K]')
        xlim([min(Temp) max(Temp)])
        legend(["Einstein Model \Theta = " + num2str(theta)],'location','southeast');
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';


        figure
        plot(f,Temp,Cvtrue)
        legend(["Einstein Model \Theta = " + num2str(theta)],'Cp Fit','location','southeast');
        title('cp = a + b*T + c*T^{2} + d*T^{3} - e/T^{2}');
        ylabel('Heat Capacity [J/(mol K)]')
        xlabel('Temperature [K]')
        box on
        xlim([min(Temp) max(Temp)])
        ylim([0 130])
        
        
        figure
        hold on
        fplot(subs(G,Ho,H0)/1000,[293,max(Temp)],'-k')
        fplot(subs(H,Ho,H0)/1000,[293,max(Temp)],'--k')
        ylabel('[kJ / mol]')
        yyaxis right
        fplot(S,[293,max(Temp)],'--r')
        ylabel('[J / (mol * K)]')
        legend('G','H','S','location','best')
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'r';
        title(['Fit Cp ; Ho = ' num2str(H0/1000 + diffH/1000) ' [kJ/mol]; So = ' num2str(S0) ' [J/(mol K)]'])
        xlabel('Temperature [K]')
        box on
        grid on
        
        figure
        hold on
        fplot(subs(dG_ein,Ho,H0)/1000,[293,max(Temp)],'-k')
        fplot(subs(dH_ein,Ho,H0)/1000,[293,max(Temp)],'--k')
        ylabel('[kJ / mol]')
        yyaxis right
        fplot(dS_ein,[293,max(Temp)],'--r')
        ylabel('[J / (mol * K)]')
        legend('G','H','S','location','best')
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'r';
        title(['Einstien Model ; \Theta = ' num2str(theta) '; Ho = ' num2str(H0/1000) ' [kJ/mol]'])
        grid on
        box on
        xlabel('Temperature [K]')
        
    end
end