function thermo_plots(dH,dS,seed,Temp,x_exp)

%% Thermo reduction plots
%  dH/ds
syms x y z T real


clrs = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
leg = {};
Ts = unique(Temp);
for i = 1:length(Ts)
    [H_o(i), S_o(i)] = get_O2_thermo(Ts(i));
    leg{i} = ['T= ', num2str(Ts(i))];
end

x_frac = unique(x_exp);
x_frac(x_frac==0) = 1e-6;
x_frac(x_frac==1) = 1-1e-6;

figure
for j = 1:length(x_frac)
    mol_x= x_frac(j);
    
    subplot(2,length(x_frac),j)
    hold on 
    for k = 1:length(Ts)
        H_redr = simplify(diff(subs(dH,[y T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
        fplot(H_redr*96.487,[0 0.5],'color',clrs(k),'linewidth',2.0,'DisplayName',['T = ', num2str(round(Ts(k)))]);
    end
    ylim([0 600])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
    title(['x = ',num2str(mol_x)])
    box on
    hold off
end
legend

for j = 1:length(x_frac)
    mol_x= x_frac(j);
    subplot(2,length(x_frac),j+length(x_frac))
    hold on
    for k = 1:length(Ts)
        S_redr = simplify(diff(subs(dS,[y T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
        fplot(S_redr*96.487*1000,[0 0.5],'color',clrs(k),'linewidth',2.0,'DisplayName',['T = ', num2str(round(Ts(k)))]);
    end
    ylim([0 250])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
    box on
    hold off
end

end