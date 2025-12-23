function compare_models(model_truth,model_cef)

    syms x y z T real

    S_truth = -diff(model_truth,T);
    H_truth = model_truth + T*S_truth;
    
    S_cef = -diff(model_cef,T);
    H_cef = model_cef + T*S_cef;

    clrs = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
    leg = {};
    Ts = linspace(500,1800,7);
    for i = 1:length(Ts)
        [H_o(i), S_o(i)] = get_O2_thermo(Ts(i));
        leg{i} = ['T= ', num2str(Ts(i))];
    end
    
    x_frac = linspace(0,1,5);
    x_frac(x_frac==0) = 1e-6;
    x_frac(x_frac==1) = 1-1e-6;
    
    figure
    for j = 1:length(x_frac)
        mol_x= x_frac(j);
        subplot(1,length(x_frac),j)
        hold on 
        for k = 1:length(Ts)
            H_redr = simplify(diff(subs(H_truth,[x T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
            H_redc = simplify(diff(subs(H_cef,[x T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
            fplot(H_redr*96.487,[0 0.5],'color',clrs(k),'linewidth',2.0);
            fplot(H_redc*96.487,[0 0.5],':','color',clrs(k),'linewidth',2.0);
        end
        ylim([0 150])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title(['x = ',num2str(mol_x)])
        box on
        hold off
    end
        
    figure
    for j = 1:length(x_frac)
        mol_x= x_frac(j);
        subplot(1,length(x_frac),j)
        hold on
        for k = 1:length(Ts)
            S_redr = simplify(diff(subs(S_truth,[x T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
            S_redc = simplify(diff(subs(S_cef,[x T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
            fplot(S_redr*96.487*1000,[0 0.5],'color',clrs(k),'linewidth',2.0);
            fplot(S_redc*96.487*1000,[0 0.5],':','color',clrs(k),'linewidth',2.0);
        end
        ylim([0 200])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off
        
    end
    legend(num2str(Ts))

end