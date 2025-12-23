function compare_models_3(model1,model2,model3)

    syms x y z T real

    S1 = -diff(model1,T);
    H1 = model1 + T*S1;
    
    S2 = -diff(model2,T);
    H2 = model2 + T*S2;

    S3 = -diff(model3,T);
    H3 = model3 + T*S3;

    clrs = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
    leg = {};
    Ts = linspace(673,1373,3);
    for i = 1:length(Ts)
        [H_o(i), S_o(i)] = get_O2_thermo(Ts(i));
        leg{i} = ['T= ', num2str(Ts(i))];
    end
    
    x_frac = linspace(0,1,3);
    x_frac(x_frac==0) = 1e-6;
    x_frac(x_frac==1) = 1-1e-6;
    
    figure
    for j = 1:length(x_frac)
        mol_x= x_frac(j);
        subplot(1,length(x_frac),j)
        hold on 
        for k = 1:length(Ts)
            H_red1 = simplify(diff(subs(H1,[x T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
            H_red2 = simplify(diff(subs(H2,[x T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
            H_red3 = simplify(diff(subs(H3,[x T],[mol_x Ts(k)]),z)+ 0.5*H_o(k)/96.4875);
            fplot(H_red1*96.487,[0 0.5],'--','color',clrs(k),'linewidth',2.0);
            fplot(H_red2*96.487,[0 0.5],':','color',clrs(k),'linewidth',1.5);
            fplot(H_red3*96.487,[0 0.5],'-.','color',clrs(k),'linewidth',1.0);
        end
        ylim([0 200])
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
            S_red1 = simplify(diff(subs(S1,[x T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
            S_red2 = simplify(diff(subs(S2,[x T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
            S_red3 = simplify(diff(subs(S3,[x T],[mol_x Ts(k)]),z) + 0.5*S_o(k)/96.4875);
            fplot(S_red1*96.487*1000,[0 0.5],'--','color',clrs(k),'linewidth',2.0);
            fplot(S_red2*96.487*1000,[0 0.5],':','color',clrs(k),'linewidth',1.5);
            fplot(S_red3*96.487*1000,[0 0.5],'-.','color',clrs(k),'linewidth',1.0);
        end
        ylim([0 250])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off
        
    end
    legend(num2str(Ts))

end