function VH_compare(model)

    syms x y z T To real
    syms G [1 8] real
    syms gendA [1 8] real
    syms gendB [1 8] real
    syms gendC [1 8] real
    syms gendD [1 8] real
    syms gendHo [1 8] real
    syms gendSo [1 8] real
    syms Ho [1 80] real
    syms So [1 80] real
    syms A [1 80] real
    syms B [1 80] real
    syms C [1 80] real
    syms D [1 80] real
    syms L [1 80] real
    
    load(model,'Gnew');
    
    S = -diff(Gnew,T);
    H = Gnew + T*S;
    
    % VH Data from NC state
    VH_data = readtable('BSF_NCstate_VH.xlsx');
    x_exp = table2array(VH_data(:,1));  
    VH_delta = table2array(VH_data(:,2));   
    VH_H = table2array(VH_data(:,3)); 
    VH_H_err = table2array(VH_data(:,4));  
    VH_S = table2array(VH_data(:,5));  
    VH_S_err = table2array(VH_data(:,6));  
  
    %% Thermo reduction plots
    %  dH/ds

    
    Ts = linspace(673,1673,6);
    xu = unique(x_exp);
    
    figure
    for n = 1:length(xu)
        xs = xu(n);
        subplot(2,length(xu),n)
        hold on
        for m = 1:length(Ts)
            ts = Ts(m);
            [Ho2, ~] = get_O2_thermo(ts);
            H_red = simplify(diff(subs(H,[x T],[xs ts]),z)+ 0.5*Ho2/96.4875);
            fplot(H_red*96.487,[0 0.5],'linewidth',2.0,'DisplayName',[num2str(ts),' K']);
        end
        errorbar(VH_delta(x_exp==xs),VH_H(x_exp==xs),VH_H_err(x_exp==xs),'o-k','DisplayName','VH')
        ylim([50 110])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title(['x = ', num2str(xs)])
        box on
        hold off
        legend('location','best')
        
    end
  
    for n = 1:length(xu)
        xs = xu(n);
        subplot(2,length(xu),n+length(xu))
        hold on
        for m = 1:length(Ts)
            ts = Ts(m);
            [~, So2] = get_O2_thermo(ts);
            S_red = simplify(diff(subs(S,[x T],[xs ts]),z) + 0.5*So2/96.4875);
            fplot(S_red*96.487*1000,[0 0.5],'linewidth',2.0,'DisplayName',[num2str(ts),' K']);
        end
        errorbar(VH_delta(x_exp==xs),VH_S(x_exp==xs),VH_S_err(x_exp==xs),'o-k','DisplayName','VH')
        ylim([0 250])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off  
    end    
end
