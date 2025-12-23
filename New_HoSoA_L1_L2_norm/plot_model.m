function [G,H,S] = plot_model(model,fit_col)

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

    load(model,'Gsol');
    load(model,'all_fits');
    load(model,'end_terms');
    load(model,'ex_terms');

    G = subs(Gsol,[end_terms ex_terms],[all_fits(end-length(end_terms)+1:end,end-(fit_col)+1)' all_fits(1:length(ex_terms),end-(fit_col)+1)']);
    S = -diff(G,T);
    H = G + T*S;

    %% Thermo reduction plots
    %  dH/ds

    figure
    subplot(2,3,1)
    hold on
    x_frac = 1e-6;

    [Ho_1, So_1] = get_O2_thermo(400+273.15);
    [Ho_2, So_2] = get_O2_thermo(575+273.15);
    [Ho_3, So_3] = get_O2_thermo(750+273.15);
    [Ho_4, So_4] = get_O2_thermo(925+273.15);
    [Ho_5, So_5] = get_O2_thermo(1100+273.15);
    [Ho_6, So_6] = get_O2_thermo(1200+273.15);


    H_red_400 = simplify(diff(subs(H,[x T],[x_frac (400+273.15)]),z)+ 0.5*Ho_1/96.4875);
    H_red_575 = simplify(diff(subs(H,[x T],[x_frac (575+273.15)]),z)+ 0.5*Ho_2/96.4875);
    H_red_750 = simplify(diff(subs(H,[x T],[x_frac (750+273.15)]),z)+ 0.5*Ho_3/96.4875);
    H_red_925 = simplify(diff(subs(H,[x T],[x_frac (925+273.15)]),z)+ 0.5*Ho_4/96.4875);
    H_red_1100 = simplify(diff(subs(H,[x T],[x_frac (1100+273.15)]),z)+ 0.5*Ho_5/96.4875);

    hold on
    fplot(H_red_1100*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_925*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_750*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_575*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_400*96.487,[0 0.5],'linewidth',2.0);
    ylim([0 250])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
    title('x = 0.0')
    box on
    hold off

    subplot(2,3,2)
    hold on
    x_frac = 0.5;

    H_red_400 = simplify(diff(subs(H,[x T],[x_frac (400+273.15)]),z)+ 0.5*Ho_1/96.4875);
    H_red_575 = simplify(diff(subs(H,[x T],[x_frac (575+273.15)]),z)+ 0.5*Ho_2/96.4875);
    H_red_750 = simplify(diff(subs(H,[x T],[x_frac (750+273.15)]),z)+ 0.5*Ho_3/96.4875);
    H_red_925 = simplify(diff(subs(H,[x T],[x_frac (925+273.15)]),z)+ 0.5*Ho_4/96.4875);
    H_red_1100 = simplify(diff(subs(H,[x T],[x_frac (1100+273.15)]),z)+ 0.5*Ho_5/96.4875);

    fplot(H_red_1100*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_925*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_750*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_575*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_400*96.487,[0 0.5],'linewidth',2.0);
    ylim([0 550])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
    title('x = 0.5')
    box on
    hold off

    subplot(2,3,3)
    hold on
    x_frac = 1.0-1e-6;

    H_red_400 = simplify(diff(subs(H,[x T],[x_frac (400+273.15)]),z)+ 0.5*Ho_1/96.4875);
    H_red_575 = simplify(diff(subs(H,[x T],[x_frac (575+273.15)]),z)+ 0.5*Ho_2/96.4875);
    H_red_750 = simplify(diff(subs(H,[x T],[x_frac (750+273.15)]),z)+ 0.5*Ho_3/96.4875);
    H_red_925 = simplify(diff(subs(H,[x T],[x_frac (925+273.15)]),z)+ 0.5*Ho_4/96.4875);
    H_red_1100 = simplify(diff(subs(H,[x T],[x_frac (1100+273.15)]),z)+ 0.5*Ho_5/96.4875);

    fplot(H_red_1100*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_925*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_750*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_575*96.487,[0 0.5],'linewidth',2.0);
    fplot(H_red_400*96.487,[0 0.5],'linewidth',2.0);
    ylim([0 550])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
    title('x = 1.0')
    box on
    hold off

    % dS/ds 

    subplot(2,3,4)
    hold on
    x_frac = 1e-6;

    S_red_400 = simplify(diff(subs(S,[x T],[x_frac (400+273.15)]),z) + 0.5*So_1/96.4875);
    S_red_575 = simplify(diff(subs(S,[x T],[x_frac (575+273.15)]),z) + 0.5*So_2/96.4875);
    S_red_750 = simplify(diff(subs(S,[x T],[x_frac (750+273.15)]),z) + 0.5*So_3/96.4875);
    S_red_925 = simplify(diff(subs(S,[x T],[x_frac (925+273.15)]),z) + 0.5*So_4/96.4875);
    S_red_1100 = simplify(diff(subs(S,[x T],[x_frac (1100+273.15)]),z) + 0.5*So_5/96.4875);
    fplot(S_red_1100*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_925*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_750*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_575*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_400*96.487*1000,[0 0.5],'linewidth',2.0);
    ylim([0 250])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
    box on
    hold off

    subplot(2,3,5)
    hold on
    x_frac = 0.5;

    S_red_400 = simplify(diff(subs(S,[x T],[x_frac (400+273.15)]),z) + 0.5*So_1/96.4875);
    S_red_575 = simplify(diff(subs(S,[x T],[x_frac (575+273.15)]),z) + 0.5*So_2/96.4875);
    S_red_750 = simplify(diff(subs(S,[x T],[x_frac (750+273.15)]),z) + 0.5*So_3/96.4875);
    S_red_925 = simplify(diff(subs(S,[x T],[x_frac (925+273.15)]),z) + 0.5*So_4/96.4875);
    S_red_1100 = simplify(diff(subs(S,[x T],[x_frac (1100+273.15)]),z) + 0.5*So_5/96.4875);
    fplot(S_red_1100*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_925*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_750*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_575*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_400*96.487*1000,[0 0.5],'linewidth',2.0);
    ylim([0 250])
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
    box on
    hold off

    subplot(2,3,6)
    hold on
    x_frac = 1-1e-6;

    S_red_400 = simplify(diff(subs(S,[x T],[x_frac (400+273.15)]),z) + 0.5*So_1/96.4875);
    S_red_575 = simplify(diff(subs(S,[x T],[x_frac (575+273.15)]),z) + 0.5*So_2/96.4875);
    S_red_750 = simplify(diff(subs(S,[x T],[x_frac (750+273.15)]),z) + 0.5*So_3/96.4875);
    S_red_925 = simplify(diff(subs(S,[x T],[x_frac (925+273.15)]),z) + 0.5*So_4/96.4875);
    S_red_1100 = simplify(diff(subs(S,[x T],[x_frac (1100+273.15)]),z) + 0.5*So_5/96.4875);
    fplot(S_red_1100*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_925*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_750*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_575*96.487*1000,[0 0.5],'linewidth',2.0);
    fplot(S_red_400*96.487*1000,[0 0.5],'linewidth',2.0);
    ylim([0 250])

    %title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
    xlabel('Extent of Reduction (\delta)')
    ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
    box on
    hold off
    %title (['Full Model Ba_xSr_1-xFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
    legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')

    %% Cp
    cp1 = double(subs(diff(H*96.487,T),[x z T],[1e-6 1e-6 800]))*1000;
    fprintf('Here is the "heat capacity" of ABO3: %.4f \n',cp1)

    cp2 = double(subs(diff(H*96.487,T),[x z],[1 0]))*1000; 
    fprintf('Here is the "heat capacity of ABpO3: %.4f \n',cp2)

    cp3 = double(subs(diff(H*96.487,T),[x z],[0 0.5]))*1000;
    fprintf('Here is the "heat capacity" of ABO2.5: %.4f \n',cp3)

    cp4 = double(subs(diff(H*96.487,T),[x z],[1 0.5]))*1000;
    fprintf('Here is the "heat capacity" of ABpO2.5: %.4f \n',cp4)

    %% Entropy Check

    xx = linspace(1e-6,1.0-1e-6,100);
    zz = linspace(1e-6,0.5-1e-6,100);
    tt = linspace(300,2500,100);
    [XX, ZZ, TT] = meshgrid(xx,zz,tt);

    XX = reshape(XX,[1,100^3]);
    ZZ = reshape(ZZ,[1,100^3]);
    TT = reshape(TT,[1,100^3]);

    testS = matlabFunction(S*96.5*1000);

    stest = testS(TT,XX,ZZ);

    fprintf('The minimum Entropy value of the mesh is %.2f [J/(mol K)] \n',min(stest))
end