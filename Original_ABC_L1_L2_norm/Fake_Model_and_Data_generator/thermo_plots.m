function thermo_plots(dH,dS,seed)

%% Thermo reduction plots
%  dH/ds
    syms x y z T real

    clrs = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];

    for y_frac = 1e-6:0.2-1e-6:1.0-1e-6
    
        figure
        subplot(2,6,1)
        hold on
        x_frac = 1e-6;

        [Ho_1, So_1] = get_O2_thermo(800+273.15);
        [Ho_2, So_2] = get_O2_thermo(1000+273.15);
        [Ho_3, So_3] = get_O2_thermo(1200+273.15);
        [Ho_4, So_4] = get_O2_thermo(1400+273.15);
        [Ho_5, So_5] = get_O2_thermo(1600+273.15);

        H_red_800 =  simplify(diff(subs(dH,[x y T],[x_frac y_frac (800 +273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);

        hold on
        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);
        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 0.0')
        box on
        hold off

        subplot(2,6,2)
        hold on
        x_frac = 0.2;

        H_red_800 = simplify(diff(subs(dH,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);


        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 0.2')
        box on
        hold off

        subplot(2,6,3)
        hold on
        x_frac = 0.4;

        H_red_800 = simplify(diff(subs(dH,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);


        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 0.4')
        box on
        hold off

        subplot(2,6,4)
        hold on
        x_frac = 0.6;

        H_red_800 = simplify(diff(subs(dH,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);



        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 0.6')
        box on
        hold off


        subplot(2,6,5)
        hold on
        x_frac = 0.8;

        H_red_800 = simplify(diff(subs(dH,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);


        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 0.8')
        box on
        hold off


        subplot(2,6,6)
        hold on
        x_frac = 1.0-1e-6;

        H_red_800 = simplify(diff(subs(dH,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*Ho_1/96.4875);
        H_red_1000 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*Ho_2/96.4875);
        H_red_1200 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*Ho_3/96.4875);
        H_red_1400 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*Ho_4/96.4875);
        H_red_1600 = simplify(diff(subs(dH,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*Ho_5/96.4875);


        fplot(H_red_1600*96.487,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(H_red_1400*96.487,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(H_red_1200*96.487,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(H_red_1000*96.487,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(H_red_800*96.487,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 600])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialH/\partial\delta (kJ \bullet (mol O)^{-1})')
        title('x = 1.0')
        box on
        hold off

        %dS/ds 

        subplot(2,6,7)
        hold on
        x_frac = 1e-6;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);


        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off

        subplot(2,6,8)
        hold on
        x_frac = 0.2;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);

        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off


        subplot(2,6,9)
        hold on
        x_frac = 0.4;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);


        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off

        subplot(2,6,10)
        hold on
        x_frac = 0.6;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);

        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off


        subplot(2,6,11)
        hold on
        x_frac = 0.8;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);

        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off

        subplot(2,6,12)
        hold on
        x_frac = 1-1e-6;

        S_red_800 = simplify(diff(subs(dS,[x y T],[x_frac y_frac  (800+273.15)]),z)+ 0.5*So_1/96.4875);
        S_red_1000 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1000+273.15)]),z)+ 0.5*So_2/96.4875);
        S_red_1200 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1200+273.15)]),z)+ 0.5*So_3/96.4875);
        S_red_1400 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1400+273.15)]),z)+ 0.5*So_4/96.4875);
        S_red_1600 = simplify(diff(subs(dS,[x y T],[x_frac y_frac (1600+273.15)]),z)+ 0.5*So_5/96.4875);



        fplot(S_red_1600*96.487*1000,[0 0.5],'Color',clrs(1),'linewidth',2.0);
        fplot(S_red_1400*96.487*1000,[0 0.5],'Color',clrs(2),'linewidth',2.0);
        fplot(S_red_1200*96.487*1000,[0 0.5],'Color',clrs(3),'linewidth',2.0);
        fplot(S_red_1000*96.487*1000,[0 0.5],'Color',clrs(4),'linewidth',2.0);
        fplot(S_red_800*96.487*1000,[0 0.5],'Color',clrs(5),'linewidth',2.0);

        ylim([0 250])
        xlim([0 0.5])

        %title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
        xlabel('Extent of Reduction (\delta)')
        ylabel(' \partialS/\partial\delta (J \bullet (mol O \bullet K)^{-1})')
        box on
        hold off
        sgtitle(['A_{x}A^{\prime}_{1-x}B_{y}B^{\prime}_{1-x}O_{3-\delta} : y = ' num2str(y_frac) ' : SEED = ' num2str(seed)]);
        legend('1600 ^\circC','1400 ^\circC','1200 ^\circC','1000 ^\circC','800 ^\circC')
    end

    
%     FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
%     for iFig = 1:length(FigList)
%       FigHandle = FigList(iFig);
%       FigName   = ['Thermo_' num2str(seed) '_' num2str(iFig)];
%       savefig(FigHandle, [FigName , '.fig']);
%       saveas(FigHandle,  [FigName , '.jpg']);
%     end
%     
%     close all
end