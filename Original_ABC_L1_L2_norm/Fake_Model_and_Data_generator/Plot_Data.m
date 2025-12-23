function Plot_Data(Press,Temp,x_exp,y_exp,dd)
    %% Exp Data Plots

%     figure
%     scatter3(x_exp,dd,Temp,100,Press)
%     zlabel('Temperature (K)');
%     ylabel('\delta');
%     xlabel('Mol Fraction (x)','Rotation',10);
%     view([-46.8640862213933 8.20408561922651]);
%     box on
%     grid on
%     cb = colorbar;
%     colormap jet
%     set(gca,'ColorScale','log')
%     cb.Label.String = 'P_{O2} (Bar)';
%     cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
%     %cb.Limits = [0.01 1];
%     
%     figure
%     scatter3(Press,dd,Temp,100,x_exp)
%     zlabel('Temperature (K)');
%     ylabel('\delta');
%     xlabel('P_{O2} (Bar)');
%     view([-46.8640862213933 8.20408561922651]);
%     box on
%     grid on
%     cb = colorbar;
%     colormap jet
%     set(gca, 'XScale', 'log')
% %     set(gca,'ColorScale','log')
% %     cb.Label.String = 'P_{O2} (Bar)';
%     cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
%     %cb.Limits = [0.01 1];
%     
%     figure
%     scatter3(Temp,dd,Press,100,x_exp)
%     zlabel('Temperature (K)');
%     ylabel('\delta');
%     xlabel('P_{O2} (Bar)');
%     view([-46.8640862213933 8.20408561922651]);
%     box on
%     grid on
%     cb = colorbar;
%     colormap jet
%     set(gca, 'ZScale', 'log')
% %     set(gca,'ColorScale','log')
% %     cb.Label.String = 'P_{O2} (Bar)';
%     cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
%     %cb.Limits = [0.01 1];
%     
%     figure
%     scatter3(Temp,Press,dd,100,x_exp)
%     zlabel('Temperature (K)');
%     ylabel('\delta');
%     xlabel('P_{O2} (Bar)');
%     view([-46.8640862213933 8.20408561922651]);
%     box on
%     grid on
%     cb = colorbar;
%     colormap jet
%     set(gca, 'YScale', 'log')
% %     set(gca,'ColorScale','log')
% %     cb.Label.String = 'P_{O2} (Bar)';
%     cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
%     %cb.Limits = [0.01 1];
%     
%     figure
%     scatter3(Temp,Press,3-dd,100,x_exp)
%     zlabel('Temperature (K)');
%     ylabel('\delta');
%     xlabel('P_{O2} (Bar)');
%     view([-46.8640862213933 8.20408561922651]);
%     box on
%     grid on
%     cb = colorbar;
%     colormap jet
%     set(gca, 'YScale', 'log')
% %     set(gca,'ColorScale','log')
% %     cb.Label.String = 'P_{O2} (Bar)';
%     cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
%     %cb.Limits = [0.01 1];
    

loopy = unique(y_exp);
for y = 1:length(loopy)
    figure
    scatter3(x_exp(y_exp==loopy(y)),Press(y_exp==loopy(y)),3-dd(y_exp==loopy(y)),100,Temp(y_exp==loopy(y)))
    zlabel('3 - \delta');
    ylabel('P_{O2} (Bar)');
    xlabel('Mol Fraction (x)');
    view([-46.8640862213933 8.20408561922651]);
    box on
    grid on
    cb = colorbar;
    colormap jet
    set(gca, 'YScale', 'log')
%     set(gca,'ColorScale','log')
    cb.Label.String = 'Temperature (K)';
    cb.Position = [0.846153846153846 0.267728706624605 0.0237922810114779 0.55403785488959];
    sgtitle(['A_{x}A^{\prime}_{1-x}B_{y}B^{\prime}_{1-x}O_{3-\delta} : y = ' num2str(loopy(y))]);
    %cb.Limits = [0.01 1];
end
    
end