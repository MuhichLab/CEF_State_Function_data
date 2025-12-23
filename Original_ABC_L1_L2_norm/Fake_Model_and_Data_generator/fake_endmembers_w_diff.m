function [fake_end_terms] = fake_endmembers_w_diff(gs,end_terms)

    syms gDiffA [1 4] real  
    syms gDiffHo [1 4] real
    syms gDiffSo [1 4] real
    
    conv = 96.487; %  1 eV/K == 96.487 kJ/(mol K)
    
    gs_j = gs*conv*1000; %convert eV values to J/mol
    
    % get endmeber terms from eistein fits usign gs vaules for H0
    [~,~,~,coeffs1,S01,diffH1] = fit_einstein(700,gs_j(1),5);
    [~,~,~,coeffs2,S02,diffH2] = fit_einstein(600,gs_j(2),5);
    [~,~,~,coeffs3,S03,diffH3] = fit_einstein(450,gs_j(5),4.5);
    [~,~,~,coeffs4,S04,diffH4] = fit_einstein(500,gs_j(6),4.5);
    
    
    fake_end_terms = zeros(1,length(end_terms));
    
    A0s_end = ismember(end_terms,gDiffA);
    H0s_end = ismember(end_terms,gDiffHo);
    S0s_end = ismember(end_terms,gDiffSo);    
    
    fake_end_terms(A0s_end) = [(coeffs3(1)-coeffs1(1)) (coeffs4(1)-coeffs2(1))];
    fake_end_terms(H0s_end) = [(diffH3-diffH1) (diffH4-diffH2)];
    fake_end_terms(S0s_end) = [(S03-S01) (S04-S02)];
    
    fake_end_terms = fake_end_terms/conv/1000;
    
end