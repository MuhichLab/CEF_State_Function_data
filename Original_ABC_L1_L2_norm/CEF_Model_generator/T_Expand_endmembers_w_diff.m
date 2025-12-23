function [G_sol,G_dft,Go,Go_dft] = T_Expand_endmembers_w_diff(G_soln,G_o,cp)

    % These are set specifically for the AA'BB'O3 perovskite with reduction
    % on B lattice only

    syms T real
    syms G [1 8]
    syms gDiff [1 4] real
    syms gDiffA [1 4] real
    syms gDiffB [1 4] real
    syms gDiffC [1 4] real
    syms gDiffD [1 4] real
    syms gDiffE [1 4] real
    syms gDiffHo [1 4] real
    syms gDiffSo [1 4] real
    

    dG = gDiffHo + gDiffSo*T + (gDiffB*T^2)/2 + (gDiffC*T^3)/6 + (gDiffD*T^4)/12 + gDiffE/(2*T) + gDiffA*T*log(T);
    %% Full model subs 
    % First we build two full solution models that have the excess terms in
    % them, tho we do not do anyting with the excess terms.
    
    % Full Solution model

    G_sol = subs(G_soln,[G1 G2 G3 G4],...
        [(G5 - gDiff1) (G6 - gDiff2) (G7 - gDiff3) (G8 - gDiff4)]);
    
    G_sol = subs(G_sol,gDiff,dG);
     
    % Full Solution Model DFT Side (non-derivative no T dependece)
    % NOTE you can not sub 0 in for T in G_soln because you get NaN 
    % if C*T*ln(T) is present. So you must build from scratch
    
    G_dft = subs(G_soln,[G1 G2 G3 G4 T],...
        [(G5 - gDiff1) (G6 - gDiff2) (G7 - gDiff3) (G8 - gDiff4) 0]);

    G_dft = subs(G_dft,gDiff,gDiffHo);

    %% Go model subs Go = G_end - T*S

    % Go subsitutions
    Go = subs(G_o,[G1 G2 G3 G4],...
        [(G5 - gDiff1) (G6 - gDiff2) (G7 - gDiff3) (G8 - gDiff4)]);
  
    Go = subs(Go,gDiff,dG);
    
    % Go_dft
    Go_dft = subs(G_o,[G1 G2 G3 G4 T],...
        [(G5 - gDiff1) (G6 - gDiff2) (G7 - gDiff3) (G8 - gDiff4) 0]);

    Go_dft = subs(Go_dft,gDiff,gDiffHo);
    
    %% Correct for number of desire heat capacity parameters
    if cp == 1
        G_sol = subs(G_sol,[gDiffB gDiffC gDiffD gDiffE],zeros(1,16));
        Go = subs(Go,[gDiffB gDiffC gDiffD gDiffE],zeros(1,16));
    elseif cp == 2
        G_sol = subs(G_sol,[gDiffC gDiffD gDiffE],zeros(1,12));
        Go = subs(Go,[gDiffC gDiffD gDiffE],zeros(1,12));       
    elseif cp == 3 
        G_sol = subs(G_sol,[gDiffD gDiffE],zeros(1,8));
        Go = subs(Go,[gDiffD gDiffE],zeros(1,8));
    elseif cp == 4 
        G_sol = subs(G_sol,gDiffE,zeros(1,4));
        Go = subs(Go,gDiffE,zeros(1,4));        
    end
    

end
