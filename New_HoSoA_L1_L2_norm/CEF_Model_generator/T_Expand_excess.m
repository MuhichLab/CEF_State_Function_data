function [Gsol,Gdft,Gex] = T_Expand_excess(G_sol,G_dft,G_ex,cp)

    % These are set specifically for the AA'BB'O3 perovskite with reduction
    % on B lattice only
    
    syms T real
    syms Ho [1 80] real
    syms So [1 80] real
    syms A [1 80] real
    syms B [1 80] real
    syms C [1 80] real
    syms D [1 80] real
    syms E [1 80] real 
    syms L [1 80] real
    
     %% Cp model
     
     cp_exp = A + T*B + T^2*C + T^3*D + E/T^2;
     
     % H, S, and G
     dH = int(cp_exp,T) + Ho;   % h is constant from integration
     dS = int(cp_exp/T,T) + So; % s is constant from integration
     dG = expand(dH - T*dS);

 %% Expand L terms in excess to be G(T) terms

    Gsol = subs(G_sol,L,dG);
    Gex = subs(G_ex,L,dG);
    
    % DFT is at T=0
    Gdft = subs(G_dft,L,Ho);

     %% Correct for number of desire heat capacity parameters
    if cp == 1
        Gsol = subs(Gsol,[B C D E],zeros(1,length(L)*4));
        Gex = subs(Gex,[B C D E],zeros(1,length(L)*4));
    elseif cp == 2
        Gsol = subs(Gsol,[C D E],zeros(1,length(L)*3));
        Gex = subs(Gex,[C D E],zeros(1,length(L)*3));       
    elseif cp == 3 
        Gsol = subs(Gsol,[D E],zeros(1,length(L)*2));
        Gex = subs(Gex,[D E],zeros(1,length(L)*2));
    elseif cp == 4 
        Gsol = subs(Gsol,E,zeros(1,length(L)));
        Gex = subs(Gex,E,zeros(1,length(L)));
    end   
    
end
