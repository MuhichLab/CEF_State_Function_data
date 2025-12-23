function [G_soln,G_ex,G_o] = CEF(xc,yc)

%{ 
    Pass-ins:
    xc == 1 if substition on A lattice - 0 otherwise
    yc == 1 if substition on B lattice - 0 otherwise

    Pass-outs:
    G_soln == Full solution model (symbolic function)
    G_ex   == Model of only excess terms (symbolic function)
    G_o    == Model of only endmebers and configurational entropy (symbolic function)

    Notes:
    This function creates a symbolic matlab expression for the CEF expression
    of an AA'BB'O3 -- 3 lattice 2 compostions on each lattice no reduction on
    A lattice

    The chemical eqution form is: 
    A(1-x)A'xB(1-y)B'yO3-z

    The paramters in the symbolic fucntions at the end of this script are 
    not sequential which may make things difficult later. We fix that with
    a different script.
%}
    
    % Constants
    m = 5/6*log(5/6) + 1/6*log(1/6); % 1/6 comes from the selction of 0.5 delta(z) as an endmember 
    r = 8.617333262E-5; % eV/K;

    syms x y z
    % Site fractions
    % these expressions are in terms of mol frac on A lattice (x) mol 
    % fraction on B lattice (y) and extent of reduction delta (z).
    
    Ya = 1-x; 
    Yap = x;             
    Ybr = 2*z*(1-y);
    Ybo = (1-2*z)*(1-y);
    Ybpr = 2*z*y;
    Ybpo = (1-2*z)*y;   
    Yo_o = 1-z/3;
    Yo_Va = z/3;    

    %% Entropy
    % Need to account if there is no x or no y
    syms R
    
    if xc~=0 && yc~=0
        S = -R*(Ya*log(Ya) + Yap*log(Yap) + Ybr*log(Ybr) + Ybo*log(Ybo)...
            + Ybpr*log(Ybpr) + Ybpo*log(Ybpo))...
            - 3*R*(Yo_o*log(Yo_o)+Yo_Va*log(Yo_Va));
        % sum r defines above allows easier conversion from ev to kJ if desired
        S = subs(S,R,r);
        
    elseif xc==0 && yc~=0
        S = -R*(Ybr*log(Ybr) + Ybo*log(Ybo)...
            + Ybpr*log(Ybpr) + Ybpo*log(Ybpo))...
            - 3*R*(Yo_o*log(Yo_o)+Yo_Va*log(Yo_Va));
        % sum r defines above allows easier conversion from ev to kJ if desired
        S = subs(S,[x R],[0 r]);
        
    elseif yc==0 && xc~=0
        S = -R*(Ya*log(Ya) + Yap*log(Yap) + Ybr*log(Ybr) + Ybo*log(Ybo))...
            - 3*R*(Yo_o*log(Yo_o)+Yo_Va*log(Yo_Va));
        % sum r defines above allows easier conversion from ev to kJ if desired
        S = subs(S,[y R],[0 r]);
        
    elseif xc==0 && yc==0
        S = -R*(Ybr*log(Ybr) + Ybo*log(Ybo))...
            - 3*R*(Yo_o*log(Yo_o)+Yo_Va*log(Yo_Va));
        % sum r defines above allows easier conversion from ev to kJ if desired
        S = subs(S,[ x y R],[0 0 r]);
    end

    %% Endmember
    % These are currently hard coded for all possbile endmember terms for a
    % AA'BB'O3 system
    
    syms G [1 8] 
    syms Go G_ABrO G_ApBrO G_ABprO G_ApBprO T M  
    
    % Here M = (5/6*log(5/6)+1/6*log(1/6)) is used to simplify the script
    % and m is plugged in for M later.
    
    % endmebers where B sites are not reduced 
    G_ABO = G1;
    G_ABVa = G1 - 3/2*Go;
    
    G_ApBO = G2;
    G_ApBVa = G2 - 3/2*Go;
    
    G_ABpO = G3;
    G_ABpVa = G3 - 3/2*Go;
    
    G_ApBpO = G4;
    G_ApBpVa = G4 - 3/2*Go;
    
    % endmebers where B sites are reduced
    G_ABrVa = G_ABrO - G_ABO + G_ABVa;
    G_ABrO = solve(((G5 - 1/6*G_ABrVa - 3*R*T*M)/(5/6))==G_ABrO,G_ABrO);
    
    G_ApBrVa = G_ApBrO - G_ApBO + G_ApBVa;
    G_ApBrO = solve(((G6 - 1/6*G_ApBrVa - 3*R*T*M)/(5/6))==G_ApBrO,G_ApBrO);
    
    G_ABprVa = G_ABprO - G_ABpO + G_ABpVa;
    G_ABprO = solve(((G7 - 1/6*G_ABprVa - 3*R*T*M)/(5/6))==G_ABprO,G_ABprO);
    
    G_ApBprVa = G_ApBprO - G_ApBpO + G_ApBpVa;
    G_ApBprO = solve(((G8 - 1/6*G_ApBprVa - 3*R*T*M)/(5/6))==G_ApBprO,G_ApBprO); 

    G_ABrVa   = G_ABrO - G_ABO + G_ABVa;  
    G_ApBrVa  = G_ApBrO - G_ApBO + G_ApBVa;  
    G_ABprVa  = G_ABprO - G_ABpO + G_ABpVa; 
    G_ApBprVa = G_ApBprO - G_ApBpO + G_ApBpVa;
    
    % 16 possible endmembers
    G_end = simplify(Ya*Ybr*Yo_o*G_ABrO...
        + Ya*Ybr*Yo_Va*G_ABrVa...
        + Ya*Ybo*Yo_o*G_ABO...
        + Ya*Ybo*Yo_Va*G_ABVa...
        + Ya*Ybpr*Yo_o*G_ABprO...
        + Ya*Ybpr*Yo_Va*G_ABprVa...
        + Ya*Ybpo*Yo_o*G_ABpO...
        + Ya*Ybpo*Yo_Va*G_ABpVa...
        + Yap*Ybr*Yo_o*G_ApBrO...
        + Yap*Ybr*Yo_Va*G_ApBrVa...
        + Yap*Ybo*Yo_o*G_ApBO...
        + Yap*Ybo*Yo_Va*G_ApBVa...
        + Yap*Ybpr*Yo_o*G_ApBprO...
        + Yap*Ybpr*Yo_Va*G_ApBprVa...
        + Yap*Ybpo*Yo_o*G_ApBpO...
        + Yap*Ybpo*Yo_Va*G_ApBpVa);
    
    % sub m and r
    G_end = subs(G_end,[M R], [m r]);

    %% Excess
    % These are currently hard coded for all possbile excess terms for a
    % AA'BB'O3 system
    
    syms L [1 80] real

    G_ex = Yap*Ya*Ybr*Yo_o*(L1 + (Ya - Yap)*L2)...     % A Lattice interactions
        + Yap*Ya*Ybr*Yo_Va*(L3 + (Ya - Yap)*L4)...
        + Yap*Ya*Ybo*Yo_o*( L5 + (Ya - Yap)*L6)...
        + Yap*Ya*Ybo*Yo_Va*(L7 + (Ya - Yap)*L8)...
        + Yap*Ya*Ybpr*Yo_o*(L9 + (Ya - Yap)*L10)...        
        + Yap*Ya*Ybpr*Yo_Va*(L11 + (Ya - Yap)*L12)...    
        + Yap*Ya*Ybpo*Yo_o*(L13 + (Ya - Yap)*L14)...    
        + Yap*Ya*Ybpo*Yo_Va*(L15 + (Ya - Yap)*L16)...     
        + Ybr*Ybo*Ya*Yo_o*(L17 + (Ybr - Ybo)*L18)...   % B Lattice interactions   
        + Ybr*Ybo*Ya*Yo_Va*(L19 + (Ybr - Ybo)*L20)...  
        + Ybr*Ybo*Yap*Yo_o*(L21 + (Ybr - Ybo)*L22)...
        + Ybr*Ybo*Yap*Yo_Va*(L23 + (Ybr - Ybo)*L24)...
        + Ybpr*Ybpo*Ya*Yo_o*(L25 + (Ybpr - Ybpo)*L26)...
        + Ybpr*Ybpo*Ya*Yo_Va*(L27 + (Ybpr - Ybpo)*L28)...        
        + Ybpr*Ybpo*Yap*Yo_o*(L29 + (Ybpr - Ybpo)*L30)...    
        + Ybpr*Ybpo*Yap*Yo_Va*(L31 + (Ybpr - Ybpo)*L32)...    
        + Ybpr*Ybr*Ya*Yo_o*(L33 + (Ybr - Ybpr)*L34)...    
        + Ybpr*Ybr*Ya*Yo_Va*(L35 + (Ybr - Ybpr)*L36)...   
        + Ybpr*Ybr*Yap*Yo_o*(L37 + (Ybr - Ybpr)*L38)...          
        + Ybpr*Ybr*Yap*Yo_Va*(L39 + (Ybr - Ybpr)*L40)...
        + Ybpo*Ybr*Ya*Yo_o*(L41 + (Ybr - Ybpo)*L42)...
        + Ybpo*Ybr*Ya*Yo_Va*(L43 + (Ybr - Ybpo)*L44)...
        + Ybpo*Ybr*Yap*Yo_o*(L45 + (Ybr - Ybpo)*L46)...        
        + Ybpo*Ybr*Yap*Yo_Va*(L47 + (Ybr - Ybpo)*L48)...    
        + Ybpo*Ybo*Ya*Yo_o*(L49 + (Ybo - Ybpo)*L50)...
        + Ybpo*Ybo*Ya*Yo_Va*(L51 + (Ybo - Ybpo)*L52)...
        + Ybpo*Ybo*Yap*Yo_o*(L53 + (Ybo - Ybpo)*L54)...        
        + Ybpo*Ybo*Yap*Yo_Va*(L55 + (Ybo - Ybpo)*L56)...         
        + Ybpr*Ybo*Ya*Yo_o*(L57 + (Ybpr - Ybo)*L58)...
        + Ybpr*Ybo*Ya*Yo_Va*(L59 + (Ybpr - Ybo)*L60)...
        + Ybpr*Ybo*Yap*Yo_o*(L61 + (Ybpr - Ybo)*L62)...        
        + Ybpr*Ybo*Yap*Yo_Va*(L63 + (Ybpr - Ybo)*L64)...     
        + Yo_o*Yo_Va*Ya*Ybr*(L65 + (Yo_o - Yo_Va)*L66)... % O Lattice interactions
        + Yo_o*Yo_Va*Ya*Ybo*(L67 + (Yo_o - Yo_Va)*L68)...    
        + Yo_o*Yo_Va*Ya*Ybpr*(L69 + (Yo_o - Yo_Va)*L70)...   
        + Yo_o*Yo_Va*Ya*Ybpo*(L71 + (Yo_o - Yo_Va)*L72)...  
        + Yo_o*Yo_Va*Yap*Ybr*(L73 + (Yo_o - Yo_Va)*L74)...
        + Yo_o*Yo_Va*Yap*Ybo*(L75 + (Yo_o - Yo_Va)*L76)...
        + Yo_o*Yo_Va*Yap*Ybpr*(L77 + (Yo_o - Yo_Va)*L78)...
        + Yo_o*Yo_Va*Yap*Ybpo*(L79 + (Yo_o - Yo_Va)*L80);  
  
    %% Complete Solution Model 
    % Put everything together
    
    G_soln = G_end - T*S + G_ex;
    if xc==0 && yc~=0
        G_soln = subs(G_soln,x,0);
        G_ex = subs(G_ex,x,0);
    elseif yc==0 && xc~=0
        G_soln = subs(G_soln,y,0);
        G_ex = subs(G_ex,y,0);
    elseif xc==0 && yc==0
        G_soln = subs(G_soln,[x y],[0 0]);
        G_ex = subs(G_ex,[x y],[0 0]);
    end
    
    %% Model no excess terms we call this the zeroth order model G_o
    
    G_o = G_end - T*S;
    if xc==0 && yc~=0
        G_o = subs(G_o,x,0);
    elseif yc==0 && xc~=0
        G_o = subs(G_o,y,0); 
    elseif xc==0 && yc==0
        G_o = subs(G_o,[x y],[0 0]);    
    end
    
end
