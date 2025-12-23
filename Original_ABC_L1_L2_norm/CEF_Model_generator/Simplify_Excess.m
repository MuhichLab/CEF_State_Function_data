function [G_soln,G_ex,subin,subout] = Simplify_Excess(G_soln,G_ex,xc,yc)

%{ 
    Pass-ins:
    G_soln == Full solution model from CEF script (symbolic function)
    G_ex   == Model of only excess terms from CEF script (symbolic function)
    xc == 1 if substition on A lattice - 0 otherwise
    yc == 1 if substition on B lattice - 0 otherwise

    Pass-outs:
    G_soln == Full solution model with parallel excess terms removed
    G_ex   == Model of only excess terms with parallel excess terms removed
    subin  == excess terms left that were parallel to other excess terms
    subout == excess terms removed that were parallel to subin excess terms

    Notes:
    It is possible for some excess terms to be linearaly identical 
    (parallel) to one another thus making the optimization of those 
    paramters impossible. We remove one of the terms in the pair that are
    linear and note which terms they are with subin and subout.

    The paramters in the symbolic fucntions at the end of this script are 
    not sequential which may make things difficult later. We fix that with
    a different script.
%}

    %% Simpplify Excess terms
    % It is possible for some excess terms to be linearaly identical 
    % (parallel to one another) thus making the optimization of those 
    % paramters impossible. We remove one of the terms in the pair that are
    % linear and note which terms they are.
    
    if xc~=0 && yc~=0

        syms x y z
        % redefine site fractions based on actual compostition
        Ya = 1-x; 
        Yap = x;             
        Ybr = 2*z*(1-y);
        Ybo = (1-2*z)*(1-y);
        Ybpr = 2*z*y;
        Ybpo = (1-2*z)*y;   
        Yo_o = 1-z/3;
        Yo_Va = z/3; 
    
    elseif xc==0 && yc~=0
        
        x = 0;
        syms y z
        % redefine site fractions based on actual compostition
        Ya = 1-x; 
        Yap = x;             
        Ybr = 2*z*(1-y);
        Ybo = (1-2*z)*(1-y);
        Ybpr = 2*z*y;
        Ybpo = (1-2*z)*y;   
        Yo_o = 1-z/3;
        Yo_Va = z/3; 
        
    elseif yc==0 && xc~=0
        
        y = 0;
        syms x z
        % redefine site fractions based on actual compostition
        Ya = 1-x; 
        Yap = x;             
        Ybr = 2*z*(1-y);
        Ybo = (1-2*z)*(1-y);
        Ybpr = 2*z*y;
        Ybpo = (1-2*z)*y;   
        Yo_o = 1-z/3;
        Yo_Va = z/3;
        
    elseif xc==0 && yc==0
        
        y = 0;
        x = 0;
        syms z
        % redefine site fractions based on actual compostition
        Ya = 1-x; 
        Yap = x;             
        Ybr = 2*z*(1-y);
        Ybo = (1-2*z)*(1-y);
        Ybpr = 2*z*y;
        Ybpo = (1-2*z)*y;   
        Yo_o = 1-z/3;
        Yo_Va = z/3; 
        
    end

    % Define all possibel excess terms based on the given compostion (xc,
    % yc)
    l1 = Yap*Ya*Ybr*Yo_o;      l2 = Yap*Ya*Ybr*Yo_o*(Ya - Yap);
    l3 = Yap*Ya*Ybr*Yo_Va;     l4 = Yap*Ya*Ybr*Yo_Va*(Ya - Yap);
    l5 = Yap*Ya*Ybo*Yo_o;      l6 = Yap*Ya*Ybo*Yo_o*(Ya - Yap);
    l7 = Yap*Ya*Ybo*Yo_Va;     l8 = Yap*Ya*Ybo*Yo_Va*(Ya - Yap);
    l9 = Yap*Ya*Ybpr*Yo_o;     l10 = Yap*Ya*Ybpr*Yo_o*(Ya - Yap);      
    l11 = Yap*Ya*Ybpr*Yo_Va;   l12 = Yap*Ya*Ybpr*Yo_Va*(Ya - Yap);   
    l13 = Yap*Ya*Ybpo*Yo_o;    l14 = Yap*Ya*Ybpo*Yo_o*(Ya - Yap);
    l15 = Yap*Ya*Ybpo*Yo_Va;   l16 = Yap*Ya*Ybpo*Yo_Va*(Ya - Yap);   
    l17 = Ybr*Ybo*Ya*Yo_o;     l18 = Ybr*Ybo*Ya*Yo_o*(Ybr - Ybo);    
    l19 = Ybr*Ybo*Ya*Yo_Va;    l20 = Ybr*Ybo*Ya*Yo_Va*(Ybr - Ybo); 
    l21 = Ybr*Ybo*Yap*Yo_o;    l22 = Ybr*Ybo*Yap*Yo_o*(Ybr - Ybo);
    l23 = Ybr*Ybo*Yap*Yo_Va;   l24 = Ybr*Ybo*Yap*Yo_Va*(Ybr - Ybo);
    l25 = Ybpr*Ybpo*Ya*Yo_o;   l26 = Ybpr*Ybpo*Ya*Yo_o*(Ybpr - Ybpo);
    l27 = Ybpr*Ybpo*Ya*Yo_Va;  l28 = Ybpr*Ybpo*Ya*Yo_Va*(Ybpr - Ybpo);      
    l29 = Ybpr*Ybpo*Yap*Yo_o;  l30 = Ybpr*Ybpo*Yap*Yo_o*(Ybpr - Ybpo);   
    l31 = Ybpr*Ybpo*Yap*Yo_Va; l32 = Ybpr*Ybpo*Yap*Yo_Va*(Ybpr - Ybpo);   
    l33 = Ybpr*Ybr*Ya*Yo_o;    l34 = Ybpr*Ybr*Ya*Yo_o*(Ybr - Ybpr);    
    l35 = Ybpr*Ybr*Ya*Yo_Va;   l36 = Ybpr*Ybr*Ya*Yo_Va*(Ybr - Ybpr);  
    l37 = Ybpr*Ybr*Yap*Yo_o;   l38 = Ybpr*Ybr*Yap*Yo_o*(Ybr - Ybpr);        
    l39 = Ybpr*Ybr*Yap*Yo_Va;  l40 = Ybpr*Ybr*Yap*Yo_Va*(Ybr - Ybpr);
    l41 = Ybpo*Ybr*Ya*Yo_o;    l42 = Ybpo*Ybr*Ya*Yo_o*(Ybr - Ybpo);
    l43 = Ybpo*Ybr*Ya*Yo_Va;   l44 = Ybpo*Ybr*Ya*Yo_Va*(Ybr - Ybpo);
    l45 = Ybpo*Ybr*Yap*Yo_o;   l46 = Ybpo*Ybr*Yap*Yo_o*(Ybr - Ybpo);     
    l47 = Ybpo*Ybr*Yap*Yo_Va;  l48 = Ybpo*Ybr*Yap*Yo_Va*(Ybr - Ybpo);   
    l49 = Ybpo*Ybo*Ya*Yo_o;    l50 = Ybpo*Ybo*Ya*Yo_o*(Ybo - Ybpo);
    l51 = Ybpo*Ybo*Ya*Yo_Va;   l52 = Ybpo*Ybo*Ya*Yo_Va*(Ybo - Ybpo);
    l53 = Ybpo*Ybo*Yap*Yo_o;   l54 = Ybpo*Ybo*Yap*Yo_o*(Ybo - Ybpo);      
    l55 = Ybpo*Ybo*Yap*Yo_Va;  l56 = Ybpo*Ybo*Yap*Yo_Va*(Ybo - Ybpo);        
    l57 = Ybpr*Ybo*Ya*Yo_o;    l58 = Ybpr*Ybo*Ya*Yo_o*(Ybpr - Ybo);
    l59 = Ybpr*Ybo*Ya*Yo_Va;   l60 = Ybpr*Ybo*Ya*Yo_Va*(Ybpr - Ybo);
    l61 = Ybpr*Ybo*Yap*Yo_o;   l62 = Ybpr*Ybo*Yap*Yo_o*(Ybpr - Ybo);     
    l63 = Ybpr*Ybo*Yap*Yo_Va;  l64 = Ybpr*Ybo*Yap*Yo_Va*(Ybpr - Ybo); 
    l65 = Yo_o*Yo_Va*Ya*Ybr;   l66 = Yo_o*Yo_Va*Ya*Ybr*(Yo_o - Yo_Va); 
    l67 = Yo_o*Yo_Va*Ya*Ybo;   l68 = Yo_o*Yo_Va*Ya*Ybo*(Yo_o - Yo_Va);  
    l69 = Yo_o*Yo_Va*Ya*Ybpr;  l70 = Yo_o*Yo_Va*Ya*Ybpr*(Yo_o - Yo_Va);  
    l71 = Yo_o*Yo_Va*Ya*Ybpo;  l72 = Yo_o*Yo_Va*Ya*Ybpo*(Yo_o - Yo_Va); 
    l73 = Yo_o*Yo_Va*Yap*Ybr;  l74 = Yo_o*Yo_Va*Yap*Ybr*(Yo_o - Yo_Va);
    l75 = Yo_o*Yo_Va*Yap*Ybo;  l76 = Yo_o*Yo_Va*Yap*Ybo*(Yo_o - Yo_Va);
    l77 = Yo_o*Yo_Va*Yap*Ybpr; l78 = Yo_o*Yo_Va*Yap*Ybpr*(Yo_o - Yo_Va);
    l79 = Yo_o*Yo_Va*Yap*Ybpo; l80 = Yo_o*Yo_Va*Yap*Ybpo*(Yo_o - Yo_Va);

    ll = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14,...
        l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26,...
        l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38,...
        l39, l40, l41, l42, l43, l44, l45, l46, l47, l48, l49, l50,...
        l51, l52, l53, l54, l55, l56, l57, l58, l59, l60, l61, l62,...
        l63, l64, l65, l66, l67, l68, l69, l70, l71, l72, l73, l74,...
        l75, l76, l77, l78, l79, l80];
    
    syms L [1 80]
    % set subout and subin as symbolic arrays for place holder (preallocate)
    subout = L;
    subin = L;
    subzero = zeros(1,80);
    lastidx = 0;
    % loop through all combinations and find parallel terms - there's
    % probably a faster way to do this.
    for u = 1:length(ll)
        for t = 1:length(ll)
            if u ~= t && u < t
                if isSymType(sym(ll(u)/ll(t)),'integer | real') && ll(u) ~= 0 && ll(t) ~= 0 
                    subin(lastidx+1) = L(u);
                    subout(lastidx+1) = L(t);
                    lastidx = lastidx + 1;
                end
            end
        end
    end
  
    if lastidx > 0
        subin = subin(1:lastidx);
        subout = subout(1:lastidx);
        subzero = subzero(1:lastidx);
        for f = 1:length(subin)
            fprintf("The parameter %s is parallel to %s and set to 0\n",subout(f),subin(f))
        end    
    end  
    
    % Remake the solution model with the removed excess terms
    
    if lastidx > 0
        if xc==0 && yc~=0
            G_soln = subs(G_soln,subout,subzero);
            G_ex = subs(G_ex,subout,subzero);
        elseif yc==0 && xc~=0
            G_soln = subs(G_soln,subout,subzero);
            G_ex = subs(G_ex,subout,subzero);
        elseif xc==0 && yc==0
            G_soln = subs(G_soln,subout,subzero);
            G_ex = subs(G_ex,subout,subzero);
        end
    end
end