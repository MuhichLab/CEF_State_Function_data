function [fake_excess_terms] = fake_excess(ex_terms)

    syms Ho [1 80] real
    syms So [1 80] real
    syms A [1 80] real
    syms B [1 80] real
    syms C [1 80] real
    syms D [1 80] real
    syms E [1 80] real

    fake_excess_terms = zeros(1,length(ex_terms));
    
    A0s_ex = ismember(ex_terms,A);
    H0s_ex = ismember(ex_terms,Ho);
    S0s_ex = ismember(ex_terms,So);
    
    pda = makedist('Lognormal','mu',0,'sigma',1e-3);
    a0 = random(pda,sum(A0s_ex),1);
    loga0 = log(a0);
    
    pds0 = makedist('Lognormal','mu',0,'sigma',1e-3);
    s0 = random(pds0,sum(S0s_ex),1);
    logs0 = log(s0);    

    pdh0 = makedist('Lognormal','mu',0,'sigma',1);
    h0 = random(pdh0,sum(H0s_ex),1);
    logh0 = log(h0);     

    fake_excess_terms(A0s_ex) = loga0;
    fake_excess_terms(H0s_ex) = logh0;
    fake_excess_terms(S0s_ex) = logs0;
    
end