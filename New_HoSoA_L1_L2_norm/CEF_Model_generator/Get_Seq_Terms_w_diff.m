function [Gsol, Gdft,ex_terms, end_terms, ex_terms_og, end_terms_og, numLs]...
    = Get_Seq_Terms_w_diff(Gsol,G_soln,Gdft,Aprime,Bprime)
    %% determine number of endmember and excess terms

    syms x y z T real
    syms G [1 8] real
    syms gDiffA [1 4] real
    syms gDiffB [1 4] real
    syms gDiffC [1 4] real
    syms gDiffD [1 4] real
    syms gDiffE [1 4] real
    syms gDiffHo [1 4] real
    syms gDiffSo [1 4] real
    syms Ho [1 80] real
    syms So [1 80] real
    syms A [1 80] real
    syms B [1 80] real
    syms C [1 80] real
    syms D [1 80] real
    syms E [1 80] real
    syms L [1 80] real

    
    tot_params = symvar(Gsol);
    tot_Ls = symvar(G_soln);
    % remove indepedent varibales T x y z   
    tot_params(tot_params==T) = [];
    tot_params(tot_params==x) = [];
    tot_params(tot_params==y) = [];
    tot_params(tot_params==z) = [];
    tot_Ls(tot_Ls==T) = [];
    tot_Ls(tot_Ls==x) = [];
    tot_Ls(tot_Ls==y) = [];
    tot_Ls(tot_Ls==z) = [];

    ex_terms_og = setdiff(tot_params,[gDiffA gDiffB gDiffC gDiffD gDiffHo gDiffSo]);
    end_terms_og = setdiff(tot_params,[A B C D Ho So]);
    L_terms = setdiff(tot_Ls,G);
    numLs = length(L_terms);

    %% excess and possibly endmeber terms are non sequential, we can fix that

    %excess terms correction
    Gsol = subs(Gsol,intersect(ex_terms_og,A),A(1:numLs));
    Gsol = subs(Gsol,intersect(ex_terms_og,B),B(1:numLs));
    Gsol = subs(Gsol,intersect(ex_terms_og,C),C(1:numLs));
    Gsol = subs(Gsol,intersect(ex_terms_og,D),D(1:numLs));
    Gsol = subs(Gsol,intersect(ex_terms_og,Ho),Ho(1:numLs));
    Gsol = subs(Gsol,intersect(ex_terms_og,So),So(1:numLs));

    Gdft = subs(Gdft,intersect(ex_terms_og,Ho),Ho(1:numLs));

    %endmeber terms correction
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffA),gDiffA(1:(Aprime+Bprime)*2));
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffB),gDiffB(1:(Aprime+Bprime)*2));
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffC),gDiffC(1:(Aprime+Bprime)*2));
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffD),gDiffD(1:(Aprime+Bprime)*2));
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffHo),gDiffHo(1:(Aprime+Bprime)*2));
    Gsol = subs(Gsol,intersect(end_terms_og,gDiffSo),gDiffSo(1:(Aprime+Bprime)*2));

    Gdft = subs(Gdft,intersect(end_terms_og,gDiffHo),gDiffHo(1:(Aprime+Bprime)*2));

    %% Get updated list of excess and endmember terms for pass ins to other functions

    tot_params = symvar(Gsol);

    % remove indepedent varibales T x y z
    tot_params(tot_params==T) = [];
    tot_params(tot_params==x) = [];
    tot_params(tot_params==y) = [];
    tot_params(tot_params==z) = [];

    ex_terms = setdiff(tot_params,[gDiffA gDiffB gDiffC gDiffD gDiffHo gDiffSo]);
    end_terms = setdiff(tot_params,[A B C D Ho So]);

end