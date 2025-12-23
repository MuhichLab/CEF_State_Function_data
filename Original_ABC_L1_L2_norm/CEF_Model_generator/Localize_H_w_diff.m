function [Gsol,Go,Gdft,Go_dft] = Localize_H_w_diff(Gsol,Go,Gdft,Go_dft,gs)

    syms gDiffHo [1 4] real
    syms G [1 8] real

    % These come from the DFT data - end point values per mol ABO3
    g1 = gs(1); 
    g2 = gs(2); 
    g3 = gs(3); 
    g4 = gs(4); 
    g5 = gs(5); 
    g6 = gs(6); 
    g7 = gs(7); 
    g8 = gs(7);

    % Full Solution Model
    Gsol = simplify(expand(subs(Gsol,[gDiffHo1 gDiffHo2 gDiffHo3 gDiffHo4 G5 G6 G7 G8],...
        [((g5 - g1) + gDiffHo1) ((g6 - g2) + gDiffHo2)...
        ((g7 - g3) + gDiffHo3) ((g8 - g4) + gDiffHo4) g5 g6 g7 g8])));

    % Go Solution Model
    Go = simplify(expand(subs(Go,[gDiffHo1 gDiffHo2 gDiffHo3 gDiffHo4 G5 G6 G7 G8],...
        [((g5 - g1) + gDiffHo1) ((g6 - g2) + gDiffHo2)...
        ((g7 - g3) + gDiffHo3) ((g8 - g4) + gDiffHo4) g5 g6 g7 g8])));

    % Gdft Solution model  (No T terms)
    Gdft = simplify(expand(subs(Gdft,[gDiffHo1 gDiffHo2 gDiffHo3 gDiffHo4 G5 G6 G7 G8],...
        [((g5 - g1) + gDiffHo1) ((g6 - g2) + gDiffHo2)...
        ((g7 - g3) + gDiffHo3) ((g8 - g4) + gDiffHo4) g5 g6 g7 g8])));
    
    % Go_dft Solution model  (No T terms)
    Go_dft = simplify(expand(subs(Go_dft,[gDiffHo1 gDiffHo2 gDiffHo3 gDiffHo4 G5 G6 G7 G8],...
        [((g5 - g1) + gDiffHo1) ((g6 - g2) + gDiffHo2)...
        ((g7 - g3) + gDiffHo3) ((g8 - g4) + gDiffHo4) g5 g6 g7 g8])));

end