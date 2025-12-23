function red_ev = get_reduction_energy(X,Z,eV)

    red_ev = [];
    xvals = unique(X);
        for i = 1:length(xvals)
            deltas = Z(X==xvals(i));
            deltas = sort(deltas);
            [J,~] = size(deltas);
            if J > 1 && deltas(1) == 0
                z = eV(X==xvals(i));
                for j = 1:J
                    redev = z(j) - z(1);
                    red_ev = [red_ev ; xvals(i) deltas(j) redev];
                end
                
            end
            
        end
    
end