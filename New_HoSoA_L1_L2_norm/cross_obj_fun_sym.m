function error = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
    x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,muhg_o,Temp,Gdft,dG_dy,d_chem,...
    lambda,T_ref,dp_terms,num_ex_params,self_solve,ex_terms,end_terms,lambda_L1,lambda_L2,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)

    % Don't forget muhg_o was passed in as -muhg_o
    
    if length(guess_excess) < num_ex_params
        og_guess = guess_excess;
        guess_excess = zeros(1,length(dp_terms));
        n = 1;
        for k = (1:length(dp_terms))
            if dp_terms(k) == 1
                guess_excess(k) = og_guess(n);
                n = n + 1;
            end
        end
    end
    

    dGdy = dG_dy;
    G_dft = Gdft;
    dG_dy = subs(dG_dy,ex_terms,guess_excess);
    dG_dy = subs(dG_dy,end_terms,guess_end);

    
    Gdft = subs(Gdft,intersect(ex_terms,Ho),guess_excess(ismember(ex_terms,Ho)));
    Gdft = subs(Gdft,intersect(end_terms,gDiffHo),guess_end(ismember(end_terms,gDiffHo)));
    
    dG_dy = matlabFunction(dG_dy);
    Gdft = matlabFunction(Gdft);

    %% EXP Error
%     % For Timing    
    if T_ref ~= 0
        if self_solve
            d0 = solve_dref(dG_dy,x_vals,d_chem,T_ref);
        else
            d0 = X0;
        end
        dd_adj = dd;
        for i = 1:length(dd)
            for k = 1:length(x_vals)
                if x_exp(i) == x_vals(k)
                    dref = d0(k);
                    break
                end
            end
            dd_adj(i) = dd_adj(i)+dref;
        end
    else
        dd_adj = dd;
    end
        
    % chem_pot error (Root Mean Sqaure Error)


        if any(symvar(dGdy)=='x') && any(symvar(dGdy)=='y')
            predict = double(dG_dy(Temp,x_exp,y_exp,dd_adj));
        elseif any(symvar(dGdy)=='x') && all(symvar(dGdy)~='y')
            predict = double(dG_dy(Temp,x_exp,dd_adj));
        elseif all(symvar(dGdy)~='x') && any(symvar(dGdy)=='y')
            predict = double(dG_dy(Temp,y_exp,dd_adj));
        elseif all(symvar(dGdy)~='x') && all(symvar(dGdy)~='y')
            predict = double(dG_dy(Temp,dd_adj));
        end

        errorEXP = sum((predict-muhg_o).^2);%./length(muhg_o));

            %% DFT error wit red ev error

%         Tot  eV Error (Root Mean Sqaure Error)
    if any(symvar(G_dft)=='x') && any(symvar(G_dft)=='y')
        predict = double(Gdft(X,Y,Z));
    elseif any(symvar(G_dft)=='x') && all(symvar(G_dft)~='y')
        predict = double(Gdft(X,Z));
    elseif all(symvar(G_dft)~='x') && any(symvar(G_dft)=='y')
        predict = double(Gdft(Y,Z));
    elseif all(symvar(G_dft)~='x') && all(symvar(G_dft)~='y')
        predict = double(Gdft(Z));
    end

    rmse_tot_ev = sum((predict-eV).^2);%./length(Y));

    % Reduciton Energy Error (Root Mean Sqaure Error)
    % red_ev output is [X ,Y, reduciton eV]
    red_ev = get_reduction_energy(X,Z,eV);

    predict2 = double(Gdft(red_ev(:,1),red_ev(:,2)) - Gdft(red_ev(:,1),0));
    rmse_red_ev = sum((predict2-red_ev(:,3)).^2);%./length(red_ev(:,2)));

    errorDFT = rmse_red_ev + rmse_tot_ev;

    % Total error to minimize

    error = lambda(2)*errorEXP + lambda(1)*errorDFT + ...
        lambda_L1*sum(abs(guess_excess)) + lambda_L2*sum(guess_excess.^2) + ...
        lambda_L1*sum(abs(guess_end)) + lambda_L2*sum(guess_end.^2);


end

