function [dp_terms,spot] = Significance_and_Drop_Logic...
    (guess_end,guess_excess,dp_terms,num_terms,T_ref,d_chem,X,Y,Z,eV,x_exp,y_exp,dd,...
    muhg_o,Temp,Gdft,dG_dy,x_vals,y_vals,X0,Y0,spot_track,exp_Weights,dft_Weights,...
    lambda,num_ex_params,ex_terms,end_terms,numLs,original_error,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)

    % Finite Diff
%     [drop_this, combos1, combos2, combos3] =...
%         drop_me_finite_diff_4th_order_para(guess_end,guess_excess,X,Z,eV,x_exp,dd,...
%         muhg_o,T_ref,d_chem,Temp,Gdft,dG_dy,x_vals,original_error,d0,spot_track,...
%         exp_Weights,dft_Weights,lambda,num_ex_params,ex_terms,end_terms,numLs,...
%         x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D);
    
    [drop_this, combos1, combos2, combos3] =...
        drop_me_finite_diff_2nd_order_para(guess_end,guess_excess,X,Y,Z,eV,x_exp,y_exp,dd,...
    muhg_o,T_ref,d_chem,Temp,Gdft,dG_dy,x_vals,y_vals,original_error,X0,Y0,spot_track,...
    exp_Weights,dft_Weights,lambda,num_ex_params,ex_terms,end_terms,numLs,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D);
    
    % Just percent change
%     [drop_this, combos1, combos2, combos3] = drop_me_2(guess_drop,...
%         X,Y,Z,x_exp,dd,-muhg_o,Temp,G_dft,dG_dy,x_vals,original_error,...
%         d0,spot_track,exp_Weights,dft_Weights,lambda);

    drop_por1 = drop_this(1:length(combos1));
    drop_por2 = drop_this(length(combos1)+1:length(combos2)+length(combos1));
    drop_por3 = drop_this(length(combos1)+length(combos2)+1:length(combos1)+length(combos2)+length(combos3));
    
    fprintf('The min of drop_this is: %0.2f\n',min(drop_this));
    fprintf('The max of drop_this is: %0.2f\n',max(drop_this));
 
    % Drop Logic
    if num_terms > 3
    
        % Find the least significant terms
        min_val1 = min(drop_por1(drop_por1~=0));
        min_val2 = min(drop_por2(drop_por2~=0));
        min_val3 = min(drop_por3(drop_por3~=0));

        spot1 = find(drop_por1 == min_val1);
        spot2 = find(drop_por2 == min_val2);
        spot3 = find(drop_por3 == min_val3);

        combos1(spot1);
        combos2(spot2,:);
        combos3(spot3,:);

%         fprintf('Combo 1 are these terms and this Value: %d       : %.3e \n',combos1(spot1),min_val1);
%         fprintf('Combo 2 are these terms and this Value: %d %d   : %.3e \n',combos2(spot2,:),min_val2);
%         fprintf('Combo 3 are these terms and this Value: %d %d %d : %.3e \n',combos3(spot3,:),min_val3);

       % Logic to compare htese terms against eachother
       if (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 < min_val3
            new_min = min(drop_por1(combos2(spot2,:)));
            spot = find(drop_por1 == new_min);
            for f = 0:num_ex_params/numLs-1
                dp_terms(spot+numLs*f) = 0;
            end
 
       elseif (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 > min_val3
            spot = setdiff(combos3(spot3,:),combos2(spot2,:));
            for f = 0:num_ex_params/numLs-1
                dp_terms(spot+numLs*f) = 0;
            end
  
       else

            if min_val2 < min_val3
                new_min = min(drop_por1(combos2(spot2,:)));
                spot = find(drop_por1 == new_min);
                for f = 0:num_ex_params/numLs-1
                    dp_terms(spot+numLs*f) = 0;
                end

            elseif min_val3 < min_val2
                new_min = min(drop_por1(combos3(spot3,:)));
                spot = find(drop_por1 == new_min);
                for f = 0:num_ex_params/numLs-1
                    dp_terms(spot+numLs*f) = 0;
                end
            end
       end

    elseif num_terms <= 3
        min_val1 = min(drop_por1(drop_por1>0));
        spot1 = find(drop_por1 == min_val1);
        spot = spot1;
        for f = 0:num_ex_params/numLs-1
            dp_terms(spot+numLs*f) = 0;
        end

    end
end