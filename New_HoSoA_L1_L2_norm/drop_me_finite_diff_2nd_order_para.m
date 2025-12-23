function [drop, combos1, combos2, combos3] =...
    drop_me_finite_diff_2nd_order_para(guess_end,guess_excess,X,Y,Z,eV,x_exp,y_exp,dd,...
    muhg_o,T_ref,d_chem,Temp,Gdft,dG_dy,x_vals,y_vals,original_error,X0,Y0,spot_track,...
    exp_Weights,dft_Weights,lambda,num_ex_params,ex_terms,end_terms,numLs,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)
    
    orig = guess_excess;
%**************************************************************************
    
    combos1 = nchoosek(1:numLs,1);
    combos2 = nchoosek(1:numLs,2);
    combos3 = nchoosek(1:numLs,3);
   
    % h is step size of parameter adjust (0.01 = 1% change in paramters)
    h = 0.005;
    
    
 %% %%%%%%%%%%%%%%%%%%%%% Look at Curvature of Ls %%%%%%%%%%%%%%%%%

%     dh = linspace(-h*3,h*3,20);
%     steps = zeros(length(combos1),length(dh));
%     test = zeros(length(combos1),1);
%     for p = 1:length(dh)
%         testh = dh(p);
%         s1 = tic();
%         for  k = 1:length(combos1)
% %             if neg1(k,1) == 0
% %                 continue
% %             end
%             guess_excess = orig;
%             for f = 0:num_ex_params/numLs-1
%                 guess_excess(combos1(k,:)+(f*numLs))=guess_excess(combos1(k,:)+(f*numLs))*(1+testh);
%             end
%             test(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Z,eV,...
%             x_exp,dd,d0,x_vals,-muhg_o, Temp,Gdft,dG_dy,d_chem,exp_Weights,...
%             dft_Weights,lambda,T_ref,ones(num_ex_params,1),num_ex_params,...
%             false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,...
%             gDiffHo,gDiffSo,Ho,So,A,B,C,D);
%         end 
% %         toc(s1)
%         steps(:,p) = test;
%     end
%     
%     steps(~any(steps,2),:) = [];
%     figure
%     plot(dh,steps)
%     legend
%     ylim([original_error-original_error*0.01 original_error+original_error*0.01])
% %     ylim([.6 .65])
%     title("This is run %d",sum(spot_track==0))    
%     
    %% pre fill plus and neg with things we don't need to calcualte
    
    plus1 = ones(length(combos1)+length(combos2)+length(combos3),1);
    count = 1;
    for k = 1:length(combos1)
        if ismember(combos1(k),spot_track)
            plus1(count,1) = 0;
            count = count + 1;
            continue
        end
        count = count + 1;
    end  
    for k = 1:length(combos2)
        if ismember(combos2(k,1),spot_track) || ismember(combos2(k,2),spot_track)
            plus1(count,1) = 0;
            count = count + 1;
            continue
        end
        count = count + 1;
    end
    for k = 1:length(combos3)
        if ismember(combos3(k,1),spot_track) || ismember(combos3(k,2),spot_track) || ismember(combos3(k,3),spot_track)
            plus1(count,1) = 0;
            count = count + 1;
            continue
        end 
        count = count + 1;
    end
    
    count = 1;
    neg1 = ones(length(combos1)+length(combos2)+length(combos3),1);
    for k = 1:length(combos1)
        if ismember(combos1(k),spot_track)
            neg1(count,1) = 0;
            count = count + 1;
            continue
        end
        count = count + 1;
    end
    for k = 1:length(combos2)
        if ismember(combos2(k,1),spot_track) || ismember(combos2(k,2),spot_track)
            neg1(count,1) = 0;
            count = count + 1;
            continue
        end
        count = count + 1;
    end   
    for k = 1:length(combos3)
        if ismember(combos3(k,1),spot_track) || ismember(combos3(k,2),spot_track) || ismember(combos3(k,3),spot_track)
            neg1(count,1) = 0;
            count = count + 1;
            continue
        end
        count = count + 1;
    end

    %% get y+1 saved as plus

    parfor k = 1:length(combos1)
        
        if plus1(k,1) == 0
            continue
        end
        
        guess_excess = orig;
        
        for f = 0:num_ex_params/numLs-1
            guess_excess(combos1(k,:)+(f*numLs))=guess_excess(combos1(k,:)+(f*numLs))*(1+h);
        end

        plus1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D);
        
    end
    
    parfor k = length(combos1)+1:length(combos2)+length(combos1)
        
        kk = k - length(combos1);
        if plus1(k,1) == 0
            continue
        end
        
        guess_excess = orig;
        
        for f = 0:num_ex_params/numLs-1
            guess_excess(combos2(kk,:)+(f*numLs))=guess_excess(combos2(kk,:)+(f*numLs))*(1+h);
        end
        
        plus1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)
    end

    parfor k = length(combos1)+length(combos2)+1:length(combos3)+length(combos2)+length(combos1)
        
        kk = k - length(combos2) - length(combos1);
        if plus1(k,1) == 0
            continue
        end
        guess_excess = orig;

        for f = 0:num_ex_params/numLs-1
            guess_excess(combos3(kk,:)+(f*numLs))=guess_excess(combos3(kk,:)+(f*numLs))*(1+h);
        end

        plus1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)
    end
    
    %% get y-1 saved as neg

    parfor k = 1:length(combos1)
        if neg1(k,1) == 0
            continue
        end
        guess_excess = orig;
        
        for f = 0:num_ex_params/numLs-1
            guess_excess(combos1(k,:)+(f*numLs))=guess_excess(combos1(k,:)+(f*numLs))*(1-h);
        end
 
        neg1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D) 
    end
    
    parfor k = length(combos1)+1:length(combos2)+length(combos1)
        
        kk = k - length(combos1)
        if neg1(k,1) == 0
            continue
        end
        
        guess_excess = orig;
        
        for f = 0:num_ex_params/numLs-1
            guess_excess(combos2(kk,:)+(f*numLs))=guess_excess(combos2(kk,:)+(f*numLs))*(1-h);
        end
        

        neg1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)
    end

    parfor  k = length(combos1)+length(combos2)+1:length(combos3)+length(combos2)+length(combos1)
        
        kk = k - length(combos2) - length(combos1);
        if neg1(k,1) == 0
            continue
        end
        guess_excess = orig;

        for f = 0:num_ex_params/numLs-1
            guess_excess(combos3(kk,:)+(f*numLs))=guess_excess(combos3(kk,:)+(f*numLs))*(1-h);
        end

 
        neg1(k,1) = cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
            x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,...
            d_chem,exp_Weights,dft_Weights,lambda,T_ref,ones(num_ex_params,1),...
            num_ex_params,false,ex_terms,end_terms,x,y,z,T,gDiffA,gDiffB,...
            gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)
    end    
    

    %% Calcualte second derivative of error
    
    guess_excess = orig;
    drop = zeros(length(combos1)+length(combos2)+length(combos3),1);
    wng = 0;
    for s =1:length(drop)
        if plus1(s)~= 0 && neg1(s)~=0
            % 4th order finit diff of second derivative
            drop(s) = (plus1(s) - 2*original_error +...
                neg1(s))/(h^2);
        elseif plus1(s) == 0 && neg1(s)==0
            continue
        else 
            wng = wng+1;
        end
    end
    %fprintf("Somethign went wrong with drop_me_finite this many times %d \n",wng)
end