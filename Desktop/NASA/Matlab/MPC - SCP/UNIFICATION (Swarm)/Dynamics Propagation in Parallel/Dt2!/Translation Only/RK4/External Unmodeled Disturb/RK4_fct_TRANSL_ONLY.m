function [ xi , Neig_i, ui ] = RK4_fct_TRANSL_ONLY( x_nom , n_step_RT, sol , ii,K , x_len, u_len, B, E, fi_ORF, N, T, TH, Dt, skip, Mref, flim, Vmax, G, Rcol, Rcomm)
% OUTPUTS are now the updated state and neighborhood for all the time instants k and all robots i
Trem = T - n_step_RT; % Remaining time steps

Neig_i=zeros(N,Trem); %And these instead correspond to the Neigh from the actual instant (included) to the penultimate (T-1) -> since don't care about T in that I cannot propagate anymore

i = K(ii);

% Up to TH the problem comprises collision avoidance, later than TH I've
% only got the other constraints


xi=x_nom(:,:,i); % as first approx  % % % % % % % % % 


    cvx_begin quiet
                          
        variable xi(x_len , Trem+1) % Where the first value is the actual one
        variable ui(u_len , Trem)
        
        %% ALTERNATIVE (1) (Classic)
% %         minimize( sum( norms(ui , 1 , 1) )*Dt );
                
        %% ALTERNATIVE (2)
        if n_step_RT~=T-1
%         minimize( sum( norms(ui , 1 , 1) )*Dt );
% %         minimize( sum( norms(ui , 1 , 1) )*Dt + sum( norms( xi - repmat([Mref(:,i) ; [0 0 0]'],1,Trem+1) , 1 , 1) )*Dt );
% % % %             minimize( sum( norms( xi - repmat([Mref(:,i) ; [0 0 0]'],1,Trem+1) , 1 , 1) )*Dt );  % Goes very fast to convergence -> not good if I want to minimize also ctrl action
        minimize( sum( norms(ui , 1 , 1) )*Dt + sum( norms( xi(4:6,:) , 1 , 1) )*Dt + sum( norms( xi(1:3,:) - repmat(Mref(:,i),1,Trem+1) , 1 , 1) )*Dt );
% %         minimize( sum( norms(ui , 1 , 1) )*Dt + sum( norms( xi(:,:) , 1 , 1) )*Dt );
        else
        minimize( sum(   norms( xi(:,end) - [ Mref(:,i) ; [0 0 0]' ] , 1 , 1 )  )*Dt );
        end                                       %%% Which is the Final Condt Constr%%%
% %         and exclude Final Cdtn Cnstrt for the final state   
        


        %% Initial Condition
        xi(:,1) == x_nom(:,n_step_RT,i); % Initial Value is the one @ the actual time instant 
        %% Final Condition
        if n_step_RT~=T-1
        xi(:,end) == [ Mref(:,i) ; [0 0 0]' ];
        end
        
        
        if Trem >= TH  % n_step_RT <= T - TH
            
            for n_step = n_step_RT:1:n_step_RT+TH-1 % So that I have TH-1 elements  
                nn_step=n_step+1;
                
                %% Saturation Constraint
                norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;  % TH elements
                
                %% Velocity Constraint
                norm( xi(4:6,n_step-n_step_RT+2) , Inf ) <= Vmax;  %%%%% I changed the norm to INF here.  % TH+1 elements
                
                %% Dynamic Constraint
                x_OE = deval( sol , n_step*Dt );
                [ Ad , Bd , ~ ] = UNI_Swarm_Transl_StateMatrices( x_nom(:,n_step,i) , x_OE , B, Dt );
                xi(:,nn_step-n_step_RT+1) == Ad * xi(:,n_step-n_step_RT+1) + Bd * ui(:,n_step-n_step_RT+1) + E * fi_ORF(:,i) ;  % TH+1 elements
                
                % Collision Avoidance
                % Robot i Neighbors Subset (Based on Nominal traj of other robots)
                TF=norms( repmat(x_nom(1:3,n_step,i),1,N)-reshape(x_nom(1:3,n_step,:),3,N) , 2,1) <= Rcomm;     % Since if I am to use xi in order to define the neighbors of i I'll generate a constraint (because I'll have to use the cvx variable cx) then of course it is also theoretically correct to use the nominal traj also for i
                ct=0; Ni_k=[];
                for kk=1:length(TF)
                    if TF(kk)==1 && kk~=i
                        ct=ct+1;
                        Ni_k(ct)=kk;
                    end
                end
                
                Neig_i(1:length(Ni_k),n_step-n_step_RT+1) = Ni_k;
                
                % Collision Avoidance Constraint
                for j = Ni_k
                    ( x_nom(:,n_step,i)-x_nom(:,n_step,j) )' * G' * G * ( xi(:,n_step-n_step_RT+1)-x_nom(:,n_step,j) ) >= Rcol * norm( G*(x_nom(:,n_step,i)-x_nom(:,n_step,j)) , 2);
                end
                
            end
            
            if n_step_RT ~= T - TH
                
            for n_step = n_step_RT+TH:skip:T-1      %skip=Dt2/Dt1  % and meanwhile the other elements remain zero or old
                nn_step=n_step + skip; %%!!!                
                %% Saturation Constraint
                norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;  % Up to Trem
                %% Velocity Constraint
                norm( xi(4:6,n_step-n_step_RT+2) , Inf ) <= Vmax;
                %% Dynamic Constraint
                x_OE = deval( sol , n_step*Dt );
                if nn_step-n_step_RT+1<=Trem+1
                [ Ad , Bd , ~ ] = UNI_Swarm_Transl_StateMatrices( x_nom(:,n_step,i) , x_OE , B, Dt );
                xi(:,nn_step-n_step_RT+1) == Ad * xi(:,n_step-n_step_RT+1) + Bd * ui(:,n_step-n_step_RT+1) + E * fi_ORF(:,i) ;  % brings me from k0 to k1 -> so involves T1! (T_T involves k_T-1 to k_T)
                end
            end
            end
                
% %             for n_step = n_step_RT+TH:T-1 % So in total the 2 for loops must iterate Trem elements
% %                 nn_step=n_step+1;
% %                 
% %                 %% Saturation Constraint
% %                 norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;  % Up to Trem
% %                 
% %                 %% Velocity Constraint
% %                 norm( xi(4:6,n_step-n_step_RT+2) , Inf ) <= Vmax;
% %                 
% %                 %% Dynamic Constraint
% %                 x_OE = deval( sol , n_step*Dt );
% %                 [ Ad , Bd , ~ ] = UNI_Swarm_Transl_StateMatrices( x_nom(:,n_step,i) , x_OE , B, Dt );
% %                 xi(:,nn_step-n_step_RT+1) == Ad * xi(:,n_step-n_step_RT+1) + Bd * ui(:,n_step-n_step_RT+1) + E * fi_ORF(:,i) ;  % brings me from k0 to k1 -> so involves T1! (T_T involves k_T-1 to k_T)
% %             end
% %             end
        
        elseif n_step_RT > T - TH  &&  n_step_RT~=T-1
            
            for n_step = n_step_RT:1:T-1
                nn_step=n_step+1;
                
                %% Saturation Constraint
                norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;
                
                %% Velocity Constraint
                norm( xi(4:6,n_step-n_step_RT+2) , Inf ) <= Vmax;
                
                %% Dynamic Constraint
                x_OE = deval( sol , n_step*Dt );
                [ Ad , Bd , ~ ] = UNI_Swarm_Transl_StateMatrices( x_nom(:,n_step,i) , x_OE , B, Dt );
                xi(:,nn_step-n_step_RT+1) == Ad * xi(:,n_step-n_step_RT+1) + Bd * ui(:,n_step-n_step_RT+1) + E * fi_ORF(:,i) ;  % brings me from k0 to k1 -> so involves T1! (T_T involves k_T-1 to k_T)
                
                                 
                % Collision Avoidance
                % Robot i Neighbors Subset (Based on Nominal traj of other robots)
                TF=norms( repmat(x_nom(1:3,n_step,i),1,N)-reshape(x_nom(1:3,n_step,:),3,N) , 2,1) <= Rcomm;     % Since if I am to use xi in order to define the neighbors of i I'll generate a constraint (because I'll have to use the cvx variable cx) then of course it is also theoretically correct to use the nominal traj also for i
                ct=0; Ni_k=[];
                for kk=1:length(TF)
                    if TF(kk)==1 && kk~=i
                        ct=ct+1;
                        Ni_k(ct)=kk;
                    end
                end
                
                Neig_i(1:length(Ni_k),n_step-n_step_RT+1) = Ni_k;
                
                % Collision Avoidance Constraint
                for j = Ni_k
                    ( x_nom(:,n_step,i)-x_nom(:,n_step,j) )' * G' * G * ( xi(:,n_step-n_step_RT+1)-x_nom(:,n_step,j) ) >= Rcol * norm( G*(x_nom(:,n_step,i)-x_nom(:,n_step,j)) , 2);
                end
                                
            end
            
            
        elseif n_step_RT==T-1
            
            %% ALTERNATIVE (1) (Classic)
                % Last gap: the only thing I can find is the control action (can't
                % propagate state since the constraints are on initial and
                % final state, so it'd be overconstrained->NaN):
                n_step = T-1;
                %% Saturation Constraint
                    norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;

                %% No Collision Avoidance, since they are already collision free at T-1 (due to prev.step) and also at T by design.     

                %% So here the only constraints are on Initial and Final state, and control action.

            
            %% ALTERNATIVE (2)
                n_step = T-1;
                nn_step=n_step+1;
                %% Saturation Constraint
                norm( ui(:,n_step-n_step_RT+1) , Inf ) <= flim;
                %% Velocity Constraint
                norm( xi(4:6,n_step-n_step_RT+2) , Inf ) <= Vmax;

                %% Dynamic Constraint
                x_OE = deval( sol , n_step*Dt );
                [ Ad , Bd , ~ ] = UNI_Swarm_Transl_StateMatrices( x_nom(:,n_step,i) , x_OE , B, Dt );
                xi(:,nn_step-n_step_RT+1) == Ad * xi(:,n_step-n_step_RT+1) + Bd * ui(:,n_step-n_step_RT+1) + E * fi_ORF(:,i) ;  % brings me from k0 to k1 -> so involves T1! (T_T involves k_T-1 to k_T)


                % Collision Avoidance
                % Robot i Neighbors Subset (Based on Nominal traj of other robots)
                TF=norms( repmat(x_nom(1:3,n_step,i),1,N)-reshape(x_nom(1:3,n_step,:),3,N) , 2,1) <= Rcomm;     % Since if I am to use xi in order to define the neighbors of i I'll generate a constraint (because I'll have to use the cvx variable cx) then of course it is also theoretically correct to use the nominal traj also for i
                ct=0; Ni_k=[];
                for kk=1:length(TF)
                    if TF(kk)==1 && kk~=i
                        ct=ct+1;
                        Ni_k(ct)=kk;
                    end
                end

                Neig_i(1:length(Ni_k),n_step-n_step_RT+1) = Ni_k;

                % Collision Avoidance Constraint
                for j = Ni_k
                    ( x_nom(:,n_step,i)-x_nom(:,n_step,j) )' * G' * G * ( xi(:,n_step-n_step_RT+1)-x_nom(:,n_step,j) ) >= Rcol * norm( G*(x_nom(:,n_step,i)-x_nom(:,n_step,j)) , 2);
                end
        end
            
    cvx_end
    
xi = full(xi);
ui = full(ui);
end
% fix last value->it's the noise