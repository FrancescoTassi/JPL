function [ dxdt ] = odefcn_RegSys_I_ExtDist ( t, x, x_Ref, Dt,tf, m, fext_ORF, Kp, Kd, flim, sol, xx0_i,i ,endColumn,timeline )
% global uu_save uu_discrete 
dxdt=zeros(length(x),1);
muE=3.986005e14; %m^3/s^2
kJ2=2.633e10*10^6; %m^2/s^2

% x(4:7) is now the quaternion that rotates from I axes to instantaneous body axes! 
% BUT I am assuming that at each time step I can measure this value!
%
% So I defined a control based on Inertial-Body rotations, without looking
% at what happens in the middle (orbit etc)
if t>=tf
    n_step_actual = tf/Dt;
else
    n_step_actual = ceil(t/Dt);
end

% % rho = x(1:3);     rho_Ref = x_Ref( 1:3 , n_step_actual );
% % rho_d = x(8:10);  rho_d_Ref = x_Ref( 4:6 , n_step_actual );
x_Ref_t=zeros(size(x_Ref,1),1);
for k=1:size(x_Ref,1)
    x_Ref_t(k) = interp1([0 (1:tf/Dt)*Dt] , [xx0_i(k) x_Ref(k,:)], t);
end
% rho = x(1:3 , count);     rho_Ref = x_Ref_t( 1:3 );
% rho_d = x(4:6 , count);   rho_d_Ref = x_Ref_t( 4:6 );

%% If Measurement Noises:
rho = x(1:3) ;      rho_Ref = x_Ref_t( 1:3 );
rho_d = x(4:6) ;    rho_d_Ref = x_Ref_t( 4:6 );

% % it_indx = find(t==timeline);
% rho_Meas = x(1:3) + meas_noise(1:3 , it_indx) ;       % THE NOISY ONE OF THE MEASURE IS ONLY USED FOR THE CONTROLLER 
% rho_d_Meas = x(4:6) + meas_noise(4:6 , it_indx) ;     % THE INTEGRATOR INSTEAD REPRESENTS THE REAL SYS-> INDEP FROM MEAS.

rho_Meas = x(1:3) -1+2*1*rand(3,1);       % THE NOISY ONE OF THE MEASURE IS ONLY USED FOR THE CONTROLLER 
rho_d_Meas = x(4:6) -.01+2*.01*rand(3,1);

%% - REGULATORs

%    Position Controller:
u = Kp * (rho_Ref - rho_Meas) + Kd * (rho_d_Ref - rho_d_Meas);  %% ACTUALLY IN REALITY I AM FEEDING TO THE SYS NOT THE rho BUT THE MEASURED VALUE, WHICH IS THE SAME FOR THE WHOLE INSTANT DT, provided that SamplingTime=Dt (like this I am saying that Ts<Dt, which is better because I have to separate the Dt_MPC from Dt_sampling, much faster)
% u = Kp * ([rho_Ref; rho_d_Ref]-[rho; rho_d]) + Kd * ([]-[dxdt(1:3); dxdt(8:10)]) ;

%    Attitude Controller:
% % e = q_Ref(2:4, n_step_actual) - q(2:4);  %%%%%%%%%%%  BUT q must be the instant quat FROM ZERO ROTATION(F_I) TO THE ACTUAL F_b
% if q(1) < 0
%     e(2:4) = - e(2:4); % Conjugate
% end
% % tau_c =   Kq*e;                          

%    Angular Speed Controller:
% % omlim=zeros(3,1);
% % for iq=1:3
% %     if abs(vf_I_i(iq)-v_I(iq))<=ang_th
% %         omlim(iq)=interp1([0 ang_th], [0 omM], abs(vf_I_i(iq)-v_I(iq)));
% %         if t<tf/4 % && abs(vf(i,iq)-v(i,iq))>ang_th   So that I can start from zero speed
% %             omlim(iq)=interp1([0 t_th], [0 omlim(iq)], t);
% %         end
% %     end
% % end
% % % for iq=1:3
% % %     omlim(iq)=interp1([0 tf/4 tf*3/4 tf], [0 omM omM 0], t); %trapezoidal speed profile (if they all arrive at final time)
% % % end

% % tau_c = tau_c + Cq * [ ( sign(vf_I_i(1)-v_I(1))*omlim(1) - x(11) ) ;
% %                       ( sign(vf_I_i(2)-v_I(2))*omlim(2) - x(12) ) ;
% %                       ( sign(vf_I_i(3)-v_I(3))*omlim(3) - x(13) )] ;
% omlim_save(:,i)=[sign(vf_I(i,1)-v_I(i,1))*omlim(1) sign(vf_I(i,2)-v_I(i,2))*omlim(2) sign(vf_I(i,3)-v_I(i,3))*omlim(3)];

%  %    Angular Acceleration Controller:
%  for iq=1:3  % not good since the velocity prof is not trapez ->try making it such. But like this I am only controlling at beginning and end (always since t_th=50, and we're always trying to accelerate probably)
%      if t(2)<t_th
%          omdlim(iq) = sign(xx(i,10+iq))*om_dM;
%      elseif abs(vf_I(i,iq)-v_I(i,iq))<=.5
%          omdlim(iq) = - sign(xx(i,10+iq))*om_dM;
%          %                             omdlim(iq)=interp1([0 ang_th], [0 om_dM], abs(vf_I(i,iq)-v_I(i,iq)));
%          %                              if t(2)<t_th % && abs(vf(i,iq)-v(i,iq))>ang_th   So that I can start from zero speed
%          %                                 omdlim(iq)=interp1([0 t_th], [0 omdlim(iq)], t(2));
%          %                              end
%      else
%          omdlim(iq) = 0;
%      end
%  end
%  tau_c= J_b * omdlim ;
%  taux_c= taux_c + tau_c(1) ;
%  tauy_c= tauy_c + tau_c(2) ;
%  tauz_c= tauz_c + tau_c(3) ;
%  omdlim_save(:,i)=omdlim;
                        
%% - Saturation 
        for ch=1:3
            if u(ch)>flim
                u(ch)=flim;
            elseif u(ch)<-flim
                u(ch)=-flim;
            end
% %             if tau_c(ch)>Tlim
% %                 tau_c(ch)=Tlim;
% %             elseif tau_c(ch)<-Tlim
% %                 tau_c(ch)=-Tlim;
% %             end
        end
 
% %         Tau_c_I=tau_c;
        
%% - Actuators Noises
% % u = u + u_noise(:,it_indx); %-.001*u+2*.001*u.*rand(3,1);    % GIVES PROBLEM TO THE ODE -> DEFINE A RND NOISE OUTSIDE ODE, FOR EACH TIME INST, AND PASS IT IN
u = u -1e-1+2*1e-1*rand(3,1);

% Tau_c_I = Tau_c_I -.01*Tau_c_I+2*.01*Tau_c_I.*rand(3,1); 



%% - REAL SYSTEM
    
    
% % Tau_c_b = quatrotate( x(4:7)' , Tau_c_I' )';

%% State Propagation  
% Angular Velocity in I
% % om_I=[x(11) x(12) x(13)]';
% % om_skw_I=[0 -x(13) x(12);
% %           x(13) 0 -x(11);
% %           -x(12) x(11) 0];
      
%% Rotate to Body Frame through the INSTANTANEOUS quaternion qI_b:
% % om_b = quatrotate( x(4:7)' , om_I' )';
% % om_skw_b=[0 -om_b(3) om_b(2);
% %           om_b(3) 0 -om_b(1);
% %           -om_b(2) om_b(1) 0];

% %% Using Axis Angle:
% theta_axang=2*atan2(norm(x(5:7)),x(4));
% Ax_axang=zeros(1,3);
% for ii=1:3
%     if theta_axang~=0
%         Ax_axang(ii)=x(4+ii)/sin(theta_axang/2);
%     end
% end
% % A_ORF_b=[theta_axang*Ax_axang(1) theta_axang*Ax_axang(2) theta_axang*Ax_axang(3)]';
% % om_b=A_ORF_b*norm(om_ORF);
% axang_ORF_b=[Ax_axang theta_axang];
% R_ORF_b = vrrotvec2mat(axang_ORF_b);
% 
% om_b=R_ORF_b'*om_ORF;
% om_skw_b=[0 -om_b(3) om_b(2);
%           om_b(3) 0 -om_b(1);
%           -om_b(2) om_b(1) 0];

% % % Rotational Equation is integrated in 	ORF Frame
% % % omd_ORF = inv(J_b) * (Tau_b - om_skw_ORF*J_b*om_ORF  - om_skw_ORF*h_b + Tau_cc_b);     % miscellaneous-ORF (so neglecting the fact that the inertia will vary)
% % omd_ORF = inv(J_b) * (Tau_b - om_skw_ORF*J_b*om_ORF  - om_skw_ORF*h_ORF + Tau_cc_ORF);   % ORF (so neglecting the fact that the inertia will vary)
% % dxdt(11:13) = omd_ORF;

% Rotational Equation is integrated in BODY Frame
% so that the inertia tensor is constant.
% % omd_b = inv(J_b) * (Tau_b - om_skw_b*J_b*om_b + Tau_c_b);               % b
% % omd_I = quatrotate( quatinv(x(4:7)'), omd_b' )';                        % I
% % dxdt(11:13) = omd_I;

% %In axis angle notation in order to obtain the reverse rotation I simply
% %swap the rotation sign (?)
% omd_b = inv(J_b) * (Tau_b - om_skw_b*J_b*om_b  - om_skw_b*h_b + Tau_cc_b);       % BODY
% % A_b_ORF=-A_ORF_b;
% % omd_ORF=A_b_ORF*norm(omd_b);
% omd_ORF=R_ORF_b*omd_b;
% dxdt(11:13) = omd_ORF;

% Alternative Quaternion Notation
% % g = [x(5) x(6) x(7)]'; q1 = [ x(4) ];                   % Quaternion from I to body (so in I)
% % dxdt(4) = - 1/2 * om_I' * g;
% % dxdt(5:7) = -1/2 * om_skw_I * g + 1/2 * q1 * om_I;  % Quaternion derivative from I to body (so in I)

%% Translational Equation
% "Gravitational" Stiffness Term - K_orb
x_OE = deval( sol , t );
r=x_OE(1);  % R0
h=x_OE(3);
i1=x_OE(4);
theta=x_OE(6);

R0_ORF = [-r 0 0];
R0_ORF_skw=[0 -R0_ORF(3) R0_ORF(2);
          R0_ORF(3) 0 -R0_ORF(1);
          -R0_ORF(2) R0_ORF(1) 0];
      
OM1 = -kJ2*sin(2*i1)*sin(theta)/h/r^3;    % ORF
OM2 = 0;
OM3 = h/r^2;
OM_skw_ORF = [0 -OM3 OM2; OM3 0 -OM1; -OM2 OM1 0];

Korb = OM_skw_ORF*(OM_skw_ORF*eye(3)) - muE * (2*norm(R0_ORF)^2*eye(3) + 3*R0_ORF_skw*(R0_ORF_skw*eye(3)))/norm(R0_ORF)^5;

f_noises = zeros(3,1) -1+2*1*rand(3,1); % 
% UNMODELED IMPULSE
if t>=500 && t<=505
    f_noises = f_noises + 1*ones(3,1);
end

% rho_i , rhod_i
dxdt(1:3) = rho_d;                                                                                                                              
dxdt(4:6) = ( 1/m * ( u + fext_ORF + f_noises - 2*m*OM_skw_ORF*rho_d - m*Korb*rho ) );   % ORF   + UNMODELED DISTURBANCES


% Save
% uu_discrete(:,n_step_actual,i)=u; %  If I want to compare it only in the Dt pts with the result of the optimization
% [~,i_time]=min(abs(timeline-t));
% uu_save(:,endColumn-1+i_time,i)=u; 

global uu_i
[~,i_time]=min(abs(timeline-t));
uu_i(:,i_time)=u;
end

