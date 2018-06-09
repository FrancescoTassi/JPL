%% MPC - SCP of ROBOTIC SWARM (DETERMINISTIC APPROACH) - Translation + Attitude - PARFOR DYN PROPAGATION
close all, clear all, clc

%% SINCE OFFLINE GOES, MEANS THE PROBL IS WHEN TREM REDUCES, SO NOT ENOUGH TIME TO REARRANGE ->either too high forces and velocities of the bots (ctrl too high), or not sufficient force (ctr too low). But changing only T doesn't matter since Trem changes
%at least the offline goes

%% System Parameters
x_len=2*3;
u_len=1*3;

r=5e-2; 
m=1;    % 1e-3 Kg

% Parameters:
% Central Landmark Robot
H=[10; 10; 10];
h=length(H(1,:));
nn=1; 
rr=10;
[X, Y, Z, Mref, l, N] = l_gen_HEX ( nn, rr, H) 

workSpace = 500;
X0=workSpace*rand(1,N);
Y0=workSpace*rand(1,N);
Z0=workSpace*rand(1,N);
for i=1:N
    if randi([0 1])
    X0(i)=-X0(i);
    end
    if randi([0 1])
    Y0(i)=-Y0(i);
    end
    if randi([0 1])
    Z0(i)=-Z0(i);
    end
end
M=[X0; Y0; Z0]; Minit=M;

%% Dynamic Parameters
Umax = 1; %1;  %.001  tf300 dt30  %.01  tf500 dt50

Vmax = 1;     % 5;  or increase mass    Have to go slower and less force   or try pump everyth up and it takes what it needs

%% Orbital Parameters
% Mean Rotation of the Orbit
% Supposing GEOSYNCHRONOUS and KEPLERIAN orbit (so that v=v_tangt):
%   OM= [0 0 sqrt(muE/rs^3)]'; %cross(x(20:22),x(23:25))/norm(x(20:22))^2;
rs=42164e3; %m  
muE=3.986005e14; %m^3/s^2
gz=muE/rs^2;

% fi_ORF = repmat(-.01+2*.01*rand(3,1),1,N);                %%%%%%%%%%%%%%!
% fi_ORF = repmat(-1+2*1*rand(3,1),1,N);
fi_ORF = repmat([ -m*gz 0 0]',1,N);    

%% State Initial Conditions
rho_d0 = zeros(3,N);
x0 = [M ; rho_d0];


%% Time Parameters
%  Continuous
t0=0;
tf = 800;   %300 By increasing this time they could be forced to take longer routes
Dt = 5;  %30        %Also try reducing Dt ->increase T

% MPC Time Horizon (It's a number of steps to be done at each optimization, not a continuous time) 
TH = 60;  %2  % Now I am solving very few steps (TH steps with coll avoid and T-TH without)-> faster
tH=TH*Dt;

% NO tf = TH*Dt1 + (T-TH)*Dt2 -> tf = T * Dt1
% T=TH+1/Dt2*(tf-TH*Dt1);
T = tf/Dt;     %   Dt1 matters, since when I reach k0>T-TH, I only use Dt1
if TH>T
    disp('Error: TH > T')
    return
end

skip = 1;

% global uu_save 
global uu_discrete  
global uu_i
uu_save=zeros(3,1,N);
uu_discrete=zeros(3,T,N);

%% Orbital Elements Integration (LVLH Origin)
vx0=0; % Radial Velocity -> keplerian orbit
h0=1;
i0=0; % Geosynchronous
theta0=pi/2; % so LVLH origin is on Y axis
x_OE_0 = [rs vx0 h0 i0 0 theta0];
% x_OE_0 = [6878 vx0 pi/4 pi/3 0 0];

sol = ode45 ( @(t,x) odefcn_OE( t, x ), [t0 tf+Dt], x_OE_0 ); % CONTINUOUS TIME!!!!!!! [0:Dt:Dt*(T-1)] 
%% State Matrices
B = [ zeros(3)
      1/m*eye(3) ]; 
E = [ zeros(3)
      1/m*eye(3) ];  
  
%% Control Gain Matrices
Kp = 1e0*eye(3);  %1e1 better position matching - 1e-1 better ctr effort matching    -> so use 1e0         % To improve the control I have to increase these, but without generating instability
Kd = 1e0*eye(3);  % REMEMBER THAT u = Kp x_e + Kd v_e(which is smaller, so need bigger Kd to counterwight)
%% When the error is {affine}=={blabla} and happens in the following oprimizations
%% could be that the control requires too much vel or force, and the optimiz becomes unfeasible
%% -> reduce Vmax/Umax but also the gains!

%% Collision Avoidance Parameters
Rcol = r;%*3;
G = [eye(3) zeros(3)];
Rcomm = 2; % Communication Radius (Dictates the Neighborhood)  

%% Non Linear System Initial Conditions
xx_save = [];
for i=1:N                                                      
    xx_save(:,1,i) = [M(:,i)' rho_d0(:,i)' ];
end

state = size(xx_save,1);

xx0 = xx_save;
%%      %%%%%%%%%%%%%%%%%%%%%%% OFFLINE %%%%%%%%%%%%%%%%%%%%%%%

%% Initial Trajectory definition with Collision Avoidance only up to TH
% This will be the trajectory followed in the first horizon, since at the
% beginning I do not have nominal trajectories.
x_nom = zeros(x_len, T, N);

x_nom0=zeros(x_len,T,N);
for ii=1:N
    x_nom0(:,:,ii) = repmat(x0(:,ii), 1, T); % Initialization considering the same nominal trajectory for each time k, equal to the initial conditions
end

tic
parfor i = 1:N
    
    [ x_nom(:,:,i), ~] = RK4_fct_TRANSL_ONLY( x_nom0 , 1 , sol , i,[1:N], x_len, u_len, B, E, fi_ORF, N, T, TH, Dt, skip, Mref, Umax, Vmax, G, Rcol, Rcomm);
        
end
toc
  
%% Plot
fig1=figure;
% figure('units','normalized','outerposition',[0 0 1 1])
plot3(X,Y,Z,'s','MarkerSize', 5), hold on, grid on 
plot3(X0,Y0,Z0,'x'), hold on, xlabel x, ylabel y, zlabel z, grid on
xlabel 'X [m]', ylabel 'Y [m]', zlabel 'Z [m]';
for ii=1:N
    plot3(x_nom(1,:,ii), x_nom(2,:,ii), x_nom(3,:,ii), '--')
end



%% SCP Iteration
x = zeros(x_len, T, N); 
u = zeros(u_len,T-1,N);
Neig = zeros(N, T-1, N); % j X k X i: j Neighbors of i, for each time k 
xi_tmp=[]; Neig_i_tmp=[]; ui_tmp=[];
x_old=zeros(size(x));
SCP_convergence=zeros(N,T,2);
eps_SCP=60; % 2 without noise
w=1;


K=1:N;

tic
while ~isempty(K)  %length(K)~=0

    parfor ii = 1:length(K)            % exact robot indices - so xi of the robots that alreeady converged will not be modified and thei last value will remain
%        Where xi=x(:,:,i)
        [ xi_tmp(:,:,ii) , Neig_i_tmp(:,:,ii) , ui_tmp(:,:,ii) ] = RK4_fct_TRANSL_ONLY( x_nom , 1 , sol , ii,K , x_len, u_len, B, E, fi_ORF, N, T, TH, Dt, skip, Mref, Umax, Vmax, G, Rcol, Rcomm);

    end
    
    for ii = 1:length(K) 
        i = K(ii);
        x(:,:,i) = xi_tmp(:,:,ii);
        Neig(:,:,i) = Neig_i_tmp(:,:,ii);
        u(:,:,i) = ui_tmp(:,:,ii);
    end
    toc
    
    flag=zeros(1,N);
    
    parfor i = 1 : N          %%%%%%%%%%% % K %
            
        % The new trajectories I found in the previous for loop will be the
        % nominal ones for the following SCP iteration (w+1)
        x_nom(:,:,i) = x(:,:,i);
        
        SCP_convergence(i,:,w) = norms( x(:,:,i)-x_old(:,:,i) , Inf , 1); % Error for every robot i at all times k, at each SCP iteraton w (ixTxw)
        
        count=0;
        if w>=2 && sum(norms( x(:,:,i)-x_old(:,:,i) , Inf , 1)<=eps_SCP)==T 
            for n_step = 1:T-1
                for j = Neig(:,n_step,i)'
                    if j==0
                        count=count+1;
                    else
                        if norm(x(1:3,n_step,i)-x(1:3,n_step,j))>Rcol
                            count=count+1;
                        end
                    end
                end
            end
        end
        
        if count==N*(T-1)
%             K(find(K==i))=[];
            flag(i)=1;
        else
            flag(i)=0;
        end
        x_old_save(:,:,i)=x_old(:,:,i); %saving the older one, even though I don't need it because if <eps than one or the other is the same
        
        x_old(:,:,i) = x(:,:,i); % So =x_nom(:,:,i) ?
    end
    
    for i=1:N
        if flag(i)==1
            K(find(K==i))=[];
        end
    end
    
    K
    w=w+1
end
toc
%% 
 w=w-1;
 x = x_old_save; % Is the Approximate solution to the Trajectory Optimization Problem [Probl3],
               % ( Approximated by the Decentralized Convex Problem [Probl5] )
                
         
% Plots
for ii=1:N
    plot3(x(1,:,ii), x(2,:,ii), x(3,:,ii))
end

% i-th Trajectory Convergence in Time

% Trajectory Convergence of a generic i-th Robot wrt SCP iterations
legstr=[];
figure
for in=1:w
    plot(SCP_convergence(end,:,in)) , hold on
    legstr{in}=strcat('w=', num2str(in));
end
grid on, title('i-th Robot Trajectory Convergence', 'Interpreter', 'latex'), xlabel ('Time Intervals', 'Interpreter', 'latex'), ylabel ('$||\epsilon_{SCP}||_1$', 'Interpreter', 'latex') 
legend(legstr)

legstr=[];
figure
for in=1:w
    plot(SCP_convergence(1,:,in)) , hold on
    legstr{in}=strcat('w=', num2str(in));
end
grid on, title('i-th Robot Trajectory Convergence', 'Interpreter', 'latex'), xlabel ('Time Intervals', 'Interpreter', 'latex'), ylabel ('$||\epsilon_{SCP}||_1$', 'Interpreter', 'latex') 
legend(legstr)

legstr=[];
figure
for in=1:w
    plot(SCP_convergence(2,:,in)) , hold on
    legstr{in}=strcat('w=', num2str(in));
end
grid on, title('i-th Robot Trajectory Convergence', 'Interpreter', 'latex'), xlabel ('Time Intervals', 'Interpreter', 'latex'), ylabel ('$||\epsilon_{SCP}||_1$', 'Interpreter', 'latex') 
legend(legstr)

tk=1:T;
% State Parameters after SCP (considering only the last point)
figure, suptitle('    OFFLINE INITIAL TRAJECTORY GENERATION')
% rho_i
ax1=subplot('311'); plot(tk,x(1,:,end), tk,x(2,:,end), tk,x(3,:,end) ,'LineWidth',2)
set(gca,'FontSize',13)
title('Translational Components $\rho_i$', 'Interpreter', 'latex'), ylabel ('Position [m]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg1=legend('$\rho_x$','$\rho_y$','$\rho_z$'); set(leg1, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
% Translational Speed rho_d_i
ax2=subplot('312'); plot(tk,x(4,:,end),tk,x(5,:,end),tk,x(6,:,end) ,'LineWidth',2)
set(gca,'FontSize',13)
title('Translational Velocity Components $\dot{\rho_i}$', 'Interpreter', 'latex'), ylabel ('Speed [m/s]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg2=legend('$\dot{\rho_x}$','$\dot{\rho_y}$','$\dot{\rho_z}$'); set(leg2, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
% Control Action u_i
ax3=subplot('313'); plot(tk,[0 u(1,:,end)],tk,[0 u(2,:,end)],tk,[0 u(3,:,end)] ,'LineWidth',2)
set(gca,'FontSize',13)
title('Control Action $u_i$', 'Interpreter', 'latex'), ylabel ('Force [N]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg3=legend('$u_x$','$u_y$','$u_z$'); set(leg3, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
xlabel ('Time Intervals', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold')
linkaxes([ax1,ax2,ax3],'x'), xlim([1,T])
%%      %%%%%%%%%%%%%%%%%%%%%%% ONLINE %%%%%%%%%%%%%%%%%%%%%%%
pause
clc
%%      MPC - Real Time
ind_save=0; n_step_RT_old=0; fi_ORF_noise=zeros(3,N);
x_Ref = x; t_save=[]; tt_save =0; 
counter = 0;

Ts = .5;



% Time Loop
% 0. Time Starts
tic
 
% 1. Upload Initial Control Sequence u to the Real System 
        % - Non Lin Real Sys with input u0

% 2. MEANWHILE(!) Update k0 to Current (Initial) Time and Begin the Prediction over the Horizon
%    As soon as I finish upload the new control sequence to the Real
%    System, and then wait for the next sample(n_stepRT) and start again
%    the calculation over the TH.

% 3. So after finishing the MPC calulation on the horizon, so after t_RUN,
%    evaluate the new state of the Real System by integrating the Non
%    Linear System Model from the actual time until actual+t_RUN

% 4. That one is going to be the new measure for the restart of the MPC
 

t_actual = toc; %    tf-t_actual>tH
k0 = floor( t_actual/Dt ); %     k0 in [0,...,T-1] Indicates the Discrete Time Instants
n_step_RT = k0 + 1; %    Real Time (Index)  -  It's the time gap number I am in  n_step_RT in [1,...,T]

while n_step_RT < T 
    if n_step_RT~=n_step_RT_old
        t_save=[t_save; t_actual];

% %     x_old = x;
% %     u_old = u;
% %     Neig_old = Neig;
% %         % Noisy Measure (Simulating the effects of the Real System as Noises added)
% %     parfor i=1:N %%%(maybe a for is quicker)
% %         x(:,n_step_RT,i) = x(:,n_step_RT,i) + [-0.5+2*0.5*rand(3,1); -.2+2*.2*rand(3,1)];
% %         fi_ORF_noise(:,i)=fi_ORF(:,1)-.2*max(fi_ORF(:,1))+2*.2*max(fi_ORF(:,1))*rand(3,1);
% %     end
% %     %% No, I should feed the obtained u to the Non linear system and after I finish
% %     %   the MPC-SCP I should take the actual value from the non lin system
    

    %% Iniitial MPC value is the last value of the Real System:
    x(:, n_step_RT, :) = xx_save([1:3 4:6], end, :)  + ( -[1*ones(3,1);.01*ones(3,1)]+2*[1*ones(3,1);.01*ones(3,1)].*rand(6,1) ) ; % ( -[.1*max(xx_save(1:3)); .01*max(xx_save(8:10))]+2*[.1*max(xx_save(1:3))*rand(3,1); .01*max(xx_save(8:10))*rand(3,1)] ) ;  % + Measurement Noise
    
    %%  SCP Iteration
    xi_tmp=[]; Neig_i_tmp=[]; ui_tmp=[];
    x_old=zeros(size(x));
    SCP_convergence=zeros(N,T,2);
    w = 1;
    
    K=1:N;
    while ~isempty(K)
        
        parfor ii = 1:length(K)           % exact robot indices - so xi of the robots that alreeady converged will not be modified and thei last value will remain
            
            [ xi_tmp(:,:,ii) , Neig_i_tmp(:,:,ii) , ui_tmp(:,:,ii) ] = RK4_fct_TRANSL_ONLY( x , n_step_RT, sol , ii,K , x_len, u_len, B, E, fi_ORF, N, T, TH, Dt, skip, Mref, Umax, Vmax, G, Rcol, Rcomm);
            
        end
        
        for ii = 1:length(K)            % Need to apply from the next instant (index:2), since I am @ 1
            i = K(ii);
            x(:,n_step_RT:end,i) = xi_tmp(:,:,ii);
            Neig(:,n_step_RT:end,i) = Neig_i_tmp(:,:,ii);
            u(:,n_step_RT:end,i) = ui_tmp(:,:,ii);
            % try creating another fct just for the final steps in which xi and ui are longer(like saying that I want to reach tgt earlier than tf, and stay there for a while)
        end
        
        flag=zeros(1,N);
        
        parfor i = 1 : N        %%%%%%%%%%% % K %
            
            % The new trajectories I found in the previous for loop will be the
            % nominal ones for the following SCP iteration (w+1)
% %             x_nom(:,:,i) = x(:,:,i);
            
            SCP_convergence(i,:,w) = norms( x(:,:,i)-x_old(:,:,i) , Inf , 1); % Error for every robot i at all times k, at each SCP iteraton w (ixTxw)
            
            count=0;
            if w>=2 && sum(norms( x(:,:,i)-x_old(:,:,i) , Inf , 1)<=eps_SCP)==T
                for n_step = 1:T-1
                    for j = Neig(:,n_step,i)'
                        if j==0
                            count=count+1;
                        else
                            if norm(x(1:3,n_step,i)-x(1:3,n_step,j))>Rcol
                                count=count+1;
                            end
                        end
                    end
                end
            end
            
            if count==N*(T-1)
                %             K(find(K==i))=[];
                flag(i)=1;
            else
                flag(i)=0;
            end
            x_old_save(:,:,i)=x_old(:,:,i); %saving the older one, even though I don't need it because if <eps than one or the other is the same
            
            x_old(:,:,i) = x(:,:,i); % So =x_nom(:,:,i) ?
        end
        
        for i=1:N
            if flag(i)==1
                K(find(K==i))=[];
            end
        end
        
        K;
        w=w+1;
    end
    w=w-1;
    x = x_old_save; 
    %%%%%%%%%%%%%%%%% end SCP 
    
    %% Meanwhile the Real System goes...
    t_run = toc - t_actual;
    if t_run<Dt
        pause(Dt-t_run);
        t_run = toc - t_actual;
    end
    
            t_actualOld = t_actual;
            t_actual = toc % considering all other propagation algorithms that follows as in the mpc
                            % since in reality I would give the new
                            % trajectory as soon as MPC finishes, but here
                            % after MPC finishes I have to propagate the
                            % Real System (that in reality would mean just
                            % to acquire the measure).
    
        
    % - REAL SYSTEM
    % a. integrate from t_actual up to t_actual+n_run*Dt in cont time using
    % the old traj  -> write x and u up to the instant ceil((t_run+t_actual)/Dt)
    % which will be where the next MPC (the one just found (new)) starts
    counter=counter+1;
    
% % %     timeline = linspace(ceil(t_actualOld/Dt)*Dt , ceil((t_actualOld+t_run)/Dt)*Dt, valuesXgap);
% % %     if count==1
% % %         timeline = linspace(t_actualOld , ceil((t_actualOld+t_run)/Dt)*Dt, valuesXgap);
% % %     end 
% % %    
% % %     for yk = 1:n_split
% % %         split_timeline=timeline((yk-1)*valuesXsplit+1:yk*valuesXsplit);
% % %         if yk~=1
% % %             split_timeline = [timeline((yk-1)*valuesXsplit)  split_timeline]; % valuesXsplit+1 values
% % %         end
        
    tspan = [ tt_save(end) , ceil((t_actualOld+t_run)/Dt)*Dt ];

    valuesXgap = length( tspan(1):Ts:tspan(end) );  % CAREFUL! NOW valuesXgap CHANGES EACH TIME
    tt_temp=zeros(1,valuesXgap,N);
    xx_temp=zeros(state,valuesXgap,N); uu_temp=zeros(3,valuesXgap,N);
    timeline=tspan(1):Ts:tspan(end);
    endColumn=size(tt_save,2);

    x0_vec=xx_save(:, end, :);
    
    uu_i=zeros(3,valuesXgap);

    parfor i = 1 : N

    x0=x0_vec(:,:,i);

    [tt_i, xx_i] = odeRK4 ( @(t_ode,x_ode) odefcn_RegSys_I_ExtDist ( t_ode, x_ode, x_Ref(:,:,i), Dt,tf, m, fi_ORF(:,i), Kp, Kd, Umax, sol, xx0([1:3 4:6],1,i) , i,endColumn,timeline ), tspan , x0 , Ts );
    tt_temp(:,:,i) = tt_i;
    xx_temp(:,:,i) = xx_i;
    uu_temp(:,:,i) = uu_i;

    end
    for i=1:N
        tt_save(:, endColumn:endColumn+valuesXgap-1, i) = tt_temp(:,:,i);
        xx_save(:, endColumn:endColumn+valuesXgap-1, i) = xx_temp(:,:,i);
        uu_save(:, endColumn:endColumn+valuesXgap-1, i) = uu_temp(:,:,i);
%         tt_save(1, valuesXgap*(count-1)+1:valuesXgap*(count-1)+1+valuesXgap-1,i) = tt_temp(:,:,i);
%         xx_save(:, valuesXgap*(count-1)+1:valuesXgap*(count-1)+2+valuesXgap-1, i) = xx_temp(:,:,i);
    % THIS IS WHAT THE REAL SYS DOES, DON'T HAVE TO ADD NOISE OF MEASURE HERE BUT IN ODE
    end
    
    
% % %     end
    
    x_Ref = x;


        
    
    else
        t_actual = toc %    tf-t_actual>tH
    end
    
    %% Actual Time Update
    k0 = floor( t_actual/Dt );%     k0 in [0,...,T-1]
    n_step_RT_old=n_step_RT;
    n_step_RT = k0 + 1 %    Real Time

end  

% %     % Finish integrating the system up to T...
    if ceil((t_actualOld+t_run)/Dt)*Dt < T  
        
        tspan = [ tt_save(end) , tf ];
        valuesXgap = length( tspan(1):Ts:tspan(end) );  % CAREFUL! NOW valuesXgap CHANGES EACH TIME
        tt_temp=zeros(1,valuesXgap,N);
        xx_temp=zeros(state,valuesXgap,N); uu_temp=zeros(3,valuesXgap,N);
        timeline=tspan(1):Ts:tspan(end);
        endColumn=size(tt_save,2);
        x0_vec=xx_save(:, end, :);
        
        uu_i=zeros(3,valuesXgap);
        parfor i = 1 : N
            x0=x0_vec(:,:,i);
            [tt_i, xx_i] = odeRK4 ( @(t_ode,x_ode) odefcn_RegSys_I_ExtDist ( t_ode, x_ode, x_Ref(:,:,i), Dt,tf, m, fi_ORF(:,i), Kp, Kd, Umax, sol, xx0([1:3 4:6],1,i) , i,endColumn,timeline ), tspan , x0 , Ts );
            tt_temp(:,:,i) = tt_i;
            xx_temp(:,:,i) = xx_i;
            uu_temp(:,:,i) = uu_i;
        end
        for i=1:N
            tt_save(:, endColumn:endColumn+valuesXgap-1, i) = tt_temp(:,:,i);
            xx_save(:, endColumn:endColumn+valuesXgap-1, i) = xx_temp(:,:,i);
            uu_save(:, endColumn:endColumn+valuesXgap-1, i) = uu_temp(:,:,i);
        end
        
    end
        
toc    
  



%%  Plots
figure(fig1), hold on
for ii=1:N
    plot3(x(1,:,ii), x(2,:,ii), x(3,:,ii),'-o')
end

tk=1:T;
% State Parameters after SCP (considering only the last point)
figure, suptitle('      AFTER MPC - SCP')
% rho_i
ax1=subplot('311'); plot(tk,x(1,:,end), tk,x(2,:,end), tk,x(3,:,end) ,'LineWidth',2)
set(gca,'FontSize',13)
title('Translational Components $\rho_i$', 'Interpreter', 'latex'), ylabel ('Position [m]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg1=legend('$\rho_x$','$\rho_y$','$\rho_z$'); set(leg1, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
% Translational Speed rho_d_i
ax2=subplot('312'); plot(tk,x(4,:,end),tk,x(5,:,end),tk,x(6,:,end) ,'LineWidth',2)
set(gca,'FontSize',13)
title('Translational Velocity Components $\dot{\rho_i}$', 'Interpreter', 'latex'), ylabel ('Speed [m/s]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg2=legend('$\dot{\rho_x}$','$\dot{\rho_y}$','$\dot{\rho_z}$'); set(leg2, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
% Control Action u_i
ax3=subplot('313'); plot(tk,[0 u(1,:,end)],tk,[0 u(2,:,end)],tk,[0 u(3,:,end)] ,'LineWidth',2)
set(gca,'FontSize',13)
title('Control Action $u_i$', 'Interpreter', 'latex'), ylabel ('Force [N]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold') 
leg3=legend('$u_x$','$u_y$','$u_z$'); set(leg3, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
xlabel ('Time Intervals', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold')
linkaxes([ax1,ax2,ax3],'x'), xlim([1,T])

% MPC VS REAL SYSTEM
figure, % suptitle('      Nonlinear Dynamic System')
% STATE (last bot)
% rho_i
ax1=subplot('211'); plot(tt_save(1,:,end),xx_save(1,:,end) , tt_save(1,:,end),xx_save(2,:,end) , tt_save(1,:,end),xx_save(3,:,end) , tk*Dt,x(1,:,end),'-.',tk*Dt,x(2,:,end),'-.',tk*Dt,x(3,:,end),'-.' ,'LineWidth',2)
set(gca,'FontSize',14)
title('Translational Components $\rho_i$', 'Interpreter', 'latex'), ylabel ('\boldmath$\rho_i$ [m]', 'Interpreter', 'latex','FontSize',18,'FontWeight','bold') 
leg1=legend('$\rho_x$','$\rho_y$','$\rho_z$' , '$\rho_{k,x}$','$\rho_{k,y}$','$\rho_{k,z}$'); set(leg1, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
xlabel ('Time [s]', 'Interpreter', 'latex','FontSize',16,'FontWeight','bold')
% Translational Speed rho_d_i
ax2=subplot('212'); plot(tt_save(1,:,end),xx_save(4,:,end) , tt_save(1,:,end),xx_save(5,:,end) , tt_save(1,:,end),xx_save(6,:,end) , tk*Dt,x(4,:,end),'-.',tk*Dt,x(5,:,end),'-.',tk*Dt,x(6,:,end),'-.' ,'LineWidth',2)
set(gca,'FontSize',14)
title('Translational Velocity Components $\dot{\rho_i}$', 'Interpreter', 'latex'), ylabel ('\boldmath$\dot{\rho}_i$ [m/s]', 'Interpreter', 'latex','FontSize',18,'FontWeight','bold') 
leg2=legend('$\dot{\rho_x}$','$\dot{\rho_y}$','$\dot{\rho_z}$' , '$\dot{\rho}_{k,x}$','$\dot{\rho}_{k,y}$','$\dot{\rho}_{k,z}$'); set(leg2, 'Interpreter', 'latex','FontSize',14,'FontWeight','bold');
xlabel ('Time [s]', 'Interpreter', 'latex','FontSize',16,'FontWeight','bold')
linkaxes([ax1,ax2],'x')


figure
% MPC VS REAL SYSTEM
% Translation Control Effort - DISCRETE (last bot)
tk_c=1:T-1;
ax1=subplot('311'); plot(tk_c,u(1,:,end) ,  tk_c,uu_discrete(1,2:end,end),'-s','MarkerSize',4)
title('MPC vs Real System - Control Efforts', 'Interpreter', 'latex'), ylabel ('Force [N]', 'Interpreter', 'latex') 
leg1=legend('$u_{k,x} - MPC_{opt}$','$u_{k,x} - Real$'); set(leg1, 'Interpreter', 'latex');
ax2=subplot('312'); plot(tk_c,u(2,:,end) ,  tk_c,uu_discrete(2,2:end,end),'-s','MarkerSize',4)
ylabel ('Force [N]', 'Interpreter', 'latex') 
leg2=legend('$u_{k,y} - MPC_{opt}$','$u_{k,y} - Real$'); set(leg2, 'Interpreter', 'latex');
ax3=subplot('313'); plot(tk_c,u(3,:,end) ,  tk_c,uu_discrete(3,2:end,end),'-s','MarkerSize',4)
ylabel ('Force [N]', 'Interpreter', 'latex') 
leg3=legend('$u_{k,z} - MPC_{opt}$','$u_{k,z} - Real$'); set(leg3, 'Interpreter', 'latex');
xlabel ('Time Intervals', 'Interpreter', 'latex')
linkaxes([ax1,ax2,ax3],'x'), xlim([1,T-1])      % could not compare discrete and cont times in same plot

% Translation Control Effort - CONTINUOUS (last bot)
figure
plot(tt_save(1,:,end),uu_save(1,:,end), tt_save(1,:,end),uu_save(2,:,end) , tt_save(1,:,end),uu_save(3,:,end) ,'LineWidth',2 );
set(gca,'FontSize',15)
title('Translational Control Efforts', 'Interpreter', 'latex'), ylabel ('$\bf{u_i\ [N]}$', 'Interpreter', 'latex')
leg1=legend('\boldmath$u_x$','\boldmath$u_y$','\boldmath$u_z$'); set(leg1, 'Interpreter', 'latex');
xlabel ('\bf{Time [s]}', 'Interpreter', 'latex')



% All bots 3D
fig2=figure;
% figure('units','normalized','outerposition',[0 0 1 1])
plot3(X,Y,Z,'s','MarkerSize', 5), hold on, grid on 
plot3(X0,Y0,Z0,'x')
xlabel ('X [m]','FontWeight','bold'), ylabel ('Y [m]','FontWeight','bold'), zlabel ('Z [m]','FontWeight','bold');
set(gca,'FontSize',20)
for ii=1:N
    plot3(x(1,:,ii), x(2,:,ii), x(3,:,ii),'-.','LineWidth',2);
    plot3(xx_save(1,:,ii),xx_save(2,:,ii),xx_save(3,:,ii),'LineWidth',2);
end


% All Bots Position Error
figure, hold on, grid
for ii=1:N
ax1=subplot('311');hold on,
plot(tt_save(1,:,ii),(xx_save(1,:,ii)-Mref(1,ii)) , 'LineWidth',2) 
set(gca,'FontSize',13)
ylabel ('\boldmath$\rho_x-\rho_x^{Ref}$ [m]', 'Interpreter', 'latex','FontSize',16,'FontWeight','bold') 
ax2=subplot('312');hold on,
plot(tt_save(1,:,ii),(xx_save(2,:,ii)-Mref(2,ii)) , 'LineWidth',2) 
set(gca,'FontSize',13)
ylabel ('\boldmath$\rho_y-\rho_y^{Ref}$ [m]', 'Interpreter', 'latex','FontSize',16,'FontWeight','bold') 
ax3=subplot('313');hold on,
plot(tt_save(1,:,ii),(xx_save(3,:,ii)-Mref(3,ii)) , 'LineWidth',2) 
set(gca,'FontSize',13)
ylabel ('\boldmath$\rho_z-\rho_z^{Ref}$ [m]', 'Interpreter', 'latex','FontSize',16,'FontWeight','bold') 
xlabel ('Time [s]', 'Interpreter', 'latex','FontSize',14,'FontWeight','bold')
linkaxes([ax1,ax2,ax3],'x')
end

