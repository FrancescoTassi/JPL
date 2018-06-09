function [ Ad , Bd , cd ] = UNI_Swarm_Transl_StateMatrices( xi_nom_k , x_OE , B, Dt )
%% State Matrices which are dependent from the Nominal (Equilibrium) State and from the Orbital Elements

muE=3.986005e14; %m^3/s^2
kJ2=2.633e10*10^6; %m^2/s^2   %%%%%%%%%%%%wrong?
% muE=3.986005e5;   %km^3/s^2
% kJ2=2.633e10; %km^2/s^2

%% State @ Equilibrium (Nominal) @ Time Instant k
rho1 = xi_nom_k(1);
rho2 = xi_nom_k(2);
rho3 = xi_nom_k(3);

%% Orbital Parameters @ Time Instant k (considered in equilibrium)
% x_OE is a 6 elements vector containing the Orbital Parameters at that time
% instant k:
r=x_OE(1);  % R0
vx=x_OE(2); % RADIAL VELOCITY
h=x_OE(3);
i=x_OE(4);
Om=x_OE(5);
theta=x_OE(6);

%% Orbital OMEGA (angular speed in LVLH Frame (x,y,z))
%  ri = R0 + rho_i (in vector form):
% r_ORF = sqrt( (r+rho1)^2+rho2^2+rho3^2 );  % If r and rho are referred to ORF -> ORF
%    so r = r_ORF = R0_ORF_x (which is a distance >0) -> R0_ORF=[-r 0 0]
R0_ORF = [-r 0 0];
R0_ORF_skw=[0 -R0_ORF(3) R0_ORF(2);
          R0_ORF(3) 0 -R0_ORF(1);
          -R0_ORF(2) R0_ORF(1) 0];
      
OM1 = -kJ2*sin(2*i)*sin(theta)/h/r^3;    % ORF
OM2 = 0;
OM3 = h/r^2;

OM_skw_ORF = [0 -OM3 OM2; OM3 0 -OM1; -OM2 OM1 0];

%% State Matrices @ Nominal Conditions
% gyro = [ rho2*((2*h*vx)/r^3 + (kJ2*sin(2*theta)*sin(i)^2)/r^5) - r*(kJ2/r^5 + muE/r^3 - kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) + (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (5*kJ2*sin(i)^2*sin(theta)^2)/r^5) + rho1*(kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - h^2/r^4 - (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + sin(i)*sin(theta)*((2*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (2*kJ2*sin(i)*sin(theta))/r^4) - (kJ2*rho3*sin(2*i)*sin(theta))/r^5;
%          rho3*((kJ2*sin(2*i)*cos(theta))/r^5 - (3*kJ2*vx*sin(2*i)*sin(theta))/(h*r^4) + (8*kJ2^2*cos(i)*cos(theta)*sin(i)^3*sin(theta)^2)/(h^2*r^6)) - rho2*(h^2/r^4 - muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (kJ2^2*sin(2*i)^2*sin(theta)^2)/(h^2*r^6)) - rho1*((2*h*vx)/r^3 + (kJ2*sin(2*theta)*sin(i)^2)/r^5) + cos(theta)*sin(i)*((2*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (2*kJ2*sin(i)*sin(theta))/r^4);
%          rho3*(kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (kJ2^2*sin(2*i)^2*sin(theta)^2)/(h^2*r^6)) + cos(i)*((2*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (2*kJ2*sin(i)*sin(theta))/r^4) - rho2*((kJ2*sin(2*i)*cos(theta))/r^5 - (3*kJ2*vx*sin(2*i)*sin(theta))/(h*r^4) + (8*kJ2^2*cos(i)*cos(theta)*sin(i)^3*sin(theta)^2)/(h^2*r^6)) - (kJ2*rho1*sin(2*i)*sin(theta))/r^5 ];
% dgyro_drhox = [ kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - rho1*((5*kJ2*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + (3*muE*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (35*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + (10*kJ2*sin(i)*sin(theta)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - r*((5*kJ2*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + (3*muE*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (35*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + (10*kJ2*sin(i)*sin(theta)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - h^2/r^4 - (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + sin(i)*sin(theta)*((2*kJ2*sin(i)*sin(theta))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (5*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2));
%                 cos(theta)*sin(i)*((2*kJ2*sin(i)*sin(theta))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (5*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - (2*h*vx)/r^3 - rho2*((5*kJ2*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + (3*muE*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (35*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + (10*kJ2*sin(i)*sin(theta)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - (kJ2*sin(2*theta)*sin(i)^2)/r^5;
%                 cos(i)*((2*kJ2*sin(i)*sin(theta))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (5*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - rho3*((5*kJ2*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + (3*muE*(2*r + 2*rho1))/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (35*kJ2*(2*r + 2*rho1)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/(2*((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + (10*kJ2*sin(i)*sin(theta)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - (kJ2*sin(2*i)*sin(theta))/r^5 ];
% dgyro_drhoy = [ (2*h*vx)/r^3 - rho1*((5*kJ2*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (35*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2) + (10*kJ2*cos(theta)*sin(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - r*((5*kJ2*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (35*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2) + (10*kJ2*cos(theta)*sin(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - sin(i)*sin(theta)*((10*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (2*kJ2*cos(theta)*sin(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) + (kJ2*sin(2*theta)*sin(i)^2)/r^5;
%                  kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - rho2*((5*kJ2*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (35*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2) + (10*kJ2*cos(theta)*sin(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) + muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - h^2/r^4 - (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - cos(theta)*sin(i)*((10*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (2*kJ2*cos(theta)*sin(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (kJ2^2*sin(2*i)^2*sin(theta)^2)/(h^2*r^6);
%                 (3*kJ2*vx*sin(2*i)*sin(theta))/(h*r^4) - rho3*((5*kJ2*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho2)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (35*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2) + (10*kJ2*cos(theta)*sin(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - (kJ2*sin(2*i)*cos(theta))/r^5 - cos(i)*((10*kJ2*rho2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (2*kJ2*cos(theta)*sin(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2)) - (8*kJ2^2*cos(i)*cos(theta)*sin(i)^3*sin(theta)^2)/(h^2*r^6) ];
% dgyro_drhoz = [ sin(i)*sin(theta)*((2*kJ2*cos(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (10*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - rho1*((5*kJ2*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + (10*kJ2*cos(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (35*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) - r*((5*kJ2*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + (10*kJ2*cos(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (35*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) - (kJ2*sin(2*i)*sin(theta))/r^5;
%                 cos(theta)*sin(i)*((2*kJ2*cos(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (10*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - rho2*((5*kJ2*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + (10*kJ2*cos(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (35*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + (kJ2*sin(2*i)*cos(theta))/r^5 - (3*kJ2*vx*sin(2*i)*sin(theta))/(h*r^4) + (8*kJ2^2*cos(i)*cos(theta)*sin(i)^3*sin(theta)^2)/(h^2*r^6);
%                 cos(i)*((2*kJ2*cos(i))/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) - (10*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2)) - rho3*((5*kJ2*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) + (3*muE*rho3)/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + (10*kJ2*cos(i)*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i)))/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (35*kJ2*rho3*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(9/2)) + kJ2/((r + rho1)^2 + rho2^2 + rho3^2)^(5/2) + muE/((r + rho1)^2 + rho2^2 + rho3^2)^(3/2) - (5*kJ2*(rho3*cos(i) + sin(i)*sin(theta)*(r + rho1) + rho2*cos(theta)*sin(i))^2)/((r + rho1)^2 + rho2^2 + rho3^2)^(7/2) - (kJ2^2*sin(2*i)^2*sin(theta)^2)/(h^2*r^6) ];

Korb = OM_skw_ORF*(OM_skw_ORF*eye(3)) - muE * (2*norm(R0_ORF)^2*eye(3) + 3*R0_ORF_skw*(R0_ORF_skw*eye(3)))/norm(R0_ORF)^5;


A = [ zeros(3)   eye(3) 
      -Korb      -2*OM_skw_ORF ];
c = zeros(6,1); 

%% Discrete Time State Matrices
Ad = expm(A*Dt);
Bd = A\( Ad-eye(length(A)) ) * B ;
cd = c * Dt;


end

