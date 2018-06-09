function [ dxdt ] = odefcn_OE ( t, x )
% Orbital Elements EOM
% State: x_oe = {r vx h i Om theta}'

muE=3.986005e14;   %m^3/s^2
kJ2=2.633e10*10^6; %m^2/s^2

dxdt=zeros(length(x),1);
r=x(1);
vx=x(2);
h=x(3);
i=x(4);
Om=x(5);
theta=x(6);

dxdt(1) = vx;
dxdt(2) = -muE/r^2+h^2/r^3-kJ2/r^4*(1-3*sin(i)^2*sin(theta)^2);
dxdt(3) = -kJ2*sin(i)^2*sin(2*theta)/r^3;
dxdt(4) = -kJ2*sin(2*i)*sin(2*theta)/(2*h*r^3);
dxdt(5) = -2*kJ2*cos(i)*sin(theta)^2/h/r^3;
dxdt(6) = h/r^2+2*kJ2*cos(i)^2*sin(theta)^2/h/r^3;


