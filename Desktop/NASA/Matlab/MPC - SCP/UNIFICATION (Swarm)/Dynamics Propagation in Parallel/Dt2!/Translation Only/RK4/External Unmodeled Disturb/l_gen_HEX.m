function [ X, Y, Z, Mref, l, ng ] = l_gen_HEX ( nn, r, H)
% HEXAGONAL grid (Final Desired Spatial Configuration)

% INPUT;
% C = [x y z] grid center (Target)
% nn = # number of subintervals
% R, the radius of the disk
% H = [x_1 ... x_h; y_1 ... y_h; z_1 ... z_h] (Obstacles)
% OUTPUT:
% l = (n*n)X(n*n+1+h) target distance matrix
% ng = Number of points

C=H(:,1); H(:,1)=[];
ng = disk_grid_count ( nn, r, C )
cg = disk_grid ( nn, r, C, ng )
% plot(cg(1,:),cg(2,:),'*'), hold on

X=cg(1,:); 
Y=cg(2,:);
% Elliptic Paraboloid
c=1;
a=2*sqrt((C(1)^2+C(2)^2)*c/C(3)); %sqrt(C(2)^2/(C(3)-C(1)^2/a^2)); % In order to center the paraboloi in Cz
b=a;
Z = c*( cg(1,:).^2/a^2 + cg(2,:).^2/b^2);
% plot3(cg(1,:),cg(2,:),Z,'*'), hold on, grid on
        
Mref=[X; Y; Z];


% Finding the Grid
l=zeros(ng);
for i=1:ng
    for j=1:ng
        l(i,j)=norm([X(i); Y(i); Z(i)]-[X(j); Y(j); Z(j)]);
    end
end

for i=1:ng
    l(i,ng+1)=norm([X(i); Y(i); Z(i)]-[C(1); C(2); C(3)]);
end

h=length(H(1,:));
if h
    for k=1:h
        for i=1:ng
            l(i,ng+1+k)=norm([X(i); Y(i); Z(i)]-[H(1,k); H(2,k); H(3,k)]);
        end
    end
end
