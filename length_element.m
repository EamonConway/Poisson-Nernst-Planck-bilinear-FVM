function [length, normal_x, normal_y] = length_element(x,y);
%% length_element - A trial function
% to calculate the length of a face after a bilinear transformation. 
t = [-1/sqrt(3),1/sqrt(3)];
t = 0;
length = zeros(1,4);
normal_x = length;
normal_y = length;
c_eta = [0.5 0 0.5 0];
c_xi  = [0 0.5 0 0.5];
a_eta= [-1 0 1 0];
a_xi = [0 1 0 -1];

    for j = 1:4
eta = c_eta(j)*(t + a_eta(j));
xi = c_xi(j)*(t + a_xi(j));

% plot(xi,eta,'.')

N1 = 0.25*(1-xi)*(1-eta);
N2 = 0.25*(1+xi)*(1-eta);
N3 = 0.25*(1+xi)*(1+eta);
N4 = 0.25*(1-xi)*(1+eta);

mid_x = x(1)*N1 + x(2)*N2 + x(3)*N3 + x(4)*N4;
mid_y = y(1)*N1 + y(2)*N2 + y(3)*N3 + y(4)*N4;
dN1 = [-0.25*(1-eta);-0.25*(1-xi)];
dN2 = [0.25*(1-eta); -0.25*(1+xi)];
dN3 = [0.25*(1+eta); 0.25*(1+xi)];
dN4 = [-0.25*(1+eta); 0.25*(1-xi)];


dx = x(1)*dN1 + x(2)*dN2 + x(3)*dN3 + x(4)*dN4;
dy = y(1)*dN1 + y(2)*dN2 + y(3)*dN3 + y(4)*dN4;

%% THis will only work forface 1at this opoint in time.
length_map = 2*sqrt((c_xi(j)*dx(1))^2 + (c_eta(j)*dx(2))^2 + (c_xi(j)*dy(1))^2 + (c_eta(j)*dy(2))^2);
length(j) = length_map;
% J = [dx(:)';dy(:)'];
J = [dx(:),dy(:)];
vec_n = [round(sin(j*pi/2));round(-cos(j*pi/2))];
normal = (J)\vec_n;
normal = normal/norm(normal);
normal_x(j) = normal(1);
normal_y(j) = normal(2);
    end
end
