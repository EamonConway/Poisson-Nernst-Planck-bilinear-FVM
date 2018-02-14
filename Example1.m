%% Example of high aspect ratio simulation of PNP equations.


N = 15; % number of nodes in vertical direction.
M = 101; % number of nodes in horizontal direction
L_x = 1; % horizontal length
L_y = 1; % vertical length

x = linspace(0,L_x,M); % x-location of nodes.
y = linspace(0,L_y,N); % y-location of nodes

%% Create Node and element structure.

[X,Y] = meshgrid(x,y); % create mesh.

Xvec = X(:);
Yvec = Y(:);

num_nodes = length(x)*length(y);
num_elements = (length(x)-1)*(length(y)-1);

nodes = [(1:num_nodes)',Xvec,Yvec];
first = (1:num_nodes-N)';
first(N:N:end) = [];
second=(N+1:num_nodes)';
second(N:N:end) = [];
third = (N+1:num_nodes)';
third(1:N:end) = [];
fourth = (1:num_nodes-N)';
fourth(1:N:end) = [];

element = [(1:num_elements)', first,second,third,fourth];

% Calculate control volume size and length of faces
Control_Volume = zeros(num_nodes,1);
DL = [(1:num_elements)',zeros(num_elements,4)];
for j = 1:element(end,1)
    
    temp_nodes = element(j,2:end);
    temp_node_coord = nodes(temp_nodes,2:end);
    temp_centre = mean(temp_node_coord);
    
    for i = 1:4
        point_1 = temp_node_coord(i,:);
        point_2 = 0.5*(temp_node_coord(i,:) + temp_node_coord(mod(i,4)+1,:));
        point_3 = temp_centre;
        point_4 = 0.5*(temp_node_coord(i,:) + temp_node_coord(4 - mod(4-i+1,4),:));
        points = [point_1;point_2;point_3;point_4];
        r_c    = mean(points(:,2));
        temp_volume = r_c*polyarea(points(:,1),points(:,2));
        Control_Volume(temp_nodes(i)) = Control_Volume(temp_nodes(i)) + temp_volume;
        face_cenre = temp_centre - point_2;
        DL(j,i+1) = norm(point_2-temp_centre);
       
    end
end

%% Defining Parameters for simulation
D = [2.83e-9,1.83e-9]; % Diffusion Coefficient
R = 8.31447215; % Gas Constant
T = 298.15; % Temperature
F = 9.6485339e4; % Faradays constant
z = [-1, 1];% Valence
epsilon_0               = 8.854187817e-12;
epsilon                 = 80*epsilon_0; % Permittivity

%% Nondimensional Parameters 
C_0   = 100;
Phi_0 = R*T/F;
T_0   = 5e-3;
L_0   = (20e-9);

%% Parameter structure.
Parameters.D = D;
Parameters.R = R;
Parameters.T = T;
Parameters.F = F;
Parameters.z = z;
Parameters.permit = epsilon;
Parameters.L_0  = L_0;
Parameters.C_0 = C_0;
Parameters.Phi_0 = Phi_0;
Parameters.T_0 = T_0;

% Nondimensional numbers
Parameters.mu  = T_0/L_0^2;
Parameters.varphi = F/(R*T)*Phi_0;
Parameters.psi = F*L_0^2*C_0/Phi_0;

Parameters.L_x  = max(nodes(:,2));
Parameters.L_y  = max(nodes(:,3));
Parameters.InitC = (100/C_0);
Parameters.N     = length(nodes); % number of nodes
Parameters.DL = DL;
InitConc = Parameters.InitC;
Parameters.V = @(t) -0.07/Phi_0*tanh(1000*t); % Voltage at working electrode.

%% Initial Conditions
y0  = zeros(3*num_nodes,1);
y0(1:3:end) =  InitConc;
y0(2:3:end) =  InitConc;
yp0 = zeros(3*num_nodes,1);
t_0 = 0; % Start time
t_end = 1; % end time
[dfdyp,dfdy] = jacobian_calc_element(nodes,element); % Calculate Jacobian structure.
% time_points = [];
options = odeset('JPattern',{dfdy,dfdyp},'InitialStep',1e-6,'AbsTol',1e-12,'RelTol',1e-6,'Stats','on'); 
% options = odeset('InitialStep',1e-6,'Events',@eventfun); 

yp0 = -bilinear_system(t_0,y0,yp0,Control_Volume,Parameters,element,nodes,BoundaryElements);

[TOUT,XOUT] = ode15i(@(t,y,yp) bilinear_system(t,y,yp,Control_Volume,Parameters,element,nodes,BoundaryElements),[t_0 t_end],y0,yp0,options);
% XOUT is nondimensional.
C_1 = XOUT(:,1:3:end)*Parameters.C_0; % Dimensional Solution to C_1
C_2 = XOUT(:,2:3:end)*Parameters.C_0; % Dimensional Solution to C_2;
Phi = XOUT(:,3:3:end)*Parameters.Phi_0; % Dimensional Solution to Phi.

%% Plot
figure
plot3(nodes(:,2),nodes(:,3),C_1(end,:),'.')
figure
plot3(nodes(:,2),nodes(:,3),C_2(end,:),'.')
figure
plot3(nodes(:,2),nodes(:,3),Phi(end,:),'.')
