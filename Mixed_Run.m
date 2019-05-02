filename = 'TrialMixed.msh'
[elementsTri,elementsQuad,BoundaryElements, nodes] = Element_Read_InTriBi(filename);


ref = BoundaryElements(:,1)==7|BoundaryElements(:,1)==8|BoundaryElements(:,1)==9;

BoundaryElements = BoundaryElements(ref,:);
BoundaryElements = nan;

%% Calculate control volume sizes.
Control_Volume = zeros(length(nodes)*3,1);
Control_Volume(1:3:end) = Volume_Calc(elementsTri,nodes)+Volume_Calc_Quad(elementsQuad,nodes);
Control_Volume(2:3:end) = Control_Volume(1:3:end);
Control_Volume(3:3:end) = Control_Volume(1:3:end);

%% Parameters for simulation
D = [2.83e-9,1.83e-9];
R = 8.31447215;
T = 298.15;
F = 9.6485339e4;
z = [-1, 1];
epsilon_0               = 8.854187817e-12;
epsilon                 = 80*epsilon_0;
sigma = -1*1.60217662e-19/(1e-9)^2;

C_0   =100;
Phi_0 = R*T/F;

L_0   = (20e-9);
T_0   = L_0^2/D(1);
Parameters.z = z;
Parameters.mu  = D*T_0/L_0^2;
Parameters.varphi = F/(R*T)*Phi_0;
Parameters.psi = Phi_0*epsilon/(F*L_0^2*C_0);
Parameters.sigma = @(t) sigma/(F*L_0*C_0)*tanh(1*t);
InitialConcentration = 100/C_0;
%% Voltage equation
Parameters.V = @(t)  -0.01/Phi_0;

%% Initial Conditions
InitConc    = zeros(length(nodes),1) + InitialConcentration;
% InitConc = nodes(:,2)
InitPhi     = 0*nodes(:,2);
y0          = zeros(length(nodes)*3,1);

%% To Solve for steady state
y0(1:3:end) = InitConc;
y0(2:3:end) = InitConc;
y0(3:3:end) = -0.0/Phi_0*nodes(:,2);

yp0             = zeros(3*length(nodes),1);
tic
yp0             = PNP_system_mixed(0,y0,yp0,elementsTri,elementsQuad,BoundaryElements, nodes, Parameters, Control_Volume);
toc
t0              = 0;
tend = 3000;

%% Define Mass Matrix (Will also include Structure of Jacobian)
diagonal = ones(3*length(nodes),1);
diagonal(3:3:end) = 0;
M = spdiags(diagonal,0,3*length(nodes),3*length(nodes));
dfdy = jacobian_calc_element_tri_quad(nodes,elementsTri,elementsQuad);

%% Solve equation System using ODE15i
options = odeset('RelTol',1e-4,'AbsTol',1e-6,'InitialStep',1e-6,'JPattern',{dfdy,M});
tic
[TOUT,XOUT] = ode15i(@(t,y,yp) PNP_system_mixed(t,y,yp,elementsTri,elementsQuad,BoundaryElements, nodes, Parameters, Control_Volume),[0 tend],y0,yp0,options);
sim_time(j) = toc;
fprintf('Simulation time = %f \n',sim_time(j))

subplot(3,1,1)
plot3(L_0*nodes(:,2),L_0*nodes(:,3),C_0*XOUT(end,1:3:end),'.')
title('Concentration of species one')
subplot(3,1,2)
plot3(L_0*nodes(:,2),L_0*nodes(:,3),C_0*XOUT(end,2:3:end),'.')
title('Concentration of species two')
subplot(3,1,3)
plot3(L_0*nodes(:,2),L_0*nodes(:,3),Phi_0*XOUT(end,3:3:end),'.')
title('Potential')