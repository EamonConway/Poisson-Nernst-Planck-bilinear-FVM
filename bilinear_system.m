function rr = bilinear_system(t,y,yp,Control_Volume,Parameters,element,nodes,BoundaryElements)
DL = Parameters.DL;


Source = source_bilinear(t,y,yp,Parameters,element,nodes);
Unsteady = unsteady_bilinear(t,y,yp,Parameters,element,nodes);
Control_Volume_L = ones(Parameters.N*3,1);

Flux = bilinear_flux_mex(t,y,Parameters,element,nodes,DL);

Control_Volume_L((1:3:end)) = Control_Volume;
Control_Volume_L((2:3:end)) = Control_Volume;
Control_Volume_L((3:3:end)) = Control_Volume;

rr = Unsteady + 1./Control_Volume_L.*Flux - Source;

%% Dirichlet Boundary Conditions
ref = nodes(nodes(:,2)==0,1);
rr(3*ref) = y(3*ref); % Ground Electrode

ref = nodes(nodes(:,2)==1,1);
rr(3*ref) = y(3*ref) - Parameters.V(t); % Working Electrode

end