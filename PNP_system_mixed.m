function res = PNP_system_mixed(t,y,yp,elementsTri,elementsQuad, BoundaryElements, nodes, Parameters, Control_Volume)
V = Parameters.V(t);
Phi = y(3:3:end);
t
%% Calculate Unsteady
% unsteady = unsteady_term(t,y,yp,elements,nodes,Parameters);
unsteady = unsteady_term(yp);
%% Calculate flux 
% tic
fluxTri = Tri_Element_Flux(y,elementsTri,nodes,Parameters);
fluxQuad = Quad_Element_Flux(y,elementsQuad,nodes,Parameters);
% toc
BoundaryFlux = flux_boundary_phi(t,BoundaryElements ,nodes, Parameters);
flux = fluxTri + fluxQuad + BoundaryFlux;
%% Calculate Source
Source = source_term(y,Parameters);

res = unsteady + 1./Control_Volume.*flux - Source;

ref = nodes(nodes(:,2)==0,1);
res(3*ref) = y(3*ref);

ref = nodes(nodes(:,2)==1,1);
res(3*ref) = y(3*ref)-V;


