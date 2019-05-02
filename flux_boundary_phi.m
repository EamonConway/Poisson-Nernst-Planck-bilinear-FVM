function flux = flux_boundary_phi(t,elements,nodes, Parameters);
flux = 0*zeros(3*length(nodes),1);
if isnan(elements)
    return % Finish calculation of flux function if no boundary elements.
end
%% Parameters
sigma = Parameters.sigma(t); % Nondimensional surface charge

%% Integration points
tau     = 0; % First Order
w       = 2;

% tau = [-1/sqrt(3),1/sqrt(3)]; % Second Order
% w = [1,1];

% tau = [-sqrt(3/5),0,sqrt(3/5)]; % Third Order
% w = [5/9,8/9,5/9];

%% Constants for Parameterisation of each face. Eg// xi = a_xi*tau + b_xi ; eta = a_eta*tau+b_eta
a_xi = [1/4,1/4];
b_xi = [1/4,3/4];

%% Flux for boundary elements
for k = 1:length(elements(:,1))

    node_ref = elements(k,2:3);
    z = nodes(node_ref,2);
    r = nodes(node_ref,3);
    
    for j =1:2 % THis is for boundaries
        
            for alpha = 1:length(tau)
                xi = a_xi(j)*tau(alpha) + b_xi(j);
                
        N1 = (1-xi);
        N2 = (xi);
    
        dN1 = -1;
        dN2 = 1;
    
        rface = r(1)*N1 + r(2)*N2;
    
        dz  = z(1)*dN1 + z(2)*dN2;
        dr  = r(1)*dN1 + r(2)*dN2;
    
        drho = sqrt((dz(1)*a_xi(j))^2 + (dr(1)*a_xi(j))^2);

        GPhi = sigma*rface*drho;
        
        flux(3*node_ref(j)) =  flux(3*node_ref(j)) + w(alpha)*GPhi;
        
            end
    end
    
end
end