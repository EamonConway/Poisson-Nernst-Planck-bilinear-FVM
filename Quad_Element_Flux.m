function flux = Quad_Element_Flux(y,element,nodes,Parameters)
%% Initialisation 
flux = 0*y;

if isnan(element)
    return % Finish calculation of flux function if no quad elements.
end

%% Parameters
mu = Parameters.mu; % mu = D*T0/L0^2
psi = Parameters.psi; % psi = perm*Phi0/(L0^2*C0*F)
varphi = Parameters.varphi; %varphi = F/(R*T)*Phi0
valence = Parameters.z;
%% Integration points and weights
tau     = 0; % First Order
w       = 2; 
% 
% tau = [-1/sqrt(3),1/sqrt(3)]; % Second Order
% w = [1,1];

% tau = [-sqrt(3/5),0,sqrt(3/5)]; % Third Order
% w = [5/9,8/9,5/9];

%% Define all node values.
C1 = y(1:3:end);
C2 = y(2:3:end);
Phi = y(3:3:end);

%% Constants for Parameterisation of each face. Eg// xi = a_xi*tau + b_xi ; eta = a_eta*tau+b_eta
a_xi = [0,-0.5,0,0.5];
b_xi = [0,0.5,0,-0.5];
a_eta = [0.5,0,-0.5,0];
b_eta = [-0.5,0,0.5,0];

for k = 1:length(element(:,1));
    
    node_ref = element(k,2:end);
    z = nodes(node_ref,2);
    r = nodes(node_ref,3);
    C1_ref = C1(node_ref);
    C2_ref = C2(node_ref);
    Phi_ref = Phi(node_ref);
    face_1 = [1/2*(z(1)+z(2)),1/2*(r(1)+r(2))];
    face_2 = [1/2*(z(2)+z(3)),1/2*(r(2)+r(3))];
    face_3 = [1/2*(z(3)+z(4)),1/2*(r(3)+r(4))];
    face_4 = [1/2*(z(1)+z(4)),1/2*(r(1)+r(4))];
    
    h_xi  = face_3-face_1;
    h_eta = face_4-face_2;
    
    for j = 1:4 %This is quadrilateral specific
        
            for alpha = 1:length(tau)
                
                eta = a_eta(j)*tau(alpha) + b_eta(j);
                xi = a_xi(j)*tau(alpha) + b_xi(j);
        
            % Shape functions
        N1 = 1/4*(1-xi)*(1-eta);
        N2 = 1/4*(1+xi)*(1-eta);
        N3 = 1/4*(1+xi)*(1+eta);
        N4 = 1/4*(1-xi)*(1+eta); 
        
            % r term for integrand
        rface = r(1)*N1 + r(2)*N2+ r(3)*N3 + r(4)*N4;

            % Derivative of shape functions
        dN1 = [-1/4*(1-eta);-1/4*(1-xi)];
        dN2 = [ 1/4*(1-eta);-1/4*(1+xi)];
        dN3 = [ 1/4*(1+eta); 1/4*(1+xi)];
        dN4 = [-1/4*(1+eta); 1/4*(1-xi)];
        
            % Derivative of mappings (will form Jacobian)
        dz  = z(1)*dN1 + z(2)*dN2 + z(3)*dN3 + z(4)*dN4;        
        dr  = r(1)*dN1 + r(2)*dN2 + r(3)*dN3 + r(4)*dN4;

            % Concentration interpolant
        C1face = C1_ref(1)*N1 + C1_ref(2)*N2 + C1_ref(3)*N3 + C1_ref(4)*N4;
        C2face = C2_ref(1)*N1 + C2_ref(2)*N2 + C2_ref(3)*N3 + C2_ref(4)*N4;
        
            % Canonical Gradient of concentrations and potential interpolant
        dC1 = C1_ref(1)*dN1 + C1_ref(2)*dN2 + C1_ref(3)*dN3 + C1_ref(4)*dN4;
        dC2 = C2_ref(1)*dN1 + C2_ref(2)*dN2 + C2_ref(3)*dN3 + C2_ref(4)*dN4;
        dPhi = Phi_ref(1)*dN1 + Phi_ref(2)*dN2 + Phi_ref(3)*dN3 + Phi_ref(4)*dN4;
            
            % Differential Length Element
        drho = sqrt((dr(1)*a_xi(j)+dr(2)*a_eta(j))^2+(dz(1)*a_xi(j)+dz(2)*a_eta(j))^2);
            
            % Jacobian of mapping XI
        J = [dz(1) dr(1);dz(2) dr(2)];
        
            % Unit Normal
        n = [0 1;-1 0]*(J'*[a_xi(j);a_eta(j)]);
        n = n/norm(n);
        
            % Upwind check
        UpdPhi = J\dPhi;
        h_1 = 1/norm(UpdPhi)*UpdPhi(:)'*h_xi(:);
        h_2 = 1/norm(UpdPhi)*UpdPhi(:)'*h_eta(:);
        h = norm(h_1) + norm(h_2);
        Pe = norm(UpdPhi)*h;
        direction = -valence(1)*UpdPhi(:)'*n(:);
        if isnan(h)
            
        elseif Pe>=2 && direction>0
            C1face = C1_ref(j);
            C2face = C2_ref(mod(j,4)+1);
        elseif h>=2 && direction<0
            C1face = C1_ref(mod(j,4)+1);
            C2face = C2_ref(j);
%             t  
        end
        
            % Integrand at tau(alpha)
        GC1 = -mu(1)*((J\(dC1 + valence(1)*varphi*dPhi*C1face))'*n)*rface*drho;
        GC2 = -mu(2)*((J\(dC2 + valence(2)*varphi*dPhi*C2face))'*n)*rface*drho;
        GPhi = psi*((J\dPhi)'*n)*rface*drho;
        
            % Gaussian Quadrature
        flux(3*node_ref(j)-2) =  flux(3*node_ref(j)-2) + w(alpha)*GC1;
        flux(3*node_ref(mod(j,4)+1)-2) = flux(3*node_ref(mod(j,4)+1)-2) -  w(alpha)*GC1;
        
        flux(3*node_ref(j)-1) =  flux(3*node_ref(j)-1) + w(alpha)*GC2;
        flux(3*node_ref(mod(j,4)+1)-1) = flux(3*node_ref(mod(j,4)+1)-1) - w(alpha)*GC2;
        
        flux(3*node_ref(j))   =  flux(3*node_ref(j)) + w(alpha)*GPhi;
        flux(3*node_ref(mod(j,4)+1)) = flux(3*node_ref(mod(j,4)+1)) - w(alpha)*GPhi;

            end
    end

end
