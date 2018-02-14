function Source = source_bilinear(~,y,~,Parameters,~,~)

psi = Parameters.psi;
z   = Parameters.z;


Source = zeros(3*Parameters.N,1);
Source((3:3:end)) = -psi*(z(1)*y((1:3:end)) + z(2)*y((2:3:end)));


end