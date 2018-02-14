function Unsteady = unsteady_bilinear(~,~,yp,Parameters,~,~,~)

Unsteady = zeros(Parameters.N*3,1);
Unsteady(1:3:end) = yp(1:3:end);
Unsteady(2:3:end) = yp(2:3:end);

Unsteady = Unsteady(:);

end