function Unsteady = unsteady_term(yp)

Unsteady = 0*yp;
Unsteady(1:3:end) = yp(1:3:end);
Unsteady(2:3:end) = yp(2:3:end);

Unsteady = Unsteady(:);

end