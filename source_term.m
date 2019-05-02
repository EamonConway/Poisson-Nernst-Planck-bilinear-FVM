function Source = source_term(y,Parameters)
Source = 0*y;
z = Parameters.z;
Source((3:3:end)) = -(z(1)*y((1:3:end)) + z(2)*y((2:3:end)));


end