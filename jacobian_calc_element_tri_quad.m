function [dfdy] = jacobian_calc_element_tri_quad(nodes,elementsTri,elementsQuad)
tic
% Quad_Length = length(
iQuad = [];
jQuad = [];
iTri  = [];
jTri  = [];
% k = i + 1;

%% Create adjacency matrix
if ~isnan(elementsQuad)
    iQuad = zeros(4*4*length(elementsQuad(:,1)),1);
    jQuad = iQuad;
for jj = 1:length(elementsQuad(:,1))

    temp_nodes  = elementsQuad(jj,2:end); 
    
    iQuad(16*(jj-1)+1:16*(jj-1)+4)      = temp_nodes(1);
    jQuad(16*(jj-1)+1:16*(jj-1)+4)      = temp_nodes;

    iQuad(16*(jj-1)+5:16*(jj-1)+8)      = temp_nodes(2);
    jQuad(16*(jj-1)+5:16*(jj-1)+8)      = temp_nodes;

    iQuad(16*(jj-1)+9:16*(jj-1)+12)      = temp_nodes(3);
    jQuad(16*(jj-1)+9:16*(jj-1)+12)      = temp_nodes;

    iQuad(16*(jj-1)+13:16*(jj-1)+16)      = temp_nodes(4);
    jQuad(16*(jj-1)+13:16*(jj-1)+16)      = temp_nodes;

end
iQuad(iQuad==0) = [];
jQuad(jQuad==0) = [];
A = unique([iQuad,jQuad],'rows');
iQuad = A(:,1);
jQuad = A(:,2);
end



if ~isnan(elementsTri)
    iTri = zeros(3*3*length(elementsTri(:,1)),1);
	jTri = iTri;
for jj = 1:length(elementsTri(:,1))

    temp_nodes  = elementsTri(jj,2:end); 
    
    iTri(9*(jj-1)+1:9*(jj-1)+length(temp_nodes))      = temp_nodes(1);
    jTri(9*(jj-1)+1:9*(jj-1)+length(temp_nodes))      = temp_nodes;

    iTri(9*(jj-1)+(4):9*(jj-1)+2*length(temp_nodes))      = temp_nodes(2);
    jTri(9*(jj-1)+(4):9*(jj-1)+2*length(temp_nodes))      = temp_nodes;

    iTri(9*(jj-1)+7:9*(jj-1)+3*length(temp_nodes))      = temp_nodes(3);
    jTri(9*(jj-1)+7:9*(jj-1)+3*length(temp_nodes))      = temp_nodes;

%     i(length(temp_nodes)^2*(jj-1)+(3*length(temp_nodes)+1):length(temp_nodes)^2*(jj-1)+416)      = temp_nodes(4);
%     j(length(temp_nodes)^2*(jj-1)+(3*length(temp_nodes)+1):length(temp_nodes)^2*(jj-1)+16)      = temp_nodes;

end
iTri(iTri==0) = [];
jTri(jTri==0) = [];
A = unique([iTri,jTri],'rows');
iTri = A(:,1);
jTri = A(:,2);
end
i =[iQuad;iTri];
j = [jQuad;jTri];

i_vec = [3*i-2;3*i-2;...
    3*i-1;3*i-1;...
    3*i;3*(1:length(nodes))';3*(1:length(nodes))'];
    
    
j_vec = [3*j-2;3*j;...
    3*j-1;3*j; ...
    3*j;3*(1:length(nodes))'-2;3*(1:length(nodes))'-1];

k = 0*i_vec + 1;
dfdy = sparse(i_vec,j_vec,k);
end