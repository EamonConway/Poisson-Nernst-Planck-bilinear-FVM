function [dfdyp,dfdy] = jacobian_calc_element(nodes,element)
tic
main_diag = ones(length(nodes)*3,1);
main_diag(3:3:end) = 0;
dfdyp = spdiags(main_diag,0,length(nodes)*3,length(nodes)*3);
i = zeros(4*4*element(end,1),1);
j = i;
% k = i + 1;

%% Create adjacency matrix
for jj = 1:length(element)

    temp_nodes  = element(jj,2:end); 
    
    i(16*(jj-1)+1:16*(jj-1)+4)      = temp_nodes(1);
    j(16*(jj-1)+1:16*(jj-1)+4)      = temp_nodes;

    i(16*(jj-1)+5:16*(jj-1)+8)      = temp_nodes(2);
    j(16*(jj-1)+5:16*(jj-1)+8)      = temp_nodes;

    i(16*(jj-1)+9:16*(jj-1)+12)      = temp_nodes(3);
    j(16*(jj-1)+9:16*(jj-1)+12)      = temp_nodes;

    i(16*(jj-1)+13:16*(jj-1)+16)      = temp_nodes(4);
    j(16*(jj-1)+13:16*(jj-1)+16)      = temp_nodes;

end

i(i==0) = [];
j(j==0) = [];
A = unique([i,j],'rows');
i = A(:,1);
j = A(:,2);

i_vec = [3*i-2;3*i-2;...
    3*i-1;3*i-1;...
    3*i;3*(1:length(nodes))';3*(1:length(nodes))'];
    
    
j_vec = [3*j-2;3*j;...
    3*j-1;3*j; ...
    3*j;3*(1:length(nodes))'-2;3*(1:length(nodes))'-1];

k = 0*i_vec + 1;
dfdy = sparse(i_vec,j_vec,k);

toc
end