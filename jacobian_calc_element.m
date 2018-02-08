function [dfdyp,dfdy] = jacobian_calc_element(nodes,element,BoundaryElements)
tic
small = 0*1e-8;
main_diag = ones(length(nodes)*3,1);
main_diag(3:3:end) = small;
dfdyp = spdiags([0;main_diag;0],0,length(nodes)*3+2,length(nodes)*3+2);
% num_nodes_left = length(nodes(nodes(:,2)==0),1);
% num_nodes_right = length(nodes(nodes(:,2)==max(nodes(:,2)),1));
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
dfdy_temp = sparse(i_vec,j_vec,k);
dfdy = speye(3*length(nodes)+2,3*length(nodes)+2);
dfdy(2:end-1,2:end-1) = dfdy_temp;
clearvars dfdy_temp


%% Boundary 4 is left electrode
% ref = BoundaryElements(:,1)==1;
ref = BoundaryElements(:,1)==1;
ref = BoundaryElements(ref,2:end);
boundary4_nodes = unique(ref(:));
% boundary4_nodes = nodes(nodes(:,2) == min(nodes(:,2)),1);
%% Boundary 5 is right electrode
% ref = BoundaryElements(:,1)==2;
ref = BoundaryElements(:,1)==4;
ref = BoundaryElements(ref,2:end);
boundary5_nodes = unique(ref(:));
% boundary5_nodes = nodes(nodes(:,2) == max(nodes(:,2)),1);

%Left Electrode
% boundary4_nodes = nodes(nodes(:,2) == 0,1);
 dfdy(1,3*boundary4_nodes - 2 + 1) = 1;
 dfdy(1,3*boundary4_nodes + 1) = 1;

 dfdy(3*boundary4_nodes-2+1,1)=1;
 dfdy(3*boundary4_nodes+1,1) = 1;

dfdy(3*boundary4_nodes - 2 + 1,end) = 1; %Depends
dfdy(3*boundary4_nodes + 1,end) = 1;

dfdy(end,3*boundary4_nodes - 2 + 1) = 1; %Depends
dfdy(end,3*boundary4_nodes + 1) = 1;
% i_L_relation    = zeros(length(dfdy),1);
% i_L_relation(3*boundary4_nodes -2 + 1) = 1;

%Right Electrode
dfdy(end,3*boundary5_nodes - 2 + 1) = 1;
dfdy(end,3*boundary5_nodes + 1) = 1;
dfdy(3*boundary5_nodes - 2 + 1,end) = 1;
dfdy(3*boundary5_nodes + 1,end) = 1;
 
dfdy(3*boundary5_nodes - 2 + 1,1) = 1;
dfdy(3*boundary5_nodes + 1,1) = 1;
dfdy(1,3*boundary5_nodes - 2 + 1) = 1;
dfdy(1,3*boundary5_nodes + 1) = 1;

main_diag = ones(length(nodes)*3,1);
main_diag(3:3:end) = 0;
dfdyp = spdiags([0;main_diag;0],0,length(nodes)*3+2,length(nodes)*3+2);

% spy(dfdy)
toc
% end