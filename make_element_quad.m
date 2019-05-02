function [elements, nodes] = make_element_quad(j)
N = 21;
L_x = 1;
L_y = 1;
lambda_d1 = 0.2;
lambda_d2 = 0.2;
total_d_length = lambda_d1+lambda_d2;
nodes_doublelayer = j*10+1;

nodes_doublelayer = j*(floor(N*lambda_d1))+1;
% nodes_doublelayer = (floor(N*lambda_d1))+1;
% nodes_doublelayer = ww*50;
middle = linspace(lambda_d1,1 - lambda_d2,ceil(N*(1-total_d_length)));
beg    = linspace(0, lambda_d1,nodes_doublelayer);
end_nodes = linspace(1 - lambda_d2,1,nodes_doublelayer);
% M = 3000
% N = 15 + (ww-1)*10;
% N = 3;

% % % % % Example of high aspect ratio
% x = linspace(0,L_x,M);
% x = MeshBoth(M,0,1,1.02);
x = unique([beg,middle,end_nodes]);
% x = linspace(0,L_x,N);
% x = MeshBoth(nodes_doublelayer,0, L_x,1.08);
% % % x = [0 rand(1,98) 1];
% % % x = sort(x);
y = linspace(0,L_y,N);

[X,Y] = meshgrid(x,y);

Xvec = X(:);
Yvec = Y(:);
% % 
num_nodes = length(x)*length(y);
num_elements = (length(x)-1)*(length(y)-1);


nodes = [(1:num_nodes)',Xvec,Yvec];
first = (1:num_nodes-N)';
first(N:N:end) = [];
second=(N+1:num_nodes)';
second(N:N:end) = [];
third = (N+1:num_nodes)';
third(1:N:end) = [];
fourth = (1:num_nodes-N)';
fourth(1:N:end) = [];

% % 
elements = [(1:num_elements)', first,second,third,fourth];
end