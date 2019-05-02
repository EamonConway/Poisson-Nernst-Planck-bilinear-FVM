function Control_Volume = Volume_Calc(elements,nodes)

    %% Initialisation
    Control_Volume = 0*nodes(:,1);
    
    if isnan(elements)
        return
    end
    
    %% Loop over elements and add area
for i = 1:length(elements(:,1));
    node_ref = elements(i,2:end);
    midpoint = mean(nodes(node_ref,2:3));
    x = nodes(node_ref,2);
    y = nodes(node_ref,3);
    for j =1:3
        face1 = [0.5*(x(j) + x(mod(j,3)+1)),0.5*(y(j) + y(mod(j,3)+1))];
        face2 = [0.5*(x(j) + x(mod(j-2,3)+1)),0.5*(y(j) + y(mod(j-2,3)+1))];
        
        vertex = [nodes(node_ref(j),2:3);face1;midpoint;face2];
        r_k = mean(vertex(:,2));
%         r_k = 1;
        Control_Volume(node_ref(j)) = Control_Volume(node_ref(j)) + r_k*polyarea(vertex(:,1),vertex(:,2));
    end
   
end