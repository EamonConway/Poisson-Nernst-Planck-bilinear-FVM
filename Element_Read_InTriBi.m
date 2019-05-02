function [elementsTri, elementsQuad, BoundaryElements, nodes] = Element_Read_InTriBi(meshname)
% open mesh and read node and element data from Gmsh 'mesh' file structure
% This file is a modified version of a file originally created by Elliot
% Carr (2012).
%
% Modified by Brody Foy, September 2013
[fid, message] = fopen(meshname, 'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    switch tline
        case '$MeshFormat'
            disp('reading $MeshFormat')
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            MF=sscanf(tline,'%g %g %g')';   % read Mesh Format
            tline = fgetl(fid);             % read $EndMeshFormat
            if ~ischar(tline) || ~strcmp(tline,'$EndMeshFormat'), break, end
        case '$Nodes'
            disp('reading $Nodes')
            tline = fgetl(fid);
            no_nodes=sscanf(tline,'%g'); nodes=[];
            %nodes = fscanf(fid, '%g %g %g %g', [4 no_nodes])'; % read nodal data
            for i=1:no_nodes
                tline = fgetl(fid);
                nodes(i,:)=sscanf(tline,'%g %g %g %g');  % read nodes
            end
            tline = fgetl(fid);            % read $EndNodes
            if ~ischar(tline) || ~strcmp(tline,'$EndNodes'), break, end
        case '$Elements'
            disp('reading $Elements')
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            no_elements=sscanf(tline,'%g'); type15=[]; type1=[]; elementsTri=[]; elementsQuad = [];
            no_type15=0; no_type1=0; no_type2=0; no_type3=0;
            for i=1:no_elements
                tline = fgetl(fid);
                type=sscanf(tline,'%g %g');
                if type(2)==15
                    no_type15=no_type15+1;
                    type15(no_type15,:)=sscanf(tline,'%g %g %g %g %g %g');  % read elements
                elseif type(2)==1
                    no_type1=no_type1+1;
                    type1(no_type1,:)=sscanf(tline,'%g %g %g %g %g %g %g');  % read elements
%                   This is for the old triangle mesh
                elseif type(2)==2
                    no_type2=no_type2+1;
                    elementsTri(no_type2,:)=sscanf(tline,'%g %g %g %g %g %g %g %g');
                 % read elements
                elseif type(2)==3
                    no_type3=no_type3+1;
                    elementsQuad(no_type3,:)=sscanf(tline,'%g %g %g %g %g %g %g %g %g');
                end

                
            end
            tline = fgetl(fid);            % read $EndElements
            if ~ischar(tline) || ~strcmp(tline,'$EndElements'), break, end
        otherwise
            disp('Unknown type encountered...'), break
    end
end
disp('Closing Meshfile...')



nodes = nodes(:,1:3);

[m,~] = size(elementsTri);
if m==0
    elementsTri = nan;
else
elementsTri = [[1:m]', elementsTri(:,(end-2):end)];
end

[m,~] = size(elementsQuad);

if m ==0
    elementsQuad = nan;
else
elementsQuad = [[1:m]', elementsQuad(:,(end-3):end)];
end

BoundaryElements = type1(:,5:7);
fclose('all');