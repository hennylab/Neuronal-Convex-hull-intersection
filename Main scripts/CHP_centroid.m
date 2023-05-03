%Intersection_realpositions_NormCellB
%Use: Calculate the coordinates for centroid of each neuron convex hull.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

folder = '';        %%%Folder that contains .txt files. OPTIONAL: also, add in this folder the matrix with cell body position

cd(folder);
filelist = dir('*.txt');

%%%%If you will use neurons position in their region, load a .mat matrix with cell body
%%%%coordinates for each neuron, listed in the same order as 'filelist'.

% load xyz_position.mat     %%%Uncomment if you will use neuron position

centroidT = [];

for i=1:length(filelist)        %%%leer los archivos    
    
    cd(folder);
    sRootName = filelist(i).name;

    folder2 = strcat(folder, sRootName);               %%%concatenar string
    formatSpec = '%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    delimiter = '\t';
    startRow = 2;
    fileID = fopen(folder2,'r');                          %%%lo deja guardado en un lugar de la memoria
    c = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');       %%%leer archivo guardado en memoria
    fclose(fileID);
    
    %%%%To extract dendritic segments starts and endings:
    
    %%%Without correcting neuron position:
    
    xyz = [c{3} c{4} c{5};c{6} c{7} c{8}];
    
%     %%%Correcting neuron position:
% 
%     xyz = [c{3}+xyz_position(i,1) c{4}+xyz_position(i,2) c{5}+xyz_position(i,3);...
%             c{6}+xyz_position(i,1) c{7}+xyz_position(i,2) c{8}+xyz_position(i,3)];    

    xyz = unique(xyz,'rows');

    Hull_centroid = centroidv2(xyz);
    
    centroidT = [centroidT; Hull_centroid ];

end

save('centroidT','centroidT')
