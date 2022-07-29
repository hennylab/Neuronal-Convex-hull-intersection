%ConvexHull2LinePoints_txt
%Use: Transform convex hulls polyhedrons edges to points. Used as input for convex
%hull intersection analysis with inhull function. 
%Use with .txt files containing Segment Points - Dendrites data obtained from Neurolucida Explorer.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

folder = '';        %%%Folder that contains .txt files. OPTIONAL: also, add in this folder matrix with cell body position

cd(folder);
fileslist = dir('*.txt');

%%%%If you will use neurons position in their region, load a .mat matrix with cell body
%%%%coordinates for each neuron, listed in the same order as 'filelist'.

% load xyz_position.mat     %%%Uncomment if you will use neuron position

n = 1000;      %%%Number of points per convex hull edge. If you require to create convex hulls for shape normalization, n = 10000.

for i=1:length(fileslist)          
    
    cd(folder);
    filename = fileslist(i).name;

    folder2 = strcat(folder, filename);               
    formatSpec = '%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    delimiter = '\t';
    startRow = 2;
    fileID = fopen(folder2,'r');                         
    c = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    
    %%%%To extract dendritic segments starts and endings:
    
    %%%Without correcting neuron position:
    
    xyz = [c{3} c{4} c{5};c{6} c{7} c{8}];
    
%     %%%Correcting neuron position:
% 
%     xyz = [c{3}+xyz_position(i,1) c{4}+xyz_position(i,2) c{5}+xyz_position(i,3);...
%             c{6}+xyz_position(i,1) c{7}+xyz_position(i,2) c{8}+xyz_position(i,3)];    

    xyz = unique(xyz,'rows');

   [K1,v1] = convhulln(xyz);
    
   Hull_1 = [];
    
    for np1 = 1: size(K1,1)
        for np2 = 1:3
            if np2 ==1          
                P1 = xyz(K1(np1,3),:);
                P2 = xyz(K1(np1,1),:);
            else
                P1 = xyz(K1(np1,np2-1),:);
                P2 = xyz(K1(np1,np2),:);
            end
            
            t = linspace(0,1,n)';
            Hull_1 = [Hull_1;(1-t)*P1 + t*P2];
        end
    end
    
    Hull_1 = unique(Hull_1,'rows');
            
    namefinal = ['Hull',num2str(i)];
    str2 = [namefinal,'= Hull_1;'];
    
    eval(str2);
    
    cd(folder);
    
    mkdir Linepoints_CHP;
    sDir3 = strcat(folder,'Linepoints_CHP');
    cd(sDir3);
    
    save(namefinal,namefinal)
    close all
    
end

