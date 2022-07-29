%Dendrites2linepoints
%Use: Transform dendritic segments to points. Used as input for convex
%hull intersection analysis with inhull function, to identify dendrites inside intersection.
%Use with .txt files containing Segment Points - Dendrites data obtained from Neurolucida Explorer.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

folder = '';    %%%Folder containing segment points - dendrites (.txt) files + xyz_position (optional)

cd(folder);
filelist = dir('*.txt');

load xyz_position.mat

for i=1:length(filelist)   
 
    cd(folder);
    filename = filelist(i).name;

    folder2 = strcat(folder, filename);              
    formatSpec = '%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    delimiter = '\t';
    startRow = 2;
    fileID = fopen(folder2,'r');                         
    c = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    
    %%%%To extract dendritic segments starts and endings:
    
    %%%Without correcting neuron position:
    
    xyz = [c{3} c{4} c{5} c{6} c{7} c{8}];
    
%     %%%Correcting neuron position:
% 
%     xyz = [c{3}+xyz_position(k,1) c{4}+xyz_position(k,2) c{5}+xyz_position(k,3) ...
%             c{6}+xyz_position(k,1) c{7}+xyz_position(k,2) c{8}+xyz_position(k,3)]; 
    
    xyz(:,7) = vecnorm(xyz(:,4:6)-xyz(:,1:3),2,2);
      
    xyz = unique(xyz,'rows');
    
    xyz_2 = [];
        
    for dt = 1:size(xyz,1)
        
        if xyz(dt,7) >= 0.001
            n = xyz(dt,7)./0.001;
        else
            n = 1;
        end
        t = linspace(0,1,n+1)';
       
        xyz_2 = [xyz_2;(1-t)*xyz(dt,1:3) + t*xyz(dt,4:6);xyz(dt,4:6)];  
    end
    
    xyz_2 = unique(xyz_2,'rows');
    
    namefinal = ['Neu',num2str(i)];
    str2 = [namefinal,'= xyz_2;'];
    
    eval(str2);
    
    cd(folder);
    
    mkdir Linepoints_dendrites;
    sDir3 = strcat(folder,'Linepoints_dendrites');
    cd(sDir3);
    
    save(namefinal,namefinal)
    close all
    
end

