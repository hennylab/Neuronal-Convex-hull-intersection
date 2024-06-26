%Intersection_realpositions_NormCellB
%Use: To obtain the intersected volume between pairs of neurons convex
%hulls, in their real positions in the region. Normalized cell body is applied.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

Intersection_realpositions_NormCellB_T = [];          %%%Name of the output matrix that contains data analyses.
dataNN = [];

folder = '';            %%%Folder containing xyz_position and centroidT files
folder2 = '';           %%%Folder containing Linepoints_CHP files

cd(folder);

load xyz_position.mat
load centroidT.mat

nCorr = nchoosek(1:15,2);     %%%Consider the total number of neurons. For our VTA DA sample, n = 15.
sCorr = size(nCorr,1);

for c = 1:sCorr
    tic           %%%OPTIONAL: to measure run time of each loop iteration
%     (also active toc, at the end of the loop)

    %%%Load Linepoints_CHP files

    cd(folder2);
    
    hull_1 = ['Hull',num2str(nCorr(c,1))];
    namehull = ['Hull',num2str(nCorr(c,1)),'.mat'];
    
    load(namehull)
    
    Hull_1 = eval(hull_1);
    
    vars = {strcat('Hull',num2str(nCorr(c,1)))};
    clear(vars{:})
    
    hull_2 = ['Hull',num2str(nCorr(c,2))];
    namehull = ['Hull',num2str(nCorr(c,2)),'.mat'];
    
    load(namehull)
    
    Hull_2 = eval(hull_2);
    
    vars = {strcat('Hull',num2str(nCorr(c,2)))};
    clear(vars{:})
        
    %%%Normalize cell body position
    
    Hull_1 = Hull_1 - centroidT(nCorr(c,1),:);
    Hull_2 = Hull_2 - centroidT(nCorr(c,2),:);
    
    Hull_1 = Hull_1 + xyz_position(nCorr(c,1),:);
    Hull_2 = Hull_2 + xyz_position(nCorr(c,2),:);
    
    %%%Calculate intersection
              
    [K1,v1] = convhulln(Hull_1);
    [K2,v2] = convhulln(Hull_2);

    inter1 = inhull(Hull_2,Hull_1,K1);
    Hull2_Inter = Hull_2(inter1,:);
    inter2 = inhull(Hull_1,Hull_2,K2);
    Hull1_Inter = Hull_1(inter2,:);
    
    inter = [Hull1_Inter;Hull2_Inter];
    
    if isempty(inter)==1 || size(inter,1) < 4
        v3 = 0;
    else
        [~,v3] = convhulln(inter,{'QJ'});
    end

    d = norm(xyz_position(nCorr(c,1),:)-xyz_position(nCorr(c,2),:));
        
    Intersection_realpositions_NormCellB_T = [Intersection_realpositions_NormCellB_T; nCorr(c,1) nCorr(c,2) ...
        xyz_position(nCorr(c,1),:) xyz_position(nCorr(c,2),:) d v1 v2 v3 v3/(v1+v2-v3)*100];
   
    cd(folder);
    save('Intersection_realpositions_NormCellB_T','Intersection_realpositions_NormCellB_T')

    clear Hull_1 Hull_2 inter
    
    c
    toc
        
end

