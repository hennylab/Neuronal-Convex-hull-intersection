%Modeling_CHPinter_Ogdata
%Use: Calculate the intersected volume of pairs of neurons convex hulls, in
%a set of distances between cell bodies and axes.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22
close all
clear
clc

folder = '';            %%%Folder containing xyz_position file
folder2 = '';           %%%Folder containing Linepoints_CHP files

cd(folder);

% load xyz_position.mat    %%%Load xyz_position if neuron position in their region was used in the ConvexHull2LinePoint script

Modeling_CHPinter_Ogdata_T = [];

timeT = [];
ttime = 0;
p = 3;              %%%Number of axes tested.
dN = 0:10:400;      %%%Sample of distances of cell bodies to test (µm).
dT = length(dN);
c1c2 = 210;         %%%Number of pairs test per distance.

for pos = 1:3
for d = dN
for c = 1:15
    for c2 = 1:15
        
        if c==c2
            continue
        end
        
        tstart = tic;  
        
    cd(folder2);

    %%%Load Linepoints_CHP files

    hull_1 = ['Hull',num2str(c)];
    namehull = ['Hull',num2str(c),'.mat'];
    
    load(namehull)
    
    Hull_1 = eval(hull_1);
    
    clear(hull_1)
    
    hull_2 = ['Hull',num2str(c2)];
    namehull = ['Hull',num2str(c2),'.mat'];
    
    load(namehull)
    
    Hull_2 = eval(hull_2);
    
    clear(hull_2)
    
    %%%Calculate intersection
    
%     %%%Substract neuron original position (if Linepoints convex hull was
%     %%%generated with the neurons position in their region)
%     
%     Hull_1 = Hull_1(:,1:3)-xyz_position(c,:);           
%     Hull_2 = Hull_2(:,1:3)-xyz_position(c2,:);          
    
    Hull_2(:,pos) = Hull_2(:,pos) + d;
  
    [K1,v1] = convhulln(Hull_1);
    [K2,v2] = convhulln(Hull_2);
    
    inter1 = inhull(Hull_2,Hull_1,K1);
    in2 = double(inter1);
    clear in    
    inT = find(in2==1);
    Hull2_Inter = Hull_2(inT,:);
    
    inter2 = inhull(Hull_1,Hull_2,K2);
    in2 = double(inter2);
    clear in    
    inT = find(in2==1);
    Hull1_Inter = Hull_1(inT,:);
    
    inter = [Hull1_Inter;Hull2_Inter];
    
    if isempty(inter)==1 || size(inter,1) < 4
        v3 = 0;
    else
        [~,v3] = convhulln(inter,{'QJ'});
    end
      
    Modeling_CHPinter_Ogdata_T = [Modeling_CHPinter_Ogdata_T; c c2 pos d v1 v2 v3 v3/(v1+v2-v3)*100]; 

    clear Hull_1 Hull_2
    
    cd(folder);

    save('Modeling_CHPinter_Ogdata_T','Modeling_CHPinter_Ogdata_T')
    
    %%%OPTINAL: to track run time of the script
    
    ctime = toc(tstart);
    ttime = ttime+ctime;
    timeT = [timeT;ctime];

    [pos d c c2 ctime ttime (mean(timeT)*p*dT*c1c2-ttime)]

    end
end
end
end
