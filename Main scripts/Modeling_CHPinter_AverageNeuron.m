%Modeling_CHPinter_AverageNeuron
%Use: Calculate the intersected volume of pairs of neurons convex hulls, in
%a set of distances between cell bodies and axes. Use a determined
%convex hull as the average neuron, according to a defined criteria (for
%example, average dendritic length, average isotropy index, etc...)
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

Modeling_CHPinter_AverageNeuron_T = [];

for nN = 1:15
    cd(folder2);
    
    hull_1 = ['Hull',num2str(nN)];
    namehull = ['Hull',num2str(nN),'.mat'];
    
    load(namehull)  
        
    Hull_1 = eval(hull_1);
    
    clear(hull_1)
        
    [~,v] = convhulln(Hull_1);
    
    namefinal2 = ['v',num2str(nN)];
    str2 = [namefinal2,'= v;'];
    eval(str2);
    
    clear Hull_1
      
end

clear v

cd(folder2);               

load 'Replace with the name of Normalized_Shape file'.mat
Hull_1m = 'Replace with the name of Normalized_Shape file';

%%%Substract neuron original position (if Linepoints convex hull was
%%%generated with the neurons position in their region)
% Hull_1m = Hull_1m - xyz_position('row corresponding to the neuron',:);

[Km,vm] = convhulln(Hull_1m);

timeT = [];
ttime = 0;
p = 3;              %%%Number of axes tested.
dN = 0:10:400;      %%%Sample of distances of cell bodies to test (µm).
dT = length(dN);
c1c2 = 210;         %%%Number of pairs test per distance.

for pos = 1:3
for d = dN
for c1 = 1:15 
for c2 = 1:15

    tstart = tic;
            
    if c1 ==c2
        continue
    end
       
    vol_1 = ['v',num2str(c1)];
    vol_2 = ['v',num2str(c2)];
                
    Hull_1 = Hull_1m;
    Hull_2 = Hull_1m;
    
    v1a = eval(vol_1);
    v2a = eval(vol_2);
    
    K1a = Km;
    K2a = Km;
        
    Hull_1_weigth = nthroot(v1a/vm,3);             
    Hull_2_weigth = nthroot(v2a/vm,3);
    
    Hull_1 = Hull_1.*Hull_1_weigth;                 
    Hull_2 = Hull_2.*Hull_2_weigth;
    Hull_2(:,pos) = Hull_2(:,pos) + d;
    
    in1 = inhull(Hull_2,Hull_1,K1a);
    Hull2_Inter = Hull_2(in1,:);
    
    in2 = inhull(Hull_1,Hull_2,K2a);
    Hull1_Inter = Hull_1(in2,:);

    inter = [Hull1_Inter;Hull2_Inter];
    
    if isempty(inter)==1 || size(inter,1) < 4
        v3a = 0;
    else
        [~,v3a] = convhulln(inter,{'QJ'});
    end
    
    Modeling_CHPinter_AverageNeuron_T = [Modeling_CHPinter_AverageNeuron_T; pos d c1 c2 v1a v2a v3a v3a/(v1a+v2a-v3a)*100];

    clear Hull_1 Hull_2 inter
    
    cd(folder);
    
    save('Modeling_CHPinter_AverageNeuron_T','Modeling_CHPinter_AverageNeuron_T')
    
    %%%OPTINAL: to track run time of the script
    
    ctime = toc(tstart);
    ttime = ttime+ctime;
    timeT = [timeT;ctime];

    [pos d c1 c2 ctime ttime (mean(timeT)*p*dT*c1c2-ttime)]
end
end
end
end
