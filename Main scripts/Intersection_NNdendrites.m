%Intersection_NNdendrites
%Use: To obtain the nearest neighbor from dendritic segments inside the
%intersected volume. Script adapted to use with neuron positioned in their
%region.
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

Intersection_NNdendrites_T = [];          %%%Name of the output matrix that contains data analyses.
dataNN = [];

folder = '';            %%%Folder containing xyz_position file
folder2 = '';           %%%Folder containing Linepoints_CHP files
folder3 = '';           %%%Folder containing Linepoints_dendrites files

cd(folder);

load xyz_position.mat

nCorr = nchoosek(1:15,2);     %%%Consider the total number of neurons. For our VTA DA sample, n = 15.
sCorr = size(nCorr,1);

for c = 1:sCorr
    tic           %%%OPTIONAL: to measure run time of each loop iteration (also active toc, at the end of the loop)

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
    
    %%%Load Linepoints_dendrites files
    
    cd(folder3);
    
    neu_1 = ['Neu',num2str(nCorr(c,1))];
    nameNeu = ['Neu',num2str(nCorr(c,1)),'.mat'];
    
    load(nameNeu)
    
    Neu_1 = eval(neu_1);
    
    vars = {strcat('Neu',num2str(nCorr(c,1)))};
    clear(vars{:})
     
    neu_2 = ['Neu',num2str(nCorr(c,2))];
    nameNeu = ['Neu',num2str(nCorr(c,2)),'.mat'];
    
    load(nameNeu)
    
    Neu_2 = eval(neu_2);
    
    vars = {strcat('Neu',num2str(nCorr(c,1)))};
    clear(vars{:})
    
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
        inter_denlen_1 = 0;
        inter_denlen_2 = 0;
        inter_denlen = 0;
        NN_1 = NaN;
        NN_mean = NaN;
        NN_median = NaN;
        
        NN_perc = NaN;

    else
        [K3,v3] = convhulln(inter,{'QJ'});
        inter_den1 = inhull(Neu_1,inter,K3);
        inter_den2 = inhull(Neu_2,inter,K3);
        
        Neu_1_inter = Neu_1(inter_den1,:);
        Neu_2_inter = Neu_2(inter_den2,:);
        
        inter_den = [Neu_1_inter;Neu_2_inter];
        inter_denlen_1 = size(Neu_1_inter,1)*0.001;
        inter_denlen_2 = size(Neu_2_inter,1)*0.001;
        inter_denlen = size(inter_den,1)*0.001;

        %%%NN

        Neu_1_inter_v2 = unique(round(Neu_1_inter,2),'rows');
        Neu_2_inter_v2 = unique(round(Neu_2_inter,2),'rows');
        
        [~,D]=knnsearch(Neu_1_inter_v2,Neu_2_inter_v2,'K',1);
        
        if isempty(D) == 1
            NN_1 = NaN;
        else
            D = D(D~=0);
            NN_1 = min(D);
        end
    end

    d = norm(xyz_position(nCorr(c,1),:)-xyz_position(nCorr(c,2),:));
    
    denlen_1 = size(Neu_1,1)*0.001;
    denlen_2 = size(Neu_2,1)*0.001;
    
    Intersection_NNdendrites_T = [Intersection_NNdendrites_T; nCorr(c,1) nCorr(c,2) ...
        xyz_position(nCorr(c,1),:) xyz_position(nCorr(c,2),:) d v1 v2 v3 v3/(v1+v2-v3)*100 denlen_1 denlen_2 inter_denlen NN_1];
   
    cd(folder);
    save('Intersection_NNdendrites_T','Intersection_NNdendrites_T')

    clear Hull_1 Hull_2 inter
        
    toc
    c
end

