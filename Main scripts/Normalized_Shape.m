%Normalized shape
%Use: To obtain the normalized shape from a sample of convex hulls
%(obtained from the CHP_normalization script).
%Version 1.0
%Author: Rafael Gatica
%Email: rigatica@uc.cl
%Release: 29/07/22

close all
clear
clc

folder = '\';            %%%Folder containing convex hull data obtained from CHP_normalization.

nHulls = 15;        %%%Number of convex hulls to analyze. For our VTA DA neurons sample, n = 15.

c1 = [];
c2 = [];
c3 = [];

for nN = 1:nHulls
                
    cd(folder);      
    
    hull_1 = ['Hull',num2str(nN)];
    namehull = ['Hull',num2str(nN),'.mat'];
    
    load(namehull)
    
    Hull_1 = eval(hull_1);
    
    c1 = [c1 Hull_1(:,2)];
    c2 = [c2 Hull_1(:,3)];
    c3 = [c3 Hull_1(:,4)];
    
    clear Hull_1
   
end

normalized_shape = [mean(c1,2) mean(c2,2) mean(c3,2)];

lF = 20;
t = linspace(0,1,lF)';

[K1,v1] = convhulln(normalized_shape);

Hull_1 = [];
    
for np1 = 1: size(K1,1)
    for np2 = 1:3
        if np2 ==1          
            P1 = normalized_shape(K1(np1,3),:);
            P2 = normalized_shape(K1(np1,1),:);
        else
            P1 = normalized_shape(K1(np1,np2-1),:);
            P2 = normalized_shape(K1(np1,np2),:);
        end
        Hull_1 = [Hull_1;(1-t)*P1 + t*P2];
    end
end

Hull_1 = unique(Hull_1,'rows');       

namefinal = 'Normalized_Shape_Xregion';

str2 = [namefinal,'= Hull_1;'];

eval(str2);

cd(folder);

save(namefinal,namefinal)
