%CHP_normalization
%Use: Transform linepoints convex hulls to a normalized point
%representation, with an equal number of points per convex hull. Used for
%normalized shape condition.
%Output: Matrix with 3D coordinates, for each tested angle.
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

load xyz_position.mat   

rF = 1;
lF = 10000;
t = linspace(0,1,lF)';
nHulls = 15;        %%%Number of convex hulls to analyze. For our VTA DA neurons sample, n = 15.

for nN = 1:nHulls
    tic

    cd(folder2);
    filelist = dir('*.mat');
      
    hull_1 = ['Hull',num2str(nN)];
    namehull = ['Hull',num2str(nN),'.mat'];
    
    load(namehull)
    
    Hull_1 = eval(hull_1);
    
    vars = {strcat('Hull',num2str(nN))};
    clear(vars{:})
         
    Hull_1 = unique(round(Hull_1-xyz_position(nN,:),rF),'rows');
    
    d_T = [];
    
    for axisn = 1:3

    if axisn == 1
        e1 = 1;
        e2 = 2;
    elseif axisn == 2
        e1 = 3;
        e2 = 2;
   elseif axisn == 3
        e1 = 1;
        e2 = 3;
    end

    [k1,~] = convhulln(Hull_1(:,[e1 e2]));

    Hull_1b = unique(Hull_1(k1,[e1 e2]),'rows');

    x1 = 1;
    y1 = 0;
    
    x2 = Hull_1b(:,1);
    y2 = Hull_1b(:,2);

    angleHull = [rad2deg(atan2(x1*y2-y1*x2,x1*x2+y1*y2)) Hull_1b(:,1:2)];

    angleHull = sortrows(angleHull);
    
    Hull1b_Int = [];
    
    if size(Hull_1b,1)==1 || isempty(k1)==1
        error('Not enough points for computation')
    else
        for nf = 1:1:size(angleHull,1)
            if nf ==1          
                P1 = angleHull(size(angleHull,1),2:3);
                P2 = angleHull(1,2:3);
            else
                P1 = angleHull(nf-1,2:3);
                P2 = angleHull(nf,2:3);
            end
            Hull1b_Int = [Hull1b_Int;(1-t)*P1 + t*P2];
        end
       
    end 
    
    Hull1b_Int = unique(round(Hull1b_Int,rF),'rows');

    x1 = 1;
    y1 = 0;

    x2 = Hull1b_Int(:,1);
    y2 = Hull1b_Int(:,2);

    angleHull = [rad2deg(atan2(x1*y2-y1*x2,x1*x2+y1*y2)) Hull1b_Int(:,1:2)];

    angMLDVT = []; 
    for ang = 0:359

     [~,index] = min(abs(angleHull(:,1) - ang));

     if isempty(index)==1
         error('Cannot find an angle')
     end

     angMLDVT = [angMLDVT; ang angleHull(index,2:3)];

    end

    if axisn == 2
        angMLDVT = angMLDVT(:,[1 3 2]);
    end

    k = 0;

        for aDV = 0:179

        if axisn == 2
            e1 = 2;
            e2 = 3;
        end

        k = k + 1;

        f1 = find(angMLDVT(:,1)==aDV);      
        f2 = find(angMLDVT(:,1)==180+aDV);

        pt1 = angMLDVT(f1,2:3);
        pt2 = angMLDVT(f2,2:3);

        angline = unique(round([(1-t)*pt1 + t*pt2],rF),'rows');

        [~,fa,~] = intersect(Hull_1(:,[e1 e2]),angline(:,1:2),'rows','stable');

        Hull_1a = Hull_1(fa,1:3);

        x1 = 1;
        y1 = 0;

        if axisn == 1 && aDV == 90
            c1 = 3; c2 = 2;
        elseif axisn == 1 
            c1 = 1; c2 = 3;
        elseif axisn == 2 && aDV == 90
            c1 = 1; c2 = 2;
        elseif axisn == 2
            c1 = 3; c2 = 1;
        elseif axisn == 3 && aDV == 90
            c1 = 2; c2 = 3;
        elseif axisn == 3
            c1 = 1; c2 = 2;
        end   

        x2 = Hull_1a(:,c1);
        y2 = Hull_1a(:,c2);

        angleT = [rad2deg(atan2(x1*y2-y1*x2,x1*x2+y1*y2)) Hull_1a(:,1:3)];

        angleT = sortrows(angleT);

        Hull1a_Int = [];
        if size(Hull_1a,1)==1 
            Hull1a_Int = Hull_1a(:,1:3);
        else
            for nf = 1:1:size(angleT,1)
                if nf ==1          
                    P1 = angleT(size(angleT,1),2:4);
                    P2 = angleT(1,2:4);
                else
                    P1 = angleT(nf-1,2:4);
                    P2 = angleT(nf,2:4);
                end
                Hull1a_Int = [Hull1a_Int;(1-t)*P1 + t*P2];
            end
        end

        Hull1a_Int = unique(round(Hull1a_Int,rF),'rows');

        if axisn == 1 && aDV == 90
            c3 = 2;
        elseif axisn == 1 
            c3 = 3;
        elseif axisn == 2 && aDV == 90
            c3 = 2;
        elseif axisn == 2
            c3 = 1;
        elseif axisn == 3
            c3 = 2;
        end   

        vector = find(Hull1a_Int(:,c3)==0);
        vectorx = max(Hull1a_Int(vector,1));

        if vectorx < 0 && size(Hull_1a,1)>2
            error('X coordinate is lower than 0')
        end

        x1 = 1;
        y1 = 0;

        x2 = Hull1a_Int(:,c1);
        y2 = Hull1a_Int(:,c2);

        angleT = [rad2deg(atan2(x1*y2-y1*x2,x1*x2+y1*y2)) Hull1a_Int(:,1:3)];

        for ang = 0:1:360
         [~,index] = min(abs(angleT(:,1) - ang));

         d_T = [d_T; ang angleT(index,2:4)];
        end

        [nN axisn k toc]
        end
    end
    
    namefinal = ['Hull',num2str(nN)];

    str2 = [namefinal,'= d_T;'];
    
    eval(str2);
    
    cd(folder);
    
    mkdir Normalization_CHP;
    folder3 = strcat(folder,'Normalization_CHP');
    cd(folder3);
    
    save(namefinal,namefinal)

    clear Hull_1 d_T
end