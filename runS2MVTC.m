clear all; 
close all; 
clc;
warning off;

addpath('datasets','Utility','anchor','measure')

ds ={'CCV'};
%ds ={'CCV','Caltech102','NUSWIDEOBJ','AwAfea','cifar10','YoutubeFace_sel'};
for dsi = 1:1:length(ds)
    dataName = ds{dsi}; 
    fprintf('\n Dataset:%s \n',dataName);
    data = dataName;
switch data
    case 'Caltech102'      
         load('Caltech102.mat'); load('Caltech102_anchor.mat');
         Beta =  1;
         Lambda =  10;
         L  = 16       ;

    case 'CCV'      
         load('CCV.mat');load('CCV_anchor.mat');
         Beta =   0.1; 
         Lambda =  0.00001;
         L =18;

    case 'NUSWIDEOBJ'   
         load('NUSWIDEOBJ.mat');   load ('NUSWIDEOBJ_anchor');
         Beta =1;
         Lambda =0.001;
         L= 16;

    case 'AwAfea'      
        load('AwAfea.mat');load('AwAfea_anchor.mat');
         Beta =0.1;
         Lambda =0.03; 
         L = 9; 

    case 'YoutubeFace_sel'      
         load('YoutubeFace_sel.mat');  load YoutubeFace_sel_anchor
         Beta =0.1;  
         Lambda = 0.005;  
         L =19;

    case 'cifar10'   
         load('cifar10.mat'); load cifar10_anchor
         Beta=1E-4;
         Lambda=1E-4;
         L =16;
end
    resultsAll = [];
V = size(X,2);
N= length(Y);

fprintf('The Nonlinear Anchor Embeeding：');
for it = 1:V
fprintf('%d \t',it);
   %Anchor{it} = X{it}(randsample(N,1000),:);
    dist = EuDist2(X{it},Anchor{it},0); 
    sigma = mean(min(dist,[],2).^0.5)*2;
    feaVec = exp(-dist/(2*sigma*sigma));
    X{it} = bsxfun(@minus, feaVec', mean(feaVec',2));
end
clear feaVec dist sigma dist Anchor it

cls_num = length(unique(Y));
total_iterations = length(L) * length(Beta) * length(Lambda);
current_iteration = 0;
for rl=1:length(L)
        for ib=1:length(Beta)
                for il=1:length(Lambda)
%% parameter setting
cls_num = length(unique(Y));n_cluster = numel(unique(Y));
V = length(X); N = size(X{1},2); 
paras.X=X;
paras.beta=Beta(ib);
paras.lambda=Lambda(il);
paras.L=L;
paras.M=cls_num;
paras.N=n_cluster;
% ---------------------------------CLUSTERING------------------------------------------------
        tic;
        pred_label= S2MVTC_function(paras);
        execution_times= toc;
% -------------------------------------- ----------------------------------------------------

        res_cluster = Clustering8Measure(Y, pred_label);
        current_iteration = current_iteration + 1;remainingTime=execution_times*(total_iterations-current_iteration);hours = floor(remainingTime / 3600);minutes = floor(mod(remainingTime, 3600) / 60);seconds = floor(mod(remainingTime, 60));
        fprintf('Iteration: %d/%d\n', current_iteration, total_iterations);
        fprintf('Estimated time remaining: %02d:%02d:%02d\n', hours, minutes, seconds);
        fprintf(['\n Beta: %.5f, Lambda: %.5f , L: %0.0f ----->\tACC:%.4f\t NMI:%.4f\t Purity:%.4f\t F-score:%.4f\t PRE:%.4f\t REC:%.4f\t AR:%.4f\t Entropy:%.4f\t ,Times = %.2f\n '], Beta(ib),Lambda(il), L(rl),res_cluster,execution_times);
       
        
        result = [Beta(ib),Lambda(il), L(rl),res_cluster,execution_times];
        resultsAll = [resultsAll; result];     

                end
        end
end

currentDateTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS'); 
filename = strcat('UnfilteredResults/', dataName, '_', currentDateTime, '_results.mat');
if ~exist('UnfilteredResults', 'dir')
    mkdir('UnfilteredResults');
end
save(filename, 'resultsAll');

clear   X Y Anchor Alpha Beta Gamma Lambda resultsAll L res_cluster execution_times Manchor
end