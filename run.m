clear all
addpath(genpath('helper_functions'));
addpath(genpath('Dependencies')); 

%-------------------------------------------------------------------

diary off; diary on;
fprintf('\nSTART TIME:    %s\n\n', datestr(now));

%-------------------------------------------------------------------

global predictionMethod gridSearchMode

gridSearchMode = 0;   % grid search mode?

predictionMethods = {'mgrnnm'};

warning off

%-------------------------------------------------------------------

global m n Sd St ds cv_setting

% The location of the folder that contains the data
path='data\';

% the different datasets
datasets={'e','ic','gpcr','nr'}%,'movielens_100k','metabolic'};

% CV parameters
m = 5;  % number of n-fold experiments (repetitions)
n = 10;%5;%10; % the 'n' in "n-fold experiment"

%-------------------------------------------------------------------

disp(['gridSearchMode = ' num2str(gridSearchMode)])
disp(' ')

similarity_types={'correlation','correlation','correlation','correlation','cosine','cosine','cosine','cosine','jaccard','jaccard','jaccard','jaccard','hamming'}%'hamming','jaccard',
%{'cosine','correlation','hamming','jaccard'}

    %
%%
   % 
    
diary off; diary on;

%-------------------------------------------------------------------
for p1=1:length(predictionMethods)
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    predictionMethod = predictionMethods{p1};

    % loop over all cross validation settings
for cvs=[1 2 3]
        disp('===========================================');
        disp(['Prediction method = ' predictionMethod])
        cv_setting = ['S' int2str(cvs)];
        switch cv_setting
            case 'S1', disp('CV Setting Used: S1 - PAIR');
            case 'S2', disp('CV Setting Used: S2 - DRUG');
            case 'S3', disp('CV Setting Used: S3 - TARGET');
        end
        disp(' ')
  
        % run chosen selection method and output CV results
        auprlist=[]; auprstdlist=[]; auclist=[]; aucstdlist=[];
        for ds=[4 3 2 1 ]
            getParameters(predictionMethod, cv_setting, ds);
            disp('-----------------------');

            fprintf('\nData Set: %s\n', datasets{ds});

            % LOAD DATA
            [Y,Sd,St,~,~]=getdata(datasets{ds},similarity_types);

            % CV experiment
            tic
            [aupr,aupr_std,auc,auc_std]=crossValidation(Y')
            auprlist=[auprlist aupr]; auprstdlist=[auprstdlist aupr_std]; auclist=[auclist auc]; aucstdlist=[aucstdlist auc_std];
            disp(' ')
            toc

            disp('-----------------------');
            diary off; diary on;
        end
        
        disp('===========================================');
        diary off; diary on;
%    save(['gs_cvsetting/' num2str(m) 'runsOf' num2str(n) 'foldcv_' char(predictionMethods(p)) '_S' num2str(cvs) '.mat' ],'auprlist','auprstdlist','auclist','aucstdlist')
    

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    diary off; diary on;
end
end
diary off;

%-------------------------------------------------------------------