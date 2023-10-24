% In this script, we provide an example of the information-theoretic
% analyses we use in the associated study
% 
% To keep computations simple, we analyze a single pair of LFP channels 
% recorded simultaneously from M1 and DLS during a day of reach-to-grasp 
% training. We include one reach feature, the maximum reaching velocity 
% during each trial. In the study we compute a null hypothesis for each 
% information quantity via trial-shuffling which we have not included for
% simplicity.
% 
% The script will compute these information quantities:
% 
% 1) Mutual Information (MI) between a neural signal and the included 
% reach feature, computed across trials for each neural signal seperately 
% at each time bin.
% 
% 2) Shared or redundant information that a pair of neural signals carry 
% about the included reach feature. Shared information is computed using 
% the partial information decomposition (PID) framework. For each time 
% point, shared information is computed for a range of temporal delays 
% between the neural signals, to explore the possibility of one signal 
% carrying movement information that is present in the other signal at a 
% later time, and vice versa.
% 
% 3) The Directed Information (DI) between a pair of neural signals. DI is 
% computed at each time point, across temporal delays, to capture the the 
% total amount of information transferred between two neural signals.
% 
% 4) The Feature-specific Information Transfer (FIT) between a pair of
% neural signals. FIT is computed at each time point, across temporal 
% delays, capturing the movement-specific information transferred between 
% two neural signals.
%
%% SETUP

% Set MATLAB to current folder (containing this script)

    cd(fileparts(matlab.desktop.editor.getActiveFilename))

% add path to code for information theoretic analysis

    path_to_info_analysis_code = fullfile(pwd,'Code');
    addpath(genpath(path_to_info_analysis_code));

% compile code
    mex('-g', '-largeArrayDims', '-outdir', path_to_info_analysis_code,...
        fullfile(pwd,'Code','MI','direct_method.c'));

%% LOAD NEURAL SIGNALS AND MOVEMENT FEATURES
    
% load neural signals with the format = channels/units x time bins x trials
    
    load(fullfile('ExampleData','NeuralSignals.mat'));

% load movement features with one value per trial and format = 1 x trials;

    load(fullfile('ExampleData','ReachFeature.mat'));

%% PREPROCESS NEURAL SIGNALS AND MOVEMENT FEATURES 

% temporal rebinning

    inputData = NeuralSignals; % input data
    signalLabel = {'M1', 'DLS'}; % label/area of origin for each neural signal provided, used for plotting
    newTimeBin = 10; % new time bin in units of the input time bin - input data must be divisable by newTimeBin  
    method = 'mean'; % options are 'sum', 'mean', and 'movmean'
    timeJump = 10; % input only used for moving mean method option, 'sum' and 'mean' use nonoverlapping bins
    NeuralSignals_binned = temporal_rebinning(inputData,newTimeBin,'mean',timeJump);

% rebin neural signals and reach values

    NeuralSignalBin = 5; % number of bins for neural signals   
    for numSignals = 1:size(NeuralSignals_binned,1)
        for numTrials = 1:size(NeuralSignals_binned,3) 
            NeuralSignals_binned(numSignals,:,numTrials) = eqpop(NeuralSignals_binned(numSignals,:,numTrials),NeuralSignalBin); % rebin each trial's neural signals with equipopulated bins
        end
    end
    
    ReachFeatureBin = 3; % number of bins for reach features    
    ReachFeature_binned = eqpop(ReachFeature,ReachFeatureBin); % rebin reach features across trials with equipopulated bins

% plot raw and and binned data

    for numSignals = 1:size(NeuralSignals,1)
        figure;
            subplot(1,2,1); hold on;
                for numTrials = 1:10
                    plot((5*numTrials)+zscore(NeuralSignals(numSignals,:,numTrials)),'k');
                end
                xlim([1 size(NeuralSignals,2)])
                yticklabels('')
                xlabel('time')
                ylabel('single trials')
                title([signalLabel{numSignals}, ' neural signals'])
            subplot(1,2,2); hold on;
                for numTrials = 1:10
                    plot((5*numTrials)+zscore(NeuralSignals_binned(numSignals,:,numTrials)),'k');
                end
                xlim([1 size(NeuralSignals_binned,2)])
                yticklabels('')
                xlabel('time')
                ylabel('single trials')
                title([signalLabel{numSignals}, ' binned neural signals'])
    end
    
%% COMPUTE MUTUAL INFORMATION

% set parameters for mutual information 

    MI_windowSize = 5; % number of consecutive time bins used to compute mutual information
    MI_timeJump = 2; % time jump used when computing mutual information over time
    MI_params.verbose = 0;
    
% compute mutual information

    MI_overTime = [];
    for numSignals = 1:size(NeuralSignals_binned,1)
        MI_tmpSignal = [];
        sTime = 1;
        for timeStep = 1:size(NeuralSignals_binned,2)/MI_timeJump
            eTime = sTime + MI_windowSize - 1;
            if eTime>size(NeuralSignals_binned,2)
                sTime = sTime + MI_timeJump;
            else
                X = NeuralSignals_binned(numSignals,sTime:eTime,:);
                Y = ReachFeature_binned;
                Y = repmat(Y,size(X,2),1);
                X = X(:);
                Y = Y(:);
                MI = information(X',Y',MI_params,{'I'});
                MI_tmpSignal = [MI_tmpSignal MI{1}(1)];
                sTime = sTime + MI_timeJump;
            end
        end
        MI_overTime = [MI_overTime; MI_tmpSignal];
    end
    
% plot mutual information

    figure;
        hold on;
        for numSignals = 1:size(MI_overTime,1)
            plot(MI_overTime(numSignals,:))
        end
        legend(signalLabel)
        title('Mutual Information')
        xlabel('time')
        ylabel('bits')
    
%% COMPUTE PARTIAL INFORMATION DECOMPOSITION

% set parameters for partial information decomposition (PID)

    PID_windowSize = 10; % number of consecutive time bins used for partial information decomposition
    PID_timeJump = 5; % time jump used for partial information decomposition over time
    PID_delay = -50:50; % delays to consider as temporal lag between neural signals
    PID_params.method = 'dr'; % estimation method: dr = Direct method
    PID_params.bias = 'naive'; % bias correction procedure: naive = baised naive estimates
    PID_params.verbose = 0;
    PIDpairs = {[1 2]}; % neural pairs to consider when computing partial information decomposition
    
% compute PID

    PIDovertime = cell(1,length(PIDpairs));
    for numPairs = 1:length(PIDpairs)
        PID_tmpSignal = [];
        sTime = 1;
        for timeStep = 1:size(NeuralSignals_binned,2)/PID_timeJump
            eTime = sTime + PID_windowSize - 1;
            if sTime+min(PID_delay)<0 || eTime+max(PID_delay)>size(NeuralSignals_binned,2)
                sTime = sTime + PID_timeJump;
            else
                PIDdelays = [];    
                for delay = PID_delay
                    X1 = NeuralSignals_binned(PIDpairs{numPairs}(1),sTime:eTime,:);
                    X2 = NeuralSignals_binned(PIDpairs{numPairs}(2),sTime+delay:eTime+delay,:);
                    Y = ReachFeature_binned;
                    Y = repmat(Y,size(X1,2),1);
                    X1 = X1(:);
                    X2 = X2(:);
                    Y = Y(:);
                    I = PID(Y',X1',X2',PID_params);
                    PIDdelays = [PIDdelays I.shared];
                end
                PID_tmpSignal = [PID_tmpSignal PIDdelays'];
                sTime = sTime + PID_timeJump;
            end
        end
        PIDovertime{numPairs} = PID_tmpSignal;
    end
    
% plot PID

    for numPairs = 1:length(PIDpairs)
        figure; 
            hold on;
            imagesc(PIDovertime{numPairs});
            plot([0.5 size(PIDovertime{numPairs},2)+.5],[round(size(PIDovertime{numPairs},1)/2) round(size(PIDovertime{numPairs},1)/2)],'color','r','LineWidth',2)
            title(['PID shared information'])
            colorbar;
            xlim([.5 size(PIDovertime{numPairs},2)+.5])
            yticks([1:5:size(PIDovertime{1,numPairs},1)])
            yticklabels([PID_delay(1):5:PID_delay(end)])
            xlabel('time')
            ylim([.5 size(PIDovertime{numPairs},1)+.5])
            ylabel('delay (time bins)')
            h1=text(2,65,[signalLabel{PIDpairs{numPairs}(1)} ' \rightarrow ' signalLabel{PIDpairs{numPairs}(2)}],'color','w','FontSize',14);
            h2=text(2,15,[signalLabel{PIDpairs{numPairs}(2)} ' \rightarrow ' signalLabel{PIDpairs{numPairs}(1)}],'color','w','FontSize',14);
            set(h1,'Rotation',90);
            set(h2,'Rotation',90);
    end

%% COMPUTE FEATURE-SPECIFIC INFORMATION TRANSFER

% set parameters for feature-specific information transfer (FIT)

    FIT_windowSize = 10; % number of consecutive time bins used for partial information decomposition
    FIT_timeJump = 5; % time jump used for partial information decomposition over time
    FIT_delay = -50:-1; % delays to consider as temporal lag between neural signals
    FIT_params.method = 'dr'; % estimation method: dr = Direct method
    FIT_params.bias = 'naive'; % bias correction procedure: naive = baised naive estimates
    FIT_params.bin_method_X = 'none';
    FIT_params.bin_method_Y = 'none';
    FIT_params.bin_method_C = 'none';
    FIT_params.verbose = 0;
    FITpairs = {[1 2]}; % neural pairs to consider when computing partial information decomposition
    
% compute FIT in both directions

    FITovertime = cell(2,length(FITpairs)); 
    for numPairs = 1:length(FITpairs)
        FIT_tmpSignal_fwd = [];
        FIT_tmpSignal_bwd = [];
        sTime = 1;
        for timeStep = 1:size(NeuralSignals_binned,2)/FIT_timeJump
            eTime = sTime + FIT_windowSize - 1;
            if sTime+min(FIT_delay)<0 || eTime>size(NeuralSignals_binned,2)
                sTime = sTime + FIT_timeJump;
            else
                FITdelays_fwd = [];    
                FITdelays_bwd = [];    
                for delay = FIT_delay
                    
                    % forward direction
                    Y = NeuralSignals_binned(FITpairs{numPairs}(2),sTime:eTime,:);
                    hY = NeuralSignals_binned(FITpairs{numPairs}(2),sTime+delay:eTime+delay,:);
                    hX = NeuralSignals_binned(FITpairs{numPairs}(1),sTime+delay:eTime+delay,:);
                    C = ReachFeature_binned;
                    C = repmat(C,size(Y,2),1);
                    Y = Y(:);
                    hY = hY(:);
                    hX = hX(:);
                    C = C(:);
                    FITdelays_fwd = [FITdelays_fwd FIT_choice(C',hX',hY',Y',FIT_params)];    

                    % backward direction
                    Y = NeuralSignals_binned(FITpairs{numPairs}(1),sTime:eTime,:);
                    hY = NeuralSignals_binned(FITpairs{numPairs}(1),sTime+delay:eTime+delay,:);
                    hX = NeuralSignals_binned(FITpairs{numPairs}(2),sTime+delay:eTime+delay,:);
                    Y = Y(:);
                    hY = hY(:);
                    hX = hX(:);
                    FITdelays_bwd = [FITdelays_bwd FIT_choice(C',hX',hY',Y',FIT_params)];    
                    
                end
                FIT_tmpSignal_fwd = [FIT_tmpSignal_fwd FITdelays_fwd'];
                FIT_tmpSignal_bwd = [FIT_tmpSignal_bwd FITdelays_bwd'];
                sTime = sTime + FIT_timeJump;
            end
        end
        FITovertime{1,numPairs} = FIT_tmpSignal_fwd;
        FITovertime{2,numPairs} = FIT_tmpSignal_bwd;
    end
    
% plot FIT in both directions

    for numPairs = 1:size(FITovertime,2)
        
        z_max = max([max(FITovertime{1,numPairs},[],'all') max(FITovertime{2,numPairs},[],'all')]);
        y_max = max([max(mean(FITovertime{1,numPairs})) max(mean(FITovertime{2,numPairs}))]);
        
        figure; 
            subplot(2,2,1); hold on;
                imagesc(FITovertime{1,numPairs},[0 z_max]);
                colorbar;
                xlim([.5 size(FITovertime{1,numPairs},2)+.5])
                ylim([.5 size(FITovertime{1,numPairs},1)+.5])
                yticks([1:5:size(FITovertime{1,numPairs},1)])
                yticklabels([FIT_delay(1):5:FIT_delay(end)])
                ylabel('M1 past (time bin)')
                xlabel('time')
                title(['FIT ' signalLabel{FITpairs{numPairs}(1)} ' \rightarrow ' signalLabel{FITpairs{numPairs}(2)} ' direction'])                                
            subplot(2,2,2); hold on;
                plot(mean(FITovertime{1,numPairs}),'k');
                xlim([1 size(FITovertime{1,numPairs},2)])
                ylabel('FIT (bits)')
                xlabel('time')
                ylim([0 y_max])
            subplot(2,2,3); hold on;
                imagesc(FITovertime{2,numPairs},[0 z_max]);
                colorbar;
                xlim([.5 size(FITovertime{2,numPairs},2)+.5])
                ylim([.5 size(FITovertime{2,numPairs},1)+.5])
                yticks([1:5:size(FITovertime{2,numPairs},1)])
                yticklabels([FIT_delay(1):5:FIT_delay(end)])
                ylabel('DLS past (time bin)')
                xlabel('time')
                title(['FIT ' signalLabel{FITpairs{numPairs}(2)} ' \rightarrow ' signalLabel{FITpairs{numPairs}(1)} ' direction'])                                
            subplot(2,2,4); hold on;
                plot(mean(FITovertime{2,numPairs}),'k');
                xlim([1 size(FITovertime{2,numPairs},2)])
                ylabel('FIT (bits)')
                xlabel('time')
                ylim([0 y_max])

    end

%% COMPUTE DIRECTED INFORMATION

% set parameters for directed information (DI)

    DI_windowSize = 10; % number of consecutive time bins used for partial information decomposition
    DI_timeJump = 5; % time jump used for partial information decomposition over time
    DI_delay = -50:-1; % delays to consider as temporal lag between neural signals
    DI_params.method = 'dr'; % estimation method: dr = Direct method
    DI_params.bias = 'naive'; % bias correction procedure: naive = baised naive estimates
    DI_params.bin_method_X = 'none';
    DI_params.bin_method_Y = 'none';
    DI_params.verbose = 0;
    DIpairs = {[1 2]}; % neural pairs to consider when computing partial information decomposition
    
% compute DI in both directions

    DIovertime = cell(2,length(DIpairs)); 
    for numPairs = 1:length(DIpairs)
        DI_tmpSignal_fwd = [];
        DI_tmpSignal_bwd = [];
        sTime = 1;
        for timeStep = 1:size(NeuralSignals_binned,2)/DI_timeJump
            eTime = sTime + DI_windowSize - 1;
            if sTime+min(DI_delay)<0 || eTime>size(NeuralSignals_binned,2)
                sTime = sTime + DI_timeJump;
            else
                DIdelays_fwd = [];    
                DIdelays_bwd = [];    
                for delay = DI_delay
                    
                    % forward direction
                    Y = NeuralSignals_binned(DIpairs{numPairs}(2),sTime:eTime,:);
                    hY = NeuralSignals_binned(DIpairs{numPairs}(2),sTime+delay:eTime+delay,:);
                    hX = NeuralSignals_binned(DIpairs{numPairs}(1),sTime+delay:eTime+delay,:);
                    Y = Y(:);
                    hY = hY(:);
                    hX = hX(:);
                    DIdelays_fwd = [DIdelays_fwd compute_DI(hX',hY',Y',DI_params)];    

                    % backward direction
                    Y = NeuralSignals_binned(DIpairs{numPairs}(1),sTime:eTime,:);
                    hY = NeuralSignals_binned(DIpairs{numPairs}(1),sTime+delay:eTime+delay,:);
                    hX = NeuralSignals_binned(DIpairs{numPairs}(2),sTime+delay:eTime+delay,:);
                    Y = Y(:);
                    hY = hY(:);
                    hX = hX(:);
                    DIdelays_bwd = [DIdelays_bwd compute_DI(hX',hY',Y',DI_params)];    

                end
                DI_tmpSignal_fwd = [DI_tmpSignal_fwd DIdelays_fwd'];
                DI_tmpSignal_bwd = [DI_tmpSignal_bwd DIdelays_bwd'];
                sTime = sTime + DI_timeJump;
            end
        end
        DIovertime{1,numPairs} = DI_tmpSignal_fwd;
        DIovertime{2,numPairs} = DI_tmpSignal_bwd;
    end
    
% plot DI

    for numPairs = 1:size(DIovertime,2)
        
        z_max = max([max(DIovertime{1,numPairs},[],'all') max(DIovertime{2,numPairs},[],'all')]);
        y_max = max([max(mean(DIovertime{1,numPairs})) max(mean(DIovertime{2,numPairs}))]);
        
        figure; 
            subplot(2,2,1); hold on;
                imagesc(DIovertime{1,numPairs},[0 z_max]);
                colorbar;
                xlim([.5 size(DIovertime{1,numPairs},2)+.5])
                ylim([.5 size(DIovertime{1,numPairs},1)+.5])
                yticks([1:5:size(DIovertime{1,numPairs},1)])
                yticklabels([DI_delay(1):5:DI_delay(end)])
                ylabel('M1 past (time bin)')
                xlabel('time')
                title(['DI ' signalLabel{DIpairs{numPairs}(1)} ' \rightarrow ' signalLabel{DIpairs{numPairs}(2)} ' direction'])                                
            subplot(2,2,2); hold on;
                plot(mean(DIovertime{1,numPairs}),'k');
                xlim([1 size(DIovertime{1,numPairs},2)])
                ylabel('DI (bits)')
                xlabel('time')
                ylim([0 y_max])
            subplot(2,2,3); hold on;
                imagesc(DIovertime{2,numPairs},[0 z_max]);
                colorbar;
                xlim([.5 size(DIovertime{2,numPairs},2)+.5])
                ylim([.5 size(DIovertime{2,numPairs},1)+.5])
                yticks([1:5:size(DIovertime{2,numPairs},1)])
                yticklabels([DI_delay(1):5:DI_delay(end)])
                ylabel('DLS past (time bin)')
                xlabel('time')
                title(['DI ' signalLabel{DIpairs{numPairs}(2)} ' \rightarrow ' signalLabel{DIpairs{numPairs}(1)} ' direction'])                                
            subplot(2,2,4); hold on;
                plot(mean(DIovertime{2,numPairs}),'k');
                xlim([1 size(DIovertime{2,numPairs},2)])
                ylabel('DI (bits)')
                xlabel('time')
                ylim([0 y_max])

    end
