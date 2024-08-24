%% script which calculates DPgrams for s007 for f2/f1 = 1.11
% there are two level paradigms
% 1. L2 = 50 or L2 = 55 dB FPL and L1 varied



clear all;



% place files to the list

Twin = 100e-3; % duration of analysis window (100 ms fixed)
Nonset = 3; % first three an. windows get out

% will be given (data in experiments)


pathname = 'Experiments/s007/r111/L2const'
FldInf = dir(pathname);



%% for L2 = 55 dB FPL
L2 = 55;
L1var = 35:5:70;

recMean = 0;


for k=1:length(L1var)
    
    recMean = 0;
    numR = 0;
    clear recPmat;
    for j=1:size(FldInf,1)
        
        L1loc = strfind(FldInf(j).name,['L1_' num2str(L1var(k))]);
        L2loc = strfind(FldInf(j).name,['L2_' num2str(L2)]);
        
        if (~isempty(L1loc))&&(~isempty(L2loc))
            load([pathname '/' FldInf(j).name]); % load recordedSignalFolderInfo(j).name
            fcutoff =300; % cutoff freq of 2order high pass butterworth filter
            [b,a] = butter(4,fcutoff/(fs/2),'high');
            nahravka = filter(b,a,recordedSignal);
            
            %recordMean = recordMean + nahravka(:,1);
            
            recP = nahravka(:,1);
            numR = numR + 1;
            recMean = recMean + recP;
            recPmat(:,numR) = recP;
            
        end
    end
    
    recMean = recMean/numR;
    
    % if data were collected, extract DP-gram by LSF technique
    if length(recMean)>1
        fdp1 = 2*F2start/f2f1-F2start; % logsweep start frequency  (for 'tone' type, only f1 is necessary)
        fdp2 = 2*F2stop/f2f1-F2stop; % logsweep end frequency    (for 'tone' type, only f1 is necessary)
        
        D = 4170; % skip onset part due to latency of the sound card
        
        octpersec = 0.5;  % oct
        nfilt = 5512;                        % window size (125ms, recommended for 0.5 oct/s)
        nstep = 51;                          % step size
        dt = 1/fs;
        recMeanS = recMean(D:end);
        recPmatS = recPmat(D:end,:);
        
        twin = nfilt/fs;
        tstep = nstep/fs;
        if fdp1<fdp2   % upward log sweep
            T = log2(fdp2/fdp1)/octpersec;     % poustelo se rychlosti 0.5 oktavy/s                               % 6 second swee
        else  % downward log sweep
            T = log2(fdp1/fdp2)/octpersec;     % poustelo se rychlosti 0.5 oktavy/s                               % 6 second swee
        end
        nsize = min(floor(T*fs), length(recMeanS));
        nfilt = floor(twin * fs + 0.5);
        nstep = floor(tstep * fs + 0.5);
        ntotal = floor((nsize-nfilt)/nstep) + 1;
        
        windowType = 2; % hann window
        Wvec = computeWin(0:nfilt-1, windowType); % han win
        
        
        params.nskip = D;%                         % ?samples to skip, default = 0 (e.g. due to system delay, propogation delay)
        
        params.windowType = 2;                      % window function (-1 = rectangular, 0 = welch, 1 = hamming, 2 = hann [default], 3 = parzen, 4 = blackman)
        params.fs = 44.1e3;                          % sampling rate (default = 44100);
        
        params.T = T;     % poustelo se rychlosti 0.5 oktavy/s                               % 6 second sweep
        params.nfilt = 5512;                        % window size (125ms, recommended for 0.5 oct/s)
        params.nstep = 200;                          % step size
        params.component(1).type = 'logsweep';      % component type (values can be 'logsweep', 'tone', 'linear', 'custom'; for 'custom' a phase function, pf needs to be defined)
        
        params.component(1).f1 = fdp1;
        params.component(1).f2 = fdp2;
        params.GainMicA = GainMicA;
        params.f2f1 = f2f1;
        % oae = nseoae(recMean, params);
        
        %freq = oae.f; % from fdp to f2
        
        Coae = 1;
        %  Poae = giveTPaWin(freq,oae.mag.*exp(sqrt(-1)*oae.phase),GainMicA,fs,Coae,fx);
        
        
        params.dFnoise = -50; % position of the lowest signal frequency in data
        [~, yp, NumOK] = nseoae541(params,recPmat(:,1));
        for k2=2:size(recPmat,2)
            [oae2, yp, NumOK] = nseoae541(params,recPmat(:,k2),yp,NumOK);
            
        end
        
        freq = oae2.f;
        
        Poae = giveTPaWin(freq,oae2.mag.*exp(sqrt(-1)*oae2.phase),GainMicA,fs,Coae,fx);
        Nfloor = abs(giveTPaWin(freq,oae2.nFloor,GainMicA,fs,Coae,fx));
        RMSerr = oae2.rmserror;
        L1o = L1var(k);
        L2o = L2;
        save([pathname '/DPgrams/dpgrL1_' num2str(L1var(k)) 'dB_L2_' num2str(L2) 'dB.mat'],'freq','f2f1','L2o','L1o','Poae','Nfloor','RMSerr');
    end
end

