%% wavelet filtering of DP grams and calculation of SL and LL components

% build not equal 0 - performs calculation of wavelet filtering
% build equal 0 - only shows the graphs using calculated variables (must stored
% in the memmory!!!)

close all

build = 1;

if build
    clear DPgrL2c55 NFgrL2c55 decL2c55CRa DPgrL2c50 NFgrL2c50 decL2c50NLa decL2c50CRa decL2c50NLph decL2c50CRph decL2c55NLa decL2c55NLph
    clear decL2c55CRph NFdecL2c55NLa decL2c55NLaI decL2c55NLphUI DPgrL1c50 NFgrL1c50 decL1c50NLa decL1c50CRa decL1c50NLph decL1c50CRph NFdecL1c50NLa
    clear decL1c55NLaI decL1c55NLphUI decL1c50NLaI decL1c50NLphUI decL2c55NLaI L1var55 L1var50 decL2c50NLaI decL2c50NLphUI
    clear DPgrL2c45 decL2c45NLa decL2c45CRa decL2c45NLph decL2c45CRph
    % parameters for wavelet analysis
    
    parW.nmm=25;
    parW.nb=6; % n. 1/3 oct bands
    parW.tt=13; % latency power law: tt*f(kHz)^-bbb
    parW.bbb=1;
    parW.p0=2e-5;
    parW.tc=-50;
    parW.tch=50;
    parW.c1=-0.3846;
    parW.c2=0.3846;
    parW.c3=1.;
    parW.c4=2.5;
    parW.c5=3;
    
    plotflag = 0;
    
    
    
    %% for L2 const and L1 varied
    
    FileNames = dir('Experiments/s039/L2const/DPgrams/');
    
    % for L1 = 50 dB and L2 is varied
    L2 = 45;
    nL1c = 1;
    for k=3:length(FileNames)
        
        
        if ~isempty(strfind(FileNames(k).name,['L2_' num2str(L2)]))
            load([FileNames(k).folder '/' FileNames(k).name])
            DPgrL2c45(:,nL1c) = fliplr(Poae);
            NFgrL2c45(:,nL1c) = fliplr(Nfloor);
            L1var45(nL1c) = L1o;
            nL1c = nL1c+1;
        end
        
    end
    
    
    L2 = 55; % for L1 =55 dB and L2 is varied
    nL1c = 1;
    for k=3:length(FileNames)
        
        
        if ~isempty(strfind(FileNames(k).name,['L2_' num2str(L2)]))
            load([FileNames(k).folder '/' FileNames(k).name])
            DPgrL2c55(:,nL1c) = fliplr(Poae);
            NFgrL2c55(:,nL1c) = fliplr(Nfloor);
            L1var55(nL1c) = L1o;
            nL1c = nL1c+1;
        end
        
    end
    
    
    f2r1 = fliplr(freq)./(2/f2f1 - 1)/1e3;
    
    
    %
    
    for k=1:size(DPgrL2c45,2)
        
        [fxL2_45, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,DPgrL2c45(:,k).',plotflag);
        
        decL2c45NLa(:,k) = DPgNLa;
        decL2c45CRa(:,k) = DPgCRa;
        decL2c45NLph(:,k) = DPgNLph;
        decL2c45CRph(:,k) = DPgCRph;
        decL2c45WL(k) = wlgraph;
        
        [fxL2_45, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,NFgrL2c45(:,k).',plotflag);
        
        NFdecL2c45NLa(:,k) = DPgNLa;  % noise floor
        
        [fxL2_55, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,DPgrL2c55(:,k).',plotflag);
        
        decL2c55NLa(:,k) = DPgNLa;
        decL2c55CRa(:,k) = DPgCRa;
        decL2c55NLph(:,k) = DPgNLph;
        decL2c55CRph(:,k) = DPgCRph;
        decL2c55WL(k) = wlgraph;
        
        contr_factor = 5;
        cm = 35;
        
        [fxL2_55, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,NFgrL2c55(:,k).',plotflag);
        
        NFdecL2c55NLa(:,k) = DPgNLa;
        
        % save s003bL2cRes.mat decL2c55CRph decL2c50CRph decL2c50NLa decL2c50NLph decL2c55NLa decL2c55CRa decL2c50CRa decL2c55NLph L1var50 L1var55 fx
    end
    
    
    
    
end

%% visualization of results for the paper

figure;

set(gcf, 'Units', 'centimeters');
% we set the position and dimension of the figure ON THE SCREEN
%
% NOTE: measurement units refer to the previous settings!
afFigurePosition = [1 1 12 10];
% [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% we link the dimension of the figure ON THE PAPER in such a way that
% it is equal to the dimension on the screen
%
% ATTENTION: if PaperPositionMode is not ’auto’ the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');

iFontSize = 10;
strFontName = 'Helvetica';
LWtick = 3.1;
LW = 1.2; % linewidth
MS = 8; % markersize

leftCornerX = 0.12;
leftCornerY = 0.57;
wEnv = 0.3;
hEnv = 0.36;
xD = 0.22;
yD = 0.08;
Fmin = 800;
Fmax = 4000;
Tdmin = -10;
Tdmax = 25;


% waterfall(fx,L1var55,decL2c55NLa')
% view([-1 45])
% subplot(2,1,2);
% waterfall(fx,L1var55,unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))
% view([-1 45])

L1i55 = L1var55(1):0.1:L1var55(end);

cycle = 2*pi;
%decL2c55NLphU = (unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))'/cycle;
decL2c55NLphU = (unwrap((decL2c55NLph')))'/cycle;
decL2c45NLphU = (unwrap((decL2c45NLph')))'/cycle;

%decL2c55NLphU(95:192,:) = decL2c55NLphU(95:192,:)-1;
%decL2c55NLphU(225:251,:) = decL2c55NLphU(225:251,:)-1;
%decL2c55NLphU(252:end,:) = decL2c55NLphU(252:end,:)-2;

%decL2c55NLphU = unwrap((unwrap((decL2c55NLph')))')/cycle;
%decL2c55NLphU = ((decL2c55NLph')*2*pi)'/cycle;
L1i55 = L1var55(1):0.1:L1var55(end);


L1i45 = L1var45(1):0.1:L1var45(end);



for k=1:size(decL2c55NLa,1)
    
    decL2c55NLaI(k,:) =  interp1(L1var55,decL2c55NLa(k,:),L1i55,'pchip');
    decL2c55NLphUI(k,:) =  interp1(L1var55,decL2c55NLphU(k,:),L1i55,'pchip');
    
end


for k=1:size(decL2c45NLa,1)
    
    decL2c45NLaI(k,:) =  interp1(L1var45,decL2c45NLa(k,:),L1i45,'pchip');
    decL2c45NLphUI(k,:) =  interp1(L1var45,decL2c45NLphU(k,:),L1i45,'pchip');
    
end



PhLim = [-2 1];

haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
[c1A,c1H] = contourf(fxL2_45,L1i45,decL2c45NLaI',30); colormap gray;

colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY 0.015 hEnv]);
text(4220,43.5,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
f2s_1 = 1350;
f2s_2 = 2800;
f2s_3 = 3600;

idxF2_1 = find(fxL2_45>=f2s_1,1,'first');
idxF2_2 = find(fxL2_45>=f2s_2,1,'first');
idxF2_3 = find(fxL2_45>=f2s_3,1,'first');

plot([fxL2_45(idxF2_1) fxL2_45(idxF2_1)],[L1i45(1) L1i45(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')

text(1000,61,'(E)','fontsize',iFontSize+3,'fontname',strFontName,'rotation',0,'horizontalalignment','center');
text(5300,61,'(F)','fontsize',iFontSize+3,'fontname',strFontName,'rotation',0,'horizontalalignment','center');
text(9000,61,'(G)','fontsize',iFontSize+3,'fontname',strFontName,'rotation',0,'horizontalalignment','center');
text(13400,61,'(H)','fontsize',iFontSize+3,'fontname',strFontName,'rotation',0,'horizontalalignment','center');
set(haR1A1,...
    'ticklength',[0.03 0.01],...
    'fontsize',iFontSize,...
    'ytick',[46,48,50,52,54,56,58,60],...
    'xtick',[1000,2000,3000,4000],...
    'ylim',[45 55],...
    'xlim',[700,3300],...
    'xticklabel','',...
    'fontname',strFontName,...
    'clim',[-20 10]);
text(2150,46.2,'{\itL}_2 = 45 dB','backgroundcolor','white','fontsize',iFontSize-2)

text(100,41,'{\itL}_1 (dB FPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );

haR1Ph1 = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
contourf(fxL2_55,L1i55,decL2c55NLaI',30); colormap gray; hold on;
plot([fxL2_55(idxF2_1) fxL2_55(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')

colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY-hEnv-yD 0.015 hEnv]);

text(2150,56,'{\itL}_2 = 55 dB','backgroundcolor','white','fontsize',iFontSize-2)

set(haR1Ph1,...
    'ticklength',[0.03 0.01],...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'ytick',[56,58,60,62,64,66,68,70],...
    'xtick',[1000,2000,3000,4000],...
    'xticklabel',{'1','2','3','4'},...
    'ylim',[55 63],...
    'xlim',[700,3300],...
    'clim',[-20 10]);

xlabel('{\itf}_2 (Hz)','fontsize',iFontSize,'fontname',strFontName );


haCa1 = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');

plot(L1var45,decL2c45NLa(idxF2_1,:),'ko-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55NLa(idxF2_1,:),'kx-','linewidth',LW,'markersize',MS);
plot(L1var45,decL2c45CRa(idxF2_1,:),'ro:','linewidth',LW-0.3,'markersize',MS);
plot(L1var55,decL2c55CRa(idxF2_1,:),'rx:','linewidth',LW-0.3,'markersize',MS);
plot(L1var45,NFdecL2c45NLa(idxF2_1,:),'ko:','linewidth',LW-0.5,'markersize',MS);
plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'kx:','linewidth',LW-0.5,'markersize',MS);


set(haCa1,...
    'ticklength',[0.02 0.01],...
    'fontsize',iFontSize,...
    'ylim',[-15,17],...
    'ytick',[-15,-10, -5, 0, 5, 10, 15],...
    'xticklabel','',...
    'fontname',strFontName);

ylabel('Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName );


haCph1 = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');


plot(L1var45,unwrap(decL2c45NLph(idxF2_1,:))/cycle,'ko-','linewidth',LW,'markersize',MS);
plot(L1var55,unwrap(decL2c55NLph(idxF2_1,:))/cycle,'kx-','linewidth',LW,'markersize',MS);
plot(L1var45,unwrap(decL2c45CRph(idxF2_1,:))/cycle,'ro:','linewidth',LW-0.3,'markersize',MS);
plot(L1var55,unwrap(decL2c55CRph(idxF2_1,:))/cycle,'rx:','linewidth',LW-0.3,'markersize',MS);

legend({['{\itL}_2 = ' num2str(45) ' dB'], ['{\itL}_2 = 55 dB']},...
    'position',[0.81 0.925 0.07 0.04],'fontsize',iFontSize-2,'orientation','vertical');

set(haCph1,...
    'ticklength',[0.02 0.01],...
    'ytick',[-0.4,-0.2,0,0.2,0.4],...
    'fontsize',iFontSize,...
    'fontname',strFontName);

xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
ylabel('Phase (cycles)','fontsize',iFontSize,'fontname',strFontName );

% ,'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_3,:),'ko-','linewidth',LW,'markersize',MS);




iResolution = 600;
strFileName = ['Figures/s039L2const'];
% %print -depsc ResultsJASADP17/ExpPokus.eps
print('-depsc', strcat(strFileName, '.eps'));
print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));


% ,'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_3,:),'ko-','linewidth',LW,'markersize',MS);
