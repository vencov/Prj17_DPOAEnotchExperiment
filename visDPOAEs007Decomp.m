%% wavelet filtering of DP grams and calculation of SL and LL components
% data for subject s001 (Fig. 1 and Fig. 2)


close all;

build = 1;

if build
    clear DPgrL2c55 NFgrL2c55 decL2c55CRa DPgrL2c50 NFgrL2c50 decL2c50NLa decL2c50CRa decL2c50NLph decL2c50CRph decL2c55NLa decL2c55NLph
    clear decL2c55CRph NFdecL2c55NLa decL2c55NLaI decL2c55NLphUI DPgrL1c50 NFgrL1c50 decL1c50NLa decL1c50CRa decL1c50NLph decL1c50CRph NFdecL1c50NLa
    clear decL1c55NLaI decL1c55NLphUI decL1c50NLaI decL1c50NLphUI decL2c55NLaI L1var55 L1var50 decL2c50NLaI decL2c50NLphUI

    FileNames = dir('Experiments/s007/r111/L2const/DPgrams/');
    
    L2 = 55; % for L2 =55 dB and L1 is varied
    nL1c = 1;
    
    for k=3:length(FileNames)
        
        
        if ~isempty(strfind(FileNames(k).name,['L2_' num2str(L2)]))
            load([FileNames(k).folder '/' FileNames(k).name])
            DPgrL2c55(:,nL1c) = Poae;
            NFgrL2c55(:,nL1c) = Nfloor;
            L1var55(nL1c) = L1o;
            nL1c = nL1c+1;
        end
    end
    
    
    f2r1 = freq./(2/f2f1 - 1)/1e3;
    
    
    % parameters for wavelet analysis
    
    parW.nmm=25;
    parW.nb=6; % n. 1/3 oct bands
    parW.tt=13; % latency power law: tt*f(kHz)^-bbb
    parW.bbb=1;
    parW.p0=2e-5;
    parW.tc=-50;
    parW.tch=50;
    parW.c1=-0.3846;  % will yield 5 ms
    parW.c2=0.3846;
    parW.c3=1.;
    parW.c4=2.5;
    parW.c5=3;
    
    plotflag = 0;
    
    for k=1:size(DPgrL2c55,2)
        
        [fx, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,DPgrL2c55(:,k).',plotflag);
        
        decL2c55Alla(:,k) = DPgAlla;
        decL2c55a(:,k) = DPgNLa;
        decL2c55CRa(:,k) = DPgCRa;
        decL2c55NLph(:,k) = DPgNLph;
        decL2c55CRph(:,k) = DPgCRph;
        decL2c55WL(k) = wlgraph;
        
        contr_factor = 5;
        cm = 35;
        
        [fx, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
            fc, DPtoctAlla, DPtoctNLa, DPtoctCRa,wlgraph] = WLTdecompFCE(parW,1000*f2r1,NFgrL2c55(:,k).',plotflag);
        NFdecL2c55NLa(:,k) = DPgNLa;
        
    end
    
end


%% visualization of results for the paper


% Fig. 1 showing L2 = 55 dB FPL and L1 varied with 5 dB step
% 
% figure;
% 
% set(gcf, 'Units', 'centimeters');
% % we set the position and dimension of the figure ON THE SCREEN
% %
% % NOTE: measurement units refer to the previous settings!
% afFigurePosition = [1 1 23 16];
% % [pos_x pos_y width_x width_y]
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % we link the dimension of the figure ON THE PAPER in such a way that
% % it is equal to the dimension on the screen
% %
% % ATTENTION: if PaperPositionMode is not â€™autoâ€™ the saved file
% % could have different dimensions from the one shown on the screen!
% set(gcf, 'PaperPositionMode', 'auto');
% 
% iFontSize = 12;
% strFontName = 'Helvetica';
% LWtick = 3.1;
% LW = 1.2; % linewidth
% MS = 8; % markersize
% 
% leftCornerX = 0.07;
% leftCornerY = 0.55;
% wEnv = 0.36;
% hEnv = 0.4;
% xD = 0.18;
% yD = 0.05;
% Fmin = 800;
% Fmax = 4000;
% Tdmin = -10;
% Tdmax = 25;
% 
% 
% L1i55 = L1var55(1):0.1:L1var55(end);
% 
% cycle = 2*pi;
% %decL2c55NLphU = (unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))'/cycle;
% % to correct unwrapping
% decL2c55NLphU = (unwrap((decL2c55NLph')))'/cycle;
% decL2c55NLphU(1:12,:) = decL2c55NLphU(1:12,:)+1;
% decL2c55NLphU(137:156,:) = decL2c55NLphU(137:156,:)-1;
% decL2c55NLphU(224:238,:) = decL2c55NLphU(224:238,:)-1;
% decL2c55NLphU(255:end,:) = decL2c55NLphU(255:end,:)-1;
% 
% 
% for k=1:size(decL2c55a,1)
%     decL2c55AllaI(k,:) =  interp1(L1var55,decL2c55Alla(k,:),L1i55,'pchip');
%     decL2c55NLaI(k,:) =  interp1(L1var55,decL2c55a(k,:),L1i55,'pchip');
%     decL2c55CRaI(k,:) =  interp1(L1var55,decL2c55CRa(k,:),L1i55,'pchip');
%     decL2c55NLphUI(k,:) =  interp1(L1var55,decL2c55NLphU(k,:),L1i55,'pchip');
%     
% end
% 
% 
% PhLim = [-2 1];
% 
% haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% [c1A,c1H] = contourf(fx,L1i55,decL2c55NLaI',30); colormap gray;
% 
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY 0.02 0.4]);
% text(4550,50,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% f2s_1 = 1700;
% f2s_2 = 2400;
% f2s_3 = 3600;
% 
% idxF2_1 = find(fx>=f2s_1,1,'first');
% idxF2_2 = find(fx>=f2s_2,1,'first');
% idxF2_3 = find(fx>=f2s_3,1,'first');
% 
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% 
% set(haR1A1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'clim',[-20 10],...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1000', '1500', '2000', '2500', '3000', '3500', '4000'});
% 
% 
% text(670,28,'{\itL}_1 (dB FPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );
% 
% haR1Ph1 = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% contourf(fx,L1i55,decL2c55NLphUI',30); colormap gray; hold on;
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% 
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY-hEnv-yD 0.02 0.4]);
% text(4550,50,'Phase (cycles)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% 
% set(haR1Ph1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1000', '1500', '2000', '2500', '3000', '3500', '4000'});
% 
% xlabel('{\itf}_2 (Hz)','fontsize',iFontSize,'fontname',strFontName );
% 
% zlim(PhLim);
% caxis(PhLim);
% 
% 
% haR1Acut = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55a(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55a(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55a(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'k+:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_2,:),'kx:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_3,:),'ko:','linewidth',LW-0.5,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_1,:),'r+:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_2,:),'rx:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_3,:),'ro:','linewidth',LW-0.3,'markersize',MS);
% 
% set(haR1Acut,...
%     'ticklength',[0.02 0.01],...
%     'ylim',[-30,10],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName);
% ylabel('Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName );
% 
% legend({['{\itf}_2 = ' num2str(f2s_1/1e3) 'kHz'], ['{\itf}_2 = ' num2str(f2s_2/1e3) 'kHz'],['{\itf}_2 = ' num2str(f2s_3/1e3) 'kHz']},...
%     'position',[0.65 0.88 0.07 0.1]);
% 
% haR1Phcut = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55NLphU(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55NLphU(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% 
% set(haR1Phcut,...
%     'ticklength',[0.02 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName);
% xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% ylabel('Phase (cycles)','fontsize',iFontSize,'fontname',strFontName );
% 
% iResolution = 600;
% strFileName = ['Figures/s007L2const55dB'];
% % %print -depsc ResultsJASADP17/ExpPokus.eps
% % print('-depsc', strcat(strFileName, '.eps'));
% print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
% print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));

%% updated figure with level map

% L_1s001 = [57 63 55 50];
% L_2s001 = [50 55 48 42];
% L_1s003 = [56, 61, 50, 55];
% L_2s003 = [50 55 43 48];
% 
% L_1s004 = [64];
% L_2s004 = [55];
% 
% L_1s039 = [53, 61];
% L_2s039 = [45 55];
% 
% L_2m = [46 50 52 56];
% L_1m = [50 56 58 62];
% 
% % Fig. 1 showing L2 = 55 dB FPL and L1 varied with 5 dB step
% 
% figure;
% 
% set(gcf, 'Units', 'centimeters');
% % we set the position and dimension of the figure ON THE SCREEN
% %
% % NOTE: measurement units refer to the previous settings!
% afFigurePosition = [1 1 23 13];
% % [pos_x pos_y width_x width_y]
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % we link the dimension of the figure ON THE PAPER in such a way that
% % it is equal to the dimension on the screen
% %
% % ATTENTION: if PaperPositionMode is not â€™autoâ€™ the saved file
% % could have different dimensions from the one shown on the screen!
% set(gcf, 'PaperPositionMode', 'auto');
% 
% iFontSize = 12;
% strFontName = 'Helvetica';
% LWtick = 3.1;
% LW = 1.2; % linewidth
% MS = 8; % markersize
% 
% leftCornerX = 0.052;
% leftCornerY = 0.56;
% wEnv = 0.24;
% hEnv = 0.38;
% xD = 0.145;
% yD = 0.06;
% Fmin = 800;
% Fmax = 4000;
% Tdmin = -10;
% Tdmax = 25;
% 
% 
% L1i55 = L1var55(1):0.1:L1var55(end);
% 
% cycle = 2*pi;
% %decL2c55NLphU = (unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))'/cycle;
% % to correct unwrapping
% decL2c55NLphU = (unwrap((decL2c55NLph')))'/cycle;
% decL2c55NLphU(1:12,:) = decL2c55NLphU(1:12,:)+1;
% decL2c55NLphU(137:156,:) = decL2c55NLphU(137:156,:)-1;
% decL2c55NLphU(224:238,:) = decL2c55NLphU(224:238,:)-1;
% decL2c55NLphU(255:end,:) = decL2c55NLphU(255:end,:)-1;
% 
% 
% for k=1:size(decL2c55a,1)
%     
%     decL2c55NLaI(k,:) =  interp1(L1var55,decL2c55a(k,:),L1i55,'pchip');
%     decL2c55CRaI(k,:) =  interp1(L1var55,decL2c55CRa(k,:),L1i55,'pchip');
%     decL2c55NLphUI(k,:) =  interp1(L1var55,decL2c55NLphU(k,:),L1i55,'pchip');
%     
% end
% 
% 
% PhLim = [-2 1];
% 
% haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% [c1A,c1H] = contourf(fx,L1i55,decL2c55NLaI',30); colormap gray;
% 
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY 0.02 hEnv]);
% text(4750,50,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% f2s_1 = 1700;
% f2s_2 = 2400;
% f2s_3 = 3600;
% 
% idxF2_1 = find(fx>=f2s_1,1,'first');
% idxF2_2 = find(fx>=f2s_2,1,'first');
% idxF2_3 = find(fx>=f2s_3,1,'first');
% 
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% text(3500,40,'A','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% set(haR1A1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'clim',[-20 10],...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});
% 
% 
% text(640,28,'{\itL}_1 (dB FPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );
% 
% haR1Ph1 = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% contourf(fx,L1i55,decL2c55NLphUI',30); colormap gray; hold on;
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% text(3500,60,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY-hEnv-yD 0.02 hEnv]);
% text(4750,50,'Phase (cycles)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% 
% set(haR1Ph1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});
% 
% xlabel('{\itf}_2 (kHz)','fontsize',iFontSize,'fontname',strFontName );
% 
% zlim(PhLim);
% caxis(PhLim);
% 
% 
% haR1Acut = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55a(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55CRa(idxF2_1,:),'r+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55a(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55a(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'k+:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_2,:),'kx:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_3,:),'ko:','linewidth',LW-0.5,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_1,:),'r+:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_2,:),'rx:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_3,:),'ro:','linewidth',LW-0.3,'markersize',MS);
% text(29,-10,'Amplitude (dB SPL)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% set(haR1Acut,...
%     'ticklength',[0.02 0.01],...
%     'ylim',[-30,10],...
%     'fontsize',iFontSize,...
%     'xticklabel','',...
%     'fontname',strFontName);
% text(36,6,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% 
% legend({['{\itf}_2 = ' num2str(f2s_1/1e3) 'kHz'], ['{\itf}_2 = ' num2str(f2s_2/1e3) 'kHz'],['{\itf}_2 = ' num2str(f2s_3/1e3) 'kHz']},...
%     'position',[0.7 leftCornerY-hEnv 0.12 0.17]);
% 
% haR1Phcut = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55NLphU(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55NLphU(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% 
% set(haR1Phcut,...
%     'ticklength',[0.02 0.01],...
%     'fontsize',iFontSize,...
%     'ylim',[-1,0.65],...
%     'fontname',strFontName);
% 
% xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% text(29,-0.2,'Phase (cycles)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% 
% 
% 
% haL1L2 = axes('position',[leftCornerX+2*wEnv+1.5*xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L_2s001,L_1s001,'k+','linewidth',LW,'markersize',MS);
% plot(L_2s003,L_1s003,'kx','linewidth',LW,'markersize',MS);
% plot(L_2s004,L_1s004,'ko','linewidth',LW,'markersize',MS);
% plot(L_2s039,L_1s039,'k*','linewidth',LW,'markersize',MS);
% plot(L_2m,L_1m,'r+-','linewidth',LW,'markersize',MS);
% 
% set(haL1L2,...
%     'ticklength',[0.02 0.01],...
%     'ylim',[47,65],...
%     'xlim',[40,60],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName);
% 
% ylabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% xlabel('{\itL}_2 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% text(41,63.2,'D','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% legend({'s001','s003','s004','s039','model'},...
%     'position',[0.905 0.59 0.08 0.17],'numcolumns',1);
% 
% 
% iResolution = 600;
% strFileName = ['Figures/s007L2const55dB'];
% % %print -depsc ResultsJASADP17/ExpPokus.eps
% % print('-depsc', strcat(strFileName, '.eps'));
% print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
% print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));
% 



%% second update/figure with level map
% 
% L_1s001 = [57 63 55 50];
% L_2s001 = [50 55 48 42];
% L_1s003 = [56, 61, 50, 55];
% L_2s003 = [50 55 43 48];
% 
% L_1s004 = [64];
% L_2s004 = [55];
% 
% L_1s039 = [53, 61];
% L_2s039 = [45 55];
% 
% L_2m = [46 50 52 56];
% L_1m = [50 56 58 62];
% 
% % Fig. 1 showing L2 = 55 dB FPL and L1 varied with 5 dB step
% 
% figure;
% 
% set(gcf, 'Units', 'centimeters');
% % we set the position and dimension of the figure ON THE SCREEN
% %
% % NOTE: measurement units refer to the previous settings!
% afFigurePosition = [1 1 23 15];
% % [pos_x pos_y width_x width_y]
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % we link the dimension of the figure ON THE PAPER in such a way that
% % it is equal to the dimension on the screen
% %
% % ATTENTION: if PaperPositionMode is not â€™autoâ€™ the saved file
% % could have different dimensions from the one shown on the screen!
% set(gcf, 'PaperPositionMode', 'auto');
% 
% iFontSize = 12;
% strFontName = 'Helvetica';
% LWtick = 3.1;
% LW = 1.2; % linewidth
% MS = 8; % markersize
% 
% leftCornerX = 0.052;
% leftCornerY = 0.64;
% wEnv = 0.24;
% hEnv = 0.3;
% xD = 0.145;
% yD = 0.03;
% Fmin = 800;
% Fmax = 4000;
% Tdmin = -10;
% Tdmax = 25;
% 
% 
% L1i55 = L1var55(1):0.1:L1var55(end);
% 
% cycle = 2*pi;
% %decL2c55NLphU = (unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))'/cycle;
% % to correct unwrapping
% decL2c55NLphU = (unwrap((decL2c55NLph')))'/cycle;
% decL2c55NLphU(1:12,:) = decL2c55NLphU(1:12,:)+1;
% decL2c55NLphU(137:156,:) = decL2c55NLphU(137:156,:)-1;
% decL2c55NLphU(224:238,:) = decL2c55NLphU(224:238,:)-1;
% decL2c55NLphU(255:end,:) = decL2c55NLphU(255:end,:)-1;
% decL2c55CRphU = (unwrap((decL2c55CRph')))'/cycle;
% 
% for k=1:size(decL2c55a,1)
%     
%     decL2c55NLaI(k,:) =  interp1(L1var55,decL2c55a(k,:),L1i55,'pchip');
%     decL2c55CRaI(k,:) =  interp1(L1var55,decL2c55CRa(k,:),L1i55,'pchip');
%     decL2c55NLphUI(k,:) =  interp1(L1var55,decL2c55NLphU(k,:),L1i55,'pchip');
%     
% end
% 
% 
% PhLim = [-2 1];
% 
% haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% [c1A,c1H] = contourf(fx,L1i55,decL2c55AllaI',30); colormap gray;
% 
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY 0.02 hEnv]);
% text(4750,50,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% f2s_1 = 1700;
% f2s_2 = 2400;
% f2s_3 = 3600;
% 
% idxF2_1 = find(fx>=f2s_1,1,'first');
% idxF2_2 = find(fx>=f2s_2,1,'first');
% idxF2_3 = find(fx>=f2s_3,1,'first');
% 
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% text(3500,40,'A','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% set(haR1A1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'clim',[-20 10],...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});
% 
% 
% text(640,28,'{\itL}_1 (dB FPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );
% 
% haR1Ph1 = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% contourf(fx,L1i55,decL2c55NLaI.',30); colormap gray; hold on;
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% text(3500,60,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY-hEnv-yD 0.02 hEnv]);
% text(4750,50,'Phase (cycles)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% 
% set(haR1Ph1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'clim',[-20 10],...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});
% 
% haR1Ph2 = axes('position',[leftCornerX leftCornerY-2*hEnv-2*yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% contourf(fx,L1i55,decL2c55CRaI',30); colormap gray; hold on;
% plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
% plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
% text(3500,60,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% colorbar('position',[leftCornerX+wEnv+0.01 leftCornerY-hEnv-yD 0.02 hEnv]);
% text(4750,50,'Phase (cycles)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
% 
% set(haR1Ph1,...
%     'ticklength',[0.03 0.01],...
%     'fontsize',iFontSize,...
%     'fontname',strFontName,...
%     'clim',[-20 10],...
%     'xtick',[1000 1500 2000 2500 3000 3500 4000],...
%     'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});
% 
% 
% xlabel('{\itf}_2 (kHz)','fontsize',iFontSize,'fontname',strFontName );
% 
% % zlim(PhLim);
% % caxis(PhLim);
% 
% 
% haR1Acut = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55a(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55a(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55a(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'k+:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_2,:),'kx:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_3,:),'ko:','linewidth',LW-0.5,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_1,:),'r+:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_2,:),'rx:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_3,:),'ro:','linewidth',LW-0.3,'markersize',MS);
% %text(29,-10,'Amplitude (dB SPL)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% set(haR1Acut,...
%     'ticklength',[0.02 0.01],...
%     'ylim',[-30,10],...
%     'fontsize',iFontSize,...
%     'xticklabel','',...
%     'fontname',strFontName);
% text(36,6,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% 
% legend({['{\itf}_2 = ' num2str(f2s_1/1e3) 'kHz'], ['{\itf}_2 = ' num2str(f2s_2/1e3) 'kHz'],['{\itf}_2 = ' num2str(f2s_3/1e3) 'kHz']},...
%     'position',[0.7 leftCornerY-hEnv 0.12 0.17]);
% 
% haR1Phcut = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55NLphU(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55NLphU(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% 
% set(haR1Phcut,...
%     'ticklength',[0.02 0.01],...
%     'fontsize',iFontSize,...
%     'ylim',[-1,0.65],...
%     'fontname',strFontName);
% 
% xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% text(29,-0.2,'Phase (cycles)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% 
% 
% 
% haR1Acut = axes('position',[leftCornerX+2*wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55CRa(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55CRa(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55CRa(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'k+:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_2,:),'kx:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,NFdecL2c55NLa(idxF2_3,:),'ko:','linewidth',LW-0.5,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_1,:),'r+:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_2,:),'rx:','linewidth',LW-0.3,'markersize',MS);
% % plot(L1var55,decL2c55CRa(idxF2_3,:),'ro:','linewidth',LW-0.3,'markersize',MS);
% text(29,-10,'Amplitude (dB SPL)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% set(haR1Acut,...
%     'ticklength',[0.02 0.01],...
%     'ylim',[-30,10],...
%     'fontsize',iFontSize,...
%     'xticklabel','',...
%     'fontname',strFontName);
% text(36,6,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
% 
% legend({['{\itf}_2 = ' num2str(f2s_1/1e3) 'kHz'], ['{\itf}_2 = ' num2str(f2s_2/1e3) 'kHz'],['{\itf}_2 = ' num2str(f2s_3/1e3) 'kHz']},...
%     'position',[0.7 leftCornerY-hEnv 0.12 0.17]);
% 
% haR1Phcut = axes('position',[leftCornerX+2*wEnv+xD leftCornerY-hEnv  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
% plot(L1var55,decL2c55NLphU(idxF2_1,:),'k+-','linewidth',LW,'markersize',MS);
% plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
% plot(L1var55(3:end),decL2c55NLphU(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
% 
% set(haR1Phcut,...
%     'ticklength',[0.02 0.01],...
%     'fontsize',iFontSize,...
%     'ylim',[-1,0.65],...
%     'fontname',strFontName);
% 
% xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
% text(29,-0.2,'Phase (cycles)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
% 
% 
% 
% 
% iResolution = 600;
% strFileName = ['Figures/s007L2const55dB'];
% % %print -depsc ResultsJASADP17/ExpPokus.eps
% % print('-depsc', strcat(strFileName, '.eps'));
% print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
% print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));
% 
% 


%% plot DPgrams in time-frequency domain yielded by wavelet filterbank

%return;

figure;

set(gcf, 'Units', 'centimeters');
% we set the position and dimension of the figure ON THE SCREEN
%
% NOTE: measurement units refer to the previous settings!
afFigurePosition = [1 1 23 6.7];
% [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% we link the dimension of the figure ON THE PAPER in such a way that
% it is equal to the dimension on the screen
%
% ATTENTION: if PaperPositionMode is not â€™autoâ€™ the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');

iFontSize = 12;
strFontName = 'Helvetica';
LWtick = 3.1;
LW = 2.2;

leftCornerX = 0.07;
leftCornerY = 0.22;
wEnv = 0.27;
hEnv = 0.75;
xD = 0.05;
yD = 0.05;
Fmin = 800;
Fmax = 4000;
Tdmin = -10;
Tdmax = 25;

L1vv = 65;
idxL = find(L1var55==L1vv)
cm1 = max(max(decL2c55WL(idxL).coewd));
L1vv = 60;
idxL = find(L1var55==L1vv)
cm2 = max(max(decL2c55WL(idxL).coewd));
L1vv = 50;
idxL = find(L1var55==L1vv)
cm3 = max(max(decL2c55WL(idxL).coewd));
cm = max([cm2 cm3 cm1]);
%decL2c55WL
contr_factor = 3;
haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add');
L1vv = 50;
idxL = find(L1var55==L1vv)
contourf((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,(-decL2c55WL(idxL).M/2+1:decL2c55WL(idxL).M/2)*decL2c55WL(idxL).tsm,1-exp(-abs(decL2c55WL(idxL).coewd(1:decL2c55WL(idxL).M,:)./(cm/contr_factor))));colormap gray;hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c1*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c4*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c2*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
axis ([Fmin Fmax Tdmin Tdmax])
text(1000,21,'A','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')

set(haR1A1,...
    'ticklength',[0.02 0.01],...
    'xticklabel',{},...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'xtick',[1000 2000 3000 4000],...
    'xticklabel',{'1', '2', '3', '4'});


%xlabel ('Frequency (dB FPL)')
xlabel ('{\itf}_2 (kHz)')
ylabel ('DPOAE delay (ms)')
%text(2000,27,'{\itf}_2/{\itf}_1 = 1.11')
text(3000,19,'{\itL}_1 = 50 dB\newline{\itL}_2 = 55 dB','backgroundcolor','white')

haR1A2 = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add');
L1vv = 60;
idxL = find(L1var55==L1vv)
contourf((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,(-decL2c55WL(idxL).M/2+1:decL2c55WL(idxL).M/2)*decL2c55WL(idxL).tsm,1-exp(-abs(decL2c55WL(idxL).coewd(1:decL2c55WL(idxL).M,:)./(cm/contr_factor))));colormap gray;hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c1*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c4*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c2*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
axis ([Fmin Fmax Tdmin Tdmax])
%xlabel ('Frequency (dB FPL)')
xlabel ('{\itf}_2 (kHz)')
% ylabel ('DPOAE delay (ms)')
text(3000,19,'{\itL}_1 = 60 dB\newline{\itL}_2 = 55 dB','backgroundcolor','white')
text(1000,21,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')

set(haR1A2,...
    'ticklength',[0.02 0.01],...
    'xticklabel',{},...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'xtick',[1000 2000 3000 4000],...
    'xticklabel',{'1', '2', '3', '4'});


haR1A3 = axes('position',[leftCornerX+2*wEnv+2*xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add');
L1vv = 65;
idxL = find(L1var55==L1vv)
contourf((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,(-decL2c55WL(idxL).M/2+1:decL2c55WL(idxL).M/2)*decL2c55WL(idxL).tsm,1-exp(-abs(decL2c55WL(idxL).coewd(1:decL2c55WL(idxL).M,:)./(cm/contr_factor))));colormap gray;hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c1*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c4*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
plot((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfwH,decL2c55WL(idxL).c2*decL2c55WL(idxL).tt*((1:decL2c55WL(idxL).NW)*decL2c55WL(idxL).dfw).^-decL2c55WL(idxL).bbb,'white', 'LineWidth', 2);hold on
axis ([Fmin Fmax Tdmin Tdmax])
xlabel ('{\itf}_2 (kHz)')
%ylabel ('DPOAE delay (ms)')
text(3000,19,'{\itL}_1 = 65 dB\newline{\itL}_2 = 55 dB','backgroundcolor','white')
text(1000,21,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
set(haR1A3,...
    'ticklength',[0.02 0.01],...
    'xticklabel',{},...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'xtick',[1000 2000 3000 4000],...
    'xticklabel',{'1', '2', '3', '4'});


iResolution = 600;
strFileName = ['Figures/s007WLdecompL2c'];
% %print -depsc ResultsJASADP17/ExpPokus.eps
% print('-depsc', strcat(strFileName, '.eps'));
print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));

%% updated figure with level map and CR component inestead of phase

L_1s001 = [57 63 55 50];
L_2s001 = [50 55 48 42];
L_1s003 = [56, 61, 50, 55];
L_2s003 = [50 55 43 48];

L_1s004 = [64];
L_2s004 = [55];

L_1s039 = [53, 61];
L_2s039 = [45 55];

L_2m = [46 50 52 56];
L_1m = [50 56 58 62];

% Fig. 1 showing L2 = 55 dB FPL and L1 varied with 5 dB step

figure;

set(gcf, 'Units', 'centimeters');
% we set the position and dimension of the figure ON THE SCREEN
%
% NOTE: measurement units refer to the previous settings!
afFigurePosition = [1 1 23 12];
% [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% we link the dimension of the figure ON THE PAPER in such a way that
% it is equal to the dimension on the screen
%
% ATTENTION: if PaperPositionMode is not â€™autoâ€™ the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');

iFontSize = 12;
strFontName = 'Helvetica';
LWtick = 3.1;
LW = 1.2; % linewidth
MS = 8; % markersize

leftCornerX = 0.052;
leftCornerY = 0.56;
wEnv = 0.24;
hEnv = 0.39;
xD = 0.145;
yD = 0.03;
Fmin = 800;
Fmax = 4000;
Tdmin = -10;
Tdmax = 25;
wEnv2 = 0.24;
hEnv2 = 0.39;


L1i55 = L1var55(1):0.1:L1var55(end);

cycle = 2*pi;
%decL2c55NLphU = (unwrap((decL2c55NLph'-decL2c55NLph(:,1)')*2*pi))'/cycle;
% to correct unwrapping
decL2c55NLphU = (unwrap((decL2c55NLph')))'/cycle;
decL2c55NLphU(1:12,:) = decL2c55NLphU(1:12,:)+1;
decL2c55NLphU(137:156,:) = decL2c55NLphU(137:156,:)-1;
decL2c55NLphU(224:238,:) = decL2c55NLphU(224:238,:)-1;
decL2c55NLphU(255:end,:) = decL2c55NLphU(255:end,:)-1;


for k=1:size(decL2c55a,1)
    
    decL2c55NLaI(k,:) =  interp1(L1var55,decL2c55a(k,:),L1i55,'pchip');
    decL2c55CRaI(k,:) =  interp1(L1var55,decL2c55CRa(k,:),L1i55,'pchip');
    decL2c55NLphUI(k,:) =  interp1(L1var55,decL2c55NLphU(k,:),L1i55,'pchip');
    
end


PhLim = [-2 1];

haR1A1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
[c1A,c1H] = contourf(fx,L1i55,decL2c55NLaI',30); colormap gray;

colorbar('position',[leftCornerX+wEnv+0.005 leftCornerY 0.02 hEnv]);
text(4730,50,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');
f2s_1 = 1600;
f2s_2 = 2020;
f2s_3 = 3600;

idxF2_1 = find(fx>=f2s_1,1,'first');
idxF2_2 = find(fx>=f2s_2,1,'first');
idxF2_3 = find(fx>=f2s_3,1,'first');

plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
%plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
text(3500,40,'A','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
set(haR1A1,...
    'ticklength',[0.03 0.01],...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'clim',[-20 10],...
    'xtick',[1000 1500 2000 2500 3000 3500 4000],...
    'xticklabel','');


text(640,28,'{\itL}_1 (dB FPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );

haR1Ph1 = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
contourf(fx,L1i55,decL2c55CRaI.',30); colormap gray; hold on;
plot([fx(idxF2_1) fx(idxF2_1)],[L1i55(1) L1i55(end)],'--','color','red')
plot([fx(idxF2_2) fx(idxF2_2)],[L1i55(1) L1i55(end)],'--','color','red')
%plot([fx(idxF2_3) fx(idxF2_3)],[L1i55(1) L1i55(end)],'--','color','red')
text(3500,60,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
colorbar('position',[leftCornerX+wEnv+0.005 leftCornerY-hEnv-yD 0.02 hEnv]);
text(4730,50,'Amplitude (dB SPL)','fontsize',iFontSize,'fontname',strFontName,'rotation',-90,'horizontalalignment','center');

set(haR1Ph1,...
    'ticklength',[0.03 0.01],...
    'fontsize',iFontSize,...
    'fontname',strFontName,...
    'xtick',[1000 1500 2000 2500 3000 3500 4000],...
    'xticklabel',{'1', '1.5', '2', '2.5','3', '3.5', '4'});

xlabel('{\itf}_2 (kHz)','fontsize',iFontSize,'fontname',strFontName );

zlim([-20 10]);
caxis([-20 10]);


haR1Acut = axes('position',[leftCornerX+wEnv+xD leftCornerY+hEnv-hEnv2  wEnv2 hEnv2],'xscale','lin','yscale','lin','nextplot','add','box','on');
plot(L1var55,decL2c55a(idxF2_1,:),'ko-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55a(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55CRa(idxF2_1,:),'ro-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55CRa(idxF2_2,:),'rx-','linewidth',LW,'markersize',MS);
%plot(L1var55,decL2c55Alla(idxF2_1,:),'bo-','linewidth',LW,'markersize',MS);
%plot(L1var55,decL2c55Alla(idxF2_2,:),'bx-','linewidth',LW,'markersize',MS);
%plot(L1var55(3:end),decL2c55a(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);
plot(L1var55,NFdecL2c55NLa(idxF2_1,:),'k+:','linewidth',LW-0.5,'markersize',MS);
plot(L1var55,NFdecL2c55NLa(idxF2_2,:),'kx:','linewidth',LW-0.5,'markersize',MS);
%plot(L1var55,NFdecL2c55NLa(idxF2_3,:),'ko:','linewidth',LW-0.5,'markersize',MS);
% plot(L1var55,decL2c55CRa(idxF2_1,:),'r+:','linewidth',LW-0.3,'markersize',MS);
% plot(L1var55,decL2c55CRa(idxF2_2,:),'rx:','linewidth',LW-0.3,'markersize',MS);
% plot(L1var55,decL2c55CRa(idxF2_3,:),'ro:','linewidth',LW-0.3,'markersize',MS);
text(29,-10,'Amplitude (dB SPL)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );
set(haR1Acut,...
    'ticklength',[0.02 0.01],...
    'ylim',[-29,10],...
    'fontsize',iFontSize,...
    'xticklabel','',...
    'fontname',strFontName);
text(36,6,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')

legend({['{\itf}_2 = ' num2str(f2s_1/1e3) 'kHz'], ['{\itf}_2 = ' num2str(f2s_2/1e3) 'kHz']},...
    'position',[0.7 leftCornerY-hEnv 0.12 0.14]);

haR1Phcut = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
plot(L1var55,decL2c55NLphU(idxF2_1,:),'ko-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55NLphU(idxF2_2,:),'kx-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55CRphU(idxF2_1,:),'ro-','linewidth',LW,'markersize',MS);
plot(L1var55,decL2c55CRphU(idxF2_2,:),'rx-','linewidth',LW,'markersize',MS);
%plot(L1var55(3:end),decL2c55NLphU(idxF2_3,3:end),'ko-','linewidth',LW,'markersize',MS);

set(haR1Phcut,...
    'ticklength',[0.02 0.01],...
    'fontsize',iFontSize,...
    'ylim',[-0.5,1.5],...
    'fontname',strFontName);

xlabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
text(29,0.5,'Phase (cycles)','rotation',90,'fontsize',iFontSize,'fontname',strFontName,'horizontalalignment','center' );



haL1L2 = axes('position',[leftCornerX+2*wEnv+1.5*xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add','box','on');
plot(L_2s001,L_1s001,'k+','linewidth',LW,'markersize',MS);
plot(L_2s003,L_1s003,'kx','linewidth',LW,'markersize',MS);
plot(L_2s004,L_1s004,'ko','linewidth',LW,'markersize',MS);
plot(L_2s039,L_1s039,'k*','linewidth',LW,'markersize',MS);
plot(L_2m,L_1m,'r+-','linewidth',LW,'markersize',MS);

set(haL1L2,...
    'ticklength',[0.02 0.01],...
    'ylim',[47,65],...
    'xlim',[40,60],...
    'fontsize',iFontSize,...
    'fontname',strFontName);

ylabel('{\itL}_1 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
xlabel('{\itL}_2 (dB FPL)','fontsize',iFontSize,'fontname',strFontName );
text(41,63.2,'D','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold')
legend({'s001','s003','s004','s039','model'},...
    'position',[0.905 0.59 0.08 0.17],'numcolumns',1);


iResolution = 600;
strFileName = ['Figures/s007L2const55dBREV'];
% %print -depsc ResultsJASADP17/ExpPokus.eps
% print('-depsc', strcat(strFileName, '.eps'));
print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));

