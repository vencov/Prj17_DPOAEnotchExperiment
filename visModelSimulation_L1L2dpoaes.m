%% script showing Fig. 4 with model simulations

% vaclav.vencovsky@gmail.com


close all
clear all;

F2 = 1600; 
F1 = 1400;

Fdp = 2*F1-F2;


% calculate the base waves
gain = 1.05;
   
phi1 = 0; % phases
phi2 = 0;

FolderRes2 = 'Simulations/';

fs = 600e3;
Tan = 19e-3; % start of analysis window
Nonset = round(Tan*fs);


load([FolderRes2 'JAJAdpL1L2mapF1_' num2str(F1) 'Hz_F2_' num2str(F2) 'Hz_gain' num2str(gain*100) 'R21.mat']);

DPf2_1600 = oaedp.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);
A2 = A2List;
A1 = A1List;


%%

 
figure;

set(gcf, 'Units', 'centimeters');
% we set the position and dimension of the figure ON THE SCREEN
%
% NOTE: measurement units refer to the previous settings!
afFigurePosition = [1 1 23 15];
% [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% we link the dimension of the figure ON THE PAPER in such a way that
% it is equal to the dimension on the screen
%
% ATTENTION: if PaperPositionMode is not ’auto’ the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');

iFontSize = 12;
strFontName = 'Helvetica';
LWtick = 2.2;
LW = 1.3;

leftCornerX = 0.06;
leftCornerY = 0.64;
yD = 0.04;
xD= 0.07;
xD2 = 0.04
wEnv = 0.26;
hEnv = 0.28;
Xlims = [20 70];
Ylims = [40 70];
Clims = [10 50];
A1i = A1(1):0.1:A1(end);
A2i = A2(1):0.1:A2(end);
[A1i,A2i] = meshgrid(A1(1):0.1:A1(end),A2(1):0.1:A2(end));

DPf2_1600i = griddata(A2,A1,(abs(DPf2_1600)),A2i,A1i,'cubic');
%c1 = contour(A1mq,A2mq,interp2(A1m,A2m,20*log10(abs(DPoaeFT2kg1(:,6:end))),A1mq,A2mq,'cubic'),30,'linewidth',LW);

haF1G1 = axes('position',[leftCornerX leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add'); 
contourf(A2i,A1i,20*log10(abs(DPf2_1600i)),30); colormap gray;
plot(A2,50*ones(size(A2)),'k-','linewidth',LW-0.2,'markersize',12,'color','r');
plot(A2,58*ones(size(A2)),'k--','linewidth',LW-0.2,'markersize',12,'color','r');
plot(48,58,'go','linewidth',LW-0.2,'markersize',8);
plot(50*ones(size(A1)),A1,'r-','linewidth',LW-0.2,'markersize',12,'color','k');
plot(56*ones(size(A1)),A1,'r--','linewidth',LW-0.2,'markersize',12,'color','k');
plot(52,58,'mo','linewidth',LW-0.2,'markersize',8);
plot(58,58,'bo','linewidth',LW-0.2,'markersize',8);

set(haF1G1,...
    'xgrid','on',...
    'ygrid','on',...
    'zgrid','on',...
    'xlim',Xlims,...
    'ylim',Ylims,...
    'clim',Clims,...
    'box','on',...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );

colorbar('position',[leftCornerX+0.06 leftCornerY+hEnv+0.01 wEnv-0.06 0.02],'orientation','horizontal');

text(12,50,'{\itL}_1 (dB SPL)','fontsize',iFontSize,'rotation',90,'fontname',strFontName );

text(15,74,'{\itL}_{\rmDP}\newline(dB re 1 a.u.)','fontsize',iFontSize-3,'rotation',0,'fontname',strFontName);

text(22,66,'A','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'backgroundcolor','white','fontweight','bold');
text(28,66,'{\itf}_2 = 1.6 kHz','fontsize',iFontSize-3,'rotation',0,'fontname',strFontName,'backgroundcolor','white');

xlabel('{\itL}_2 (dB SPL)','fontsize',iFontSize,'fontname',strFontName)

haF1A1 = axes('position',[leftCornerX+wEnv+xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add');

idxA2_56 = find(A2>=56,1,'first');
idxA2_50 = find(A2>=50,1,'first');
idxA1_58 = find(A1>=58,1,'first');
idxA1_50 = find(A1>=50,1,'first');
LW = 1.2;
MS = 8;
plot(A1,20*log10(abs(DPf2_1600(:,idxA2_50))),'k-','linewidth',LW,'markersize',MS);
plot(A1,20*log10(abs(DPf2_1600(:,idxA2_56))),'k--','linewidth',LW,'markersize',MS);

plot(A2,20*log10(abs(DPf2_1600(idxA1_50,:))),'r-','linewidth',LW,'markersize',MS);
plot(A2,20*log10(abs(DPf2_1600(idxA1_58,:))),'r--','linewidth',LW,'markersize',MS);

plot(48,20*log10(abs(DPf2_1600(idxA1_58,15))),'go','linewidth',LW-0.2,'markersize',8);
plot(52,20*log10(abs(DPf2_1600(idxA1_58,17))),'mo','linewidth',LW-0.2,'markersize',8);
plot(58,20*log10(abs(DPf2_1600(idxA1_58,20))),'bo','linewidth',LW-0.2,'markersize',8);


cycle = 2*pi;

text(22,35,'Amplitude (dB re 1 a.u.)','fontsize',iFontSize,'rotation',90,'fontname',strFontName,'horizontalalignment','center' );
text(32,47,'B','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'fontweight','bold');
xlabel('{\itL}_1 or {\itL}_2 (dB SPL)','fontsize',iFontSize,'fontname',strFontName)
XlimsA1 = [30 70];
set(haF1A1,...
    'box','on',...
    'xlim',XlimsA1,...
    'ticklength',[0.02 0.01],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );


haF1Ph1 = axes('position',[leftCornerX+2*wEnv+2*xD leftCornerY  wEnv hEnv],'xscale','lin','yscale','lin','nextplot','add');
plot(A1,unwrap(angle(DPf2_1600(:,idxA2_50)))/cycle,'k-','linewidth',LW,'markersize',MS);
plot(A1,unwrap(angle(DPf2_1600(:,idxA2_56)))/cycle,'k--','linewidth',LW,'markersize',MS);

plot(A2,unwrap(angle(DPf2_1600(idxA1_50,:)))/cycle,'r-','linewidth',LW,'markersize',MS);
plot(A2,unwrap(angle(DPf2_1600(idxA1_58,:)))/cycle,'r--','linewidth',LW,'markersize',MS);
plot(56,unwrap(angle(DPf2_1600(find(A1==56,1,'first'),idxA2_50)))/cycle,'k*','linewidth',LW,'markersize',MS);
plot(62,unwrap(angle(DPf2_1600(find(A1==62,1,'first'),idxA2_56)))/cycle,'k*','linewidth',LW,'markersize',MS);
plot(46,unwrap(angle(DPf2_1600(idxA1_50,find(A2==46,1,'first'))))/cycle,'r*','linewidth',LW,'markersize',MS);
plot(52,unwrap(angle(DPf2_1600(idxA1_58,find(A2==52,1,'first'))))/cycle,'r*','linewidth',LW,'markersize',MS);

xlabel('{\itL}_1 or {\itL}_2 (dB SPL)','fontsize',iFontSize,'fontname',strFontName)
text(22,0,'Phase (cycles)','fontsize',iFontSize,'rotation',90,'fontname',strFontName,'horizontalalignment','center' );
set(haF1Ph1,...
    'box','on',...
    'ylim',[-1.2 1.2],...
    'ticklength',[0.02 0.01],...
    'xlim',XlimsA1,...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );
text(32,1,'C','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'fontweight','bold');


leg1 = legend('{\itL}_2 = 50 dB, {\itL}_1 varied','{\itL}_2 = 56 dB,  {\itL}_1 varied', ...
    '{\itL}_1 = 50 dB, {\itL}_2 varied', '{\itL}_1 = 58 dB , {\itL}_2 varied');

set(leg1,'Interpreter','tex','fontname',strFontName,'fontsize',iFontSize-4,'position',...
    [leftCornerX+2*wEnv+0.6*xD leftCornerY+0.278 0.15 0.08],'NumColumns',2);

  

load([FolderRes2 'JASAELdpL1L2mapF1_1400Hz_F2_1600Hz_gain105R21.mat'])
hEnv2 = 0.22;
fs = 600e3;
Tan = 19e-3; % start of analysis window
Nonset = round(Tan*fs);
RR1{1,1} = RR1{1,1}.*exp(-sqrt(-1)*2*pi*F1*Nonset*1/fs);
RR1{1,2} = RR1{1,2}.*exp(-sqrt(-1)*2*pi*F1*Nonset*1/fs);
RR1{1,3} = RR1{1,3}.*exp(-sqrt(-1)*2*pi*F1*Nonset*1/fs);
RR2{1,1} = RR2{1,1}.*exp(-sqrt(-1)*2*pi*F2*Nonset*1/fs);
RR2{1,2} = RR2{1,2}.*exp(-sqrt(-1)*2*pi*F2*Nonset*1/fs);
RR2{1,3} = RR2{1,3}.*exp(-sqrt(-1)*2*pi*F2*Nonset*1/fs);
RRdp{1,1} = RRdp{1,1}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);
RRdp{1,2} = RRdp{1,2}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);
RRdp{1,3} = RRdp{1,3}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);


UUrdp{1,1} = UUrdp{1,1}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);
UUrdp{1,2} = UUrdp{1,2}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);
UUrdp{1,3} = UUrdp{1,3}.*exp(-sqrt(-1)*2*pi*Fdp*Nonset*1/fs);

haTWa = axes('position',[leftCornerX leftCornerY-hEnv-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 
plot(x,20*log10(abs(RR1{1,1})/0.001),'g','linewidth',LW)
plot(x,20*log10(abs(RR2{1,1})/0.001),'g--','linewidth',LW)
plot(x,20*log10(abs(RRdp{1,1})/0.001),'g:','linewidth',LW)
plot(x,20*log10(abs(RR1{1,2})/0.001),'m','linewidth',LW)
plot(x,20*log10(abs(RR2{1,2})/0.001),'m--','linewidth',LW)
plot(x,20*log10(abs(RRdp{1,2})/0.001),'m:','linewidth',LW)
plot(x,20*log10(abs(RR1{1,3})/0.001),'b','linewidth',LW)
plot(x,20*log10(abs(RR2{1,3})/0.001),'b--','linewidth',LW)
plot(x,20*log10(abs(RRdp{1,3})/0.001),'b:','linewidth',LW)


text(1.8,27,'{\it f_{\rmDP}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')
text(1.6,39,'{\it f_{\rm1}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')
text(1.345,25,'{\it f_{\rm2}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')
plot([1.43 1.43],[10, 100],'-','color',[.5 .5 .5])
%zlabel('DPOAE amplitude (a.u.)','fontsize',iFontSize+4,'fontname',strFontName)
%view([-17 59])
Xlims = [1,2]
Ylims = [15,40]
ylabel('Ampl. (dB re)','fontsize',iFontSize,'fontname',strFontName)
set(haTWa,...   
    'xlim',Xlims,...
    'ylim',Ylims,...
    'yscale','lin',...
    'xticklabel','',...
    'ticklength',[0.02 0.01],...
    'box','on',...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );
text(1.05,37,'D','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'fontweight','bold');

haTWph = axes('position',[leftCornerX leftCornerY-hEnv-hEnv2-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 
cycle = 2*pi
plot(x,unwrap(angle(RR1{1,1}))/cycle+1,'g','linewidth',LW)
plot(x,unwrap(angle(RR2{1,1}))/cycle+1,'g--','linewidth',LW)
plot(x,unwrap(angle(RRdp{1,1}))/cycle,'g:','linewidth',LW)
plot(x,unwrap(angle(RR1{1,2}))/cycle+1,'m-','linewidth',LW)
plot(x,unwrap(angle(RR2{1,2}))/cycle+1,'m--','linewidth',LW)
plot(x,unwrap(angle(RRdp{1,2}))/cycle,'m:','linewidth',LW)
plot(x,unwrap(angle(RR1{1,3}))/cycle+1,'b','linewidth',LW)
plot(x,unwrap(angle(RR2{1,3}))/cycle+1,'b--','linewidth',LW)
plot(x,unwrap(angle(RRdp{1,3}))/cycle,'b:','linewidth',LW)

text(1.55,0.5,'{\it f_{\rmDP}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')
text(1.61,-1,'{\it f_{\rm1}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')
text(1.38,-1,'{\it f_{\rm2}}','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top')

ylabel('Phase (cycles)','fontsize',iFontSize,'fontname',strFontName)
xlabel('{\it x} (cm)','fontsize',iFontSize,'fontname',strFontName)
YlimsPh = [-3,0.99]
set(haTWph,...   
    'xlim',Xlims,...
    'ylim',YlimsPh,...
    'yscale','lin',...
    'box','on',...
    'ticklength',[0.02 0.01],...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );

haUnla = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 
text(1.78,88,'$\hat{U}_{\omega_{\rm DP}}^{\rm NL}(x)$','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal',...
    'verticalalignment','top','interpreter','LaTeX')
%text(0.82,2.8e5,'$\hat{U}_{n_{DP}}^N$','fontsize',iFontSize,'fontname',strFontName,'fontweight','normal','rotation',0,'horizontalalignment','center','interpreter','latex')
plot(x,20*log10(abs(UUrdp{1,1})),'g','linewidth',LW)
plot(x,20*log10(abs(UUrdp{1,2})),'m','linewidth',LW)
plot(x,20*log10(abs(UUrdp{1,3})),'b','linewidth',LW)
plot([1.43 1.43],[60, 100],'-','color',[.5 .5 .5])
text(1.05,86,'E','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'fontweight','bold');
%zlabel('DPOAE amplitude (a.u.)','fontsize',iFontSize+4,'fontname',strFontName)
%view([-17 59])
Xlims = [1,2]
Ylims = [60,90]
ylabel('Ampl. (dB)','fontsize',iFontSize,'fontname',strFontName)
set(haUnla,...   
    'xlim',Xlims,...
    'ylim',Ylims,...
    'yscale','lin',...
    'xticklabel','',...
    'ticklength',[0.02 0.01],...
    'box','on',...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );


haUnlph = axes('position',[leftCornerX+wEnv+xD leftCornerY-hEnv-hEnv2-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 
cycle = 2*pi
plot(x,unwrap(angle(UUrdp{1,1}))/cycle,'g','linewidth',LW)
plot(x,unwrap(angle(UUrdp{1,2}))/cycle,'m','linewidth',LW)
plot(x,unwrap(angle(UUrdp{1,3}))/cycle,'b','linewidth',LW)
ylabel('Phase (cycles)','fontsize',iFontSize,'fontname',strFontName)
xlabel('{\it x} (cm)','fontsize',iFontSize,'fontname',strFontName)
YlimsPh = [-3,0.99]
set(haUnlph,...   
    'xlim',Xlims,...
    'ylim',YlimsPh,...
    'yscale','lin',...
    'box','on',...
    'ticklength',[0.02 0.01],...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );


haDPa = axes('position',[leftCornerX+2*wEnv+2*xD leftCornerY-hEnv-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 

addpath('BoxModel')
Y = FreqDomainA_MeTal(Fdp,0,0,0,1.05); % freq domain solution for base wave (TW for linear model at Fdp frequency)
rmpath('BoxModel')

text(1.05,17.6,'F','fontsize',iFontSize+2,'rotation',0,'fontname',strFontName,'fontweight','bold');
plot(x,(fliplr(abs(cumsum(fliplr(UUrdp{1,1}).*flipud(Y).')))),'g','linewidth',LW)
plot(x,(fliplr(abs(cumsum(fliplr(UUrdp{1,2}).*flipud(Y).')))),'m','linewidth',LW)
plot(x,(fliplr(abs(cumsum(fliplr(UUrdp{1,3}).*flipud(Y).')))),'b','linewidth',LW)
plot([1.43 1.43],[0, 100],'-','color',[.5 .5 .5])

Xlims = [1,2]
Ylims = [0,20]
ylabel('Ampl. (a.u.)','fontsize',iFontSize,'fontname',strFontName)
set(haDPa,...   
    'xlim',Xlims,...
    'ylim',Ylims,...
    'ticklength',[0.02 0.01],...
    'yscale','lin',...
    'xticklabel','',...
    'box','on',...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );
text(1.48,17.5,'$\int_L^x \hat{U}_{\omega_{\rm DP}}^{\rm NL}(s)\xi^{(2)}_{\omega_{\rm DP}}(s)ds$','fontsize',iFontSize-2,'fontname',strFontName,'fontweight','normal','interpreter','latex')

haDPph = axes('position',[leftCornerX+2*wEnv+2*xD leftCornerY-hEnv-hEnv2-yD  wEnv hEnv2],'xscale','lin','yscale','lin','nextplot','add'); 
cycle = 2*pi;


plot(x,(fliplr(unwrap(angle(cumsum(fliplr(UUrdp{1,1}).*flipud(Y).')))/cycle))-17,'g','linewidth',LW)
plot(x,(fliplr(unwrap(angle(cumsum(fliplr(UUrdp{1,2}).*flipud(Y).')))/cycle))-18,'m','linewidth',LW)
plot(x,(fliplr(unwrap(angle(cumsum(fliplr(UUrdp{1,3}).*flipud(Y).')))/cycle))-16,'b','linewidth',LW)
ylabel('Phase (cycles)','fontsize',iFontSize,'fontname',strFontName)
xlabel('{\it x} (cm)','fontsize',iFontSize,'fontname',strFontName)
YlimsPh = [-3,0.99]
set(haDPph,...   
    'xlim',Xlims,...    
    'ylim',YlimsPh,...    
    'yscale','lin',...
    'box','on',...
    'xtick',[1,1.2,1.4,1.6,1.8,2],...
    'ticklength',[0.02 0.01],...
    'fontsize',iFontSize,...
    'fontname',strFontName...
   );


iResolution = 600;  
strFileName = ['../Results/Figures/MOH20DPOAE/model'];
print('-depsc', strcat(strFileName, '.eps'));
print('-dpng', sprintf('-r%d', iResolution), strcat(strFileName, '.png'));
print('-djpeg', sprintf('-r%d', iResolution), strcat(strFileName, '.jpg'));

