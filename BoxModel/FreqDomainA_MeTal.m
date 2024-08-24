function [Y sigma_ow Pd0]=FreqDomainA_MeTal(F,Inp,plotflag,compare,gain);
%function [Y sigma_ow]=FreqDomainA_MeTal(F,Inp,plotflag,compare,gain);
%
% Cochlear model solution in frequency domain
% ME + Cochlea
%
% F: Frequency in Hz
% Inp: Level in dB
% plotflag: 1 - plot figure, 0 - no plot
% compare: 1 - compare with result without ME, 0 - no comparison
% gain: gain of cochlear amplifier
%
% Y: BM displacement
% sigma_ow: OW displacement
% Pd0: pressure difference near stapes (input to the cochlea)
%
% Vencovsky, Vetesnik 2017

 if nargin<5, gain = 1.05; end
 if nargin<4,  compare=1; end % 1 - compare TW without ME and with ME
 if nargin <3, plotflag=1; end
 if nargin <2, Inp=10; end
 if nargin <1, F=4000; end

  N=800;
  rebuild_flag=1; 
  
  [x,Gs,G,M,stiff,DampSp,undamp,bigamma,wm2,Qbm,Qs,mME,kME,hME,Gme,Sty,...
    gammaAir,Ve,GammaMi,Sow,W,Pea] = AlldataFme(N,rebuild_flag,gain);
  
  
  Gs=Gs(:);
  %U=zeros(N,1);U(IND-20:IND+20)=AM;
%  U=Gs*AM;
  I=sqrt(-1);omega=2*pi*F;omega_2=omega*omega;
  AM=db2inputME(Inp);%*omega^2;
  Bh=DampSp./omega;
  Bk=diag(stiff./omega_2);
  %-----------------------------
  % Vetesnik Gummer 2012
  %Gme=Gs*G(1,:);
%   fac=1;
%   Ke=fac*3.2400e+009; %dyn/cm^5
%   Re=fac*4e5;   %dyn-s/cm^5;
%   Me=fac*17.6;   %4.4 g/cm^4; 
%   As=0.03;
%   
  %M3=As*(Ke/omega_2+I*Re/omega-Me)-Gs(1);%
  %Gmew=(1/M3)*Gme/100;

  % Vencovsky Vetesnik 2017
  % middle-ear model according to Talmadge
  Zme = mME + Sow*Qs/W(1) - I*hME/omega-kME/omega_2; % middle ear
  
  %MhMe = Gs/(mOW+Gs(0));
  %Gmew=MhMe*G(1,:);
  %GME = Sow*Gs*Qbm/(W(1)*Zme);
  
  %undamp(300:310) = 0;
  Th=I*omega*bigamma;Tk=wm2;TM=diag(undamp./(omega_2-Th-wm2));
  
  %A= M + G - GME - I*Bh - Bk + TM;
  
  A=G+M-Sow*Gs*Qbm/(W(1)*Zme) - I*Bh-Bk+TM;
  
  Input = Gs*Sow*Gme*AM/(Zme*omega_2);
  
 % Input = Gs;
  Y=A\Input;
    
  sigma_ow = (-Sow*Qbm*Y/W(1)-Sow*Gme*AM/omega_2)/Zme;
  %IsoIntY = Y/sigma_ow;
  
  Pd0 = Qbm*Y*omega_2/W+Gs(1)*sigma_ow*omega_2/W;
  
 if compare==1 % compare with the model without ME
  Aw= M + G - I*Bh - Bk + TM;
  Uw = Gs(:)*AM;
  Yw = -Aw\Uw;
 end
  
 if plotflag==1
    
     if compare==0
       cycle=2*pi;
       Ya=abs(Y); Yph=unwrap(angle(Y))/cycle;
       
        tm=diag(TM);
      TMa=abs(tm);TMph=unwrap(angle(tm))/pi;
      
      Ap=G+M-Sow*Gs*Qbm/(W(1)*Zme) - I*Bh-Bk;
  
      Input = Gs*Sow*Gme*AM/(Zme*omega_2);
  
 % Input = Gs;
      Yp=Ap\Input;
      Ypa=abs(Yp); Ypph=unwrap(angle(Yp))/cycle;
       
%        figure(1)
%        clf
%        subplot(2,1,1)
%        plot(x,20*log10(Ya),'r')   
%        plot(x,20*log10(Ya),'r')   
%        
%   %       ylim([-100 -30])
%        xlim([0 3.5])
%        ylabel('Amplitude')
%        title(['Input-> frequency: ',num2str(F),' [Hz];',' amplitude: ', num2str(Inp), ' [dB]']) 
%        grid on
%        subplot(2,1,2)
%        plot(x,Yph,'r')
%        
       
       figure(1)
       clf
       subplot(2,1,1)
         plot(x,Ya/max(Ya),'r',x,TMa/max(TMa),'b',x,Ypa/max(Ypa),'k')
         legend('BM displacement-active','TM feedback','BM displacement-passive')
         ylabel('Amplitude')
         title(['Input-> frequency: ',num2str(F),' [Hz];',' amplitude: ', num2str(Inp), ' [dB]']) 
         grid on
       subplot(2,1,2)
        plot(x,Yph,'r',x,TMph,'b')
        xlabel('Distance from stapes [cm]')
        ylabel('Phase[\pi]')
        grid on
       
       
    %   xlabel('Distance from stapes [cm]')
  %     ylabel('Phase(cycle)')
    %   grid on
     else
       cycle=2*pi;
       Ya=abs(Y); Yph=unwrap(angle(Y))/cycle;
       Ywa=abs(Yw); Ywph=unwrap(angle(Yw))/cycle;
       
        
     
       
       figure(1)
       clf
       subplot(2,1,1)
       hold on;
       plot(x,20*log10(Ya),'r')         
       plot(x,20*log10(Ywa),'g--');
       hold off;
  %       ylim([-100 -30])
       xlim([0 3.5])
       ylabel('Amplitude')
       title(['Input-> frequency: ',num2str(F),' [Hz];',' amplitude: ', num2str(Inp), ' [dB]']) 
       grid on
       subplot(2,1,2)
       hold on;
       plot(x,Yph,'r')
       plot(x,Ywph,'g--')
       hold off;
       xlabel('Distance from stapes [cm]')
       ylabel('Phase(cycle)')
       grid on
     
     end
 end
   if nargout<1,Y=[];end