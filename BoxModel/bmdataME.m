function bmdataME(N); 


 if nargin <1, N=800; end  

   H = 0.1;  %height (cm) 
   L = 3.5;  %length (cm)
   W = 1;%0.029; % width of the BM (cm)
   Sow = H*W;%0.032; % (cm^2)
   rho = 1;  %water density (g/cm^3)
%    H = H/4;
   [G,Gs,x,Sh]= cgsgreenf(N, L, H);%Green's function
   Sh=2e-4*Sh;

  %Gs=Sow*Gs(:);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------- MASS -------------------------------------------------
  %m0 = 0.05;
  m0 = 1e-4; % g/cm
  M = diag(m0*ones(size(x)));
%---------------------TONOPIC MAP FOR MAN --------------------------
  cf =GetHumanCF(x,L,1,0);
%-------------------- STIFFNESS --------------------------------------
  stiff=(H/3+m0)*((2*pi*cf).^2);
  %stiff=(m0)*((2*pi*cf).^2);
%-------------------  DAMPING ----------------------------------------
  h0= 6000;
  fact=10;%????
  damp0=sqrt(stiff);damp0=h0*(damp0./damp0(1));
  damp=damp0+fact*damp0(N)*exp(x-L);
  Ce=80*damp(N);	% increased damping near appex (VV)			
  damp = damp + Ce*exp((x/L-1)/0.055);  % VV

  ShSp=sparse(diag(damp)) + sparse(Sh);
 
 % ME PARAMETERS
 
%  mME = 0.059; % stapes mass (g)
%  wME = 1500*2*pi; % Hz*(2*pi) middle ear frequency 
%  kME = (wME^2)*mME; % stiffness of ME
%  gammaME = 500; % s^-1 middle ear damping constant
% mME = 7.08; % Mass of the oval window [g] % for Sow = 1 and Sty = 1;
% hME = 100300; %20*500*mME; % damping of the ME [g/s]
% kME = 1.3102e+09; % stiffness of the ME [g/s^2]
mME = 0.0531; % Mass of the oval window [g]
hME = 1298; %20*500*mME; % damping of the ME [g/s]
kME = 6.8130e+06; % stiffness of the ME [g/s^2]
%kME = 6.8130e+05/7; % stiffness of the ME [g/s^2]

 Pea = 1e6; % adiabatic pressue in the ear canal [dyn/m^2] (assumed equal to athmospheric pressure)
 
 
 %Sty = 1; %0.49; % cm^2 effective area of tympanic membrane
 Sty = 0.49;
 gammaAir = 1.4; % ratio of specific heats of air
 Ve = 0.5; % cm^3 volume of tympanic cavity
 GammaMi = 1.4; % level ratio of incus
 Gme = Sty/Sow*GammaMi; %21.4; % mechanical gain of ossicles 
 %-------------TIME DOMAIN--------------------
 Qbm=G(1,:);
 Qs=Gs(1); 
   
 plotflag = 0;
                                    % Unit dimension of Sh is Kg/(BETA*sec) 

%%%%%%%%%%%%%  PLOTS %%%%%%%%%%%%%%%%%%%
if plotflag==1,
	DX=L/N;     % dx length in meters. As G is in Kg/BETA, G/dx is in Kg/BETA^2  
	Gfactor=1/(DX); 
	figure;
	subplot(211), 
	plot(x, G(1:20:N, :)*Gfactor,'k-'),
	title('Samples of BM-BM Green''s function [g/cm^2]', 'Interpreter', 'tex'),
	subplot(212),
	plot(x, Gs, 'k-'),
	title('Stapes-BM fluid coupling  [g/cm]')
	xlabel('Fractional distance from stapes');

    figure;
    plot(x,stiff*1e-3/1e-2);
    figure;
    plot(x,damp*1e-3/1e-2)
end
  

 save bmdataME.mat x N M stiff damp0 ShSp G Gs Qbm Qs mME kME hME Gme Sty gammaAir ...
      Ve GammaMi Sow W Pea