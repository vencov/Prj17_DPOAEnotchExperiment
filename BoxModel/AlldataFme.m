function [x,Gs,G,M,stiff,ShSp,undamp,bigamma,wm2,Qbm,Qs,mME,kME,hME,Gme,Sty,...
    gammaAir,Ve,GammaMi,Sow,W,Pea] = AlldataFme(N,rebuild_flag,AmF);
 if nargin <3,AmF=1.3;end
 if nargin <2,rebuild_flag=0;end
 if nargin <1,N=800;end
 if rebuild_flag==1
    delete bmdataME.mat
    delete TMDATA.MAT
 end
%________________________BM data_______________________________
 if exist('bmdataME.mat','file')==2,
	load bmdataME.mat x N M stiff damp0 ShSp G Gs Qbm Qs mME kME hME Gme Sty gammaAir ...
      Ve GammaMi Sow W Pea
 else,
	disp(' File bmdataME.mat does not exist. It will be created!');
    bmdataME(N); 
    disp('	Now bmdataME.mat does exist');
    load bmdataME.mat x N M stiff damp0 ShSp G Gs Qbm Qs mME kME hME Gme Sty gammaAir ...
      Ve GammaMi Sow W Pea
 end
%________________________TM data________________________________
if exist('TMDATA.MAT','file')==2,
	load TMDATA.MAT bigamma wm2 -mat
else
    disp('File TMDATA.MAT does not exist. It will be created!');
    tmdata(x);    
    disp('Now TMDATA.MAT does exist.');
    load TMDATA.MAT bigamma wm2 -mat
end

  shapeOn = hann(100);
  
  shapeAll =[shapeOn(1:50); ones(N-100,1); shapeOn(51:end)];

  Xzero = 1.0; % point where shaping function is 0
Sc = 6; % scale x axis for tanh - defines the shape of tanh function
%ScAmp = 0.1; % scale amplitude of tanh function
ScAmp = 0.12;
ScAmp = 0.15;

%UmpBas = 1.04+ScAmp*(tanh((Sc*(Xzero-x)/max(x)))); % g2
%UmpBas = 1.02+ScAmp*(tanh((Sc*(Xzero-x)/max(x)))); % g3
%UmpBas = 1.02+ScAmp*(tanh((Sc*(Xzero-x)/max(x)))); % g4
UmpBas = 1.02+ScAmp*(tanh((Sc*(Xzero-x)/max(x)))); % g2

  
  %undamp=1.34*AmF*damp0(:).*bigamma.*shapeAll;
undamp=1.34*AmF*damp0(:).*bigamma.*UmpBas';


 
   



