function TruePa = giveTrueSPL(F,Rec,gain,fs,Coae,fxC)
%function TruePa = giveTrueSPL(F,Rec)
% function which takes the calibration curve of the microphone and
% calculates pressure from the signal measured by the microphone and stored
% in Matlab
%
% F - frequency in Hz
% Rec - measured value

% deal with calibration

% load microphone data 

%load Calibration/MicRefSensit.mat Hoaemicsens;
load Calibration/CalibMicData.mat Hoaemicsens;

fxI = 200:1:fs/2; %

fx = (0:length(Hoaemicsens)-1)*fs/length(Hoaemicsens); % frequency axis
% 
% Fmin = 100; % frequency range for estimation
% Fmax = 14e3; % frequency range for estimation
% 
% FminId = find(fx>=Fmin,1,'first');
% FmaxId = find(fx>=Fmax,1,'first');
% 
% fxI = fx(FminId:FmaxId);


%CalCurve = smooth(real(Hoaemicsens),50,'loess')+sqrt(-1)*smooth(imag(Hoaemicsens),50,'loess');
%CalCurveI = interp1(fx,real(CalCurve),fxI,'pchip') + sqrt(-1)*interp1(fx,imag(CalCurve),fxI,'pchip');
CalCurveI = interp1(fx,Hoaemicsens,fxI,'pchip'); % + sqrt(-1)*interp1(fx,imag(CalCurve),fxI,'pchip');
if length(Coae)>1
    %CoaeI = interp1(fxC,real(Coae),fxI,'pchip') + sqrt(-1)*interp1(fxC,imag(Coae),fxI,'pchip');
    CoaeI = interp1(fxC,Coae,fxI,'pchip');
else
    CoaeI = ones(size(fxI));
end


TruePa = zeros(size(F));

for k=1:length(F)
    [~, idxF] = min(abs(fxI-F(k))); % find index corresponding to the frequency
    TruePa(k) = CoaeI(idxF).*(Rec(k)/(gain*CalCurveI(idxF)));
end
