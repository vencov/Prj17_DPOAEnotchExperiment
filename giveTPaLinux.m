function TruePa = giveTruePascals(F,Rec,gain)
%function TruePa = giveTrueSPL(F,Rec)
% function which takes the calibration curve of the microphone and
% calculates pressure from the signal measured by the microphone and stored
% in Matlab
%

TruePa = Rec/(0.003*(10^(gain/20)));
