function input = db2inputME(dB);

%	rescales dB SPL into corresponding input values for human cochlea  


AMo=3.5e-6;
%AMo = 134;
AMo = 424;
input = sqrt(2)*AMo*10.^(dB/20);
