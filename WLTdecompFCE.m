function [fx, DPgAlla, DPgNLa, DPgCRa, DPgCR1a, DPgAllreca, DPgAllph, DPgNLph, DPgCRph, DPgCR1ph, DPgAllrecph, ...
    fc, DPtoctAlla, DPtoctNLa, DPtoctCRa, wlgraph] = WLTdecompFCE(params,freqf2,DPgram,plotflag)

%plotflag = 1;

nmm=params.nmm;
nb=params.nb; % n. 1/3 oct bands
tt=params.tt; % latency power law: tt*f(kHz)^-bbb
bbb=params.bbb;
p0=params.p0;
tc=params.tc;
tch=params.tch;
c1=params.c1;
c2=params.c2;
c3=params.c3;
c4=params.c4;
c5=params.c5;
fDo=freqf2;
Ao=size(fDo);
NDo=Ao(2);
DPgramsm=movmean(DPgram,10);
fDomax=fDo(end);
fDomin=fDo(1);
dfD=10;
ND=round((fDomax-fDomin)/dfD);
DPmm=zeros(ND,1);
fD=fDomin+dfD*(1:ND);
for jjj=1:ND
for iii=1:NDo-1
DPmm(jjj) = DPmm(jjj) + DPgramsm(iii)*(fDo(iii)<fD(jjj))*(fDo(iii+1)>fD(jjj));
end
end
M=2048;  
%M=1024;
%M=2*2048;
Fs=dfD*M; % dummy sampling frequency
dfDk=dfD/1000;
ts = 1/Fs; % Sampling time
N=(M/2);
NW=N;
tsm=1000*ts;
NZ=ceil(fD(1)/dfD)-1;
fmin=(NZ+1)*dfD;
for kk=1:nb+1
fi(kk)=round(fmin/dfD*(2^((kk-1)/3)));
end
fc(1)=fmin*(2^(1/3)-1)/0.23;
for kk=1:nb-1
fc(kk+1)=fc(kk)*2^(1/3);  % center frequencies of 1/3 oct bands
end 
segn1 = zeros(M,1);
segn2 = zeros(M,1);
segn3 = zeros(M,1);
segn4 = zeros(M,1);
tf_coewfil1= zeros(M,NW);
tf_coewfil2= zeros(M,NW);
tf_coewfil3= zeros(M,NW);
tf_coewfil4= zeros(M,NW);
coew=zeros(M,NW);
coewd=zeros(M,NW);
wbtf=zeros(M,NW);
coewfil1= zeros(M,NW);
coewfil2= zeros(M,NW);
coewfil3= zeros(M,NW);
coewfil4= zeros(M,NW);
segnw1tf= zeros(M,NW);
segnw2tf= zeros(M,NW);
segnw3tf= zeros(M,NW);
segnw4tf= zeros(M,NW);
wavb = zeros(M,NW);  
for j=1:NZ
       f(j) = j*dfD;
       signal(j)=0;
end
    for j=NZ+1:NZ+ND  
       f(j)= j*dfD; 
       signal(j) = DPmm(j-NZ)/p0;
    end 
    for j=NZ+ND+1:M
        f(j)= j*dfD;
    end
    for j=NZ+ND+1:N
        signal(j)=0;
    end
Y=signal;
for j=1:N
Ys(j) = conj(Y(mod(N-j+1,N)+1));
end
for i=N+1:M
Y(i) = Ys(i-N);
end
      
    alpha=2*pi*dfD;
    a3=0.075*alpha;
    beta=4;
    dfw= double(alpha/2/pi/1000); 
    dfwH=dfw*1000;
for k= 1:NW   
for j=1:M
wavb(j,k)=(1/(1+(a3*k*(j-M/2)*ts)^beta))*cos(k*alpha*(j-M/2)*ts); % wavelet basis functions
end
end
for k=1:NW 
wbtf(:,k)=fft(wavb(:,k)); %wavelet FTs 
end
for k=1:NW
coewtf(:,k)=(k^0.5)*wbtf(:,k).*Y(:); 
end
for k=1:NW
coew(:,k)=ifft(coewtf(:,k));% wavelet coefficients(t,f)
end
% wavelet filtering
for k=1:NW
    for j=1:M
        coewfil1(j,k)=coew(j,k).*((j-M/2)*tsm>c1*tt*(k*dfw)^-bbb).*((j-M/2)*tsm<c2*tt*(k*dfw)^-bbb);
        coewfil2(j,k)=coew(j,k).*((j-M/2)*tsm>c2*tt*(k*dfw)^-bbb).*((j-M/2)*tsm<c3*tt*(k*dfw)^-bbb);
        coewfil3(j,k)=coew(j,k);
        coewfil4(j,k)=coew(j,k).*((j-M/2)*tsm>c2*tt*(k*dfw)^-bbb).*((j-M/2)*tsm<c4*tt*(k*dfw)^-bbb);        
    end
end
for k=1:NW 
tf_coewfil1(:,k)= fft(coewfil1(:,k)); % filtered wavelet coeffs FTs
tf_coewfil2(:,k)= fft(coewfil2(:,k));
tf_coewfil3(:,k)= fft(coewfil3(:,k));
tf_coewfil4(:,k)= fft(coewfil4(:,k));
end
%%%%%%%  filtered signal reconstuction  %%%%%%
for j=1:NW
    segnw1tf(:,j)= (j^0.5)*(tf_coewfil1(:,j).* wbtf(:,j)); 
    segnw2tf(:,j)= (j^0.5)*(tf_coewfil2(:,j).* wbtf(:,j)); 
    segnw3tf(:,j)= (j^0.5)*(tf_coewfil3(:,j).* wbtf(:,j)); 
    segnw4tf(:,j)= (j^0.5)*(tf_coewfil4(:,j).* wbtf(:,j));     
end 
for k=1:NW
    segn1(:,1) = segn1(:,1) + segnw1tf(:,k); % distortion (ZL) component
    segn2(:,1) = segn2(:,1) + segnw2tf(:,k); % basal first reflection (SL) component
    segn3(:,1) = segn3(:,1) + segnw3tf(:,k); % whole signal reconstructed (for normalization purposes)
    segn4(:,1) = segn4(:,1) + segnw4tf(:,k); % total first reflection (LL) component
end
 %%%%%%%%%%%%%%%%%%%%%%% signal in 1/3 oct bands
for kk=1:nb
Ym(kk)=sum(abs(Y(1,fi(kk):fi(kk+1)-1)).^2);  
segn1m(kk)=sum(abs(segn1(fi(kk):fi(kk+1)-1,1)).^2); 
segn2m(kk)=sum(abs(segn2(fi(kk):fi(kk+1)-1,1)).^2);
segn4m(kk)=sum(abs(segn4(fi(kk):fi(kk+1)-1,1)).^2);
ndd(kk)=fi(kk+1)-fi(kk);
end
 Ym = 10*log10(Ym./ndd);
 segn1m = 10*log10(segn1m./ndd);
 segn2m = 10*log10(segn2m./ndd);
 segn4m = 10*log10(segn4m./ndd);
for k=1:NW
span=2*round(3*NW/k)+1;
coewd(:,k)=sqrt(smooth(coew(:,k).^2,span));
cmk(k)=max(coewd(:,k));
end
cm=max(cmk);
c=['b','r','g'];
imed=round((fDomax+fDomin)/dfD/2)
fac=(20*log10(abs(segn3(imed))/abs(Y(imed))));


% data out

fx = (NZ+1:NZ+ND)*dfD; % frequency vector for detailed DP grams
DPgAlla = 20*log10(abs(Y(NZ+1:NZ+ND))); % DPOAE
DPgAllreca = 20*log10(abs(segn3(NZ+1:NZ+ND)))-fac; % DPOAE reconstructed
DPgNLa = 20*log10(abs(segn1(NZ+1:NZ+ND)))-fac; % NL component
DPgCRa = 20*log10(abs(segn4(NZ+1:NZ+ND)))-fac; % CR component
DPgCR1a = 20*log10(abs(segn2(NZ+1:NZ+ND)))-fac; % basal reflection (SL)

DPgAllph = unwrap(angle(Y(NZ+1:NZ+ND)))/(2*pi); % phases 
DPgAllrecph = unwrap(angle(segn3(NZ+1:NZ+ND)))/(2*pi);
DPgNLph = unwrap(angle(segn1(NZ+1:NZ+ND)))/(2*pi);
DPgCRph = unwrap(angle(segn4(NZ+1:NZ+ND)))/(2*pi);
DPgNLph = angle(segn1(NZ+1:NZ+ND));
DPgCRph = angle(segn4(NZ+1:NZ+ND));
DPgCR1ph = unwrap(angle(segn2(NZ+1:NZ+ND)))/(2*pi);

DPtoctAlla = Ym; % amplitude in 1/3 oct bands
DPtoctNLa = segn1m-fac;
DPtoctCRa = segn4m-fac;


wlgraph.coewd = coewd;
wlgraph.M = M;
wlgraph.NW = NW;
wlgraph.dfwH = dfwH;
wlgraph.tsm = tsm;
wlgraph.c1 = c1;
wlgraph.c4 = c4;
wlgraph.c2 = c2;
wlgraph.tt = tt;
wlgraph.dfw = dfw;
wlgraph.bbb = bbb;


if plotflag==1
figure(1)
subplot(2,1,1)
plot((NZ+1:NZ+ND)*dfD,20*log10(abs(Y(NZ+1:NZ+ND))),c(1)); hold on
plot((NZ+1:NZ+ND)*dfD,20*log10(abs(segn1(NZ+1:NZ+ND)))-fac,c(2)); hold on 
plot((NZ+1:NZ+ND)*dfD,20*log10(abs(segn4(NZ+1:NZ+ND)))-fac,c(3)); hold on %ampiezza 
plot((NZ+1:NZ+ND)*dfD,20*log10(abs(segn3(NZ+1:NZ+ND)))-fac,'m--'); hold on %ampiezza 
xlabel ('DPOAE Frequency (Hz)')
ylabel ('DPOAE level (dB SPL)')
legend('Total signal','ZL signal', 'LL signal','reconst')

%figure(2)
subplot(2,1,2)
plot((NZ+1:NZ+ND)*dfD,unwrap(angle(Y(NZ+1:NZ+ND)))/(2*pi),c(1)); hold on
plot((NZ+1:NZ+ND)*dfD,unwrap(angle(segn1(NZ+1:NZ+ND)))/(2*pi),c(2)); hold on
plot((NZ+1:NZ+ND)*dfD,unwrap(angle(segn4(NZ+1:NZ+ND)))/(2*pi),c(3)); hold on
plot((NZ+1:NZ+ND)*dfD,unwrap(angle(segn3(NZ+1:NZ+ND)))/(2*pi),'m--'); hold on
xlabel ('DPOAE Frequency (Hz)')
ylabel ('DPOAE phase (cycles)')
legend('Total signal','ZL signal', 'LL signal','reconst')

figure(2)
plot(fc,Ym,c(1)); hold on
plot(fc,segn1m-fac,c(2)); hold on
plot(fc,segn4m-fac,c(3)); hold on
xlabel ('DPOAE Center Frequency (Hz)')
ylabel ('mean DPOAE level (dB SPL)')
legend('Total signal','ZL signal','LL signal')


contr_factor=4;
% Increase contr_factor to enhance visibility in Fig.4 of low-intensity (relative to
% the ZL region) features. The ZL region will look 'over-exposed'

% coewd
% M
% NW
% dfwH
% tsm
% c1
% tt
% c4
% dfw
% bbb
% c2



figure; 
contourf((1:NW)*dfwH,(-M/2+1:M/2)*tsm,1-exp(-abs(coewd(1:M,:)./(cm/contr_factor))));colormap gray;hold on 
plot((1:NW)*dfwH,c1*tt*((1:NW)*dfw).^-bbb,'white', 'LineWidth', 2);hold on
plot((1:NW)*dfwH,c4*tt*((1:NW)*dfw).^-bbb,'white', 'LineWidth', 2);hold on
plot((1:NW)*dfwH,c2*tt*((1:NW)*dfw).^-bbb,'white', 'LineWidth', 2);hold on
axis ([0 4000 -25 25])
xlabel ('DPOAE Frequency (Hz)')
ylabel ('DPOAE delay (ms)')


end