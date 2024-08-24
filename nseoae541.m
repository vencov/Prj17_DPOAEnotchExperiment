function [oae, yp, numOK] = nseoae541(params,yn,yp,numOK)
% NSEOAE - least squares fit (lsf) signal estimate
%   Origanlly developed by Carrick Talmadge
%   Long, G. R., & Talmadge, C. L. (1997). Spontaneous otoacoustic emission frequency is modulated by heartbeat. 
%   The Journal of the Acoustical Society of America, 102, 2831?2848.
%
%   This is the version adapted by Vaclav Vencovsky to include artifact
%   rejection
%   NSEOAE541(PARAMS, yn, yp, numOK) estimates the amplitude and phase of N components of
%   the signal contained in yn.
%
%   PARAMS is a structure containing N components to be estimated.
%   EX.
%        params.T = log2(F2stop/F2start)/.5;     % poustelo se rychlosti 0.5 oktavy/s                               % 6 second sweep
%        params.nfilt = 5512;                        % window size
%        params.nstep = 51;                          % step size
%        params.component(1).type = 'linsweep';      % component type            (values can be 'logsweep', 'tone', 'linear', 'custom'; for 'custom' a phase function, pf needs to be defined)
%        params.component(1).f1 = 2*F2start/1.2-F2start;              % logsweep start frequency  (for 'tone' type, only f1 is necessary)
%        params.component(1).f2 = 2*F2stop/1.2-F2stop;              % logsweep end frequency    (for 'tone' type, only f1 is necessary)
%        y=recordedSignal;
%   Optional component parameters for 'custom' type:
%       params.component(1).phi = @(t)Phi0.*(exp(gamma*t)-1); % define a custom phase function, phi
%
%   Alternatively, you can just set phi to the pre-computed phi values:
%       params.component(1).phi = myphivector;                % length = T*fs x 1
%
%   Optional PARAMS:
%                                                   % This (nskip) is a very important parameter, as can cause misalignment 
%                                                   % between LSF model and actual signal, especially for long filter windows    
%       params.nskip = 6200;                         % ?samples to skip, default = 0 (e.g. due to system delay, propogation delay)
%
%   for components with predefined latency, an optional function can be specified
%       
%         params.component(n).delay = '1./(0.035+0.00006.*x+0.000000002.*x.^2)';    % latency function (specified in ms) as a function of frequency, x
%
%   an optional 'fit' can be used (only if a delay function is applied):
%      params.fit = -5:1:5;    % fitting values, in ms, over which to evaluate the fit
%
%   August 2011 - Simon Henin <shenin@gc.cuny.edu>
%   Copyright, Simon Henin 2011-2018
%   Version 1.3
%
%   12-2012 - added linear sweeps
%   01-14-2013 - added verified delay estimates
%   01-14-2013 - added fitting procedure
%   01-30-2013 - added versioning
%   05-08-2013 - find minimum rms_error keeping track of index (reduce
%   overhead) TODO: find more elegant solution to solving this problem
%   05-12-2013 - Added amplitude-varying components
%   07-18-2018 - Added custom phase functions


y = zeros(size(yn));

fit = 0;
res = [];
maximize_by_component = [];

%% basic error checks
if ~isfield(params, 'component'),
    error('No components specified in the model.');
end
if ~isfield(params, 'T'),
    error('T (duration) must be specified.');
end
if ~isfield(params, 'nfilt'),
    error('nfilt (filter size) must be specified.');
end
if ~isfield(params, 'nstep'),
    error('nstep (step size) must be specified.');
end

%% set optional params
if isfield(params, 'fs'),
    fs = params.fs;
end
if isfield(params, 'nskip'),
    nskip = params.nskip;
end
if isfield(params, 'p'),
    p = params.p;
end
if isfield(params, 'gain'),
    gain = params.gain;
end
if isfield(params, 'delay'),
    del = params.delay;
end
if isfield(params, 'windowType'),
    windowType = params.windowType;
end
if isfield(params, 'maximize_by_component'),
    maximize_by_component = params.maximize_by_component;
end


%% set required params
T = params.T;
nfilt = params.nfilt;
nstep = params.nstep;
dt = 1/fs;

%time vector
t = 0:dt:T;

%% initialize the phase function for each component
nPhi = size(params.component,2);    %number of components*2 (sin & cos)
nExtra = 0;                         % number of extra (e.g. amplitude varying comp.)
Phi = cell(nPhi);
PhiF = cell(nPhi);

delay_t = zeros(nPhi, length(t));
delay_phi = zeros(nPhi, length(t));
fit = 0;
niter = 1; % no fitting = only 1 iteration
idx = 1;
GainMicA = params.GainMicA;

for i=1:nPhi,
    switch params.component(i).type,
        case {'custom'},
            % this is a custom sweep with predefined phase function
            if isfield(params.component(i), 'phi')&~isempty(params.component(i).phi),
                if isa(params.component(i).phi,'function_handle'),
                    Phi{i} = params.component(i).phi;
                else
                    % phi is a vector (e.g. pre-computed phase function)
                    % generate a hacked function handle to just get samples by index
                    pf = params.component(i).phi;
                    if iscolumn(pf), pf = pf'; end
                    Phi{i} = @(x) pf(round(x*fs)+1); % t is time (in seconds), so multiply by fs to get indices
                end
                
            else
                error('Phase function is undefined');
            end
        case {'logsweep'},
        
          
            f1 = params.component(i).f1;
            f2 = params.component(i).f2;
            
            % extra components?
            if isfield(params.component(i), 'extra')
            if (params.component(i).extra),
                nExtra = nExtra+2;                       %two extra components (half sine, full sine)
                params.component(i).idx = [idx:idx+2]; 
                idx = idx+3;
            else
                params.component(i).idx = [idx]; %1(*2) component total
                idx = idx+1;
            end
            end
            
            % logsweep phase function setup
            gamma = log(f2/f1)/T;

            % The version results is slightly out of phase            
             Phi0(i) = 2*pi*f1/gamma;
             Phi{i} = @(t) Phi0(i) .* (exp(gamma.*t)-1);
            
            %Phi0(i) = (2*pi) * f1 * T/log(f2/f1);
            %Phi{i} = @(t) Phi0(i) * (f2/f1).^(t/T);
            PhiF{i} = @(t) f1.*exp(gamma .* t);
            
            % setup fit (fit contains array of delay offsets, e.g. -10:1:10)
            if isfield(params.component(i), 'fit'),
                fit = params.component(i).fit;
            end
            
            
            % Optional delay function (e.g. to extract reflection source
            % component from DPOAE). Delay function is specified as a
            % function of frequency.
            if isfield(params.component(i), 'delay'),
                if ~isempty(params.component(i).delay),
                    niter = length(fit);
                    for n=1:niter,
                        x = (f1)*(exp(gamma .* (t) )-1) + f1; % sweep frequency
                        delay_t(i,:,n) = ( feval(@(x) eval(params.component(i).delay), x) + fit(n) ) ./ 1000; % evaluate delay function at sweep frequency (in ms)
                        delay_phi(i,:,n) = computePhiTau(delay_t(i,:,n), x);                        
                    end
                end
            end
        case {'linsweep'},
            f1 = params.component(i).f1;
            f2 = params.component(i).f2;
            w = 2*pi*f1;
            wdot = pi*(f2-f1)/T;
            beta = (f2-f1)/T;
            Phi{i} = @(t) (w + wdot .* t) .* t;
            PhiF{i} = @(t) f1 + beta .* t;
            % Optional delay function (e.g. to extract reflection source
            % component from DPOAE). Delay function is specified as a
            % function of frequency.
            if isfield(params.component(i), 'delay'),
                if ~isempty(params.component(i).delay),
                    niter = length(fit);
                    for n=1:niter,
                        x = f1+beta.*t; % sweep frequency
                        delay_t(i,:,n) = ( feval(@(x) eval(params.component(i).delay), x) + fit(n) ) ./ 1000; % evaluate delay function at sweep frequency (in ms)
                        delay_phi(i,:,n) = computePhiTau(delay_t(i,:,n), x);
                    end
                end
            end
            
        case {'tone'},
            f = params.component(i).f1;
            w = 2*pi*f;
            Phi{i} = @(t) w .* t;
            PhiF{i} = @(t) f;           % does not change with time, t, but used to keep format the same.
    end
end

if size(y, 1) > size(y, 2),
    y = y';
    yn = yn';
end
y = y(nskip+1:end);
yn = yn(nskip+1:end);

%% set up analysis loop
twin = nfilt/fs;
tstep = nstep/fs;
nsize = min(floor(T*fs), length(y));
nfilt = floor(twin * fs + 0.5);
nstep = floor(tstep * fs + 0.5);
ntotal = floor((nsize-nfilt)/nstep) + 1;

%% initialize outputs
fvec = zeros((nPhi+nExtra),ntotal);      % used to keep track of frequency
tcoef = zeros(1,ntotal);        % used to keep track of current sample window time
coef = zeros((nPhi+nExtra)*2, ntotal);   % coefficients for lsf problem
del = zeros((nPhi+nExtra),ntotal);       % used to keep track of estimated delay (if applicable)
phicorr = zeros((nPhi+nExtra), ntotal);
rmserror = zeros((nPhi+nExtra), ntotal);

t1 = 0;
t2 = twin;
toff = twin/2;

% Window function  2=hann window
Wvec = computeWindow(0:nfilt-1, windowType);

v1 = cos(pi*((0.5+(0:nfilt-1))/nfilt));         % half period cosine
v2 = cos(2*pi*((0.5+(0:nfilt-1))/nfilt));       % full period cosine



fx = (0:nfilt-1)*fs/nfilt; % frequency axis in one frame

if nargin<3
%    yp = zeros(nfilt,ntotal);
    numOK = zeros(ntotal,1);
end
%remove later
nntau = []; dur = []; signal_estimate = [];
nFloor = zeros(1,ntotal);
nntau = zeros(1,ntotal);
signal_estimate = zeros(nfilt,ntotal);

% initiliaze temporary storage variables (used in iterative fitting)
    coefs = zeros((nPhi+nExtra)*2, niter);
    delay = zeros((nPhi+nExtra), niter);
    Pdelay = zeros((nPhi+nExtra),niter);
    Tdelay = zeros((nPhi+nExtra),niter);
    residue = zeros((nPhi+nExtra), niter);
    p_delay = zeros((nPhi+nExtra),ntotal);
    t_delay = zeros((nPhi+nExtra),ntotal);
    %infoN = zeros(ntotal,1);
for istep = 1:ntotal;
    % index paramters
    tcoef(istep) = t1 + toff;
%     tfilt = t2 - t1;
    tvec = (0:(nfilt-1))*dt + t1;
%     ttvec = (1:(nfilt))*dt + t1;
    ntau = floor(tcoef(istep)*fs);
    %ntau = floor(t1+10);
    
    % data indices
    i1 = floor(t1 * fs + 0.5) +1;
    i2 = i1 + nfilt - 1;
    
    Yn = abs(fft(Wvec .* yn(i1:i2))/nfilt); % fft for new signal to estimate background noise
    pfFreq = feval(PhiF{i}, (tvec)); % instantaneous frequency of fdp in the frame
    IdNoise1 = find(fx>=(pfFreq(1)*(1 - 1/(2-params.f2f1) + 1)),1,'first')+3; % find the lowest freq edge
    IdNoise2 = find(fx>=pfFreq(1),1,'first')-2; % find the highest freq edge;
    %IdNoise2 = find(fx>=pfFreq(1)+params.dFnoise,1,'first'); % find the lowest freq edge
    NoiseFloor = max(Yn(IdNoise1:IdNoise2)); % estimate noise in the signal
    NoiseFloordB = 20*log10(abs(giveTPaWin(mean(fx(IdNoise1:IdNoise2)),NoiseFloor,GainMicA,fs,1,fx))/(sqrt(2)*2e-5));
    
    % artifact rejection criteria
    if NoiseFloordB<0 % for smaller noise than 0 dB SPL include this frame
       if nargin<3
           y(i1:i2) = yn(i1:i2);
           yp(:,istep) = y(i1:i2);
       else
           y(i1:i2) = yn(i1:i2)/(numOK(istep)+1) + numOK(istep)*yp(:,istep)'/(numOK(istep)+1);
           yp(:,istep) = y(i1:i2);
       end
       numOK(istep) = numOK(istep) + 1;
       % estimate noise floor
       Y = abs(fft(Wvec .* y(i1:i2))/nfilt); % fft for new signal to estimate background noise
       nFloor(istep) = mean(abs(Y(IdNoise1:IdNoise2))); % estimate noise in the signal
    if sum(fit),
        rms_error = zeros(1,length(fit));
    end
    
   
    
    min_index = 0;    
    min_value = 0;
    for n=1:niter,
        
        % Build the model, including each component
        %F = zeros(nPhi*2, nfilt);
        F = zeros((nPhi+nExtra)*2, nfilt);
        ii = 1;
        for i=1:nPhi,
            Pdelay(i, n) = delay_phi(i, ntau, n);
            Tdelay(i, n) = delay_t(i, ntau, n);
            
            %fvec(i,istep) = feval(PhiF{i}, median(ttvec));  % frequency is the middle of the window %%--> SH, 7/18/18. No longer using function handle defined frequency. INstead differentiate the phase (see below)
            %pf = feval(Phi{i}, (tvec - delayfunc(i,ntau) )); % phase function with optional delay at frequency in middle of this window (e.g. ntau)
            pf = feval(Phi{i}, (tvec))+delay_phi(i, i1:i2, n); % phase function with optional delay at frequency in middle of this window (e.g. ntau)
            fvec(i,istep) = median( diff(pf*fs)./(2*pi) ); % differntiate the phase function to get frequency, then take the middle value
            
            F((2*ii)-1, :) = cos(pf);
            F((2*ii), :) = sin(pf);
            ii = ii+1;
%             if isfield(params.component(i), 'extra'),
%                 if params.component(i).extra,
%                     F((2*ii)-1, :) = v1.*cos(pf);
%                     F((2*ii), :) = v1.*sin(pf);
%                     ii = ii+1;
%                     F((2*ii)-1, :) = v2.*cos(pf);
%                     F((2*ii), :) = v2.*sin(pf);
%                     ii = ii+1;
%                 end
%             end
            
            %phicorr(i,istep) = -(pi/2)*(delayfunc(i,ntau)+delayfunc(i, ntau-1))*(feval(PhiF{i}, (ntau + 1)/fs) - feval(PhiF{i}, (ntau - 1)/fs));
        end
        
        
        % LSF
        
%         % compute covariance matrix ( A*A') of windowed components
%         FT = F';                         % transpose for faster processing
%         A = zeros((nPhi+nExtra)*2, (nPhi+nExtra)*2);
%         tmp = zeros(1,(nPhi+nExtra)*2);
%         nComp = (nPhi+nExtra)*2;
%         for i = 1:nComp
%             tmp = Wvec .* F(i, :);
%             A(i,:) = tmp * FT;
%         end
%        
%         
%         b =  F * (Wvec .* y(i1:i2))';
%         
%         %keyboard;
%         % solve the matrix equation A*x = b
%         coefs(:,n) = A \ b;  % x = inv(A)*b;
        
        % simplified syntax for the same code as above
        Fw = bsxfun(@times, F, Wvec);
        b =  F * (Wvec .* y(i1:i2))';
        coefs(:,n) = (Fw*F') \ b;
        
        
        est = F' * coefs(:,n);
        
        residual = Wvec .* (y(i1:i2) - est');
        %signal_estimate = [signal_estimate (y(i1:i2) - est')'];
        signal_estimate(:,istep) = (y(i1:i2) - est')';
        % estimate residual for each component using an fft. Look in
        % clostest bin
        nfft = length(residual);
        f = fs*(0:(nfft/2))/nfft;
        Y = fft(residual, nfft);
        Y = Y*2/nfft; % scale it based on fft bin size, multiply by 2, since we will only take half
        Y = abs(Y(1:floor(length(Y)/2)+1))./sqrt(2); % rms
        for i=1:nPhi, 
          [~, f_idx] = min(abs(f - fvec(i,istep)));
           residue(i,n) = Y(f_idx);
        end
        
        % fitting to max OR min rms error
%         if maximize_by_component,
%             total = hypot(coefs((2*maximize_by_component)-1,n), coefs((2*maximize_by_component),n));
%             if n==1,
%                 fit_value = total;
%                 fit_index = n;
%             else
%                 if total > fit_value,
%                     fit_value = total;
%                     fit_index = n;
%                 end
%             end
%         else
            rms_error(n) = sqrt(mean(residual.^2));
            m = sqrt(mean(residual.^2));
%             if n==1,
                fit_value = m;
                fit_index = n;
%             else
%                 if  m < fit_value,
%                     fit_value = m;
%                     fit_index = n;
%                 end
%             end
%         end
    end
    
    rmserror(:, istep) = residue(:,fit_index);
%     rmserror(:, istep) = fit_value;
    coef(:,istep) = coefs(:,fit_index);
    p_delay(:,istep) = Pdelay(:,fit_index);
    t_delay(:,istep) = Tdelay(:,fit_index);
    
    %nntau = [nntau ntau];
    nntau(istep) = ntau;
    else % above noise floor
        if nargin<3
            yp(:,istep) = y(i1:i2)*0;
            numOK(istep) = numOK(istep);
        else
            yp(:,istep) = yp(:,istep);
            numOK(istep) = numOK(istep);
        end
       % infoN(istep) = 1;
    end
        
    % increment the frame
    t1 = t1 + tstep;
    t2 = t2 + tstep;
end


%% OUTPUT RESULTS
for i=1:nPhi,
    %% convert LSF coefficients into useful data
    oae(i).time         = tcoef;
    oae(i).f            = fvec(i,:);
    oae(i).coefs        = [coef((2*i)-1,:); coef((2*i),:)];
    oae(i).mag          = hypot(coef((2*i)-1,:), coef((2*(i)),:));
    oae(i).phase        = -atan2(coef((2*i),:),coef((2*i)-1,:)); % radians
    oae(i).total_phase  = -atan2(coef((2*i),:),coef((2*i)-1,:))+p_delay(i,:); % radians
    oae(i).delay        = t_delay(i,:);
    oae(i).rmserror     = rmserror(i,:);
    oae(i).nFloor     = nFloor;
    
end
oae(1).signal_estimate = signal_estimate;
end

function Phi_tau = computePhiTau(delay, freq)
%% computePhiTau
ndata = length(delay);
ndatap = ndata -1;
Phi_tau = zeros(1, ndata);
sum = (-pi/2) * delay(1) * (freq(2)-freq(1));
Phi_tau(1) = sum;
for (n = 2:ndatap),
    sum = sum - pi/2 * (delay(n) + delay(n-1)) * (freq(n+1)-freq(n-1));
    Phi_tau(n) = sum;
end
Phi_tau(n+1) = sum;
end


function wind = computeWindow(freq, windowIndex)
    
    df = freq(2)-freq(1);

    switch windowIndex
    case -1
%	'rect'
	wind = ones(1, length(freq));
    case 0
%	'welch'
	wind = wwelch(freq, freq(1)-df, freq(length(freq))+df);
    case 1
%	'hamming'
	wind = whamming(freq, freq(1)-df, freq(length(freq))+df);
    case 2
%	'hann'
	wind = whann(freq, freq(1)-df, freq(length(freq))+df);
    case 3
%	'parzen'
	wind = wparzen(freq, freq(1)-df, freq(length(freq))+df);
    case 4
%	'blackman'
	wind = wblackman(freq, freq(1)-df, freq(length(freq))+df);
    case 5
%	'gauss'
	wind = wgauss(freq, freq(1)-df, freq(length(freq))+df, windowAux);
    otherwise
        error('internal error, undefined type %d', windowIndex);
    end
    % normalize window
%     wnorm = sqrt(sum(wind.^2));
%     wind = wind ./ wnorm;
%save ans
end

function wx = whamming(x, x1, x2)
    k = 2*pi/(x2-x1);
    phi = -k*x1;
    wx = 0.54 - 0.46 * cos(k*x + phi);
end

function wx = wwelch(x, x1, x2)
   y = (x-x1)/(x2-x1);
   wx = 4*y.*(1-y);
end

function wx = whann(x, x1, x2)
    k = 2*pi/(x2-x1);
    phi = -k*x1;
    wx = 0.5 * (1 - cos(k*x + phi));
end

function wx = wparzen(x, x1, x2)
   y = (x-x1)/(x2-x1);
   wx = 1-abs(2*y-1);
end

function wx = wblackman(x, x1, x2)
    Phi = 2*pi*(x-x1)/(x2-x1);
    wx = 0.42 - 0.5 * cos(Phi) + 0.08 * cos(2*Phi);
end

function wx = wgauss(x, x1, x2, auxdata)
    y = ((x-x1)/(x2-x1) - 0.5)*auxdata;
    wx = exp(-0.5*y.*y);
end


