

function wind = computeWin(freq, windowIndex)
    
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
