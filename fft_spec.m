function [X, Pxx, f] = fft_spec(x, window, Fs, pad);
%
% calculation of spectum of a TS with many trials
%
%  [Pxx, f] = coh_fft(x, window, Fs, pad);
%  
% Input:
%   x, y: T x N respectively, T - length of TS, N - size of ensemble
%   pad:  padding for fft window, which CAN'T increase the resolution of spec.
%   Fs:   sampling frequency
%   window: window used, e.g. hanning.
% Output:
%    Pxx: spectra
%    f:   freq

%
%  Hualou Liang, 03/30/99, FAU
%


[T, N] = size(x);  
if nargin<4,
  pad = T;
end

if nargin<3,
  Fs = 200;
end

if nargin<2,
  window = boxcar(T);  % T x 1
end

% window data here
% window = hanning(T);
x = x.*window(:,ones(1, N));

X = fft(x, pad);
nfft = size(X, 1)/2;

Pxx = sum(X.*conj(X), 2)/N;
% Pxx = Pxx(1:T);
Pxx = Pxx(1:nfft);

% coh = abs(Pxy(1:nfft)).^2 ./ (Pxx(1:nfft).*Pyy(1:nfft));

% f = [1:nfft]*Fs/size(X, 1);
f = [0:nfft-1]*Fs/size(X, 1);

