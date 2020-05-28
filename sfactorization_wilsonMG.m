function [Snew,Hnew,Znew] = sfactorization_wilsonMG(S,freq,fs)
%
% It factorizes spectral matrix to transfer function and error covariance matrix by Wilson algorithm
% Input: 
%   S:    spectral density matrix for 2 channels(channels x channels x frequency)
%   freq: array of frequency at which spectral densities are evaluated
%   fs:   sampling frequency.
% Output: %
%   Snew: improved spectral density matrix (channels x channels x frequency)
%   Hnew: transfer function(channels x channels x frequency)
%   Znew: Error covariance matrix(channels x channels)
%
%  Modified by Xiajing Gong, Drexel University, 10/2010
%
m = size(S,1);
N=length(freq)-1;
Sarr = zeros([size(S, 1) size(S, 2) 2*N]);
Sarr(:,:,1:size(S, 3)) = S;
I = repmat(eye(size(S,1),size(S,2)),[1 1 (size(S,3)-1)]); % added by KB dec 16th 2019

Sarr(:,:,2*N+2-(2:length(freq))) = mtimesx(S(:,:,2:length(freq)), 'T', I);

%perform ifft to obtain gammas
gam = real(ifft(Sarr, [], 3)*fs);
gam0 = gam(:,:,1);
h = chol(gam0);
psi = repmat(h, [1 1 size(Sarr, 3)]);
I = repmat(eye(size(Sarr,1),size(Sarr,2)),[1 1 size(Sarr,3)]); % added by KB dec 16th 2019


mtimesx('SPEEDOMP');
for iter = 1:25
    psiinv = NDInv(psi);
%     psiinv = mmx('backslash', psi, I);  %inaccurate when iter is large
    g = mtimesx(mtimesx(psiinv, Sarr), psiinv, 'C')+I;
    gp = PlusOperator(g, m, fs, freq);
    psiold = psi;
    psi = mtimesx(psi, gp);
    psierr = max(squeeze(sum(abs(psi-psiold), 1)), [], 1);
    psierrf = mean(psierr);
    if psierrf<1E-12
        break;
    end
end
gamtmp = real(ifft(psi, [], 3));
A0 = squeeze(gamtmp(:,:,1));
A0inv = inv(A0);
Znew = A0*A0.'*fs;
psi = psi(:,:,1:length(freq));
Snew = mtimesx(psi, psi, 'C');
Hnew = mtimesx(psi, repmat(A0inv, [1 1 length(freq)]));
% Serr = Sarr(:,:,1:length(freq))-mtimesx(mtimesx(Hnew, repmat(Znew, [1 1 length(freq)])), Hnew, 'C')/fs;
% Serrnorm = max(squeeze(sum(abs(Serr), 1)), [], 1);
% for k = 1:length(freq)
%     Serrnorm(k) = Serrnorm(k)/norm(Sarr(:,:,k));
% end