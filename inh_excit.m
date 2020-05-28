function [h] = inh_excit(HH)
% input: 
%    HH - 2 x 2 x Freq from Wilson factoriztion
% output:
%    h12- 
% Hneg, H21 shows large neg values
% Hpos, H21 shows large pos values
%%% Seth. For Hnegfact, either simply take its real part after ifft or remove the
%%% comment to obtain real-valued impulse response % large neg values
%VAR coecients Ak may then be recovered from H() by a matrix inversion and inverse Fourier transform
% Xiajing Gong @ drexel U 2014/5

fLn = size(HH,3);
for m = 1:fLn
   HHtemp(:,:,m) = inv(HH(:,:,m)); % Makes inv matrix of HH
end

%% 
 HHtemp = repmat(eye(size(HH,1)),[1 1 fLn])-HHtemp; %identity matrix - inv matrix
%%
for m = 1:size(HH,1)
    for n = 1:size(HH,1)
        newH = [squeeze(HHtemp(m,n,:));conj(flipud(squeeze(HHtemp(m,n,2:end))))];
        newH(1) = real(newH(1));
        h(m,n,:) = ifft(newH);
        clear newH
    end
end