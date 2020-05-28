function [h12] = inh_excit1(HH)
%% input: HH is vector of partial Transfer Function from Wilson factoriztion, e.g squeeze(H(1,2,:))

% Hneg, H21 shows large neg values
% Hpos, H21 shows large pos values
%%% For Hnegfact, either simply take its real part after ifft or remove the
%%% comment to obtain real-valued impulse response % large neg values
%newH = [HH;conj(flipud(HH(2:end)))];
%newH = [HH;conj(flipud(HH(2:end)))];
newH = [HH;conj(flipud(HH(2:end)))];
newH(1) = real(newH(1));
h12 = ifft(newH);

%
Xo = [fliplr(conj(HH(2:end)));HH];
Xo = real(ifft(ifftshift(Xo)));