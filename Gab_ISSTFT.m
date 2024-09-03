function [ S_rec ] = Gab_ISSTFT(S, L, M)
% signal reconstruction from the STFT

S_rec = sum( S,1) * sqrt(2*pi) * L/M;

end