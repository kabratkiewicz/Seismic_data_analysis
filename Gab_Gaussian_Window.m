function w = Gab_Gaussian_Window(k, sigma)
%% Input args:
% k - time samples of the window
% sigma - standard deviation of the Gaussian window


if ~exist('sigma', 'var')
 sigma = 10;
end
if ~exist('k', 'var')
 error('k range not defined');
end

%%%%%%%%%%%
A     = 1/(sqrt(2*pi)*sigma);
B     = -1 / (2*sigma^2);
k2    = k.^2;
w     = A * exp(B * k2);
w = w(:);


end

