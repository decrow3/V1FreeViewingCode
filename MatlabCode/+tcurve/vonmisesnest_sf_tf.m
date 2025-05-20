function v = vonmisesnest_sf_tf(dir, sf, tf, params)
% nested von mises function that can return an orientation or
% direction-tuned curve
%
% INPUT:
% x - radians
% params - {mu, kappa, A, b, mu_sf, sig_sf, mu_tf, sig_tf, lambda}
%   mu - central tendency   
%   kappa - dispersion
%   A - amplitude
%   b - baseline
%   mu_sf 
%   sig_sf
%   mu_tf
%   sig_tf

%   lambda - amount of direction

k = params(2);
mu = params(1);
A = params(3);
b = params(4);
mu_sf = params(5);
sig_sf = params(6);
mu_tf = params(7);
sig_tf = params(8);

if numel(params)==8
    lambda = 1;
else
    lambda = params(9);
end

von_mises = exp(k*cos(dir - mu) - k) + A*lambda*exp(k*cos(dir - mu - pi) - k);


% Tuning to SF, gaussian distribution as a function of log(sf).
SF_ = (1 / (sig_sf * sqrt(2 * pi))) * exp(-((log(sf) - mu_sf).^2) / (2 * sig_sf.^2));

% Tuning to TF, gaussian distribution as a function of log(tf), weighting
% absorped by A
TF_ = (1 / (sig_tf * sqrt(2 * pi))) * exp(-((log(tf) - mu_tf).^2) / (2 * sig_tf.^2));

% Full(?) tuning
v = b + A.*von_mises.*SF_.*TF_;
