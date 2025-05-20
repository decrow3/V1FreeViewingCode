function fit = fit_tuningcurve_sf_tf(theta, sf, tf, count, useBootstrapping)
% fit = fit_tuningcurve(th, R)
% fit orientation / direction tuning curve using maximum likelihood with
% poisson distributed spike counts

% Assumes directions are in degrees
if nargin < 3
    useBootstrapping = false;
end

% get tuning curves for min and max firing rate
thetas = unique(theta);
tuningCurve = arrayfun(@(x) mean(count(theta==x)), thetas);
tuningCurveSE = arrayfun(@(x) std(count(theta==x))./sqrt(sum(theta==x)), thetas);

if range(theta) > 100 % convert to radians
    theta = theta / 180 * pi;
end

sfs=unique(sf);
tfs=unique(tf);

% initialize min and max parameters
minFR = min(tuningCurve);
maxFR = max(tuningCurve);

% utility functions for wrapping on unit circle or half circle
wrappi = @(x) mod(x/pi, 1)*pi;
wrap2pi = @(x) mod(x/2/pi, 1)*2*pi;

%%
opts = optimset('MaxFunEval', 10e3, 'Display', 'off');

% if you want to try bounds
LB = [-inf 0 0         0    realmin, 0.25, realmin, 0.25, 0]; % lower bound
UB = [inf 25 maxFR maxFR/2 16, 16, 16, 16, 1];
%[mu0, kappa0, maxFR, minFR, mu_sf0, sig_sf0, mu_tf0, sig_tf0,lambda0];

nTrials = numel(count);

% LB = [];
% UB = [];
if useBootstrapping
    nBoot = 100; %2000;
    pBoot = zeros(nBoot, 5);
    pBoot0 = zeros(nBoot, 5);
    
    for i = 1:nBoot
        
        % sample trials
        inds = randi(nTrials, nTrials, 1);
        
        % get initialization
        r0 = sum(count(inds).*exp(1i*(theta(inds))))/sum(count(inds)); %resultant length 2pi
        r1 = sum(count(inds).*exp(1i*wrappi(theta(inds))*2))/sum(count(inds))/2; %resultant length pi
        
        th0 = wrap2pi(angle(r0)); % circular mean
        th1 = wrappi(angle(r1));
        
        cvar0 = abs(r0);
        cvar1 = abs(r1);
        
        rrat = abs(r1)/abs(r0); % ratio gives which params to start with
        
        % set parameters
        lambda0 = 1 - exp(-rrat);
        if rrat > 1
            mu0 = th1;
            kappa0 = cvar1/.3*5;
        else
            mu0 = th0;
            kappa0 = cvar0/.3*5;
        end
    
        mu_sf0=4;
        sig_sf0=4;
        mu_tf0=4;
        sig_tf0=4;

        params0 = [mu0, kappa0, maxFR, minFR, mu_sf0, sig_sf0, mu_tf0, sig_tf0,lambda0];
        
        pBoot0(i,:) = params0;
        
        % fit
        fun = @(params) tcurve.neglogli_poissGLM(tcurve.vonmisesnest_sf_tf(theta(inds),sf(inds),tf(inds), params), count(inds));
        pBoot(i,:) = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
    end
    
    paramsSD = prctile(pBoot, [16 84]); % 68% confidence intervals
    
    %build objective function
    fun = @(params) tcurve.neglogli_poissGLM(tcurve.vonmisesnest_sf_tf(theta, sf, tf, params), count);
    
    % optimization
    [phat, fval] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);

else
    
    % get initialization
    r0 = sum(count.*exp(1i*(theta)))/sum(count); %resultant length 2pi
    r1 = sum(count.*exp(1i*wrappi(theta)*2))/sum(count)/2; %resultant length pi
    
    th0 = wrap2pi(angle(r0)); % circular mean
    th1 = wrappi(angle(r1));
    
    cvar0 = abs(r0);
    cvar1 = abs(r1);
    
    rrat = abs(r1)/abs(r0); % ratio gives which params to start with
    
    % set parameters
    lambda0 = 1 - exp(-rrat);
    if rrat > 1
        mu0 = th0;
        kappa0 = cvar0/.3*5;
    else
        mu0 = th1;
        kappa0 = cvar1/.3*5;
    end

    mu_sf0=4;
    sig_sf0=4;
    mu_tf0=4;
    sig_tf0=4;

    params0 = [mu0, kappa0, maxFR, minFR, mu_sf0, sig_sf0, mu_tf0, sig_tf0, lambda0];
    
    %build objective function
    fun = @(params) tcurve.neglogli_poissGLM(tcurve.vonmisesnest_sf_tf(theta, sf, tf, params), count);
    
    % optimization
    try
        [phat, fval, ~, ~, ~, ~, H] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
        
        % error bars
        paramsSD = sqrt(diag(inv(H)))';
    catch
        phat = nan(9,1);
        paramsSD = nan(9,1);
    end
end

fval = tcurve.neglogli_poissGLM(tcurve.vonmisesnest_sf_tf(theta, sf, tf, phat), count);
fval0 = tcurve.neglogli_poissGLM(mean(count)*ones(size(count)), count);

% fit output
fit.paramsML = phat;
fit.fvalue = fval;
fit.fvalueNull = fval0;

D=2*(fval0-fval);
df=5-1;
fit.llr = D;
fit.llrpval=1-chi2cdf(D, df);

fit.paramsSD = paramsSD;
fit.tuningFun = @(th, sf_, tf_) tcurve.vonmisesnest_sf_tf(th/180*pi, sf_, tf_, phat);
fit.vonmises = @(th, sf_, tf_, params) tcurve.vonmisesnest_sf_tf(th/180*pi, sf_, tf_, params);
if exist('pBoot', 'var')
    fit.pBoot = pBoot;
    fit.pBoot0 = pBoot0;
    
    thetaPref = pBoot(:,1)/pi*180;
    halfWidth = k2hw(pBoot(:,2))/pi*180;
    
    fit.thetaPref = phat(1)/pi*180;
    fit.halfWidth = k2hw(phat(2))/pi*180;
    
    cibnd = [16 84];
    fit.thetaPrefSD = abs(prctile(thetaPref, cibnd) - fit.thetaPref);
    fit.halfWidthSD = abs(prctile(halfWidth, cibnd) - fit.halfWidth);
    fit.minFR = phat(4);
    fit.minFRSD = abs(prctile(pBoot(:,4), cibnd) - phat(4));
    fit.maxFR = phat(3);
    fit.maxFRSD = abs(prctile(pBoot(:,3), cibnd) - phat(3));
    fit.log_mu_sf= phat(5); 
    fit.log_mu_sfSD= abs(prctile(pBoot(:,5), cibnd) - phat(5));
    fit.sig_sf= phat(6);
    fit.sig_sfSD= abs(prctile(pBoot(:,6), cibnd) - phat(6));
    fit.log_mu_tf= phat(7);
    fit.log_mu_tfSD= abs(prctile(pBoot(:,7), cibnd) - phat(7));
    fit.sig_tf= phat(8);
    fit.sig_tfSD= abs(prctile(pBoot(:,8), cibnd) - phat(8));
    fit.lambda = phat(9);
    fit.lambdaSD = paramsSD(9);
else
    fit.thetaPref = phat(1)/pi*180;
    fit.thetaPrefSD = paramsSD(1)/pi*180;
    fit.halfWidth = k2hw(phat(2))/pi*180;
    fit.halfWidthSD = abs(k2hw(phat(3)+paramsSD(3))/pi*180 - fit.halfWidth);
    fit.minFR = phat(4);
    fit.minFRSD = paramsSD(4);
    fit.maxFR = phat(3);
    fit.maxFRSD = paramsSD(3);
    fit.log_mu_sf= phat(5); 
    fit.log_mu_sfSD= paramsSD(5);
    fit.sig_sf= phat(6);
    fit.sig_sfSD= paramsSD(6);
    fit.log_mu_tf= phat(7);
    fit.log_mu_tfSD= paramsSD(7);
    fit.sig_tf= phat(8);
    fit.sig_tfSD= paramsSD(8);
    fit.lambda = phat(9);
    fit.lambdaSD = paramsSD(9);
end

fit.thetas = thetas;
fit.sfs= sfs;
fit.tfs=tfs;
fit.tuningCurve = tuningCurve;
fit.tuningCurveSE = tuningCurveSE;
fit.numTrials = numel(count);

function hw = k2hw(k)
% von Mises k to half-width at half-maximum (in rad.)
% hw = acos(log(0.5)./k+1);
hw = acos(log(.5 + .5*exp(2*k))./k -1);

