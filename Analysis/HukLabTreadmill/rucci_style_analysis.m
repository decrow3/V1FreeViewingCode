
%%
% import os
% import numpy as np
% import matplotlib.pyplot as plt
% from huk_treadmill import load_file, get_run_epochs, get_good_eye_epochs, get_epochs, compute_Q, get_eye_displacements

%% Load example file
% fname = '/Users/jake/Downloads/Eyesrunning_Rocky20240802_reD.mat'
% fdir = '/Users/jake/Dropbox/Datasets/HuklabTreadmill/gratings/'
% flist = os.listdir(fdir)
% %flist = [f for f in flist if f.endswith('.mat')]
% 
% subject = 'brie'
% 
%flist = [f for f in flist if subject in f]
% 
x_segments_run = [];
x_segments_stat = [];
% 
% gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/V1simult/';
% subject = 'npx_V1_fov';
% cd('/home/huklab/Documents/NPX_pilot/V1Locomotion/Code/')
% D = load_subject(subject,gratingpath);
% cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')

% get run and eye epochs
    [t_run_starts, t_run_ends] = get_run_epochs(D,false,551,2); %plotfig=0,medfilt_window=551,speed_threshold=2)\
    %%
    [t_good_starts, t_good_ends] = get_good_eye_epochs(D, false);%plot=
    %%
    [t_starts, t_ends] = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, true);%union=true
    %%
    xsegments_run_ = get_eye_displacements(D, t_starts, t_ends, 512, 512-512/4);% Twin = 512, noverlap=512/4);
    %%
    [t_starts, t_ends] = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, false);
    %%
    [xsegments_stat_] = get_eye_displacements(D, t_starts, t_ends, 512, 512-512/4);
    x_segments_run=[x_segments_run; xsegments_run_];
    x_segments_stat=[x_segments_stat; xsegments_stat_];
    


%%
xs_run = x_segments_run';%concatenate(x_segments_run)
xs_stat = x_segments_stat';%concatenate(x_segments_stat)

disp(size(xs_run))
disp(size(xs_stat))


%%
%matplotlib ipympl
figure()
plot(D.treadTime, D.treadSpeed)
for i = 1:(length(t_run_starts))
%     fill_betweenx([-1, nanmax(D.treadSpeed)], t_run_starts(i), t_run_ends(i), color='r', alpha=.5)
end
% show()

%%
%matplotlib inline
bins = -15:1/60:15;
prun = zeros(length(bins)-1, size(xs_run,2));
pstat = zeros(length(bins)-1, size(xs_stat,2));
for i = 1:size(xs_run,2)
    prun(:,i) = histcounts(xs_run(:,i), bins);
    prun(:,i)  = gaussfilt(1:length(prun(:,1)),prun(:,i), 5);
end
for i = 1:size(xs_stat,2)
    pstat(:,i) = histcounts(xs_stat(:,i), bins);
    pstat(:,i) = gaussfilt(1:length(pstat(:,1)),pstat(:,i), 5);
end


%%
figure
subplot(1,2,1)
imagesc(log10(prun+1e-3))
%), aspect='auto', extent=[0, xs_run.shape(2), bins(1), bins(end)])
title('Running')
xlabel('Time (ms)')
ylabel('Position (deg)')

subplot(1,2,2)
imagesc(log10(pstat+1e-3));%, aspect='auto', extent=[0, xs_stat.shape(2), bins(1), bins(end)])
title('Stationary')
xlabel('Time (ms)')
ylabel('Position (deg)')



%% Loop over segments and calculate the eye trajectory displacement power spectrum

omega_x_values = 2.^linspace(-7, 3, 50);
omega_t_values = 2.^linspace(-7, 7, 30);
% omega_t_values = concatenate((append(flip(-omega_t_values), 0), omega_t_values))
omega_t_values = [flip(-omega_t_values) 0 omega_t_values];

fps = 1000;

Qrun = 0;

nseg = size(xs_run,1);
for i = 1:nseg
    xseg = xs_run(i,:);
    Qrun = Qrun + compute_Q(xseg, omega_x_values, omega_t_values, 1/fps);
end
Qrun = Qrun/nseg;

Qstat = 0;
nseg = size(xs_stat,1);
for i = 1:nseg
    xseg = xs_stat(i,:);
    Qstat = Qstat +compute_Q(xseg, omega_x_values, omega_t_values, 1/fps);
end
Qstat = Qstat/nseg;


%%

%% run main analysis

% get support for the frequency analysis
omega_x_values = 2.^linspace(-7, 3, 50);
omega_t_values = 2.^linspace(-7, 7, 30);
% omega_t_values = concatenate((append(flip(-omega_t_values), 0), omega_t_values))
omega_t_values = [flip(-omega_t_values) 0 omega_t_values];
[fx, ft] = meshgrid(omega_x_values, omega_t_values);
scatter(fx, ft, hypot(fx,ft)/1000);
% _ = plot(fx, ft, '.k', 'markersize', .01)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')
title('Frequency Support')


%% plot the same thing but for the power spectrum
%matplotlib inline
lambda = @(x) 10*log10(x + 1e-6);


% vmin = min(fun(Qstat).min(),lambda(Qrun).min());
% vmax = max(fun(Qstat).max(),lambda(Qrun).max());

vmin = min(min(lambda(Qstat(:))),min(lambda(Qrun(:))));
vmax = max(max(lambda(Qstat(:))),max(lambda(Qrun(:))));

figure(2)
subplot(1,3,1)
contourf(xfun(omega_x_values), xfun(omega_t_values), lambda(Qrun), linspace(vmin, vmax, 25));%, linspace(vmin, vmax, 25), vmin=vmin, vmax=vmax)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Q (Running)')
%%
set(gca,'XTickLabel',xinv(get(gca,'XTick')))
set(gca,'YTickLabel',xinv(get(gca,'YTick')))

subplot(1,3,2)

contourf(xfun(omega_x_values), xfun(omega_t_values),lambda(Qstat), linspace(vmin, vmax, 25));%, vmin=vmin, vmax=vmax)
%imagesc(xfun(omega_x_values), xfun(omega_t_values),lambda(Qstat))
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Q (Stationary)')
set(gca,'XTickLabel', xinv(get(gca,'XTick')))
set(gca,'YTickLabel', xinv(get(gca,'YTick')))

% plot difference
subplot(1,3,3)
contourf(xfun(omega_x_values), xfun(omega_t_values),lambda(Qrun) -lambda(Qstat), 25); %, vmin=.5, vmax=1.5, linspace(.5, 1.5, 25))
%imagesc(xfun(omega_x_values), xfun(omega_t_values),lambda(Qrun) -lambda(Qstat))
set(gca,'XTickLabel',xinv(get(gca,'XTick')))
set(gca,'YTickLabel',xinv(get(gca,'YTick')))

xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Relative Gain')
colorbar()

%subplots_adjust(hspace=.5, wspace=.5)
%% plot Q at different TF slices
figure(5)
tf0 = 2;
hold off
cmap=colormap;
for tf = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
    [~,iitf] = min(abs(omega_t_values - tf - tf0));
    h = plot(omega_x_values,lambda(Qrun(iitf,:)),'Color',cmap(tf*20,:));%, label='%.1f Hz' % tf);
    hold on
    % draw line usin same color
    plot(omega_x_values,lambda(Qstat(iitf,:)),'--','Color',cmap(tf*20,:));%, '--', color=h(1).get_color());
    % ylim(-20, 0)
end

%
% for tf in [-0.5, -1.0, -2.0, -4.0, -8.0, -12.0]:
%     iitf = argmin(abs(omega_t_values - tf - tf0))
%     h = plot(omega_x_values,lambda(Qrun[iitf,:]), label='%.1f Hz' % tf)
%     % draw line usin same color
%     plot(omega_x_values,lambda(Qstat[iitf,:]), '--', color=h(1).get_color())
%     % ylim(-20, 0)
legend('0.5','','1.0','','2.0','', '4.0','', '8.0','', '12.0','')

set(gca,'XScale','log')

for sf = [1.0, 2.0, 4.0]
    plot([sf sf],ylim, 'k--')
end
legend()
xlim([0.1 10])
xlabel('Spatial Frequency (cyc/deg)')
%title(subject)

%%
% save(['/home/huklab/Documents/RunnyEyeAnalysis/' subject '_stats.mat'],'omega_x_values','omega_t_values','Qstat','Qrun','xs_run','xs_stat','prun','pstat')

%%
%{
%% get hypothetical neural tuning

%from huk_treadmill import gaussian_derivative_RF

R_fovea = gaussian_derivative_RF(omega_x_values, omega_t_values, 
                                 omega_x_0=2, omega_y_0=3, b=.5, d=0)
R_periph = gaussian_derivative_RF(omega_x_values, omega_t_values, 
                                  omega_x_0=.5, omega_y_0=4, b=1, d=0)

% plot
figure(figsize=(9,3))
subplot(1,2,1)
contourf((omega_x_values), (omega_t_values), R_fovea, 15)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Fovea')
xlim((0,8))
ylim((0, 15))

subplot(1,2,2)
contourf((omega_x_values), (omega_t_values), R_periph, 15)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Periphery')
xlim((0,8))
ylim((0, 15))


stim_sfs = [1, 2, 4]
speeds = [1.0, 2.0]
figure(figsize=(9,9))

G_Run_fovea = zeros((length(stim_sfs), length(speeds)))
G_Run_periph = zeros((length(stim_sfs), length(speeds)))
G_Stat_fovea = zeros((length(stim_sfs), length(speeds)))
G_Stat_periph = zeros((length(stim_sfs), length(speeds)))

for istim = 1:lengthgth(stim_sfs)
    
    stim_sf = stim_sfs[istim]
    iisf = argmin(abs(omega_x_values - stim_sf))


    for i=speeds
        subplot(length(stim_sfs), length(speeds), istim*length(speeds)+i+1)
        
        tf = stim_sf * speed
        % iitf = argmin(abs(omega_t_values - tf))

        plot(omega_t_values + tf,lambda(Qrun[:,iisf]), 'k', label='Running')
        plot(omega_t_values + tf,lambda(Qstat[:,iisf]), 'k--', label='Stationary')
        xlabel('Temporal Frequency (cmaximumycles/s)')
        ylabel('dB')
        legend()
        title('Stimulus SF: %.1f cyc/deg, Speed: %.1f cyc/s' % (stim_sf, speed))
        axvline(tf, color='k', linestyle='--')
        xlim(0, 32)

        % new axis that has same x values, but new y values
        ax = gca().twinx()
        ax.plot(omega_t_values,lambda(R_fovea[:,iisf]), 'r', label='Fovea')
        ax.plot(omega_t_values,lambda(R_periph[:,iisf]), 'b', label='Periphery')
        % ax.set_ylim(0, 1)
        ax.set_xlim(0, 32)

        G_Run_fovea[istim,i] = trapz(R_fovea[:,iisf] * Qrun[:,iisf], omega_t_values)
        G_Run_periph[istim,i] = trapz(R_periph[:,iisf] * Qrun[:,iisf], omega_t_values)
        G_Stat_fovea[istim,i] = trapz(R_fovea[:,iisf] * Qstat[:,iisf], omega_t_values)
        G_Stat_periph[istim,i] = trapz(R_periph[:,iisf] * Qstat[:,iisf], omega_t_values)
    end
end
subplots_adjust(hspace=.5, wspace=.5)

figure()
plot(stim_sfs, G_Run_periph*1e3, 'b')
plot(stim_sfs, G_Stat_periph*1e3, 'r')
xlabel('Spatial Frequency (cyc/deg)')

figure()
plot(G_Run_periph / G_Stat_periph, 'bo-', label='Periphery')
plot(G_Run_fovea / G_Stat_fovea, 'ro-', label='Fovea')
axhline(1, color='k', linestyle='--')
%%
vmin = minimum(minimum(G_Run_fovea.min(), G_Run_periph.min()), minimum(G_Stat_fovea.min(), G_Stat_periph.min()))
vmax = maximum(maximum(G_Run_fovea.max(), G_Run_periph.max()), maximum(G_Stat_fovea.max(), G_Stat_periph.max()))

figure(figsize=(9,3))
subplot(2,3,1)
imshow(G_Run_fovea, vmin=vmin, vmax=vmax)
title('Running Fovea')
ylabel('Stimulus SF')
xlabel('Speed')

subplot(2,3,2)
imshow(G_Stat_fovea, vmin=vmin, vmax=vmax)
title('Stationary Fovea')
ylabel('Stimulus SF')
xlabel('Speed')

subplot(2,3,4)
imshow(G_Run_periph, vmin=vmin, vmax=vmax)
title('Running Periphery')
ylabel('Stimulus SF')
xlabel('Speed')

subplot(2,3,5)
imshow(G_Stat_periph, vmin=vmin, vmax=vmax)
title('Stationary Periphery')
ylabel('Stimulus SF')
xlabel('Speed')

subplot(2,3,3)
imshow(G_Run_fovea / G_Stat_fovea, vmin=.8, vmax=1.2)
title('Running / Stationary Fovea')
ylabel('Stimulus SF')
xlabel('Speed')

subplot(2,3,6)
imshow(G_Run_periph / G_Stat_periph, vmin=.8, vmax=1.2)
title('Running / Stationary Periphery')
ylabel('Stimulus SF')

subplots_adjust(hspace=.5, wspace=.5)



%% Do again using the histogram method

t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, valid=true)
Prun = get_Q_histogram(D, t_starts, t_ends, 
            omega_x_values = omega_x_values,
            omega_t_values = omega_t_values,
            Twin=512)

t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, valid=false)
Pstat = get_Q_histogram(D, t_starts, t_ends,
            omega_x_values = omega_x_values,
            omega_t_values = omega_t_values,
            Twin=512)


%% visualize the distribution of eye trajectories in each condition
fun = lambda x: 10*log10(x + 1e-6)
%matplotlib inline

vmin = minimum(fun(Prun['Htraj']).min(),lambda(Pstat['Htraj']).min())
vmax = maximum(fun(Prun['Htraj']).max(),lambda(Pstat['Htraj']).max())
figure(figsize=(9,3))
subplot(1,3,1)
contourf(Prun['x'], Prun['t'],lambda(Prun['Htraj']), 20, vmin=vmin, vmax=vmax)
xlabel('Position (deg)')
ylabel('Time (ms)')
title('Running')

subplot(1,3,2)
contourf(Pstat['x'], Pstat['t'],lambda(Pstat['Htraj']), 20, vmin=vmin, vmax=vmax)
xlabel('Position (deg)')
ylabel('Time (ms)')
title('Stationary')

% plot difference
subplot(1,3,3)
contourf(Pstat['x'], Pstat['t'],lambda(Pstat['Htraj']) -lambda(Prun['Htraj']), 20)
xlabel('Position (deg)')
ylabel('Time (ms)')
title('Difference')

subplots_adjust(hspace=.5, wspace=.5)

%% plot the same thing but for the power spectrum
fun = lambda x: 10*log10(x + 1e-6)
xfun = lambda x: log2(x + 1)
xinv = lambda x: 2.^x - 1
%matplotlib inline

iix = (Prun['fx'] >= 0) & (Prun['fx'] <= 32)
iiy = (Prun['ft'] >= 0) & (Prun['ft'] <= 128)

vmin = minimum(fun(Prun['Pspec'][iiy,:][:,iix]).min(),lambda(Pstat['Pspec'][iiy,:][:,iix]).min())
vmax = maximum(fun(Prun['Pspec'][iiy,:][:,iix]).max(),lambda(Pstat['Pspec'][iiy,:][:,iix]).max())

figure(figsize=(9,3))
subplot(1,3,1)
contourf(xfun(Prun['fx'](iix)), xfun(Prun['ft'][iiy]),lambda(Prun['Pspec'][iiy,:][:,iix]), 40, vmin=vmin, vmax=vmax)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Q (Running)')
set(gca,'XTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'XTick')])
set(gca,'YTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'YTick')])

subplot(1,3,2)

contourf(xfun(Pstat['fx'](iix)), xfun(Pstat['ft'][iiy]),lambda(Pstat['Pspec'][iiy,:][:,iix]), 40, vmin=vmin, vmax=vmax)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Q (Stationary)')
set(gca,'XTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'XTick')])
set(gca,'YTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'YTick')])

% plot difference
subplot(1,3,3)
% contourf(Pstat['fx'](iix), Pstat['ft'][iiy],lambda(Prun['Pspec'][iiy,:][:,iix])-fun(Pstat['Pspec'][iiy,:][:,iix]), 40)
contourf(xfun(Prun['fx'](iix)), xfun(Prun['ft'][iiy]), Prun['Pspec'][iiy,:][:,iix]/Pstat['Pspec'][iiy,:][:,iix], 40)
% set x and y tick labels to be the inverse of xfun
set(gca,'XTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'XTick')])
set(gca,'YTickLabel')(['%.1f' % (xinv(x)) for x in get(gca,'YTick')])

xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Relative Gain')
colorbar()

subplots_adjust(hspace=.5, wspace=.5)


%%
fx, ft = meshgrid(Prun['fx'], Prun['ft'])
fx[fx<.005] = .005
velocity = ft/fx

subplot(1,2,1)
contourf((Prun['fx'](iix)), (Prun['ft'][iiy]), ft[iiy,:][:,iix]/fx[iiy,:][:,iix], 40)
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Speed')
colorbar()

speed = 2
subplot(1,2,2)
dxdt = max(diff(velocity, axis=0), axis=0)
delta_approx = abs(velocity - speed) < dxdt[None,:]
% delta_approx /= delta_approx.sum(0)[None,:]
contourf((Prun['fx'](iix)), (Prun['ft'][iiy]), delta_approx[iiy,:][:,iix], 40)

spatial_freq_axis = Prun['fx'](iix)
temporal_freq_axis = Prun['ft'][iiy]
Pstim = delta_approx[iiy,:][:,iix]
Qrun = Prun['Pspec'][iiy,:][:,iix]
Qstat = Pstat['Pspec'][iiy,:][:,iix]
Qrun = Qrun / Qrun.max()
Qstat = Qstat / Qstat.max()

% convolve Pstim and Qrun along axis 0
Rrun = zeros_like(Qrun)
Rstat = zeros_like(Qstat)

for i = 1:length(spatial_freq_axis)
    Rrun[:,i] = convolve(Qrun[:,i], Pstim[:,i], mode='same')
    Rstat[:,i] = convolve(Qstat[:,i], Pstim[:,i], mode='same')
    % Rrun[:,i] = convolve(Pstim[:,i], Qrun[:,i], mode='same')
    % Rstat[:,i] = convolve(Pstim[:,i], Qstat[:,i], mode='same')
end
figure()
subplot(1,2,1)
contourf(spatial_freq_axis, temporal_freq_axis, Rrun, 40)
xlim([0,10])
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')
title('Running')
colorbar()

subplot(1,2,2)
contourf((spatial_freq_axis), (temporal_freq_axis), Rstat, 40)
xlim([0,10])
xlabel('SF (cyc/deg)')
ylabel('TF (cyc/s)')

title('Stationary')
%%

from scipy.signal import convolve

stim_sfs = array([1, 2, 4])

fun = lambda x: 10*log10(x + 1e-6)

for i = 1:(length(stim_sfs))

    sf = stim_sfs(i)
    ix_sf = argmin(abs(spatial_freq_axis - sf))

    h = plot(temporal_freq_axis,fun(Rrun[:,ix_sf]), '-', label='Running SF = %d' % sf)
    % plot a line that is the same color, but dashed
    plot(temporal_freq_axis,lambda(Rstat[:,ix_sf]), '--', color=h(1).get_color(), label='Stationary SF = %d' % sf)
    plot(temporal_freq_axis, Pstim[:,ix_sf], ':', color=h(1).get_color(), label='Stimulus SF = %d' % sf)
    
    xlabel('TF (cyc/s)')
end

legend()


%%
% title('SF = %d' % sf)
% convolve

% plot()
% plot(Pstat['Pspec'][:,ix_sf].mean(axis=1))


% Example usage
sf = log2(maximum(omega_x_values,0)+1e-6)
tf = log2(maximum(omega_t_values,0)+1e-6)
sf_, tf_ = meshgrid(sf, tf)

% fovea_params = {'sf_star': 1.0, 'sigma_s': 1.4, 'zeta_s': 0.25}
% periphery_params = {'sf_star': 0.5, 'sigma_s': 1.4, 'zeta_s': 0.15}
fovea_params = {'sf_star': 1, 'sigma_s': 1.5, 'zeta_s': 0.25, 
                  'tf_star': 4, 'sigma_t': 1.7, 'zeta_t': 0.18, 'Q':1}
periphery_params = {'sf_star': .5, 'sigma_s': 1.5, 'zeta_s': .29, 
                      'tf_star': 4, 'sigma_t': 2.8, 'zeta_t': 0.23, 'Q':1}

R_s_fovea = compute_R_s(sf, .^fovea_params)
R_s_periph = compute_R_s(sf, .^periphery_params)

plot(2.^sf, R_s_fovea, label='Fovea')
plot(2.^sf, R_s_periph, label='Periphery')
xlabel('Spatial Frequency (cyc/deg)')
ylabel('R_s')
legend()

R_fovea = compute_R_st(sf_, tf_, .^fovea_params)
R_periph = compute_R_st(sf_, tf_, .^periphery_params)

figure(figsize=(10, 5))
subplot(121)
contourf(2.^sf_, 2.^tf_, R_fovea, 30)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')
title('Fovea')
xlim([0, 10])
ylim([0, 30])

subplot(122)
contourf(2.^sf_, 2.^tf_, R_periph, 30)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')
title('Periphery')
xlim([0, 10])
ylim([0, 30])

%%


% %%
import numpy as np
import matplotlib.pyplot as plt

% Function to compute Q(k, omega)
function Q_k_omega=compute_Q(u_t, k_values, omega_values, dt)
%     """
%     Computes Q(k, omega) from a measured xi(t).
%     
%     Parameters:
%     - u_t: numpy array, the measured eye trajectory u(t).
%     - k_values: numpy array, the spatial frequency values k.
%     - omega_values: numpy array, the temporal frequency values omega.
%     - dt: float, the time step between measurements of u(t).
%     
%     Returns:
%     - Q_k_omega: 2D numpy array, Q(k, omega) for each k and omega.
%     
    
    % Time vector
    t = arange(length(u_t)) * dt
    
    % Initialize Q(k, omega)
    Q_k_omega = zeros((length(k_values), length(omega_values)), dtype=complex128)
    
    % Iterate over each k and omega
    for k = (k_values)
        for omega =(omega_values)
            % Compute the modulated signal
            modulated_signal = exp(-2j * pi * k * u_t) * exp(-2j * pi * omega * t)
            
            % Integrate over time (Fourier transform)
            integrated_value = trapz(modulated_signal, t)
            
            % Compute the magnitude squared
            Q_k_omega[i, j] = abs(integrated_value) .^ 2
        end
    end
end

    
 

% Function to compute Q(k, omega)
function Q_k_omega = compute_Q_vec(u_t, k_values, omega_values, dt)
%     """
%     Computes Q(k, omega) from a measured xi(t).
%     
%     Parameters:
%     - u_t: numpy array, the measured eye trajectory u(t).
%     - k_values: numpy array, the spatial frequency values k.
%     - omega_values: numpy array, the temporal frequency values omega.
%     - dt: float, the time step between measurements of u(t).
%     
%     Returns:
%     - Q_k_omega: 2D numpy array, Q(k, omega) for each k and omega.
%     """
    
    % Time vector
    t = arange(length(u_t)) * dt
    
    % Initialize Q(k, omega)
    Q_k_omega = zeros((length(k_values), length(omega_values)), dtype=complex128)
    
    % Iterate over each k and omega
    term1 = exp(-2j * pi * k_values[None,:,None] * u_t[:,None,None])
    term2 = exp(-2j * pi * omega_values[None,None,:] * t[:,None,None])
    modulated_signal = term1 * term2
    
            
    % Integrate over time (Fourier transform)
    integrated_value = trapz(modulated_signal, t[:,None,None], axis=0)
            
    
    Q_k_omega = abs(integrated_value) .^ 2
    
end

% Function to simulate Brownian motion for xi(t)
function u_t=simulate_brownian_motion(n_steps, dt, delta)
%     """
%     Simulates 1D Brownian motion.
%     
%     Parameters:
%     - n_steps: int, the number of time steps.
%     - dt: float, the time step size.
%     - delta: float, the standard deviation of the increments.
%     
%     Returns:
%     - u_t: numpy array, the simulated Brownian motion trajectory.
%     """
    u_t = zeros(n_steps)
    increments = random.normal(0, delta * sqrt(dt), size=n_steps-1)
    u_t[1:] = cumsum(increments)
end

% Parameters for simulation
n_steps = 512  % Number of time steps
dt = 0.001  % Time step between samples in seconds
delta = 5.0  % Standard deviation of the increments in the Brownian motion
n_fix = 500  % Number of simulations to average over

Q_k_omega_real = 0
% functionine ranges for spatial and temporal frequencies
k_values = 2.^(linspace(-5,2,50)) % Spatial frequencies
omega_values = 2.^(linspace(-5,5,50))%linspace(0.1, 32, 100)  % Temporal frequencies

for i = 1:(n_fix)
    % Simulate the trajectory u(t)
    u_t = simulate_brownian_motion(n_steps, dt, delta)

    % Compute the spectral power Q(k, omega)
    Q_k_omega = compute_Q_vec(u_t, k_values, omega_values, dt)

    % Convert to real values using the magnitude for plotting
    Q_k_omega_real += abs(Q_k_omega)
end

% Average the spectral power over all simulations
Q_k_omega_real = Q_k_omega_real/n_fix

%% Plot the spectral power as a heatmap
figure(figsize=(6, 4))
contourf(k_values, omega_values, 10*log10(Q_k_omega_real+1e-5).T, 100)
colorbar(label='dB')
xlabel('Spatial Frequency k')
ylabel('Temporal Frequency ω')
title('Spectral Power Q(k, ω) for Brownian Motion')
show()

% %% V1 derivative filter

import matplotlib.pyplot as plt
import numpy as np

% functionine the function parameters
omega = linspace(0.1, 30, 500)  % Frequency range
omega_0 = 4.0  % Reference frequency
b = .5  % Exponent parameter

% Calculate the function
r_omega = (omega / omega_0) * exp(-0.5 * (omega / omega_0).^2)
r_omega_b = r_omega .^ b

% Plot the function
figure(figsize=(8, 6))
plot(omega, r_omega_b, label=r'$r_{\omega}(\omega; \omega_0, b)$')
xlabel(r'$\omega$', fontsize=14)
ylabel(r'$r_{\omega}(\omega; \omega_0, b)$', fontsize=14)
title(r'Plot of $r_{\omega}(\omega; \omega_0, b)$', fontsize=16)
legend()
grid(true)
show()


% %%

import matplotlib.pyplot as plt
import numpy as np


function r_omega_b=gaussian_derivative_freq(omega_x, omega_t, omega_x_0, omega_t_0, b)
%     """
%     Compute the 2D Gaussian derivative function in the frequency domain.
%     
%     Parameters:
%     - omega_x: numpy array, the spatial frequency values in x.
%     - omega_t: numpy array, the temporal frequency values in t.
%     - omega_x_0: float, the reference frequency in x.
%     - omega_t_0: float, the reference frequency in t.
%     - b: float, the exponent parameter.
%     
%     Returns:
%     - r_omega_b: 2D numpy array, the 2D Gaussian derivative function.
%     """
    % Create a meshgrid for the 2D frequency space
    omega_x, omega_t = meshgrid(omega_x, omega_t)

    % Calculate the 2D Gaussian derivative function
    r_omega = (omega_x / omega_x_0) * (omega_t / omega_t_0) * exp(
        -0.5 * (omega_x / omega_x_0).^2 - 0.5 * (omega_t / omega_t_0).^2)
    r_omega_b = r_omega .^ b
    r_omega_b[isnan(r_omega_b)] = 0
end


omega_x_0 = 1  % Reference frequency for x
omega_t_0 = 4.0  % Reference frequency for t
b = .5  % Exponent parameter
r_omega_b = gaussian_derivative_freq(omega_x_values, omega_t_values, omega_x_0, omega_t_0, b)

% Plot the function
figure(figsize=(8, 6))
contourf(omega_x_values, omega_t_values, r_omega_b, 10, cmap='viridis')
colorbar(label=r'$r_{\omega}(\omega_x, \omega_t; \omega_{x_0}, \omega_{y_0}, b)$')
xlabel(r'$\omega_x$', fontsize=14)
ylabel(r'$\omega_t$', fontsize=14)
title(r'2D Gaussian Derivative in Frequency Domain', fontsize=16)
% xlim(0, 8)
ylim(-10, 30)
grid(true)

show()


% %%
import matplotlib.pyplot as plt
import numpy as np

function  r_omega_b=gaussian_derivative_RF(omega_x, omega_y, omega_x_0=1.0, omega_y_0=1.0, b=1.0, d=0)
%     """
%     2D Gaussian derivative function in the frequency domain.
%     
%     Parameters:
%     - omega_x: numpy array, the spatial frequency values in x.
%     - omega_y: numpy array, the spatial frequency values in y.
%     - omega_x_0: float, the reference frequency in x.
%     - omega_y_0: float, the reference frequency in y.
%     - b: float, the exponent parameter.
%     - d: float, the direction selectivity parameter.
% 
%     Returns:
%     - r_omega_b: 2D numpy array, the 2D Gaussian derivative function.
%     """
    
    omega_x, omega_y = meshgrid(omega_x, omega_y)

    % Calculate 
    r_omega_x = (omega_x / omega_x_0)  * exp(
        -0.5 * (omega_x / omega_x_0).^2 - 0.5 * (omega_y / omega_y_0).^2)
    r_omega_y = (omega_y / omega_y_0)
    
    % Calculate the 2D Gaussian derivative function
    r_omega = (omega_x / omega_x_0) * (omega_y / omega_y_0) * exp(
        -0.5 * (omega_x / omega_x_0).^2 - 0.5 * (omega_y / omega_y_0).^2)
    r_omega = abs(r_omega)
    d_selec = (1+d*sign(omega_x*omega_y))/2
    r_omega_b = d_selec*r_omega .^ b
end


% functionine the function parameters
omega_x = linspace(-5, 5, 500)
omega_y = linspace(-35, 35, 500)
omega_x_0 = 2.5  % Reference frequency for x
omega_y_0 = 4.0  % Reference frequency for y
b = 1  % Exponent parameter
d = 0 % direction selectivity

r_omega_b = gaussian_derivative_RF(omega_x, omega_y, omega_x_0, omega_y_0, b, d)
% Plot the function
figure(figsize=(8, 6))
contourf(omega_x, omega_y, r_omega_b, 50, cmap='viridis')
colorbar()
xlabel(r'$\omega_x$', fontsize=14)
ylabel(r'$\omega_y$', fontsize=14)
title(r'2D Gaussian Derivative in Frequency Domain', fontsize=16)
% grid(true)
show()

 %}

function [t_run_starts, t_run_ends] = get_run_epochs (D, plotfig, medfilt_window, speed_threshold)
%     '''
%     Get run epochs from treadmill data
%     '''
%     from scipy.interpolate import interp1d
%     from scipy.signal import savgol_filter, medfilt

    x_orig = D.treadSpeed;
    t_orig = D.treadTime;
    
    % interpolate nans in t_orig
%      iix = ~isnan(t_orig);
%      t_orig = interp1(find(iix), t_orig(iix),1:length(t_orig));
    % t_orig = t_orig(iix),

%     % x_orig_diff = diff(x_orig, prepend=0)
%     t_orig_diff = [0 diff(t_orig)]; % don't use
%     dt = median(t_orig_diff)

    % dx = x_orig %x_orig_diff 
    % breaks = where(isnan(dx))(1)
    % breaks = where(dx < -10)(1)

    % x_new = x_orig.copy()

    % for b in breaks:
    %     if isclose(x_orig[b], 0, atol=.1):
    %         x_new[b:] += x_orig[b-1]

    % calculate the derivative
    % dx_new = savgol_filter(x_new.(), 251, 1, deriv=1, delta=dt)
    % median filter speed

    dx_new = x_orig;
    medspeed = medfilt1(abs(dx_new), medfilt_window);

    % find run epochs
    run_epochs = (medspeed>speed_threshold);                      

    run_starts = find(diff(run_epochs)==1);
    run_ends = find(diff(run_epochs)==-1);

    if run_starts(1) > run_ends(1)
        % concatenate 0 to the fron of run_starts
        run_starts = [0; run_starts];
    end
    if run_ends(end) < run_starts(end)
        % concatenate last index to the end of run_ends
        run_ends = [run_ends; length(run_epochs)];
    end

    % convert to time 
    t_run_starts = t_orig(run_starts);
    t_run_ends = t_orig(run_ends);

    good_ix = ~(isnan(t_run_starts) | isnan(t_run_ends));
    t_run_starts = t_run_starts(good_ix);
    t_run_ends = t_run_ends(good_ix);

    % loop over and combine epochs that are too close together
    refractory = 5; % seconds
    i = 1;
    while i < length(t_run_starts)
        if t_run_starts(i+1) - t_run_ends(i) < refractory
            t_run_ends(i) = [];
            t_run_starts(i+1) = [];
        else
            i = i+ 1;
        end
    end

    % plot
    if plotfig==1
        figure()
        plot(t_orig, x_orig, '-')
        plot(t_orig, dx_new, '-')
        plot(t_orig, medspeed, '-')
        plot(t_orig, 5*run_epochs, '-')
        xlabel('Time (s)')
        ylabel('Velocity (cm/s)')

        for i = 1:(length(t_run_starts))
            % plot vertical fill between run epochs
            %fill_betweenx([0,max(x_orig)], t_run_starts(i), t_run_ends(i), color='lightblue', alpha=0.5)
        end
    end
end
    




function [t_good_starts, t_good_ends]=get_good_eye_epochs(D, plotfig)
%     '''
%     Get good eye epochs from eye data
%     '''
%     from scipy.interpolate import interp1d
    % figure()
    x_orig = D.eyePos(:,1);%interp1(D.eyeTime, D.eyePos(:,1), D.eyeTime);
    y_orig = D.eyePos(:,2);%interp1(D.eyeTime, D.eyePos(:,2), D.eyeTime);
    x_orig = x_orig - nanmedian(x_orig);
    y_orig = y_orig - nanmedian(y_orig);

    dx = [0; diff(x_orig)];
    dt = median([0; diff(D.eyeTime)]);
    % plot(D.eyeTime(),x_orig)
    % plot(D.eyeTime(), y_orig)


    numdiff = @(x) filter([1;-1],1,x);
    velpx = hypot(numdiff(x_orig), numdiff(y_orig));
    
    %figure(1); clf
    bad = filter(ones(10,1), 1, velpx>12)>0;

    bad_samples = (abs(x_orig) > 15)| (abs(y_orig) > 15)|(isnan(x_orig))|bad;
% 
%     bad_samples = logical_or(abs(x_orig) > 10, abs(y_orig) > 8)
%     bad_samples = logical_or(bad_samples, isnan(x_orig)).astype(float64())
    
    % find bad sample starts and stops
%     bad_starts = where(diff(bad_samples)==1)[0]
%     bad_ends = where(diff(bad_samples)==-1)[0]
% 
%     if bad_starts[0] > bad_ends[0]:
%         % concatenate 0 to the fron of run_starts
%         bad_starts = concatenate(([0], bad_starts))
%     end
%     if bad_ends[-1] < bad_starts[-1]:
%         % concatenate last index to the end of run_ends
%         bad_ends = concatenate((bad_ends, [length(bad_samples)-1]))
%     end
%                  
    bad_starts = find(diff(bad_samples)==1);
    bad_ends = find(diff(bad_samples)==-1);

    if bad_starts(1) > bad_ends(1)
        % concatenate 0 to the fron of bad_starts
        bad_starts = [0; bad_starts];
    end
    if bad_ends(end) < bad_starts(end)
        % concatenate last index to the end of bad_ends
        bad_ends = [bad_ends; length(bad_epochs)];
    end

    % offset the starts and stops by a 5 samples
    bad_starts = bad_starts - 50;
    bad_ends = bad_ends + 50;

    % make sure none are less than zero or greater than the length of the signal
    bad_starts = max(ones(length(bad_starts),1), bad_starts);
    bad_ends = min(length(bad_samples)*ones(length(bad_ends),1), bad_ends);
    

    % combine bad epochs that are too close together
    refractory = 5; %refractory = (.5/dt); % samples ()
    i = 1;
    while i < length(bad_starts)
        if bad_starts(i+1) - bad_ends(i) < refractory
            bad_ends(i) =[];
            bad_starts(i+1) = [];
        else
            i = i+ 1;
        end
    end

    % plot
    if plotfig==1
        figure()
        plot(D.eyeTime, x_orig)
        plot(D.eyeTime, y_orig)
    end


    % replace x_orig and y_orig with nans for bad epochs
    for i = 1:length(bad_starts)
        x_orig(bad_starts(i):bad_ends(i)) = nan;
        y_orig(bad_starts(i):bad_ends(i)) = nan;
    end

    % find good epochs by finding non-nan values
    good_epochs = ~isnan(x_orig);
    % convert to float64
    good_epochs = double(good_epochs);
    good_starts = find(diff(good_epochs)==1);
    good_ends = find(diff(good_epochs)==-1);

    if good_starts(1) > good_ends(1)
        % concatenate 0 to the fron of run_starts
        good_starts = [1; good_starts];
    end
    if good_ends(end) < good_starts(end)
        % concatenate last index to the end of run_ends
        good_ends = [good_ends; length(good_ends)];
    end

    % plot
    if plotfig==1
        for i = 1:length(good_starts)
            % plot vertical fill between run epochs
            %fill_betweenx([-20,20], D.eyeTime()[good_starts(i)], D.eyeTime()[good_ends(i)], color='lightblue', alpha=0.5)
        end
    end
    % return good epoch starts and stops as times
    t_good_starts = D.eyeTime(good_starts);
    t_good_ends = D.eyeTime(good_ends);

end

function [t_starts, t_ends]=get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, union_)%union=true
%     '''
%     Find the instersection of run and eye epochs
%     union=true -> returns run epochs where eye is valid
%     union=false -> returns stationary epochs where eye is valid
%     '''
    t_0 = min(t_run_starts(1), t_good_starts(1));
    t_1 = max(t_run_ends(end), t_good_ends(end));
    t = t_0:1e-3:t_1;%arange(t_0, t_1, 1e-3);
    run_epochs = zeros(size(t));
    good_epochs = zeros(size(t));
    for i = 1:length(t_run_starts)
        run_epochs = run_epochs | (t >= t_run_starts(i)) & (t <= t_run_ends(i));
    end
    for i = 1:length(t_good_starts)
        good_epochs = good_epochs | (t >= t_good_starts(i)) & (t <= t_good_ends(i));
    end
    % if valid, return only epochs where both run and eye are valid
    if union_
        epochs = run_epochs & good_epochs;
    else
        epochs = (~run_epochs) & good_epochs;
    end

%     % find starts and stops
%     starts = where(diff(epochs)==1)(1)
%     ends = where(diff(epochs)==-1)(1)
% 
%     if starts(1) > ends(1):
%         % concatenate 0 to the fron of run_starts
%         starts = concatenate(([0], starts));
%     end
%     if ends(end) < starts(end):
%         % concatenate last index to the end of run_ends
%         ends = concatenate((ends, [length(epochs)-1]));
%     end

    starts = find(diff(epochs)==1);
    ends = find(diff(epochs)==-1);

    if starts(1) > ends(1)
        % concatenate 0 to the fron of starts
        starts = [1 starts];
    end
    
%     if ends(end) < starts(end)
%         % concatenate last index to the end of ends
%         ends = [ends; length(epochs)];
%     end
    if ends(end) < starts(end)
        %better to just remove last start
        starts=starts(1:end-1);
    end
    assert(length(starts)==length(ends));

    % convert to time
    t_starts = t(starts);
    t_ends = t(ends);

end    

%% Function to compute Q(k, omega)
function Q_k_omega= compute_Q(u_t, k_values, omega_values, dt)
%     """
%     Computes Q(k, omega) from a measured xi(t).
%     
%     Parameters:
%     - u_t: numpy array, the measured eye trajectory u(t).
%     - k_values: numpy array, the spatial frequency values k.
%     - omega_values: numpy array, the temporal frequency values omega.
%     - dt: float, the time step between measurements of u(t).
%     
%     Returns:
%     - Q_k_omega: 2D numpy array, Q(k, omega) for each k and omega.
%     """
    
    % Time vector
    t = (1:length(u_t)) * dt;
    
    % Initialize Q(k, omega)
    Q_k_omega = zeros(length(k_values), length(omega_values));
    
    % Iterate over each k and omega
    term1 = exp(-2j * pi * reshape(k_values,1,1,[]) .* reshape(u_t,[],1,1));
    term2 = exp(-2j * pi * reshape(omega_values,1,[],1) .* reshape(t,[],1,1));
    modulated_signal = term1 .* term2;
    
    % Integrate over time (Fourier transform)
    integrated_value = squeeze(trapz(t,modulated_signal, 1));
%     integrated_value=squeeze(sum(modulated_signal,1));
    
    Q_k_omega = abs(integrated_value) * 2;
end



function xsegments=get_eye_displacements(D, t_starts, t_ends, Twin, noverlap)%Twin = 512, noverlap=None)
    
%     '''
%     get segments of eye position data in delta x, delta t
% 
%     '''
    % D is the data dictionary
    % t_starts, t_ends are the start and end times of the epochs
    % Twin (integer) in samples
    % noverlap = number of samples to overlap

%     xpos = D.eyePos[:,0].()
    xpos = D.eyePos;
    xsegments = [];

    if noverlap == 0
        noverlap = Twin / 2;
    end
    % loop over epochs
    for i = 1:length(t_starts)
        t0 = t_starts(i);
        t1 = t_ends(i);
        iix = (D.eyeTime >= t0) & (D.eyeTime <= t1);

        xsegment = xpos(iix);
        Ns = length(xsegment);

        step = Twin - noverlap;
        
        % loop over segments stepping by noverlap
        for j = 1: step:Ns-Twin
            xseg = xsegment(j:j+Twin);

            % subtract first sample to calculate displacement
            xseg = xseg - xseg(1);

            xsegments=[xsegments xseg];
        end
    end
end



function y=xfun(x)
    y = zeros(size(x),'like',x);
    iix = x >=0;
    y(iix) = log2(x(iix) + 1);
    y(~iix) = -log2(-x(~iix) + 1);
end

function y=xinv(x)
    y = zeros(size(x),'like',x);
    iix = x >=0;
    y(iix) = 2.^x(iix) - 1;
    y(~iix) = -2.^(-x(~iix)) + 1;
end

