
#%%
import os
import numpy as np
import matplotlib.pyplot as plt
from huk_treadmill import load_file, get_run_epochs, get_good_eye_epochs, get_epochs, compute_Q, get_eye_displacements

#%% Load example file
# fname = '/Users/jake/Downloads/Eyesrunning_Rocky20240802_reD.mat'
fdir = '/Users/jake/Dropbox/Datasets/HuklabTreadmill/gratings/'
flist = os.listdir(fdir)
flist = [f for f in flist if f.endswith('.mat')]

subject = 'brie'

flist = [f for f in flist if subject in f]

x_segments_run = []
x_segments_stat = []

for i_sess in range(len(flist)):
    fname = os.path.join(fdir, flist[i_sess])
    print("Loading", fname)

    D = load_file(fname)

    fps = int(1/np.median(np.diff(D['eyeTime'].flatten())))

    if fps < 999:
        print('Warning: low frame rate detected: %d' % fps)
        continue
    
    # get run and eye epochs
    t_run_starts, t_run_ends = get_run_epochs(D, speed_threshold=2.0, plot=False)
    t_good_starts, t_good_ends = get_good_eye_epochs(D, plot=False)

    t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, union=True)
    xsegments_run_ = get_eye_displacements(D, t_starts, t_ends, Twin = 512, noverlap=512//4)

    t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, union=False)
    xsegments_stat_ = get_eye_displacements(D, t_starts, t_ends, Twin = 512, noverlap=512//4)

    x_segments_run.append(xsegments_run_)
    x_segments_stat.append(xsegments_stat_)

#%%
xs_run = np.concatenate(x_segments_run)
xs_stat = np.concatenate(x_segments_stat)

print(xs_run.shape, xs_stat.shape)

#%%
%matplotlib ipympl
plt.figure()
plt.plot(D['treadTime'], D['treadSpeed'])
for i in range(len(t_run_starts)):
    plt.fill_betweenx([-1, np.nanmax(D['treadSpeed'])], t_run_starts[i], t_run_ends[i], color='r', alpha=.5)
plt.show()

#%%
%matplotlib inline
bins = np.arange(-15, 15, 1/60)
prun = np.zeros((len(bins)-1, xs_run.shape[1]))
pstat = np.zeros((len(bins)-1, xs_stat.shape[1]))
for i in range(xs_run.shape[1]):
    prun[:,i] = np.histogram(xs_run[:,i], bins=bins)[0]
    pstat[:,i] = np.histogram(xs_stat[:,i], bins=bins)[0]

from scipy.ndimage import gaussian_filter
prun = gaussian_filter(prun, 5)
pstat = gaussian_filter(pstat, 5)

plt.figure(figsize=(9,3))
plt.subplot(1,2,1)
plt.imshow(np.log10(prun+1e-3), aspect='auto', extent=[0, xs_run.shape[1], bins[0], bins[-1]])
plt.title('Running')
plt.xlabel('Time (ms)')
plt.ylabel('Position (deg)')

plt.subplot(1,2,2)
plt.imshow(np.log10(pstat+1e-3), aspect='auto', extent=[0, xs_stat.shape[1], bins[0], bins[-1]])
plt.title('Stationary')
plt.xlabel('Time (ms)')
plt.ylabel('Position (deg)')




#%% Loop over segments and calculate the eye trajectory displacement power spectrum

omega_x_values = 2**np.linspace(-7, 3, 50)
omega_t_values = 2**np.linspace(-7, 7, 30)
omega_t_values = np.concatenate((np.append(np.flip(-omega_t_values), 0), omega_t_values))

fps = 1000

Qrun = 0

nseg = xs_run.shape[0]
for i in range(nseg):
    xseg = xs_run[i,:]
    Qrun += compute_Q(xseg, omega_x_values, omega_t_values, 1/fps)

Qrun /= nseg

Qstat = 0
nseg = xs_stat.shape[0]
for i in range(nseg):
    xseg = xs_stat[i,:]
    Qstat += compute_Q(xseg, omega_x_values, omega_t_values, 1/fps)

Qstat /= nseg






#%% run main analysis

# get support for the frequency analysis
omega_x_values = 2**np.linspace(-7, 3, 50)
omega_t_values = 2**np.linspace(-7, 7, 30)
omega_t_values = np.concatenate((np.append(np.flip(-omega_t_values), 0), omega_t_values))
fx, ft = np.meshgrid(omega_x_values, omega_t_values)
_ = plt.scatter(fx, ft, s=np.hypot(fx,ft)/1000)
# _ = plt.plot(fx, ft, '.k', 'markersize', .01)
plt.xlabel('Spatial Frequency (cyc/deg)')
plt.ylabel('Temporal Frequency (cyc/s)')
plt.title('Frequency Support')


#%% plot the same thing but for the power spectrum
%matplotlib inline
fun = lambda x: 10*np.log10(x + 1e-6)
def xfun(x):
    y = np.zeros_like(x)
    iix = x >=0
    y[iix] = np.log2(x[iix] + 1)
    y[~iix] = -np.log2(-x[~iix] + 1)
    return y

def xinv(x):
    y = np.zeros_like(x)
    iix = x >=0
    y[iix] = 2**x[iix] - 1
    y[~iix] = -2**(-x[~iix]) + 1
    return y

vmin = np.minimum(fun(Qstat).min(), fun(Qrun).min())
vmax = np.maximum(fun(Qstat).max(), fun(Qrun).max())

plt.figure(figsize=(9,3))
plt.subplot(1,3,1)
plt.contourf(xfun(omega_x_values), xfun(omega_t_values), fun(Qrun), levels=np.linspace(vmin, vmax, 25), vmin=vmin, vmax=vmax)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Q (Running)')
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

plt.subplot(1,3,2)

plt.contourf(xfun(omega_x_values), xfun(omega_t_values), fun(Qstat), levels=np.linspace(vmin, vmax, 25), vmin=vmin, vmax=vmax)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Q (Stationary)')
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

# plot difference
plt.subplot(1,3,3)
plt.contourf(xfun(omega_x_values), xfun(omega_t_values), fun(Qrun) - fun(Qstat), levels=25) #, vmin=.5, vmax=1.5, levels=np.linspace(.5, 1.5, 25))
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Relative Gain')
plt.colorbar()

plt.subplots_adjust(hspace=.5, wspace=.5)

#%% plot Q at different TF slices

tf0 = 2
for tf in [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]:
    iitf = np.argmin(np.abs(omega_t_values - tf - tf0))
    h = plt.plot(omega_x_values, fun(Qrun[iitf,:]), label='%.1f Hz' % tf)
    # draw line usin same color
    plt.plot(omega_x_values, fun(Qstat[iitf,:]), '--', color=h[0].get_color())
    # plt.ylim(-20, 0)

# for tf in [-0.5, -1.0, -2.0, -4.0, -8.0, -12.0]:
#     iitf = np.argmin(np.abs(omega_t_values - tf - tf0))
#     h = plt.plot(omega_x_values, fun(Qrun[iitf,:]), label='%.1f Hz' % tf)
#     # draw line usin same color
#     plt.plot(omega_x_values, fun(Qstat[iitf,:]), '--', color=h[0].get_color())
#     # plt.ylim(-20, 0)


plt.xscale('log')

for sf in [1.0, 2.0, 4.0]:
    plt.axvline(sf, color='k', linestyle='--')

plt.legend()
plt.xlabel('Spatial Frequency (cyc/deg)')


#%% get hypothetical neural tuning

from huk_treadmill import gaussian_derivative_RF

R_fovea = gaussian_derivative_RF(omega_x_values, omega_t_values, 
                                 omega_x_0=2, omega_y_0=3, b=.5, d=0)
R_periph = gaussian_derivative_RF(omega_x_values, omega_t_values, 
                                  omega_x_0=.5, omega_y_0=4, b=1, d=0)

# plot
plt.figure(figsize=(9,3))
plt.subplot(1,2,1)
plt.contourf((omega_x_values), (omega_t_values), R_fovea, levels=15)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Fovea')
plt.xlim((0,8))
plt.ylim((0, 15))

plt.subplot(1,2,2)
plt.contourf((omega_x_values), (omega_t_values), R_periph, levels=15)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Periphery')
plt.xlim((0,8))
plt.ylim((0, 15))


stim_sfs = [1, 2, 4]
speeds = [1.0, 2.0]
plt.figure(figsize=(9,9))

G_Run_fovea = np.zeros((len(stim_sfs), len(speeds)))
G_Run_periph = np.zeros((len(stim_sfs), len(speeds)))
G_Stat_fovea = np.zeros((len(stim_sfs), len(speeds)))
G_Stat_periph = np.zeros((len(stim_sfs), len(speeds)))

for istim in range(len(stim_sfs)):
    
    stim_sf = stim_sfs[istim]
    iisf = np.argmin(np.abs(omega_x_values - stim_sf))


    for i,speed in enumerate(speeds):
        plt.subplot(len(stim_sfs), len(speeds), istim*len(speeds)+i+1)
        
        tf = stim_sf * speed
        # iitf = np.argmin(np.abs(omega_t_values - tf))

        plt.plot(omega_t_values + tf, fun(Qrun[:,iisf]), 'k', label='Running')
        plt.plot(omega_t_values + tf, fun(Qstat[:,iisf]), 'k--', label='Stationary')
        plt.xlabel('Temporal Frequency (cycles/s)')
        plt.ylabel('dB')
        plt.legend()
        plt.title('Stimulus SF: %.1f cyc/deg, Speed: %.1f cyc/s' % (stim_sf, speed))
        plt.axvline(tf, color='k', linestyle='--')
        plt.xlim(0, 32)

        # new axis that has same x values, but new y values
        ax = plt.gca().twinx()
        ax.plot(omega_t_values, fun(R_fovea[:,iisf]), 'r', label='Fovea')
        ax.plot(omega_t_values, fun(R_periph[:,iisf]), 'b', label='Periphery')
        # ax.set_ylim(0, 1)
        ax.set_xlim(0, 32)

        G_Run_fovea[istim,i] = np.trapz(R_fovea[:,iisf] * Qrun[:,iisf], omega_t_values)
        G_Run_periph[istim,i] = np.trapz(R_periph[:,iisf] * Qrun[:,iisf], omega_t_values)
        G_Stat_fovea[istim,i] = np.trapz(R_fovea[:,iisf] * Qstat[:,iisf], omega_t_values)
        G_Stat_periph[istim,i] = np.trapz(R_periph[:,iisf] * Qstat[:,iisf], omega_t_values)

plt.subplots_adjust(hspace=.5, wspace=.5)

plt.figure()
plt.plot(stim_sfs, G_Run_periph*1e3, 'b')
plt.plot(stim_sfs, G_Stat_periph*1e3, 'r')
plt.xlabel('Spatial Frequency (cyc/deg)')

plt.figure()
_ = plt.plot(G_Run_periph / G_Stat_periph, 'bo-', label='Periphery')
_ = plt.plot(G_Run_fovea / G_Stat_fovea, 'ro-', label='Fovea')
plt.axhline(1, color='k', linestyle='--')
#%%
vmin = np.minimum(np.minimum(G_Run_fovea.min(), G_Run_periph.min()), np.minimum(G_Stat_fovea.min(), G_Stat_periph.min()))
vmax = np.maximum(np.maximum(G_Run_fovea.max(), G_Run_periph.max()), np.maximum(G_Stat_fovea.max(), G_Stat_periph.max()))

plt.figure(figsize=(9,3))
plt.subplot(2,3,1)
plt.imshow(G_Run_fovea, vmin=vmin, vmax=vmax)
plt.title('Running Fovea')
plt.ylabel('Stimulus SF')
plt.xlabel('Speed')

plt.subplot(2,3,2)
plt.imshow(G_Stat_fovea, vmin=vmin, vmax=vmax)
plt.title('Stationary Fovea')
plt.ylabel('Stimulus SF')
plt.xlabel('Speed')

plt.subplot(2,3,4)
plt.imshow(G_Run_periph, vmin=vmin, vmax=vmax)
plt.title('Running Periphery')
plt.ylabel('Stimulus SF')
plt.xlabel('Speed')

plt.subplot(2,3,5)
plt.imshow(G_Stat_periph, vmin=vmin, vmax=vmax)
plt.title('Stationary Periphery')
plt.ylabel('Stimulus SF')
plt.xlabel('Speed')

plt.subplot(2,3,3)
plt.imshow(G_Run_fovea / G_Stat_fovea, vmin=.8, vmax=1.2)
plt.title('Running / Stationary Fovea')
plt.ylabel('Stimulus SF')
plt.xlabel('Speed')

plt.subplot(2,3,6)
plt.imshow(G_Run_periph / G_Stat_periph, vmin=.8, vmax=1.2)
plt.title('Running / Stationary Periphery')
plt.ylabel('Stimulus SF')

plt.subplots_adjust(hspace=.5, wspace=.5)



#%% Do again using the histogram method

t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, valid=True)
Prun = get_Q_histogram(D, t_starts, t_ends, 
            omega_x_values = omega_x_values,
            omega_t_values = omega_t_values,
            Twin=512)

t_starts, t_ends = get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, valid=False)
Pstat = get_Q_histogram(D, t_starts, t_ends,
            omega_x_values = omega_x_values,
            omega_t_values = omega_t_values,
            Twin=512)


#%% visualize the distribution of eye trajectories in each condition
fun = lambda x: 10*np.log10(x + 1e-6)
%matplotlib inline

vmin = np.minimum(fun(Prun['Htraj']).min(), fun(Pstat['Htraj']).min())
vmax = np.maximum(fun(Prun['Htraj']).max(), fun(Pstat['Htraj']).max())
plt.figure(figsize=(9,3))
plt.subplot(1,3,1)
plt.contourf(Prun['x'], Prun['t'], fun(Prun['Htraj']), levels=20, vmin=vmin, vmax=vmax)
plt.xlabel('Position (deg)')
plt.ylabel('Time (ms)')
plt.title('Running')

plt.subplot(1,3,2)
plt.contourf(Pstat['x'], Pstat['t'], fun(Pstat['Htraj']), levels=20, vmin=vmin, vmax=vmax)
plt.xlabel('Position (deg)')
plt.ylabel('Time (ms)')
plt.title('Stationary')

# plot difference
plt.subplot(1,3,3)
plt.contourf(Pstat['x'], Pstat['t'], fun(Pstat['Htraj']) - fun(Prun['Htraj']), levels=20)
plt.xlabel('Position (deg)')
plt.ylabel('Time (ms)')
plt.title('Difference')

plt.subplots_adjust(hspace=.5, wspace=.5)

#%% plot the same thing but for the power spectrum
fun = lambda x: 10*np.log10(x + 1e-6)
xfun = lambda x: np.log2(x + 1)
xinv = lambda x: 2**x - 1
%matplotlib inline

iix = (Prun['fx'] >= 0) & (Prun['fx'] <= 32)
iiy = (Prun['ft'] >= 0) & (Prun['ft'] <= 128)

vmin = np.minimum(fun(Prun['Pspec'][iiy,:][:,iix]).min(), fun(Pstat['Pspec'][iiy,:][:,iix]).min())
vmax = np.maximum(fun(Prun['Pspec'][iiy,:][:,iix]).max(), fun(Pstat['Pspec'][iiy,:][:,iix]).max())

plt.figure(figsize=(9,3))
plt.subplot(1,3,1)
plt.contourf(xfun(Prun['fx'][iix]), xfun(Prun['ft'][iiy]), fun(Prun['Pspec'][iiy,:][:,iix]), levels=40, vmin=vmin, vmax=vmax)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Q (Running)')
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

plt.subplot(1,3,2)

plt.contourf(xfun(Pstat['fx'][iix]), xfun(Pstat['ft'][iiy]), fun(Pstat['Pspec'][iiy,:][:,iix]), levels=40, vmin=vmin, vmax=vmax)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Q (Stationary)')
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

# plot difference
plt.subplot(1,3,3)
# plt.contourf(Pstat['fx'][iix], Pstat['ft'][iiy], fun(Prun['Pspec'][iiy,:][:,iix])-fun(Pstat['Pspec'][iiy,:][:,iix]), levels=40)
plt.contourf(xfun(Prun['fx'][iix]), xfun(Prun['ft'][iiy]), Prun['Pspec'][iiy,:][:,iix]/Pstat['Pspec'][iiy,:][:,iix], levels=40)
# set x and y tick labels to be the inverse of xfun
plt.gca().set_xticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_xticks()])
plt.gca().set_yticklabels(['%.1f' % (xinv(x)) for x in plt.gca().get_yticks()])

plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Relative Gain')
plt.colorbar()

plt.subplots_adjust(hspace=.5, wspace=.5)


#%%
fx, ft = np.meshgrid(Prun['fx'], Prun['ft'])
fx[fx<.005] = .005
velocity = ft/fx

plt.subplot(1,2,1)
plt.contourf((Prun['fx'][iix]), (Prun['ft'][iiy]), ft[iiy,:][:,iix]/fx[iiy,:][:,iix], levels=40)
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Speed')
plt.colorbar()

speed = 2
plt.subplot(1,2,2)
dxdt = np.max(np.diff(velocity, axis=0), axis=0)
delta_approx = np.abs(velocity - speed) < dxdt[None,:]
# delta_approx /= delta_approx.sum(0)[None,:]
plt.contourf((Prun['fx'][iix]), (Prun['ft'][iiy]), delta_approx[iiy,:][:,iix], levels=40)

spatial_freq_axis = Prun['fx'][iix]
temporal_freq_axis = Prun['ft'][iiy]
Pstim = delta_approx[iiy,:][:,iix]
Qrun = Prun['Pspec'][iiy,:][:,iix]
Qstat = Pstat['Pspec'][iiy,:][:,iix]
Qrun = Qrun / Qrun.max()
Qstat = Qstat / Qstat.max()

# convolve Pstim and Qrun along axis 0
Rrun = np.zeros_like(Qrun)
Rstat = np.zeros_like(Qstat)

for i in range(len(spatial_freq_axis)):
    Rrun[:,i] = np.convolve(Qrun[:,i], Pstim[:,i], mode='same')
    Rstat[:,i] = np.convolve(Qstat[:,i], Pstim[:,i], mode='same')
    # Rrun[:,i] = np.convolve(Pstim[:,i], Qrun[:,i], mode='same')
    # Rstat[:,i] = np.convolve(Pstim[:,i], Qstat[:,i], mode='same')

plt.figure()
plt.subplot(1,2,1)
plt.contourf(spatial_freq_axis, temporal_freq_axis, Rrun, levels=40)
plt.xlim([0,10])
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')
plt.title('Running')
plt.colorbar()

plt.subplot(1,2,2)
plt.contourf((spatial_freq_axis), (temporal_freq_axis), Rstat, levels=40)
plt.xlim([0,10])
plt.xlabel('SF (cyc/deg)')
plt.ylabel('TF (cyc/s)')

plt.title('Stationary')
#%%

from scipy.signal import convolve

stim_sfs = np.array([1, 2, 4])

fun = lambda x: 10*np.log10(x + 1e-6)

for i in range(len(stim_sfs)):

    sf = stim_sfs[i]
    ix_sf = np.argmin(np.abs(spatial_freq_axis - sf))

    h = plt.plot(temporal_freq_axis,fun(Rrun[:,ix_sf]), '-', label='Running SF = %d' % sf)
    # plot a line that is the same color, but dashed
    plt.plot(temporal_freq_axis, fun(Rstat[:,ix_sf]), '--', color=h[0].get_color(), label='Stationary SF = %d' % sf)
    plt.plot(temporal_freq_axis, Pstim[:,ix_sf], ':', color=h[0].get_color(), label='Stimulus SF = %d' % sf)
    
    plt.xlabel('TF (cyc/s)')
    

plt.legend()


#%%
# plt.title('SF = %d' % sf)
# convolve

# plt.plot()
# plt.plot(Pstat['Pspec'][:,ix_sf].mean(axis=1))


# Example usage
sf = np.log2(np.maximum(omega_x_values,0)+1e-6)
tf = np.log2(np.maximum(omega_t_values,0)+1e-6)
sf_, tf_ = np.meshgrid(sf, tf)

# fovea_params = {'sf_star': 1.0, 'sigma_s': 1.4, 'zeta_s': 0.25}
# periphery_params = {'sf_star': 0.5, 'sigma_s': 1.4, 'zeta_s': 0.15}
fovea_params = {'sf_star': 1, 'sigma_s': 1.5, 'zeta_s': 0.25, 
                  'tf_star': 4, 'sigma_t': 1.7, 'zeta_t': 0.18, 'Q':1}
periphery_params = {'sf_star': .5, 'sigma_s': 1.5, 'zeta_s': .29, 
                      'tf_star': 4, 'sigma_t': 2.8, 'zeta_t': 0.23, 'Q':1}

R_s_fovea = compute_R_s(sf, **fovea_params)
R_s_periph = compute_R_s(sf, **periphery_params)

plt.plot(2**sf, R_s_fovea, label='Fovea')
plt.plot(2**sf, R_s_periph, label='Periphery')
plt.xlabel('Spatial Frequency (cyc/deg)')
plt.ylabel('R_s')
plt.legend()

R_fovea = compute_R_st(sf_, tf_, **fovea_params)
R_periph = compute_R_st(sf_, tf_, **periphery_params)

plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.contourf(2**sf_, 2**tf_, R_fovea, levels=30)
plt.xlabel('Spatial Frequency (cyc/deg)')
plt.ylabel('Temporal Frequency (cyc/s)')
plt.title('Fovea')
plt.xlim([0, 10])
plt.ylim([0, 30])

plt.subplot(122)
plt.contourf(2**sf_, 2**tf_, R_periph, levels=30)
plt.xlabel('Spatial Frequency (cyc/deg)')
plt.ylabel('Temporal Frequency (cyc/s)')
plt.title('Periphery')
plt.xlim([0, 10])
plt.ylim([0, 30])

#%%


# %%
import numpy as np
import matplotlib.pyplot as plt

# Function to compute Q(k, omega)
def compute_Q(u_t, k_values, omega_values, dt):
    """
    Computes Q(k, omega) from a measured xi(t).
    
    Parameters:
    - u_t: numpy array, the measured eye trajectory u(t).
    - k_values: numpy array, the spatial frequency values k.
    - omega_values: numpy array, the temporal frequency values omega.
    - dt: float, the time step between measurements of u(t).
    
    Returns:
    - Q_k_omega: 2D numpy array, Q(k, omega) for each k and omega.
    """
    
    # Time vector
    t = np.arange(len(u_t)) * dt
    
    # Initialize Q(k, omega)
    Q_k_omega = np.zeros((len(k_values), len(omega_values)), dtype=np.complex128)
    
    # Iterate over each k and omega
    for i, k in enumerate(k_values):
        for j, omega in enumerate(omega_values):
            # Compute the modulated signal
            modulated_signal = np.exp(-2j * np.pi * k * u_t) * np.exp(-2j * np.pi * omega * t)
            
            # Integrate over time (Fourier transform)
            integrated_value = np.trapz(modulated_signal, t)
            
            # Compute the magnitude squared
            Q_k_omega[i, j] = np.abs(integrated_value) ** 2
    
    return Q_k_omega

# Function to compute Q(k, omega)
def compute_Q_vec(u_t, k_values, omega_values, dt):
    """
    Computes Q(k, omega) from a measured xi(t).
    
    Parameters:
    - u_t: numpy array, the measured eye trajectory u(t).
    - k_values: numpy array, the spatial frequency values k.
    - omega_values: numpy array, the temporal frequency values omega.
    - dt: float, the time step between measurements of u(t).
    
    Returns:
    - Q_k_omega: 2D numpy array, Q(k, omega) for each k and omega.
    """
    
    # Time vector
    t = np.arange(len(u_t)) * dt
    
    # Initialize Q(k, omega)
    Q_k_omega = np.zeros((len(k_values), len(omega_values)), dtype=np.complex128)
    
    # Iterate over each k and omega
    term1 = np.exp(-2j * np.pi * k_values[None,:,None] * u_t[:,None,None])
    term2 = np.exp(-2j * np.pi * omega_values[None,None,:] * t[:,None,None])
    modulated_signal = term1 * term2
    
            
    # Integrate over time (Fourier transform)
    integrated_value = np.trapz(modulated_signal, t[:,None,None], axis=0)
            
    
    Q_k_omega = np.abs(integrated_value) ** 2
    
    return Q_k_omega

# Function to simulate Brownian motion for xi(t)
def simulate_brownian_motion(n_steps, dt, delta):
    """
    Simulates 1D Brownian motion.
    
    Parameters:
    - n_steps: int, the number of time steps.
    - dt: float, the time step size.
    - delta: float, the standard deviation of the increments.
    
    Returns:
    - u_t: numpy array, the simulated Brownian motion trajectory.
    """
    u_t = np.zeros(n_steps)
    increments = np.random.normal(0, delta * np.sqrt(dt), size=n_steps-1)
    u_t[1:] = np.cumsum(increments)
    return u_t

# Parameters for simulation
n_steps = 512  # Number of time steps
dt = 0.001  # Time step between samples in seconds
delta = 5.0  # Standard deviation of the increments in the Brownian motion
n_fix = 500  # Number of simulations to average over

Q_k_omega_real = 0
# Define ranges for spatial and temporal frequencies
k_values = 2**(np.linspace(-5,2,50)) # Spatial frequencies
omega_values = 2**(np.linspace(-5,5,50))#np.linspace(0.1, 32, 100)  # Temporal frequencies

for i in range(n_fix):
    # Simulate the trajectory u(t)
    u_t = simulate_brownian_motion(n_steps, dt, delta)

    # Compute the spectral power Q(k, omega)
    Q_k_omega = compute_Q_vec(u_t, k_values, omega_values, dt)

    # Convert to real values using the magnitude for plotting
    Q_k_omega_real += np.abs(Q_k_omega)

# Average the spectral power over all simulations
Q_k_omega_real /= n_fix

#%% Plot the spectral power as a heatmap
plt.figure(figsize=(6, 4))
plt.contourf(k_values, omega_values, 10*np.log10(Q_k_omega_real+1e-5).T, levels=100)
plt.colorbar(label='dB')
plt.xlabel('Spatial Frequency k')
plt.ylabel('Temporal Frequency ω')
plt.title('Spectral Power Q(k, ω) for Brownian Motion')
plt.show()

# %% V1 derivative filter

import matplotlib.pyplot as plt
import numpy as np

# Define the function parameters
omega = np.linspace(0.1, 30, 500)  # Frequency range
omega_0 = 4.0  # Reference frequency
b = .5  # Exponent parameter

# Calculate the function
r_omega = (omega / omega_0) * np.exp(-0.5 * (omega / omega_0)**2)
r_omega_b = r_omega ** b

# Plot the function
plt.figure(figsize=(8, 6))
plt.plot(omega, r_omega_b, label=r'$r_{\omega}(\omega; \omega_0, b)$')
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$r_{\omega}(\omega; \omega_0, b)$', fontsize=14)
plt.title(r'Plot of $r_{\omega}(\omega; \omega_0, b)$', fontsize=16)
plt.legend()
plt.grid(True)
plt.show()


# %%

import matplotlib.pyplot as plt
import numpy as np


def gaussian_derivative_freq(omega_x, omega_t, omega_x_0, omega_t_0, b):
    """
    Compute the 2D Gaussian derivative function in the frequency domain.
    
    Parameters:
    - omega_x: numpy array, the spatial frequency values in x.
    - omega_t: numpy array, the temporal frequency values in t.
    - omega_x_0: float, the reference frequency in x.
    - omega_t_0: float, the reference frequency in t.
    - b: float, the exponent parameter.
    
    Returns:
    - r_omega_b: 2D numpy array, the 2D Gaussian derivative function.
    """
    # Create a meshgrid for the 2D frequency space
    omega_x, omega_t = np.meshgrid(omega_x, omega_t)

    # Calculate the 2D Gaussian derivative function
    r_omega = (omega_x / omega_x_0) * (omega_t / omega_t_0) * np.exp(
        -0.5 * (omega_x / omega_x_0)**2 - 0.5 * (omega_t / omega_t_0)**2)
    r_omega_b = r_omega ** b
    r_omega_b[np.isnan(r_omega_b)] = 0
    return r_omega_b


omega_x_0 = 1  # Reference frequency for x
omega_t_0 = 4.0  # Reference frequency for t
b = .5  # Exponent parameter
r_omega_b = gaussian_derivative_freq(omega_x_values, omega_t_values, omega_x_0, omega_t_0, b)

# Plot the function
plt.figure(figsize=(8, 6))
plt.contourf(omega_x_values, omega_t_values, r_omega_b, levels=10, cmap='viridis')
plt.colorbar(label=r'$r_{\omega}(\omega_x, \omega_t; \omega_{x_0}, \omega_{y_0}, b)$')
plt.xlabel(r'$\omega_x$', fontsize=14)
plt.ylabel(r'$\omega_t$', fontsize=14)
plt.title(r'2D Gaussian Derivative in Frequency Domain', fontsize=16)
# plt.xlim(0, 8)
plt.ylim(-10, 30)
plt.grid(True)

plt.show()


# %%
import matplotlib.pyplot as plt
import numpy as np

def gaussian_derivative_RF(omega_x, omega_y, omega_x_0=1.0, omega_y_0=1.0, b=1.0, d=0):
    """
    2D Gaussian derivative function in the frequency domain.
    
    Parameters:
    - omega_x: numpy array, the spatial frequency values in x.
    - omega_y: numpy array, the spatial frequency values in y.
    - omega_x_0: float, the reference frequency in x.
    - omega_y_0: float, the reference frequency in y.
    - b: float, the exponent parameter.
    - d: float, the direction selectivity parameter.

    Returns:
    - r_omega_b: 2D numpy array, the 2D Gaussian derivative function.
    """
    
    omega_x, omega_y = np.meshgrid(omega_x, omega_y)

    # Calculate 
    r_omega_x = (omega_x / omega_x_0)  * np.exp(
        -0.5 * (omega_x / omega_x_0)**2 - 0.5 * (omega_y / omega_y_0)**2)
    r_omega_y = (omega_y / omega_y_0)
    
    # Calculate the 2D Gaussian derivative function
    r_omega = (omega_x / omega_x_0) * (omega_y / omega_y_0) * np.exp(
        -0.5 * (omega_x / omega_x_0)**2 - 0.5 * (omega_y / omega_y_0)**2)
    r_omega = np.abs(r_omega)
    d_selec = (1+d*np.sign(omega_x*omega_y))/2
    r_omega_b = d_selec*r_omega ** b
    return r_omega_b


# Define the function parameters
omega_x = np.linspace(-5, 5, 500)
omega_y = np.linspace(-35, 35, 500)
omega_x_0 = 2.5  # Reference frequency for x
omega_y_0 = 4.0  # Reference frequency for y
b = 1  # Exponent parameter
d = 0 # direction selectivity

r_omega_b = gaussian_derivative_RF(omega_x, omega_y, omega_x_0, omega_y_0, b, d)
# Plot the function
plt.figure(figsize=(8, 6))
plt.contourf(omega_x, omega_y, r_omega_b, levels=50, cmap='viridis')
plt.colorbar()
plt.xlabel(r'$\omega_x$', fontsize=14)
plt.ylabel(r'$\omega_y$', fontsize=14)
plt.title(r'2D Gaussian Derivative in Frequency Domain', fontsize=16)
# plt.grid(True)
plt.show()

# %%
