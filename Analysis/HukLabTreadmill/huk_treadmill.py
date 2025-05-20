
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

# read matlab files (hd5 format)
def load_file(fname, subtract_t0=True):
    D = {}
    t_start = np.inf

    with h5.File(fname, 'r') as f:

        if 'D' in f.keys():
            for key in f['D'].keys():

                D[key] = f['D'][key][:].T
                # print(key, D[key].shape)
                # find experiment start time. If a key has "time" in it, check if it is the minimum time
                if 'Time' in key:
                    t_start = np.minimum(t_start, D[key][0,0])
        else:
            for key in f.keys():
                if key == '#refs#':
                    continue
                
                try:
                    D[key] = f[key][:].T
                    # print(key, D[key].shape)
                    # find experiment start time. If a key has "time" in it, check if it is the minimum time
                    if 'Time' in key:
                        t_start = np.minimum(t_start, np.nanmin(D[key][:]))

                except:
                    print('Could not load key', key)
                    continue
                
    if subtract_t0:
        # loop over keys with "Time" in the name and subtract t_start
        for key in D.keys():
            if 'Time' in key:
                try:
                    D[key] -= t_start
                except:
                    print('Could not subtract t_start from', key)

    print("Eye time sampled at: %d" %(1/np.median(np.diff(D['eyeTime'].flatten()))))
    return D


def get_run_epochs(D, plot=False, medfilt_window=551, speed_threshold=2):
    '''
    Get run epochs from treadmill data
    '''
    from scipy.interpolate import interp1d
    from scipy.signal import savgol_filter, medfilt

    x_orig = D['treadSpeed'].flatten()
    t_orig = D['treadTime'].flatten()
    
    # interpolate nans in t_orig
    iix = ~np.isnan(t_orig)
    t_orig = interp1d(np.arange(len(t_orig))[iix], t_orig[iix], kind='linear', bounds_error=False)(np.arange(len(t_orig)))

    # x_orig_diff = np.diff(x_orig, prepend=0)
    t_orig_diff = np.diff(t_orig, prepend=0)
    dt = np.median(t_orig_diff)

    # dx = x_orig #x_orig_diff 
    # breaks = np.where(np.isnan(dx))[0]
    # breaks = np.where(dx < -10)[0]

    # x_new = x_orig.copy()

    # for b in breaks:
    #     if np.isclose(x_orig[b], 0, atol=.1):
    #         x_new[b:] += x_orig[b-1]

    # calculate the derivative
    # dx_new = savgol_filter(x_new.flatten(), 251, 1, deriv=1, delta=dt)
    # median filter speed

    dx_new = x_orig
    medspeed = medfilt(np.abs(dx_new), medfilt_window)

    # find run epochs
    run_epochs = (medspeed>speed_threshold).astype(np.float64())                       

    run_starts = np.where(np.diff(run_epochs)==1)[0]
    run_ends = np.where(np.diff(run_epochs)==-1)[0]

    if run_starts[0] > run_ends[0]:
        # concatenate 0 to the fron of run_starts
        run_starts = np.concatenate(([0], run_starts))
    if run_ends[-1] < run_starts[-1]:
        # concatenate last index to the end of run_ends
        run_ends = np.concatenate((run_ends, [len(run_epochs)]))

    # convert to time 
    t_run_starts = t_orig[run_starts]
    t_run_ends = t_orig[run_ends]

    good_ix = ~np.logical_or(np.isnan(t_run_starts), np.isnan(t_run_ends))
    t_run_starts = t_run_starts[good_ix]
    t_run_ends = t_run_ends[good_ix]

    # loop over and combine epochs that are too close together
    refractory = 5 # seconds
    i = 0
    while i < len(t_run_starts)-1:
        if t_run_starts[i+1] - t_run_ends[i] < refractory:
            t_run_ends = np.delete(t_run_ends, i)
            t_run_starts = np.delete(t_run_starts, i+1)
        else:
            i += 1

    # plot
    if plot:
        plt.figure()
        plt.plot(t_orig, x_orig, '-')
        plt.plot(t_orig, dx_new, '-')
        plt.plot(t_orig, medspeed, '-')
        plt.plot(t_orig, 5*run_epochs, '-')
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (cm/s)')

        for i in range(len(t_run_starts)):
            # plot vertical fill between run epochs
            plt.fill_betweenx([0,np.max(x_orig)], t_run_starts[i], t_run_ends[i], color='lightblue', alpha=0.5)
    
    return t_run_starts, t_run_ends


def get_good_eye_epochs(D, plot=False):
    '''
    Get good eye epochs from eye data
    '''
    from scipy.interpolate import interp1d
    # plt.figure()
    x_orig = interp1d(D['eyeTime'].flatten(), D['eyePos'][:,0].flatten(), kind='linear')(D['eyeTime'].flatten())
    y_orig = interp1d(D['eyeTime'].flatten(), D['eyePos'][:,1].flatten(), kind='linear')(D['eyeTime'].flatten())
    x_orig -= np.nanmedian(x_orig)
    y_orig -= np.nanmedian(y_orig)

    dx = np.diff(x_orig.flatten(), prepend=0)
    dt = np.median(np.diff(D['eyeTime'].flatten(), prepend=0))
    # plt.plot(D['eyeTime'].flatten(),x_orig)
    # plt.plot(D['eyeTime'].flatten(), y_orig)
    bad_samples = np.logical_or(np.abs(x_orig) > 10, np.abs(y_orig) > 8)
    bad_samples = np.logical_or(bad_samples, np.isnan(x_orig)).astype(np.float64())
    
    # find bad sample starts and stops
    bad_starts = np.where(np.diff(bad_samples)==1)[0]
    bad_ends = np.where(np.diff(bad_samples)==-1)[0]

    if bad_starts[0] > bad_ends[0]:
        # concatenate 0 to the fron of run_starts
        bad_starts = np.concatenate(([0], bad_starts))
    if bad_ends[-1] < bad_starts[-1]:
        # concatenate last index to the end of run_ends
        bad_ends = np.concatenate((bad_ends, [len(bad_samples)-1]))
                                
    # offset the starts and stops by a 5 samples
    bad_starts -= 50
    bad_ends += 50

    # make sure none are less than zero or greater than the length of the signal
    bad_starts = np.maximum(0, bad_starts)
    bad_ends = np.minimum(len(bad_samples), bad_ends)

    # combine bad epochs that are too close together
    refractory = int(.5/dt) # samples ()
    i = 0
    while i < len(bad_starts)-1:
        if bad_starts[i+1] - bad_ends[i] < refractory:
            bad_ends = np.delete(bad_ends, i)
            bad_starts = np.delete(bad_starts, i+1)
        else:
            i += 1

    # plot
    if plot:
        plt.figure()
        plt.plot(D['eyeTime'].flatten(), x_orig)
        plt.plot(D['eyeTime'].flatten(), y_orig)


    # replace x_orig and y_orig with nans for bad epochs
    for i in range(len(bad_starts)):
        x_orig[bad_starts[i]:bad_ends[i]] = np.nan
        y_orig[bad_starts[i]:bad_ends[i]] = np.nan

    # find good epochs by finding non-nan values
    good_epochs = ~np.isnan(x_orig)
    # convert to float64
    good_epochs = good_epochs.astype(np.float64)
    good_starts = np.where(np.diff(good_epochs)==1)[0]
    good_ends = np.where(np.diff(good_epochs)==-1)[0]

    if good_starts[0] > good_ends[0]:
        # concatenate 0 to the fron of run_starts
        good_starts = np.concatenate(([0], good_starts))
    if good_ends[-1] < good_starts[-1]:
        # concatenate last index to the end of run_ends
        good_ends = np.concatenate((good_ends, [len(good_epochs)-1]))

    # plot
    if plot:
        for i in range(len(good_starts)):
            # plot vertical fill between run epochs
            plt.fill_betweenx([-20,20], D['eyeTime'].flatten()[good_starts[i]], D['eyeTime'].flatten()[good_ends[i]], color='lightblue', alpha=0.5)

    # return good epoch starts and stops as times
    t_good_starts = D['eyeTime'].flatten()[good_starts]
    t_good_ends = D['eyeTime'].flatten()[good_ends]

    return t_good_starts, t_good_ends


def get_epochs(t_run_starts, t_run_ends, t_good_starts, t_good_ends, union=True):
    '''
    Find the instersection of run and eye epochs
    union=True -> returns run epochs where eye is valid
    union=False -> returns stationary epochs where eye is valid
    '''
    t_0 = np.minimum(t_run_starts[0], t_good_starts[0])
    t_1 = np.maximum(t_run_ends[-1], t_good_ends[-1])
    t = np.arange(t_0, t_1, 1e-3)
    run_epochs = np.zeros_like(t)
    good_epochs = np.zeros_like(t)
    for i in range(len(t_run_starts)):
        run_epochs += (t >= t_run_starts[i]) & (t <= t_run_ends[i])
    
    for i in range(len(t_good_starts)):
        good_epochs += (t >= t_good_starts[i]) & (t <= t_good_ends[i])
    
    # if valid, return only epochs where both run and eye are valid
    if union:
        epochs = run_epochs * good_epochs
    else:
        epochs = (1-run_epochs) * good_epochs

    # find starts and stops
    starts = np.where(np.diff(epochs)==1)[0]
    ends = np.where(np.diff(epochs)==-1)[0]

    if starts[0] > ends[0]:
        # concatenate 0 to the fron of run_starts
        starts = np.concatenate(([0], starts))
    if ends[-1] < starts[-1]:
        # concatenate last index to the end of run_ends
        ends = np.concatenate((ends, [len(epochs)-1]))

    # convert to time
    t_starts = t[starts]
    t_ends = t[ends]

    return t_starts, t_ends


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
    term1 = np.exp(-2j * np.pi * k_values[None,None,:] * u_t[:,None,None])
    term2 = np.exp(-2j * np.pi * omega_values[None,:,None] * t[:,None,None])
    modulated_signal = term1 * term2
    
    # Integrate over time (Fourier transform)
    integrated_value = np.trapz(modulated_signal, t[:,None,None], axis=0)
            
    
    Q_k_omega = np.abs(integrated_value) ** 2
    
    return Q_k_omega



def get_eye_displacements(D, t_starts, t_ends, 
        Twin = 512,
        noverlap=None):
    
    '''
    get segments of eye position data in delta x, delta t

    '''
    # D is the data dictionary
    # t_starts, t_ends are the start and end times of the epochs
    # Twin (integer) in samples
    # noverlap = number of samples to overlap

    xpos = D['eyePos'][:,0].flatten()
    xsegments = []

    if noverlap is None:
        noverlap = Twin // 2

    # loop over epochs
    for i in range(len(t_starts)):
        t0 = t_starts[i]
        t1 = t_ends[i]
        iix = (D['eyeTime'].flatten() >= t0) & (D['eyeTime'].flatten() <= t1)

        xsegment = xpos[iix]
        Ns = len(xsegment)

        step = Twin - noverlap
        
        # loop over segments stepping by noverlap
        for j in range(0, Ns-Twin, step):
            xseg = xsegment[j:j+Twin]

            # subtract first sample to calculate displacement
            xseg = xseg - xseg[0]

            xsegments.append(xseg)

    return np.stack(xsegments)


# simulate neural tuning

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