import numpy as np
from scipy.signal import butter, filtfilt

# Butterworth bandpass filter
def filter_all(x, order, lowcut, highcut, fs):
    # input [trial, channel, length]
    # output [trial, channel, length]
    b,a = butter(N=order, Wn=[lowcut,highcut], fs=fs, btype='bandpass')
    output = np.zeros(x.shape, dtype = np.float32)
    for trial_idx in range(x.shape[0]):
        for channel_idx in range(x.shape[1]):
            output[trial_idx, channel_idx] = filtfilt(b, a, x[trial_idx, channel_idx])
    return output