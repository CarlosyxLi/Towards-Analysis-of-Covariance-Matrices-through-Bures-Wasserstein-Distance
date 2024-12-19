import scipy.io
import os
import numpy as np

def import_data(subject, start_time=0.5, end_time=3): 
    
    ##### How to use #####
    # subject - aa,al,av,aw,ay
    # start_time=0.5 - cut segment start at 0.5s after the cue
    # end_time=3 - cut segment end at 3s after the cue
    
    # E.g. If subject="al", start_time=1, end_time=2, then output has the following shape
    # data_train [224,118,1000], data_test [56,118,1000], label_train [224], label_test [56]
    # It is advised to have end_time < 5, since the shortest event is about 5.6s long. 
    
    data_dir = "D:/Projects/BCI Competition Data/BCI III IVa"
    data = scipy.io.loadmat(os.path.join(data_dir, "data_set_IVa_"+subject+".mat"))
    cnt = data['cnt'] # uncut signals
    mrk = data['mrk']
    pos = mrk[0][0][0][0] # positions of the cue in the EEG signals
    y = mrk[0][0][1][0] # labels of both train and test set
    
    pos_train = pos[~np.isnan(y)]
    pos_test = pos[np.isnan(y)]
    N_train = len(pos_train) # size of train set
    N_test = len(pos_test) # size of test set
    
    p = cnt.shape[1] # number of channels
    
    sampling_rate = 1000 # sampling rate in Hz
    start_pos = int(sampling_rate*start_time) # start position of cut segment
    end_pos = int(sampling_rate*end_time) # end position of cut segment
    l = end_pos-start_pos # signal length
    
    # initiate an [N,p,l] array to restore cut signals
    data_train = np.zeros((N_train, p, l), dtype = np.float32)
    data_test = np.zeros((N_test, p, l), dtype = np.float32)
    
    # cut train set
    for i in range(N_train):
        for j in range(p):
            data_train[i,j] = cnt[(pos_train[i]+start_pos):(pos_train[i]+start_pos+l),j]
    
    # cut test set
    for i in range(N_test):
        for j in range(p):
            data_test[i,j] = cnt[(pos_test[i]+start_pos):(pos_test[i]+start_pos+l),j]
    
    label = scipy.io.loadmat(os.path.join(data_dir, "true_labels_"+subject+".mat"))
    label = label['true_y'][0]
    
    label_train = label[~np.isnan(y)]
    label_test = label[np.isnan(y)]
    
    return data_train, data_test, label_train, label_test