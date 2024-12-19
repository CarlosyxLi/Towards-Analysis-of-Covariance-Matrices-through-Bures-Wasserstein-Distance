import scipy.io
import os
import numpy as np

def import_data(file_name):
    
    if "T" in file_name: # e.g., B0101T
        # training data
        data_dir = "D:/Projects/BCI Competition Data/BCI IV IIb/Training Data"
        data = scipy.io.loadmat(os.path.join(data_dir, file_name+".mat"))
        data_class1 = data['class1_data']
        data_class2 = data['class2_data']
        data_output = np.concatenate((data_class1, data_class2), axis=2)
        data_output = np.transpose(data_output, (2,0,1)) # trial*channel*length
        #labels
        label_output = np.repeat([1,2], [data_class1.shape[2],data_class2.shape[2]])
        
    elif "E" in file_name: # e.g., B0104E
        # testing data
        data_dir = "D:/Projects/BCI Competition Data/BCI IV IIb/Testing Data"
        data = scipy.io.loadmat(os.path.join(data_dir, file_name+".mat"))
        data = data['test_data']
        data_output = np.transpose(data, (2,0,1)) # trial*channel*length
        # labels
        label_dir = "D:/Projects/BCI Competition Data/BCI IV IIb/Labels"
        label = scipy.io.loadmat(os.path.join(label_dir, file_name+".mat"))
        label_output = label['classlabel'][:,0]
        
    return data_output, label_output