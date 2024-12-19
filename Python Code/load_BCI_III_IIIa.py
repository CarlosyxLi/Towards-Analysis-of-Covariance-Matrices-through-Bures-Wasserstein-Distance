import scipy.io
import os
import numpy as np

def import_data(subject): # k3b, k6b, l1b
    
    # training data
    data_dir = "D:/Projects/BCI Competition Data/BCI III IIIa/Training Data"
    data = scipy.io.loadmat(os.path.join(data_dir, subject+".mat"))
    data_class1 = data['class1_data']
    data_class2 = data['class2_data']
    data_class3 = data['class3_data']
    data_class4 = data['class4_data']
    data_train = np.concatenate((data_class1, data_class2, data_class3, data_class4), axis=2)
    data_train = np.transpose(data_train, (2,0,1)) # trial*channel*length
    # training labels
    label_train = np.repeat([1,2,3,4],[data_class1.shape[2],data_class2.shape[2],data_class3.shape[2],data_class4.shape[2]])
    # testing data
    data_dir = "D:/Projects/BCI Competition Data/BCI III IIIa/Testing Data"
    data = scipy.io.loadmat(os.path.join(data_dir, subject+"_test.mat"))
    data = data['test_data']
    data_test = np.transpose(data, (2,0,1)) # trial*channel*length
    # testing labels
    label_dir = "D:/Projects/BCI Competition Data/BCI III IIIa/Labels"
    label = scipy.io.loadmat(os.path.join(label_dir, subject+"_testlabel.mat"))
    label_test = label['label'][:,0]

    return data_train, data_test, label_train, label_test