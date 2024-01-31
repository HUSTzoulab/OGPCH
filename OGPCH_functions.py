import csv
import os
import glob 
import sklearn
import sklearn.manifold
import scipy
import scipy.io
import numpy as np
import pandas as pd 
import time
import matplotlib.pyplot as plt
from scGeneFit.functions import *
from sklearn.metrics import cohen_kappa_score
from sklearn.model_selection import train_test_split
import random
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
import seaborn as sns
import scipy.cluster.hierarchy as sch

def load_data_FLAT(dataset):
    if dataset == 'Mouse_E6_75' or dataset == 'PBMC3K' or dataset == 'CBMC8K' or dataset == 'zeisel':
        a = scipy.io.loadmat('data_files/'+dataset+'_data.mat')
        data= a['data']
    elif dataset == 'Mouse_E9_5':
        a = scipy.io.loadmat('data_files/'+dataset+'_data.mat')
        data= a['data'].toarray()
    else:
        print("currently available options are only 'Mouse_E6_75', 'CBMC8K' and 'PBMC3K'")
    
    N,d=data.shape      
    b = scipy.io.loadmat('data_files/'+dataset+'_names.mat')
    if dataset == 'Mouse_E6_75' or dataset == 'Mouse_E9_5' or dataset == 'zeisel':
        a = scipy.io.loadmat('data_files/'+dataset+'_labels2.mat')
        l_aux = a['labels2']
        names=[b['names'][i][1][0] for i in range(N)]
    else :
        a = scipy.io.loadmat('data_files/'+dataset+'_labels1.mat')
        l_aux = a['labels1']
        names=[b['names'][i][0][0] for i in range(N)]
    labels = np.array([i for [i] in l_aux])
    
    return [data, labels, names]
    
def load_data_HIE(dataset):
    if dataset == 'Mouse_E6_75' or dataset == 'zeisel':
        a = scipy.io.loadmat('data_files/'+dataset+'_data.mat')
        data= a['data']
    elif dataset == 'Mouse_E9_5':
        a = scipy.io.loadmat('data_files/'+dataset+'_data.mat')
        data= a['data'].toarray()
    else:
        print("currently available options are only  'Mouse_E6_75' and 'Mouse_E9_5'")
    
    N,d=data.shape
    a = scipy.io.loadmat('data_files/'+dataset+'_labels1.mat')
    l_aux = a['labels1']
    l_0=[l_aux[i][0] for i in range(l_aux.shape[0])]
    a = scipy.io.loadmat('data_files/'+dataset+'_labels2.mat')
    l_aux = a['labels2']
    l_1=[l_aux[i][0] for i in range(l_aux.shape[0])]
    labels=np.array([l_0, l_1])

    a = scipy.io.loadmat('data_files/'+dataset+'_names.mat')
    names0=[a['names'][i][0][0] for i in range(N)]
    names1=[a['names'][i][1][0] for i in range(N)]
    return [data, labels, [names0,names1]]

def get_scGeneFit_markers(dataset,method):
    if method == 'HIE':
        [data, labels, names] = load_data_HIE(dataset)
    elif method == 'FLAT':
        [data, labels, names] = load_data_FLAT(dataset)
    else:
        print('method error')
    N,d=data.shape

    methods='centers'
    redundancy=0.25

    markers=[]
    for i in range(20):
        num_markers=i*5+5
        if method == 'FLAT':
            markers.append(get_markers(data, labels, num_markers, method=methods, redundancy=redundancy))
        elif method == 'HIE':
            markers.append(get_markers_hierarchy(data, labels, num_markers, method=methods, redundancy=redundancy))
    pd.DataFrame(markers).to_csv('./markers/'+dataset+'/markers_'+dataset+'_'+method+'_scGeneFit.csv', index=None, header=None)
    
    return markers

def performance(X_train, y_train, X_test, y_test, clf):
    clf.fit(X_train, y_train)
    return clf.score(X_test, y_test)

def one_vs_all(dataset,num_bins=20):
    [data, labels, names] = load_data_FLAT(dataset)
    data_by_label = {}
    unique_labels = list(set(labels))
    number_classes = len(unique_labels)
    [N, d] = data.shape
    for lab in unique_labels: 
        X = [data[x, :] for x in range(len(labels)) if labels[x] == lab]
        data_by_label[lab] = X
    markers = [[] for i in range(number_classes)]
    for i in range(number_classes):
        markers[i] = [0 for j in range(i+1)]
        bins = data.max() / num_bins * range(num_bins + 1)
        for idx in range(len(markers[i])):
            c = unique_labels[idx]
            current_class = np.array(data_by_label[c])
            others = np.concatenate([data_by_label[lab]
                                     for lab in unique_labels if lab != c])
            big_dist = 0
            for gene in range(d):
                if gene not in markers[i][0:idx+1]:
                    [h1, b1] = np.histogram(current_class[:, gene], bins)
                    h1 = np.array(h1).reshape(1, -1) / current_class.shape[0]
                    [h2, b2] = np.histogram(others[:, gene], bins)
                    h2 = np.array(h2).reshape(1, -1) / others.shape[0]
                    dist = -sklearn.metrics.pairwise.additive_chi2_kernel(h1, h2)
                    if dist > big_dist:
                        markers[i][idx] = gene
                        big_dist = dist
                    
    pd.DataFrame(markers).to_csv('./markers/'+dataset+'/markers_'+dataset+'_ova.csv')
    return markers
    
def load_markers(dataset,method,way):
    if way == 'OGPCH' or way == 'scGeneFit':
        markers_dataframe = pd.read_csv('markers/'+dataset+'/markers_'+dataset+'_'+method+'_'+way+'.csv', header = None)
        markers_dataframe = np.array(markers_dataframe)
        markers = []
    
        markers = markers_dataframe.tolist()
        markers = [markers[i][0:i*5+5] for i in range(len(markers))]
        for marker in markers:
            for j in range(len(marker)):
                marker[j] = int(marker[j])
    elif way == 'ova':
        markers_dataframe = pd.read_csv('markers/'+dataset+'/markers_'+dataset+'_'+way+'.csv')
        markers_dataframe = np.array(markers_dataframe)
        markers = []
        m = markers_dataframe.tolist()
        for i in range(len(m)):   
            markers.append(m[i][1:i+2])
        for i in range(len(markers)):
            for j in range(len(markers[i])):
                markers[i][j] = int(markers[i][j])
    return markers

def get_accuracy(dataset,method,model,way):
    markers = load_markers(dataset,method,way)
    if method == 'FLAT':
        [data, labels, names] = load_data_FLAT(dataset)
        if model == 'kmeans':
            clf=NearestCentroid()
        elif model == 'knn':
            clf= KNeighborsClassifier()
        else:
            print('model error')
        accuracy = []
        for i in range(len(markers)):
            x_train, x_test, y_train, y_test = train_test_split(data[:,markers[i]], labels, random_state=1, train_size=0.8)
            accuracy.append(performance(x_train, y_train, x_test, y_test, clf))
        pd.DataFrame(accuracy).to_csv('./accuracy/'+dataset+'/accuracy_'+dataset+'_'+method+'_'+model+'_'+way+'.csv')
   
    elif method == 'HIE':
        [data, labels, names] = load_data_HIE(dataset)
        if model == 'kmeans':
            clf0=NearestCentroid()
            clf1=NearestCentroid()
        elif model == 'knn':
            clf0= KNeighborsClassifier()
            clf1= KNeighborsClassifier()
        else:
            print('model error')
        accuracy0, accuracy1 = [],[]
        for i in range(len(markers)):        
            x_train, x_test, y_train, y_test = train_test_split(data[:,markers[i]], labels[0], random_state=1, train_size=0.8)
            accuracy0.append(performance(x_train, y_train, x_test, y_test, clf0))
            x_train, x_test, y_train, y_test = train_test_split(data[:,markers[i]], labels[1], random_state=1, train_size=0.8)
            accuracy1.append(performance(x_train, y_train, x_test, y_test, clf1))
        pd.DataFrame(accuracy0).to_csv('./accuracy/'+dataset+'/accuracy_'+dataset+'_'+method+'_high_'+model+'_'+way+'.csv')
        pd.DataFrame(accuracy1).to_csv('./accuracy/'+dataset+'/accuracy_'+dataset+'_'+method+'_low_'+model+'_'+way+'.csv')
    else:
        print('method error')

def tuning_accuracy_mu(dataset,method,model):
    markers = pd.read_csv('./markers/markers_tuning/tuning_markers_'+dataset+'_'+method+'_OGPCH.csv', header = None)
    markers = np.array(markers).tolist()
    if method == 'FLAT':
        [data, labels, names] = load_data_FLAT(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
        else:
            print('model error')
            
        accuracy = []
        for i in range(len(markers)):
            x_train, x_test, y_train, y_test = train_test_split(data[:,markers[i]], labels, random_state=1, train_size=0.8)
            accuracy.append(performance(x_train, y_train, x_test, y_test, clf))
        
        pd.DataFrame(accuracy).to_csv('./accuracy/accuracy_tuning/accuracy_parameters_'+dataset+'_'+method+'_'+model+'.csv')
        
    elif method == 'HIE':
        [data, labels, names] = load_data_HIE(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
            clf1 = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
            clf1 = KNeighborsClassifier()
        else:
            print('model error')
        
        accuracy,accuracy1 = [], []
        for i in range(len(markers)):
            x_train, x_test, y_train, y_test = train_test_split(data[:,markers[i]], labels[0], random_state=1, train_size=0.8)
            accuracy.append(performance(x_train, y_train, x_test, y_test, clf))
            x_train1, x_test1, y_train1, y_test1 = train_test_split(data[:,markers[i]], labels[1], random_state=1, train_size=0.8)
            accuracy1.append(performance(x_train1, y_train1, x_test1, y_test1, clf1))

        pd.DataFrame(accuracy).to_csv('./accuracy/accuracy_tuning/accuracy_parameters_'+dataset+'_'+method+'_high_'+model+'.csv')
        pd.DataFrame(accuracy1).to_csv('./accuracy/accuracy_tuning/accuracy_parameters_'+dataset+'_'+method+'_low_'+model+'.csv')
    else:
        print('method error')    
    
def accuracy_line_chart(dataset, model):
    d1 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_FLAT_'+ model+'_OGPCH.csv')
    y1 = list(d1['0'])
    x = range(5,50,5)
    if dataset=='Mouse_E9_5':
        x = range(5,105,5)
        
    fig = plt.figure(dpi=600,figsize=(6,3.8))
    plt.xlabel("# of markers")
    plt.ylabel("Accuracy")
    
    if dataset=='CBMC8K' or dataset=='PBMC3K':
        d2 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_FLAT_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d4 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_FLAT_'+ model+'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        plt.plot(x, y1, '-o',color='lawngreen')+ plt.plot(x, y2, '-o',color='skyblue')+  plt.plot(x, y4, '-o',color='lightgray')#  
        plt.legend(['OGPCH_FLAT', 'scGeneFit', 'one vs all('+str(len(list(d4['0'])))+')'], loc =4)
        plt.title('Accuracy Tendency of '+dataset)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Accuracy_Compare.png', dpi = 600)
    if dataset=='zeisel' :
        d2 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_HIE_low_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d3 = pd.read_csv('./accuracy/'+dataset+'/accuracy_'+dataset+'_HIE_low_'+model+'_OGPCH.csv')
        y3 = list(d3['0'])
        d4 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_HIE_low_'+ model +'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        plt.plot(x, y3, '-o',color='red') +  plt.plot(x, y4, '-o',color='lightgray')+ plt.plot(x, y2, '-o',color='skyblue')#    
        plt.legend(['OGPCH_HIE', 'one vs all('+str(len(list(d4['0'])))+')', 'scGeneFit' ], loc =4) 
        plt.title('Accuracy Tendency of '+dataset)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Accuracy_Compare.png', dpi = 600)
    if dataset=='Mouse_E6_75' or dataset=='Mouse_E9_5':
        d2 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_HIE_low_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d3 = pd.read_csv('./accuracy/'+dataset+'/accuracy_'+dataset+'_HIE_low_'+model+'_OGPCH.csv')
        y3 = list(d3['0'])
        d4 = pd.read_csv('./accuracy/'+dataset+'/accuracy_' + dataset +'_HIE_low_'+ model +'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        if dataset=='Mouse_E9_5':
            y4 = [y4_ for i in range(20)]
        plt.plot(x, y3, '-o',color='red') +plt.plot(x, y1, '-o',color='lawngreen')+  plt.plot(x, y4, '-o',color='lightgray')+ plt.plot(x, y2, '-o',color='skyblue')#    
        plt.legend(['OGPCH_HIE', 'OGPCH_FLAT', 'one vs all('+str(len(list(d4['0'])))+')', 'scGeneFit' ], loc =4) 
        plt.title('Accuracy Tendency of '+dataset)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Accuracy_Compare.png', dpi = 600)

def kappa(dataset,method,model,way):
    if method == 'FLAT':
        markers = load_markers(dataset,method,way)
        [data, labels, names]= load_data_FLAT(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
        else:
            print('model error')
        pred = []
        for i in range(len(markers)):
            clf.fit(data[:,markers[i]], labels)    
            pred.append(clf.predict(data[:,markers[i]]))
        kappa = []
        for i in range(len(markers)):
            kappa.append(cohen_kappa_score(labels, pred[i]))
        pd.DataFrame(kappa).to_csv('./kappa/'+dataset+'/kappa_'+dataset+'_'+method+'_'+model+'_'+way+'.csv')
 
    elif method == 'HIE':
        markers = load_markers(dataset,method,way)
        [data, labels, names]= load_data_HIE(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
            clf1 = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
            clf1 = KNeighborsClassifier()
        else:
            print('model error')
        pred,pred1= [],[]
        for i in range(len(markers)):
            clf.fit(data[:,markers[i]], labels[0])
            clf1.fit(data[:,markers[i]], labels[1])    
            pred.append(clf.predict(data[:,markers[i]]))
            pred1.append(clf1.predict(data[:,markers[i]]))
        kappa,kappa1 = [],[]
        for i in range(len(markers)):
            kappa.append(cohen_kappa_score(labels[0], pred[i]))
            kappa1.append(cohen_kappa_score(labels[1], pred1[i]))
        pd.DataFrame(kappa).to_csv('./kappa/'+dataset+'/kappa_'+dataset+'_'+method+'_high_'+model+'_'+way+'.csv')
        pd.DataFrame(kappa1).to_csv('./kappa/'+dataset+'/kappa_'+dataset+'_'+method+'_low_'+model+'_'+way+'.csv')
    else:
        print('method error')       

def tuning_kappa_mu(dataset,method,model):
    markers = pd.read_csv('./markers/markers_tuning/tuning_markers_'+dataset+'_'+method+'_OGPCH.csv', header = None)
    markers = np.array(markers).tolist()
    if method == 'FLAT':
        [data, labels, names] = load_data_FLAT(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
        else:
            print('model error')
        pred = []
        for i in range(len(markers)):
            clf.fit(data[:,markers[i]], labels)    
            pred.append(clf.predict(data[:,markers[i]]))
        kappa = []
        for i in range(len(markers)):
            kappa.append(cohen_kappa_score(labels, pred[i]))
        
        pd.DataFrame(kappa).to_csv('./kappa/kappa_tuning/kappa_parameters_'+dataset+'_'+method+'_'+model+'.csv')
        
    elif method == 'HIE':
        [data, labels, names] = load_data_HIE(dataset)
        if model == 'kmeans':
            clf = NearestCentroid()
            clf1 = NearestCentroid()
        elif model == 'knn':
            clf = KNeighborsClassifier()
            clf1 = KNeighborsClassifier()
        else:
            print('model error')
        pred,pred1= [],[]
        for i in range(len(markers)):
            clf.fit(data[:,markers[i]], labels[0])
            clf1.fit(data[:,markers[i]], labels[1])    
            pred.append(clf.predict(data[:,markers[i]]))
            pred1.append(clf1.predict(data[:,markers[i]]))
        kappa,kappa1 = [],[]
        for i in range(len(markers)):
            kappa.append(cohen_kappa_score(labels[0], pred[i]))
            kappa1.append(cohen_kappa_score(labels[1], pred1[i]))
            
        pd.DataFrame(kappa).to_csv('./kappa/kappa_tuning/kappa_parameters_'+dataset+'_'+method+'_high_'+model+'.csv')    
        pd.DataFrame(kappa1).to_csv('./kappa/kappa_tuning/kappa_parameters_'+dataset+'_'+method+'_low_'+model+'.csv')
    else:
        print('method error')

def tuning_kappa_sam(dataset,method,model):

    if method == 'HIE':
        A = [6000,8000,10000,12000,14000]
    else:
        A = [3000,4000,5000,6000,7000]
    
    for a in A:
        markers = pd.read_csv('./markers/markers_sam/markers_'+str(a)+'_sam100_'+dataset+'_'+method+'_OGPCH.csv', header = None)
        markers = np.array(markers).tolist()
        if method == 'FLAT':
            [data, labels, names] = load_data_FLAT(dataset)
            if model == 'kmeans':
                clf = NearestCentroid()
            elif model == 'knn':
                clf = KNeighborsClassifier()
            else:
                print('model error')
            pred = []
            for i in range(len(markers)):
                clf.fit(data[:,markers[i]], labels)    
                pred.append(clf.predict(data[:,markers[i]]))
            kappa = []
            for i in range(len(markers)):
                kappa.append(cohen_kappa_score(labels, pred[i]))
            
            pd.DataFrame(kappa).to_csv('./kappa/kappa_tuning/kappa_'+str(a)+'_sampling_'+dataset+'_'+method+'_'+model+'.csv')
            
        elif method == 'HIE':
            [data, labels, names] = load_data_HIE(dataset)
            if model == 'kmeans':
                clf1 = NearestCentroid()
            elif model == 'knn':
                clf1 = KNeighborsClassifier()
            else:
                print('model error')
            pred,pred1= [],[]
            for i in range(len(markers)):
                clf1.fit(data[:,markers[i]], labels[1])    
                pred1.append(clf1.predict(data[:,markers[i]]))
            kappa,kappa1 = [],[]
            for i in range(len(markers)):
                kappa1.append(cohen_kappa_score(labels[1], pred1[i]))
                    
            pd.DataFrame(kappa1).to_csv('./kappa/kappa_tuning/kappa_'+str(a)+'_sampling_'+dataset+'_'+method+'_'+model+'.csv')
        else:
            print('method error')

def robust_voilin(dataset):
    if dataset=='zeisel' or dataset=='Mouse_E6_75' or dataset=='Mouse_E9_5':
        method = 'HIE'
        A = [6000,8000,10000,12000,14000]
    if dataset=='CBMC8K' or dataset=='PBMC3K':
        method = 'FLAT'
        A = [3000,4000,5000,6000,7000]
    b = []
    for i in A:
        a = pd.read_csv('kappa/kappa_tuning/kappa_'+str(i)+'_sampling_'+dataset+'_'+method+'_knn.csv')
        a = np.array(a).tolist()
        a = [_[1] for _ in a]
        b.append(a)
    b = np.array(b)
    b = b.flatten()
    x = [[A[j] for i in range(100)] for j in range(5)]
    x = np.array(x)
    x = x.flatten()
    
    
    plt.figure(figsize=(8,10))
    plt.tick_params(labelsize=20)
    plt.xlabel('# of constraints',fontsize=24)
    plt.ylabel('kappa',fontsize=24)
    plt.ylim(0,1)
    sns.violinplot(x=x, y=b, color="#ff7f0e")
    plt.savefig('pictures/'+dataset+'/'+ dataset+'_violin.png',dpi=600,bbox_inches = 'tight')

    
def kappa_line_chart(dataset, model):
    d1 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_FLAT_'+ model+'_OGPCH.csv')
    y1 = list(d1['0'])
    x = range(5,50,5)
    if dataset=='Mouse_E9_5':
        x = range(5,105,5)
    
    fig = plt.figure(dpi=600,figsize=(6,3.8))
    plt.xlabel("# of markers")
    plt.ylabel("Kappa")
    
    if dataset=='CBMC8K' or dataset=='PBMC3K':
        d2 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_FLAT_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d4 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_FLAT_'+ model +'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        plt.plot(x, y1, '-o',color='lawngreen')+ plt.plot(x, y2, '-o',color='skyblue')+  plt.plot(x, y4, '-o',color='lightgray')#
        plt.legend(['OGPCH_FLAT', 'scGeneFit', 'one vs all('+str(len(list(d4['0'])))+')'], loc =4)
        plt.title('Kappa Tendency of '+dataset)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Kappa_Compare.png', dpi = 600)
    
    if dataset=='zeisel' :
        d2 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_HIE_low_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d3 = pd.read_csv('./kappa/'+dataset+'/kappa_'+dataset+'_HIE_low_'+model+'_OGPCH.csv')
        y3 = list(d3['0'])
        d4 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_HIE_low_'+ model +'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        plt.plot(x, y3, '-o',color='red') +  plt.plot(x, y4, '-o',color='lightgray')+ plt.plot(x, y2, '-o',color='skyblue')#    
        plt.legend(['OGPCH_HIE', 'one vs all('+str(len(list(d4['0'])))+')', 'scGeneFit' ], loc =4)#
        plt.title('Kappa Tendency of '+dataset )
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Kappa_Compare.png', dpi = 600)

    if dataset=='Mouse_E6_75' or dataset=='Mouse_E9_5':
        d2 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_HIE_low_'+ model+'_scGeneFit.csv')
        y2 = list(d2['0'])
        d3 = pd.read_csv('./kappa/'+dataset+'/kappa_'+dataset+'_HIE_low_'+model+'_OGPCH.csv')
        y3 = list(d3['0'])
        d4 = pd.read_csv('./kappa/'+dataset+'/kappa_' + dataset +'_HIE_low_'+ model +'_ova.csv')
        y4_ = list(d4['0'])[-1]
        y4 = [y4_ for i in range(9)]
        if dataset=='Mouse_E9_5':
            y4 = [y4_ for i in range(20)]
        plt.plot(x, y3, '-o',color='red') +plt.plot(x, y1, '-o',color='lawngreen')+  plt.plot(x, y4, '-o',color='lightgray')+ plt.plot(x, y2, '-o',color='skyblue')# 
        plt.legend(['OGPCH_HIE', 'OGPCH_FLAT', 'one vs all('+str(len(list(d4['0'])))+')', 'scGeneFit' ], loc =4)#
        plt.title('Kappa Tendency of '+dataset )
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_'+ model +'_Kappa_Compare.png', dpi = 600)

def Compare_25_plot(dataset):
    if dataset == 'Mouse_E9_5':
        a = scipy.io.loadmat('data_files/'+dataset+'_data_sam.mat')
        data= a['data']
        N,d=data.shape
        b = scipy.io.loadmat('data_files/'+dataset+'_names_sam.mat')
        a = scipy.io.loadmat('data_files/'+dataset+'_labels2_sam.mat')
        l_aux = a['labels2']
        names=[b['names2'][i][0][0] for i in range(N)]
        labels = np.array([i for [i] in l_aux])

        markers0 = load_markers(dataset,'FLAT','scGeneFit')
        markers0 = markers0[4]
        markers1 = load_markers(dataset,'FLAT','ova')
        markers1 = markers1[-1]

        X_original  = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data)
        X_embedded0 = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data[:, markers0])
        X_embedded1 = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data[:, markers1])

        cmap = plt.cm.jet
        unique_names = list(set(names))
        num_labels = len(unique_names)
        colors = [cmap(int(i * 256 / num_labels)) for i in range(num_labels)]
        aux = [colors[unique_names.index(dataset)] for dataset in names]

        fig = plt.figure(dpi=600,figsize=(18,4))

        ax = fig.add_subplot(161)
        plt.ylabel('t-SNE2')
        for g in unique_names:
            i = [s for s in range(len(names)) if names[s] == g]
            ax.scatter(X_original[i, 0], X_original[i, 1],
                       c=[aux[i[0]]], s=5, label=names[i[0]])
        ax.set_title('Original data')

        ax1 = fig.add_subplot(162)
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax1.scatter(X_embedded1[i, 0], X_embedded1[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax1.set_title('one vs all ' +  str(len(markers1)) + ' markers')

        ax0 = fig.add_subplot(163)
        plt.xlabel('t-SNE1')
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax0.scatter(X_embedded0[i, 0], X_embedded0[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax0.set_title('scGeneFit 25 markers')

        markers2 = load_markers(dataset,'FLAT','OGPCH')
        markers2 = markers2[4]
        markers3_ = load_markers(dataset,'HIE','OGPCH')
        markers3 = markers3_[8]
        markers4 = markers3_[-1]

        X_embedded2 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers2])
        X_embedded3 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers3])
        X_embedded4 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers4])

        ax2 = fig.add_subplot(164)
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax2.scatter(X_embedded2[i, 0], X_embedded2[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax2.set_title('OGPCH FLAT 25 markers')

        ax3 = fig.add_subplot(165)
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax3.scatter(X_embedded3[i, 0], X_embedded3[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax3.set_title('OGPCH HIE 45 markers')

        ax4 = fig.add_subplot(166)
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax4.scatter(X_embedded4[i, 0], X_embedded4[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax4.set_title('OGPCH HIE 100 markers')

        plt.tight_layout()
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(right=0.7)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_compare_25_sam.png', dpi = 600)
    
    else:
        if dataset == 'PBMC3K' or dataset == 'CBMC8K':
            [data, labels, names]= load_data_FLAT(dataset)
            N,d=data.shape
        else:
            [data, labels, names]= load_data_HIE(dataset)
            if dataset == 'Mouse_E6_75':
                names = names[1]
            else:
                names = names[0]
            N,d=data.shape
        
        markers0 = load_markers(dataset,'FLAT','scGeneFit')
        markers0 = markers0[4]
        markers1 = load_markers(dataset,'FLAT','ova')
        markers1 = markers1[-1]
        
        X_original  = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data)
        X_embedded0 = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data[:, markers0])
        X_embedded1 = sklearn.manifold.TSNE(n_components=2, perplexity=45).fit_transform(data[:, markers1])
     
        cmap = plt.cm.jet
        unique_names = list(set(names))
        num_labels = len(unique_names)
        colors = [cmap(int(i * 256 / num_labels)) for i in range(num_labels)]
        aux = [colors[unique_names.index(dataset)] for dataset in names]

        fig = plt.figure(dpi=600,figsize=(18,4))
        
        ax = fig.add_subplot(161)
        plt.ylabel('t-SNE2')
        for g in unique_names:
            i = [s for s in range(len(names)) if names[s] == g]
            ax.scatter(X_original[i, 0], X_original[i, 1],
                       c=[aux[i[0]]], s=5, label=names[i[0]])
        ax.set_title('Original data')
        
        ax1 = fig.add_subplot(162)
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax1.scatter(X_embedded1[i, 0], X_embedded1[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax1.set_title('one vs all ' +  str(len(markers1)) + ' markers')
        
        ax0 = fig.add_subplot(163)
        plt.xlabel('t-SNE1')
        for g in np.unique(names):
            i = [s for s in range(len(names)) if names[s] == g]
            ax0.scatter(X_embedded0[i, 0], X_embedded0[i, 1],
                        c=[aux[i[0]]], s=5, label=names[i[0]])
        ax0.set_title('scGeneFit 25 markers')
        
        if dataset == 'PBMC3K' or dataset == 'CBMC8K':
            markers2_ = load_markers(dataset,'FLAT','OGPCH')
            markers2 = markers2_[4]
            X_embedded2 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers2])
            ax2 = fig.add_subplot(164)
            for g in np.unique(names):
                i = [s for s in range(len(names)) if names[s] == g]
                ax2.scatter(X_embedded2[i, 0], X_embedded2[i, 1],
                            c=[aux[i[0]]], s=5, label=names[i[0]])
            ax2.set_title('OGPCH 25 markers')
            
            if dataset == 'PBMC3K':
                markers3 = markers2_[8]
                X_embedded3 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers3])
                ax3 = fig.add_subplot(165)
                for g in np.unique(names):
                    i = [s for s in range(len(names)) if names[s] == g]
                    ax3.scatter(X_embedded3[i, 0], X_embedded3[i, 1],
                                c=[aux[i[0]]], s=5, label=names[i[0]])
                ax3.set_title('OGPCH 45 markers')
        
        elif dataset == 'Mouse_E6_75' or dataset == 'zeisel':
            markers2 = load_markers(dataset,'FLAT','OGPCH')
            markers2 = markers2[4]
            markers3_ = load_markers(dataset,'HIE','OGPCH')
            markers3 = markers3_[4]
            markers4 = markers3_[8]
            
            
            X_embedded2 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers2])
            X_embedded3 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers3])
            X_embedded4 = sklearn.manifold.TSNE(n_components=2, perplexity=40).fit_transform(data[:, markers4])
            
            ax2 = fig.add_subplot(164)
            for g in np.unique(names):
                i = [s for s in range(len(names)) if names[s] == g]
                ax2.scatter(X_embedded2[i, 0], X_embedded2[i, 1],
                            c=[aux[i[0]]], s=5, label=names[i[0]])
            ax2.set_title('OGPCH FLAT 25 markers')
            
            if dataset == 'Mouse_E6_75':
                ax3 = fig.add_subplot(165)
                for g in np.unique(names):
                    i = [s for s in range(len(names)) if names[s] == g]
                    ax3.scatter(X_embedded3[i, 0], X_embedded3[i, 1],
                                c=[aux[i[0]]], s=5, label=names[i[0]])
                ax3.set_title('OGPCH HIE 25 markers')
            
            if dataset == 'zeisel':
                ax4 = fig.add_subplot(165)
                for g in np.unique(names):
                    i = [s for s in range(len(names)) if names[s] == g]
                    ax4.scatter(X_embedded4[i, 0], X_embedded4[i, 1],
                                c=[aux[i[0]]], s=5, label=names[i[0]])
                ax4.set_title('OGPCH HIE 45 markers')

        else: 
            print("currently available options are only 'Mouse_E6_75', 'Mouse_E9_5' and 'PBMC3K'")

        plt.tight_layout()
        plt.legend(bbox_to_anchor=(1, 1))
        plt.subplots_adjust(right=0.7)
        plt.savefig('pictures/'+dataset+'/'+ dataset + '_compare_25.png', dpi = 600)
    return fig

def tree_plot(dataset):
    cell_name = pd.read_csv('cell_type/'+dataset+'_cell_type.csv',header=None)
    cell_name = np.array(cell_name).tolist()
    cell_name = [_[0] for _ in cell_name]
    
    OGPCH_HIE = pd.read_csv('markers/'+dataset+'/markers_'+dataset+'_HIE_OGPCH.csv',header=None)
    OGPCH_HIE = np.array(OGPCH_HIE).tolist()
    if dataset == 'Mouse_E9_5':
        m1 = OGPCH_HIE[-1]
    if dataset == 'Mouse_E6_75':
        m1 = OGPCH_HIE[4][:25]
    else:
        m1 = OGPCH_HIE[8][:45]
    
    a = pd.read_csv('data_files/'+dataset+'_mean.csv',header=None)
    a = np.array(a)
    a = a.T
    A = a[m1]
    A = A.T

    linkage_matrix = sch.linkage(A, method='ward')
    fig, ax = plt.subplots(figsize=(14,18))
    plt.xlabel('Distance')
    plt.title('Dendrogram')
    plt.xticks(size=14)
    plt.yticks(size=14)
    if dataset == 'Mouse_E9_5':
        cutoff = linkage_matrix[-3][2]
        dendrogram = sch.dendrogram(linkage_matrix, color_threshold=cutoff, orientation='left', labels=cell_name)
        for i in range(12):
            ax.get_yticklabels()[i].set_color("darkorange")
        for i in range(12,44):
            ax.get_yticklabels()[i].set_color("purple")

        for i in range(20,25):
            ax.get_yticklabels()[i].set_color("firebrick")
        for i in range(18,20):
            ax.get_yticklabels()[i].set_color("royalblue")
        ax.get_yticklabels()[25].set_color("darkorange")
        ax.get_yticklabels()[31].set_color("darkorange")
        ax.get_yticklabels()[32].set_color("firebrick")
        ax.get_yticklabels()[42].set_color("firebrick")
        ax.get_yticklabels()[14].set_color("slategray")
        ax.get_yticklabels()[16].set_color("darkblue")
    if dataset == 'Mouse_E6_75':
        cutoff = linkage_matrix[-2][2]
        dendrogram = sch.dendrogram(linkage_matrix, color_threshold=cutoff, orientation='left', labels=cell_name)
        for i in range(3):
            ax.get_yticklabels()[i].set_color("darkorange")
        for i in range(3,4):
            ax.get_yticklabels()[i].set_color("royalblue")
        for i in range(4,6):
            ax.get_yticklabels()[i].set_color("green")   
    if dataset == 'zeisel':        
        cutoff = linkage_matrix[-5][2]
        dendrogram = sch.dendrogram(linkage_matrix, color_threshold=cutoff, orientation='left', labels=cell_name) # 
        for i in range(6):
            ax.get_yticklabels()[i].set_color("darkorange")
        for i in range(6,10):
            ax.get_yticklabels()[i].set_color("forestgreen")
        for i in range(10,14):
            ax.get_yticklabels()[i].set_color("mediumspringgreen")
        for i in range(14,18):
            ax.get_yticklabels()[i].set_color("red")
        for i in range(18,33):
            ax.get_yticklabels()[i].set_color("purple")
        ax.get_yticklabels()[34].set_color("purple")
        for i in range(35,40):
            ax.get_yticklabels()[i].set_color("peru")
        for i in range(40,48):
            ax.get_yticklabels()[i].set_color("deeppink")
    
    plt.savefig('pictures/'+dataset+'/'+dataset+'_tree.png', dpi = 600,bbox_inches='tight')
    plt.show()

