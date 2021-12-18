import h5py
import scipy
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.metrics import accuracy_score
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score as ari
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.io import loadmat




def load_gsva(arg):
    from sklearn.preprocessing import LabelEncoder
    x = pd.read_csv(arg+"/Pathway_score.csv")
    y = pd.read_csv(arg+"/Pathway_metadata.csv")
    # print(x.shape)
    # print(y.shape)
    x = x.iloc[1:,1:]
    y = y.iloc[:,0]
    labelencoder = LabelEncoder()
    print(y.unique())
    y_num = labelencoder.fit_transform(y)
    # print(y_num)
    x = x.to_numpy()
    x = np.transpose(x)
    return x, y_num
 





# clustering datasets of sizes having samples in the range
# of thousands and incorporate the K-means clustering cost into the 
# deep dictionary learning (DDL) framework.
seed = 78638924
np.random.seed(seed)

def findclusterlabels(x, y, n_clusters,out_dir):
    x = np.transpose(x)

    # Initialization of parameters
    k1=512
    k2=256
    k3=128
    num_input_features = x.shape[0]
    num_samples = x.shape[1]
    n = 11 #number of iterations
    mu = 0.5

    d1 = np.random.randn(num_input_features,k1)
    d2 = np.random.randn(k1,k2)
    d3 = np.random.randn(k2,k3)
    z = np.random.randn(k3,num_samples)

    # H initialization
    xlabels = KMeans(n_clusters = args.n_clusters, init="k-means++").fit_predict(np.transpose(x))
    print("Length of xlabels ",xlabels.shape)


    onehotencoder = preprocessing.OneHotEncoder()
    h = onehotencoder.fit_transform(np.transpose(xlabels).reshape(-1,1)).toarray()
    print("H shape ",h.shape)
    hh = h.sum(axis=1)
    print(hh)
    print("Sum over cols of h ",hh.sum(axis=0))
    h = np.transpose(h)
    print("KMeans Clustering Result")
    print('nmi = %.4f, ari = %.4f' % (nmi(y, xlabels), ari(y, xlabels)))
    print("Confusion Matrix  = ",confusion_matrix(y,xlabels))
    cost_list=[]
    cost_list1=[]
    cost_list2=[]
    d1_list = []
    d2_list = []
    d3_list = []
    z_list = []
    
    z[z<0]=0

    bestnmi = 0
    bestari = 0
    for i in range(n):
        print("Calculating d1 d2 d3")
        #Solving d1,d2,d3 via pseudo-inverse
        P = np.dot(d3,z)
        Q = np.dot(d2,P)

        P[P<0]=0
        Q[Q<0]=0
        
        #Saving old d1,d2,d3 values
        d1_pre = d1
        d2_pre = d2
        d3_pre = d3
        z_pre = z

        d1 = np.dot(x,(np.linalg.pinv(Q)))  #d1 shape num_input_features*k1

        d2 = np.dot((np.linalg.pinv(d1)),np.dot(x,np.linalg.pinv(P))) # d2 shape k1*k2

        W = np.linalg.pinv(np.dot(d1,d2))
        d3 = np.dot(W, np.dot(x, np.linalg.pinv(z))) #d3 shape k2*k3

        #Calculating error of d1,d2,d3
        d1_error = np.linalg.norm((d1 - d1_pre), 'fro')
        d2_error = np.linalg.norm((d2 - d2_pre), 'fro')
        d3_error = np.linalg.norm((d3 - d3_pre), 'fro')

        print("Solving for z...")
        #Solving z via sylvester equation of form ax+xb=c
        hterm = np.dot(np.transpose(h),np.dot(np.linalg.inv(np.dot(h,np.transpose(h))),h))
        #print("hterm shape",hterm.shape)
        #print("Calculating a")
        a = np.dot((np.transpose(np.dot(d1,np.dot(d2,d3)))),(np.dot(d1,np.dot(d2,d3))))
        #print("Calculating b")
        b = mu*(np.dot(np.transpose(np.identity(hterm.shape[0])-np.transpose(hterm)),(np.identity(hterm.shape[0])-hterm)))
        #print("Calculating c")
        c = np.dot((np.transpose(np.dot(d1,np.dot(d2,d3)))),x)
        print("A ",a.shape)
        print("B ",b.shape)
        print("C shape", c.shape)
        print("Calculating z")
        z = scipy.linalg.solve_sylvester(a, b, c)

        z[z<0]=0

        z_error = np.linalg.norm((z - z_pre), 'fro')

        # print("Calculating h by kmeans")
        z_pred = KMeans(n_clusters = args.n_clusters, init="k-means++").fit_predict(np.transpose(z))

        if(i==9):
            filename = str(out_dir)+'/DDLKlabels_nclusters_'+str(n_clusters)+'.csv'
            np.savetxt(filename, z_pred, delimiter=',', fmt='%d')
        nmi_cur = nmi(y, z_pred)
        ari_cur = ari(y, z_pred)
        print('nmi = %.4f, ari = %.4f' % (nmi_cur, ari_cur))
        print("Confusion Matrix ",confusion_matrix(y,z_pred))

        if bestnmi < nmi_cur:
            bestnmi = nmi_cur
        if bestari < ari_cur:
            bestari = ari_cur

        h_pre = h
        h = onehotencoder.fit_transform(np.transpose(z_pred).reshape(-1,1)).toarray()

        hh = h.sum(axis=1)
        print("Sum over cols of h ",hh.sum(axis=0))

        hhh = h-np.transpose(h_pre)

        h = np.transpose(h)

        print("Calculating cost")
        #Cost Calculation
        z_cost = z
        P = np.dot(d3,z_cost)
        P[P<0]=0
        Q = np.dot(d2,P)
        Q[Q<0]=0
        cost1 = (x - (np.dot(d1,Q)))

        result1 = np.sum(cost1**2)
        cost1f = np.linalg.norm(cost1, 'fro')
        cost1f = np.square(cost1f)

        hterm = np.dot(np.transpose(h),np.dot(np.linalg.inv(np.dot(h,np.transpose(h))),h))
        cost2 = (z - (np.dot(z,hterm)))
        cost2f = np.linalg.norm(cost2, 'fro')
        cost2f = np.square(cost2f)

        costf = cost1f + (mu*cost2f)

        if (i>0):
            cost_list.append(costf)
            cost_list1.append(cost1f)
            cost_list2.append(cost2f)
            d1_list.append(d1_error)
            d2_list.append(d2_error)
            d3_list.append(d3_error)
            z_list.append(z_error)

        print("Result 1 : ", result1)
        print("Cost 1 : ", cost1f)
        print("Cost 2 : ", cost2f)
        print("Final Cost ", costf)
        print("End of %d iteration",i)
        print("==================================")
    print("cost list 1 : ",cost_list1)
    print("cost list 2 : ",cost_list2)
    print("final cost list : ",cost_list)
    print("d1 list : ",d1_list)
    print("d2 list : ",d2_list)
    print("d3 list : ",d3_list)
    print("Z list : ",z_list)
    print('bestnmi = %.4f',bestnmi, 'bestari = %.4f', bestari)

    return 

if __name__ == "__main__":
    # setting the hyper parameters
    import argparse
    parser = argparse.ArgumentParser(description='train')
    parser.add_argument('--n_clusters', default=10, type=int)
    parser.add_argument('--out_dir', default="", type=str)
    args = parser.parse_args()
    print(args)

    n_clusters = args.n_clusters
    out_dir = args.out_dir

    # load dataset
    x, y = load_gsva(out_dir)
    print("X shape ",x.shape)
    print("Y shape ",y.shape) 

    import time
    start_time = time.time()

    findclusterlabels(x, y, n_clusters,out_dir)

    print("--- %s seconds ---" % (time.time() - start_time))
    
