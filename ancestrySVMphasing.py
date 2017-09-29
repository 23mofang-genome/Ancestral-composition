#!/usr/bin/python
# -*- coding:utf-8 -*-


import sys, os
import numpy as np
import commands
from sklearn.svm import SVC


def strclacstr(arr0, arr1):
    result=[]
    for k0 in range(arr0.shape[0]):
        if arr0[k0]==arr1[k0]:
            result.append(1)
        else:
            result.append(0)
    return calcweight(result)


def calcweight(calcarray):
    weight=0
    substr=0
    flag=True
    for i in calcarray:
        if i:
            substr=substr+1
        else:
            if substr:
                weight=weight+substr*(substr+1)/2
                if flag:
                    flag=False
            substr=0
    if flag:
        weight=weight+substr*(substr+1)/2
    return weight


def my_kernel(X,Y):
    #begin=time.clock()
    ilist=[]
    forlen=X.shape[0]
    yforlen=Y.shape[0]
    for i in range(forlen):
        jlist=[]
        for j in range(yforlen):
            oldflag=strclacstr(X[i,:],Y[j,:])
            jlist.append(oldflag)
        ilist.append(jlist)
    #print "核函数输入X",X.shape
    #print "核函数输入Y",Y.shape
    #print "核函数输出",np.array(ilist).shape
    # # #print(ilist)
    #end=time.clock()
    #print("用时:"+str(end-begin))
    return np.array(ilist)


def get_y_features(filename, chro):
    y_train = []
    y_test = []
    dic_p = {'CHS':'1', 'CHB':'2', 'JPT':'3', 'KPGP':'4'}
    file_features_train = 'chro_'+chro+'_features_train.txt'
    file_features_test = 'chro_'+chro+'_features_test.txt'
    f1 = open(file_features_train, 'w')
    f2 = open(file_features_test, 'w')
    
    line_y = commands.getoutput("cat " + filename + "|grep '^#CHROM' ")
    li_y = line_y.strip().split('\t')
    

    for i in range(9,len(li_y)):
        if i%2 == 0:
            y_train.extend([dic_p[li_y[i].split('_')[0]], dic_p[li_y[i].split('_')[0]]])
        else:
            y_test.extend([dic_p[li_y[i].split('_')[0]], dic_p[li_y[i].split('_')[0]]])


    str_train = ','.join(y_train)
    str_test = ','.join(y_test)
    f1.write(str_train+'\n')
    f2.write(str_test+'\n')
    f1.close()
    f2.close()

       
def get_x_features(filename, chro):
    file_features_train = 'chro_'+chro+'_features_train.txt'
    file_features_test = 'chro_'+chro+'_features_test.txt'
    f1 = open(file_features_train, 'a')
    f2 = open(file_features_test, 'a')
    lines_x = commands.getoutput("cat " + filename + "|grep '^"+chro+"[^0-9].*'").split("\n")
    for line in lines_x:
        train = []
        test = []
        li = line.strip().split('\t')
        for i in range(9, len(li)):
            if i%2 == 0:
                train.append(li[i].replace('|', ','))
            else:
                test.append(li[i].replace('|', ','))
        str_train = ','.join(train)
        str_test = ','.join(test)
        f1.write(str_train+'\n')
        f2.write(str_test+'\n')
    f1.close()
    f2.close()

    trainSet = np.loadtxt(file_features_train, delimiter=',')
    trainSet = np.transpose(trainSet)
    testSet = np.loadtxt(file_features_test, delimiter=',')
    testSet = np.transpose(testSet)
    return trainSet, testSet


def svm(trainSet, testSet, gap, len_window, chro):
    f = open('predict_temp_result_'+chro+'.txt', 'w')
    y_train = np.ndarray.tolist(trainSet[:, 0])
    y_test = np.ndarray.tolist(testSet[:, 0])
    f.write(str(y_test)[1:-1]+'\n')
    a, length = np.shape(trainSet)
    for i in range(1,length-len_window+gap,gap):
        x_train = trainSet[:, i:i+len_window]
        x_test = testSet[:, i:i+len_window]
        x_train = np.ndarray.tolist(x_train)
        x_test = np.ndarray.tolist(x_test)
        
        clf = SVC(kernel="linear") #线性核函数
        # clf = SVC(kernel="rbf") #径向基函数

        # clf = SVC(kernel=my_kernel) #string核函数
        clf.fit(x_train, y_train)
        answer = clf.predict(x_test)
        f.write(str(np.ndarray.tolist(answer))[1:-1]+'\n')
        show_accuracy(np.array(answer), np.array(y_test),'chro1')
    f.close()

    fd = np.transpose(np.loadtxt('predict_temp_result_'+chro+'.txt', delimiter=', '))
    np.savetxt('predict_final_result_chro_'+chro+'.txt', fd, fmt='%d', delimiter=' ')
    os.system("rm "+'predict_temp_result_'+chro+'.txt')

   
def show_accuracy(a, b, tip):
    acc = a.ravel() == b.ravel()
    print tip + '正确率：', np.mean(acc)



def main():
    '''chro:染色体号；gap:窗口每次滑动距离；len_window:窗口大小'''
    """usage example: python ancestrySVMphasing.py -chro 1 -gap 100 -win 100"""
    length = len(sys.argv)
    if length < 7:
        print """usage example: 
                 python ancestrySVMphasing.py -chro 1 -gap 100 -win 100
                 chro:染色体号；gap:窗口每次滑动距离；len_window:窗口大小
        """
        os._exit(0)
    else:
        for i in range(length):
            if sys.argv[i] == "-chro":
                chro = sys.argv[i+1]
            if sys.argv[i] == "-gap":
                gap = int(sys.argv[i+1])
            if sys.argv[i] == "-win":
                len_window = int(sys.argv[i+1])
        get_y_features('1000genome23andmeLoci.CHB_CHS_JPT_KPGP.RM.8W.vcf', chro)
        trainSet, testSet = get_x_features('1000genome23andmeLoci.CHB_CHS_JPT_KPGP.RM.8W.1.Phasing.vcf', chro)
        svm(trainSet, testSet, gap, len_window, chro)


if __name__ == '__main__':
    main()





