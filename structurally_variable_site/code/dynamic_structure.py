#-*- coding:utf-8 -*-
import sys
import numpy as np
from scipy.spatial.distance import pdist
from scipy import stats

def read_shape(shape_file):
    """
    shape_file       -- A file of shape data
    ENST00000341426 3235    7.714   NULL    NULL    NULL    NULL    NULL
    """
    data = {}
    f = open(shape_file, 'r')
    i = 0
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        sent1 = sent[0].split('.')
        data[sent1[0]] = sent[3:]
    plen = i
    return data
    f.close()

def calc_dist(list_of_values1, list_of_values2, min_num=10):
    """
    list_of_values1, list_of_values2       -- Two list of float values
    
    """
    value1 = []
    value2 = []
    length1 = len(list_of_values1)
    Len1 = 0
    for i in range(length1):
        if (list_of_values1[i] != 'NULL') and (list_of_values2[i] != 'NULL'):
            value1.append(float(list_of_values1[i]))
            value2.append(float(list_of_values2[i]))
            Len1 += 1 
    x = np.array(value1)
    y = np.array(value2)
    if Len1 > min_num:
        dist1 = np.linalg.norm( x - y )
        if (sum(x) != 0) and (sum(y) != 0):
            dist2 = 1 - np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))
        else:
            dist2 = dist1
    else:
        dist1 = -1
        dist2 = -1
    
    return dist1, dist2
   

def calc_shape_test(list_of_values1, list_of_values2, min_num=10):
    """
    list_of_values1, list_of_values2       -- Two list of float values
    
    """
    score_cutoff = 0.546
    value1 = []
    value2 = []
    length1 = len(list_of_values1)
    Len1 = 0
    test_estimate = 0
    test_num = 0
    for i in range(length1):
        if (list_of_values1[i] != 'NULL') and (list_of_values2[i] != 'NULL'):
            value1.append(abs(float(list_of_values1[i]) - float(list_of_values2[i])))
            if (abs(float(list_of_values1[i]) - float(list_of_values2[i])) >= score_cutoff):
                test_estimate = test_estimate + 1
            test_num = test_num + 1
            Len1 += 1 
    x = np.array(value1)
    if Len1 >= min_num:
        l1_dist = x.mean()
    else:
        test_estimate = -1
        test_num = -1
        l1_dist = -1
    return test_estimate, test_num, l1_dist



def slide_window(Shape_value1, Shape_value2, wlen=20):
    """
    Shape_value1, Shape_value2         -- two lists of SHAPE scores
    wlen                  -- the length of window
    
    """
    len1 = len(Shape_value1)
    len2 = len(Shape_value2)

    test_est = []
    test_p = []
    ll_dist1 = []
    ll_dist2 = []
    total_dist = []
    if len1 != len2:
        return test_est, test_p, ll_dist1, ll_dist2
    for i in range(len1 - wlen):
        x = Shape_value1[i:(i+wlen)]
        y = Shape_value2[i:(i+wlen)]
        test_value, test_pvalue, d1 = calc_shape_test(x, y)
        d2, d3 = calc_dist(x,y)
        #if(test_pvalue < 0.05):
        test_est.append(test_value)
        test_p.append(test_pvalue)
        ll_dist1.append(d1)
        ll_dist2.append(d2)
    return test_est, test_p, ll_dist1, ll_dist2



if __name__ == '__main__':
    cell1 = sys.argv[1]
    cell2 = sys.argv[2]
    outwin_file = sys.argv[3]
    path1 = "/150T/zhangqf2/huangwz/total_smart_icshape/new_smartSHAPE_0412/"
    shape_file1 = path1 + cell1 + "_smartSHAPE.out"
    shape_file2 = path1 + cell2 + "_smartSHAPE.out"
    data1 = read_shape(shape_file1)
    data2 = read_shape(shape_file2) 
    
    num = 0
    fw = open(outwin_file, 'w')
    for ke in data1.keys():
        if (ke in data2.keys()):
            num += 1
    print("total_trx\t%s\t%s\t%d\n"%(cell1, cell2, num))
    num = 0
    for ke in data1.keys():
        num = num + 1
        #if num > 10:
        #    break
        if (ke in data2.keys()):
            value1 = data1[ke]
            value2 = data2[ke]
            test_est, test_p, ll_dist1, ll_dist2 = slide_window(value1, value2)
            for i in range(len(test_est)):
                if (test_est[i] != -1) and (test_p[i] != -1) and (ll_dist1[i] != -1) and (ll_dist2[i] != -1):
                    outstr = ke + "\t" + str(i) + "\t" + str(round(test_est[i],5)) + "\t" + str(round(test_p[i],5)) + "\t" + str(round(ll_dist1[i],5)) + "\t" + str(round(ll_dist2[i],5)) +"\n"
                    fw.write(outstr)
    fw.close()
