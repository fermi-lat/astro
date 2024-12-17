#!/usr/bin/env python
import sys,os
import matplotlib
matplotlib.use('Agg')

import scipy as sp
import datetime
import matplotlib as mpl
from matplotlib import pyplot as plt
import argparse


IGRF_old='test_IGRF_14_output.txt'
IGRF_new='test_IGRF_14_output.txt'
igrf_14_file = 'igrf_14.dat' # this file has been created running the code igrf_14.f (f77 igrf14.f -o igrf14.o and then ./igrf14.o > igrf_14.dat (all inside the contaier))

def readFileFortran(input_file_name,min_year=2009,max_year=2019):
    print('Reading fortran file %s' % input_file_name)

    my_dictionary={
        'date':[],
        'bEast':[],
        'bNorth':[],
        'bDown':[],
        }
        
    for l in open(input_file_name,'r').readlines():
        values = l.split()
        if len(values)!=16: raise BaseException
        year    = int(values[0])
        month   = int(values[1])
        if year>max_year or year<min_year: continue
        date   = datetime.date(year,month,1)
        my_dictionary['date'].append(date)
        my_dictionary['bNorth'].append(float(values[8])*1e-5) # nT -> G
        my_dictionary['bEast'].append(float(values[10])*1e-5) # nT -> G
        my_dictionary['bDown'].append(float(values[12])*1e-5) # nT -> G
        pass
    for k in my_dictionary.keys(): my_dictionary[k]=sp.array(my_dictionary[k])
    return my_dictionary

def readFile(input_file_name,min_year=2009,max_year=2019):

    my_dictionary={
        'date':[],
        'met':[],
        'l':[],
        'b':[],
        'lambda':[],
        'bEast':[],
        'bNorth':[],
        'bDown':[],
        'R':[],
        'verticalRigidityCutoff':[]
        }
    print('Reading file %s' % input_file_name)        
    for l in open(input_file_name,'r').readlines():
        values = l.split()
        if len(values)!=11: raise BaseException
        month   = int(values[0])
        year    = int(values[1])
        if year>max_year or year<min_year: continue
        date   = datetime.date(year,month,1)

        my_dictionary['date'].append(date)
        my_dictionary['met'].append(float(values[2]))
        my_dictionary['l'].append(float(values[3]))
        my_dictionary['b'].append(float(values[4]))
        my_dictionary['lambda'].append(float(values[5]))
        my_dictionary['bEast'].append(float(values[6]))
        my_dictionary['bNorth'].append(float(values[7]))
        my_dictionary['bDown'].append(float(values[8]))
        my_dictionary['R'].append(float(values[9]))
        my_dictionary['verticalRigidityCutoff'].append(float(values[10]))
        pass
    for k in my_dictionary.keys(): my_dictionary[k]=sp.array(my_dictionary[k])
    return my_dictionary


def compare(file1,file2):
    variable_to_plot=['l','b','lambda','bEast','bNorth','bDown','R','verticalRigidityCutoff']
    variable_to_plot=['bEast','bNorth','bDown']
    min_year=2019
    max_year=2029

    try:        igrf1 = readFile(file1,min_year,max_year)
    except:     igrf1 = readFileFortran(file1,min_year,max_year)

    try:        igrf2 = readFile(file2,min_year,max_year)
    except:     igrf2 = readFileFortran(file2,min_year,max_year)


    fig, axs = plt.subplots(len(variable_to_plot),2,sharex='col',figsize=(15,10))
    plt.subplots_adjust(hspace=0.001)
    to_print='bEast'
    print(igrf1['date'][0],igrf1[to_print][0])
    print(igrf2['date'][0],igrf2[to_print][0])
    for iy,v in enumerate(variable_to_plot):

        if v in igrf1.keys(): axs[iy,0].plot(igrf1['date'],igrf1[v],'r',label='%s %s' % (v,file1.split('/')[-1]))
        if v in igrf2.keys(): axs[iy,0].plot(igrf2['date'],igrf2[v],'g',label='%s %s' % (v,file2.split('/')[-1]))

        if v in igrf1.keys() and v in igrf2.keys():
            print(v,len(igrf1[v]),len(igrf2[v]),len((igrf1[v]-igrf2[v])/igrf1[v]))
            axs[iy,1].plot(igrf1['date'],(igrf1[v]-igrf2[v])/igrf1[v],'b',label='%s' % v)

        axs[iy,0].set_xlabel('date')
        axs[iy,0].set_ylabel(v)

        axs[iy,1].set_xlabel('date')
        axs[iy,1].set_ylabel('(I1-I2)/I1')
        axs[iy,1].set_ylim(-1e-2,1e-2)
        pass
    axs[0, 0].set_title("I1=%s (red) \n I2=%s (green)" % (file1.split('/')[-1],file2.split('/')[-1]), size="large")
    #plt.show()
    plt.savefig('IGRF_comparison.png')
    print('figure saved in IGRF_comparison.png')
    
    

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='comparison_IGRFs',
                        description='Comparebtwo sets of output for IGRF moel',
                        epilog='Fermi Large Area Telescope')
    parser.add_argument('-f1','--file1',required=True)
    parser.add_argument('-f2','--file2',required=True)
    args = parser.parse_args()
    compare(args.file1,args.file2)

    
