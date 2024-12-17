#!/usr/bin/env python
"""
Simple script to generate igrf/dgrf files from update igrf coefficient file.
This script simply parse the coefficient file and extracts the update numbers. 

@author N. Omodei <nicola.omodei@stanford.edu>
"""
import numpy as np

def dump_dictionary(GRF_dictionary):
    base_name='igrf_data'
    for k in GRF_dictionary.keys():
        if '19' in k: continue
        output_path='%s/%s' % (base_name,k)
        with open(output_path,'w') as out_file:
            out_file.write('    %s\n' % k.replace('.dat',''))
            out_file.write(' 13  6371.2 %.1f\n' % float(k.replace('.dat','').split('grf')[-1].replace('s','')))
            for x in GRF_dictionary[k]:
                out_file.write('%s\n' % x)
        print('file %s saved' % output_path)



def parse_file(filename):
    GRF_dictionary={}
    with open(filename) as my_file:
        lines = my_file.readlines()
        i=0
        while lines[i][0]=='#': i+=1
        i-=1
        #print(lines[i])
        nyears = len(lines[i].split())-3
        #print(nyears)
        i+=1
        #print(lines[i])
        types=lines[i].split()[3:]
        i+=1
        years=lines[i].split()[3:]
        i+=1   
        file_names = []
        #print(years)
        #print(types)
        for t,y in zip(types,years): 
            try: 
                y='%d' % float(y)
                t=t.lower()
            except:
                y='%ds' % float(y.split('-')[0])
                t='igrf'
            output_file_name='%s%s.dat' % (t,y)
            file_names.append(output_file_name)
            #print (output_file_name)
            GRF_dictionary[output_file_name]=[]
        #print('-----------------')
        n=len(lines)
        #n=5
        while i < n:
            for j,l in enumerate(lines[i].split()[3:]):
                output_file_name = file_names[j]                
                GRF_dictionary[output_file_name].append(float(l))
                #print(i,j,output_file_name,l)
            i+=1
        #print (GRF_dictionary)
        dump_dictionary(GRF_dictionary)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='generate_igrf_from_coeff',
                        description='Extract igrf/dgrf files from update igrf coefficient file',
                        epilog='Fermi Large Area Telescope')
    parser.add_argument('-f','--file')
    args = parser.parse_args()
    parse_file(args.file)

