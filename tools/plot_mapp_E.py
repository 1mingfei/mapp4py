#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import sys
import os
'''
Program header:
Mingfei Zhang
mingfei@umich.edu
a code to plot E against MC simulation steps
python2 code
'''
def get_start_line(filename):
    with open(filename,'r')as fin:
        i = 0
        for line in fin:
            if '|     step    |      T      |' in line:
                start = i
                break
            i = i + 1
    return start+2

def get_end_line(filename):
    with open(filename,'r')as fin:
        i = 0
        for line in fin:
            if 'time elapsed:' in line:
                start = i
                break
            i = i + 1
    return start-2

def fetch_step_energy(filename):
    start, end = get_start_line(filename), get_end_line(filename)
    data = np.zeros((end-start+1,3))
    with open(filename,'r') as fin:
        lines = fin.readlines()
        for i in range(start, end+1):
            data[i-start][0] = float(lines[i].split()[1])
            data[i-start][1] = float(lines[i].split()[5])
            data[i-start][2] = float(lines[i].split()[19])/float(lines[i].split()[21])
    print data
    fig1, ax1 = plt.subplots()
    ax1.plot(data[:,0], data[:,1])
    ax1.set(xlabel='timestep', ylabel='energy (eV)')
    ax1.grid()
    fig1.savefig("plots/step_E.pdf")

    fig2, ax2 = plt.subplots()
    ax2.plot(data[:,0], data[:,2])
    ax2.set(xlabel='timestep', ylabel='element number ratio')
    ax2.grid()
    fig2.savefig("plots/step_ele_ratio.pdf")
    return

if __name__ == "__main__":
    if (len(sys.argv) != 2):
        print("usage:\nplot_mapp_E.py <mapp output file>\n")
    else:
        if (not os.path.isdir('plots')):
           os.mkdir('plots')
           inFile = sys.argv[1]
           print inFile
           fetch_step_energy(inFile)
