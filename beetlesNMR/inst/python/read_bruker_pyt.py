
import os
import sys
import zipfile
import tarfile
import shutil
import traceback

#import nmrglue as ng
import pandas as pd
import numpy as np


# hacky way to suppress deprecation warnings
import warnings

warnings.filterwarnings("ignore")
import nmrglue as ng

def read_data(infile):
  importedFiles = []
  if os.path.isdir(infile):
            spectraDirs = os.listdir(infile)
  for spec in spectraDirs:
    if os.path.isdir(os.path.join(infile,spec)):
      importedFiles.append(readBruker(os.path.join(infile,spec)))
    # TODO: write code...
    
  final_file = fileList2DataFrame(importedFiles)
  return(final_file)
  
  
  
    
def readBruker(path):
    dic, data = ng.bruker.read(path)
    udic = ng.bruker.guess_udic(dic, data)
    uc = ng.fileiobase.uc_from_udic(udic)
    ppm_scale = uc.ppm_scale()

    with open(path+'/pdata/1/title') as fl:
        name = fl.readline().rstrip()

    name.replace(' ','_')
    dic, data = ng.bruker.read_pdata(path+'/pdata/1/')

    ppm_scale = np.linspace(max(ppm_scale), min(ppm_scale), len(data))
    return((name, ppm_scale, data))

def writeCsv(fileList, outFile):
    dataDF = fileList2DataFrame(fileList)
    dataDF.to_csv(outFile, sep='\t', index=True, header=True)

def fileList2DataFrame(fileList):
    
    fileList_ = align_scales(fileList)
    
    data = {x[0]:x[2] for x in fileList_}
    ppm = fileList_[0][1]
    dataDF = pd.DataFrame.from_dict(data, orient='columns')
    dataDF.index = ppm
    return(dataDF)

def align_scales(list_of_spectra):
    """
    Aligning all spectra scales in ppm taking into account offsets after processing in Topspin
    """
    
    ex_ppm = list_of_spectra[0][1]
    ppm_per_pt = (max(ex_ppm) - min(ex_ppm))/ex_ppm.shape[0]
    
    min_ppm, max_ppm = 666, -666
    
    for spec in list_of_spectra:
        if spec[1][0] > max_ppm:
            max_ppm = spec[1][0]
        if spec[1][-1] < min_ppm:
            min_ppm = spec[1][-1]
            
    n_pts = int(np.round((max_ppm - min_ppm)/ppm_per_pt))
    new_ppm = np.linspace(max_ppm, min_ppm, n_pts)
    
    print(new_ppm[0], new_ppm[-1])
    
    output = []
    for i in range(len(list_of_spectra)):
        output.append((list_of_spectra[i][0],
                       new_ppm,
                       extend_spec(list_of_spectra[i][2], list_of_spectra[i][1], new_ppm)))
    
    return output
    
def extend_spec(data, ppm_old, ppm_new):
    """
    Pad spectra with zeroes to fill a new ppm range
    """

    idx_max = np.argmin(abs(ppm_new - ppm_old[0]))
    
    pad_start = idx_max
    pad_end = (ppm_new.shape[0] - ppm_old.shape[0] - pad_start)
    
    if pad_end<0:
        pad_start = pad_start-1
        pad_end = 0

    return(np.concatenate((np.zeros(pad_start), data, np.zeros(pad_end)), axis=0))

