##taking a stab at multiprocessing the pypeit step
from multiprocessing import Pool, Manager
import multiprocessing
#Imports to get the codes to actually run
import sys
import io
import os
import glob
from pykoa.koa import Koa 
from astropy.table import Table,Column
#imports needed for the pypeit run code stuff 
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import simple_norm


#Input variables 
instrument = 'LRIS'
instr = 'lris'  
outdir = '/k2drpdata/LUNDQUIST/LRIS_Performance_Check/outputLR/'
datadir = '/k2drpdata/LUNDQUIST/LRIS_Performance_Check/data/KOA'
#outdir = './outputLR/'
#datadir = '/Users/egold/data/KOA'
datetimerange = '2018-01-01 00:00:00/2018-02-01 23:59:59'
targname = 'G191B2B'

side='red'
if side=='red' and datetimerange >= '2022-01-01':
    detector='keck_lris_red_markiv'
    filenamestart='LR'
elif side =='red' and datetimerange <= '2009-01-01':
    detector='keck_lris_red_orig'
    filenamestart='LR'
elif side=='red' and datetimerange >= '2009-01-01':
    detector='keck_lris_red'
    filenamestart='LR'
else:
    detector='keck_lris_blue'
    filenamestart='LB'


def pypeit_code(inp):
    print(inp)
    date_dir=inp[0]
    print(date_dir)
    if True:
        count_arc = 0
        count_flat = 0
        total = 0

        for name in glob.glob(datadir + '/' + date_dir +  '/*.fits'):
            f = fits.open(name)

            if (f[0].header ["MERCURY"] == 'on' or f[0].header["NEON"] == 'on' or f[0].header["ARGON"] == 'on' or f[0].header["CADMIUM"] == 'on' or f[0].header["ZINC"] == 'on' or f[0].header["KRYPTON"] == 'on' or f[0].header["XENON"] == 'on' or f[0].header["FEARGON"] == 'on') and f[0].header["TRAPDOOR"] == "closed":
                print("arcs found" , name)
                count_arc = count_arc +1



    #flats
            if f[0].header ["HALOGEN"] == 'on' or f[0].header ["FLAMP1"] == 'on' or f[0].header ["FLAMP2"]  == 'on':
                print("flats found" , name)
                count_flat = count_flat +1
            total = total +1


            if count_arc >= 1 and count_flat >= 3:
                print("Count Arc and Count Flat", count_arc , count_flat)
                break

    #Running Pypeit
        try:
            print("Now Running Pypeit")
            if count_flat >=3 and count_arc >= 1:
                os.system("pypeit_setup -s detector -r  "+ datadir +"/" + date_dir +"/filenamestart -c all")
                calib = "[baseprocess]\n        use_biasimage = False\n[reduce]\n       [[skysub]]\n            bspline_spacing=0.6\n               local_maskwidth=2.0\n               sky_sigrej=3.0\n"
                files=glob.glob('*/*.pypeit')

                for dataset in files:
                    print("dataset", dataset)
                    with open(dataset, "r") as file:
                        data = file.read()
                        data = data.replace("# User-defined execution parameters", str(calib))
                    with open(dataset, "w") as file:
                        file.write(data)

                    os.system('run_pypeit ' +  dataset +  ' -o')
            else:
                print("Not enough flats or arcs to run pypeit")
        except Exception as e:
            print(e)

#Running Fluxing
        try:
            print("Now Running Fluxing")

            files=glob.glob("Science/spec1d**.fits")
            if os.path.isdir('sens_outputs')==False:
                os.system('mkdir sens_outputs')
            for file in files:
                print(file)
                filename=os.path.basename(file).split('/')[-1]
                os.system('pypeit_sensfunc ' +  file  + ' -o ' + ' ./sens_outputs/'  + filename + '_sens')
            os.system('rm -r ' + 'QA/*')
            os.system('rm -r ' + 'Calibrations/*')
            os.system('rm -r ' + 'Science/*')
            os.system('rm -r ' + 'keck_lris_*/*')

        except Exception as e:
            print(e)
            print("Fluxing did not work")






if True:

    #makes directory for files to be downloaded into 
    try:
        os.mkdir(outdir)
    except:
        print(instrument + " Directory exists already", flush=True)
          
    #finds the files looking for with variables inputed
    param = dict()
    param['instrument'] = instr
    param['datetime'] = datetimerange
    param['target'] = targname
    Koa.query_criteria (param, \
        outdir + instr + "_parameters.tbl", overwrite=True )
    rec = Table.read(outdir + instr + "_parameters.tbl", format='ipac')
    print (rec)

#downlaods file to KOA directory 
    Koa.download (outdir + instr + "_parameters.tbl", \
        'ipac', \
        datadir,\
                  calibfile=1 )
masks=['long_1.0','long=1.5','long=0.7']    

date_lt = []
for mask in masks:
    #code to orgainze files into different directories for pypeit to run through them 
    from astropy.io import fits
    
    
    #Organize Science Files 
    for name in glob.glob(datadir + '/lev0/*.fits'):
        f = fits.open(name)
        try:
            if f[0].header["SLITNAME"] == mask:
                if os.path.exists(datadir + "/" + f[0].header["DATE-OBS"]):
                    os.system("cp " + name + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
                    pass
                else:
                    os.system("mkdir " + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
                    os.system("cp " + name + " " + datadir + "/" +  f[0].header["DATE-OBS"])
                date_lt.append(f[0].header["DATE-OBS"])
        except Exception as e:   
            print('Error with filename: '+name)
            print(e)
    
    
#    date_lt = []
    for name in glob.glob(datadir + '/calib/*.fits'):
        f = fits.open(name)
        try:
            if f[0].header["SLITNAME"] == mask:
                if os.path.exists(datadir + "/" + f[0].header["DATE-OBS"]):
                    os.system("cp " + name + " " + datadir + "/" + f[0].header["DATE-OBS"] )
                    pass
                else:
                    os.system("mkdir " + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
                    os.system("cp " + name + " " + datadir + "/" + f[0].header["DATE-OBS"])
#                date_lt.append(f[0].header["DATE-OBS"])        

        except Exception as e:
            print('Error with filename: '+name)
            print(e)

    date_lt = list(set(date_lt))

    print(date_lt)
    count_arc = 0 
    count_flat = 0 
    total = 0
    
    #arcs
#    pool = Pool(multiprocessing.cpu_count(), maxtasksperchild=1)
    pool = Pool(1, maxtasksperchild=1)

    jobs=[]
    for date_dir in date_lt:
        print('apply_async')
        jobs.append(pool.apply_async(pypeit_code, [[date_dir]]))  # add waiting files to pool
    

    pool.close()  # close pool
    pool.join()  # join pool

    for job in jobs:  # wait unitl all science images are finished before exiting
        try:
            job.get()
        except IOError as e:
            print(e)


if True:
    print('Cleaning up KOA')
    os.system('rm -r ' + datadir + '/*')
    os.system('rm -r ' + outdir + '/*')

