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
import glob
import os


#Input variables 
instrument = 'LRIS'
instr = 'lris'  
outdir = './outputLR/'
datetimerange = '2017-02-25 00:00:00/2017-03-25 23:59:59'
targname = 'G191B2B '

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
    '/Users/egold/data/KOA' ,\
              calibfile=1 )


#code to orgainze files into different directories for pypeit to run through them 
os.system ("ls dnload_dir2/*")
file = "dnload_dir2/lev0/LB.20170319.19187.fits"          
from astropy.io import fits

f = fits.open(file)
f[0].header ["SLITNAME"]

#Organize Science Files 
for name in glob.glob('dnload_dir2/lev0/*.fits'):
    f = fits.open(name)
    if f[0].header["SLITNAME"] == 'long_1.0':
        if os.path.exists(f[0].header["DATE-OBS"]):
         #   print("path exists")
            os.system("cp " + name + " " + f[0].header["DATE-OBS"] )
        else:
          #  print("path doesnt exist making directory " + f[0].header["DATE-OBS"])
            os.system("mkdir " + " " + f[0].header["DATE-OBS"] )
            os.system("cp " + name + " " + f[0].header["DATE-OBS"])
      #  print(name,f[0].header["DATE-OBS"])
       
# Organize Calib Files 

for file in glob.glob('dnload_dir2/calib/*.fits'):
    from astropy.io import fits
    f = fits.open(file) 
    f[0].header ["SLITNAME"]

for name in glob.glob('dnload_dir2/calib/*.fits'):
    f = fits.open(name)
    if f[0].header["SLITNAME"] == 'long_1.0':
        if os.path.exists(f[0].header["DATE-OBS"]):
            #print("path exists")
            os.system("cp " + name + " " + f[0].header["DATE-OBS"] )
        else:
           # print("path doesnt exist making directory " + f[0].header["DATE-OBS"])
            os.system("mkdir " + " " + f[0].header["DATE-OBS"] )
            os.system("cp " + name + " " + f[0].header["DATE-OBS"])
       # print(name,f[0].header["DATE-OBS"])
        
 #Find files with Arcs and Flats 
#print("PWD: ")
os.system ("pwd")
#print(" ")

for name in glob.glob ("Users/egold/data/blue/keck_lris_blue_B/dnload_dir2/calib/*.fits"):
    from astropy.io import fits
    f = fits.open(file) 
#f[0].header ["MERCURY" or "NEON" or "ARGON" or "CADMIUM" or "ZINC" or "KRYPTON" or "XENON" or "FEAERGON"]

count_arc = 0 
count_flat = 0 
total = 0



file = "2017-03-19/LB.20170319.04357.fits"
from astropy.io import fits
f = fits.open(file) 
f[0].header ["MERCURY" or "NEON" or "ARGON" or "CADMIUM" or "ZINC" or "KRYPTON" or "XENON" or "FEAERGON"]


#arcs
for name in glob.glob('2017-03-19/*.fits'):
    f = fits.open(name)
    
    if f[0].header ["MERCURY"] == 'on' or f[0].header["NEON"] == 'on' or f[0].header["ARGON"] == 'on' or f[0].header["CADMIUM"] == 'on' or f[0].header["ZINC"] == 'on' or f[0].header["KRYPTON"] == 'on' or f[0].header["XENON"] == 'on' or f[0].header["FEARGON"] == 'on':
        print("arcs found" , name)
        count_arc = count_arc +1 
        
    
#flats 
    if f[0].header ["HALOGEN"] == 'on' or f[0].header ["FLAMP1"] == 'on' or f[0].header ["FLAMP2"]  == 'on':
        print("flats found" , name)
        count_flat = count_flat +1
    total = total +1
    
    
    if count_arc >= 3 and count_flat >= 3: 
        print("Count Arc and Count Flat", count_arc , count_flat)
        break 
    
if count_flat and count_arc >= 3: 
        print("")
        os.system(run_pypeit keck_lris_blue_B.pypeit -o)
        
    
else: 
    print("not enough flats/arcs")


