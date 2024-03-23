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
outdir = './outputLR/'
datadir = '/Users/egold/data/KOA'
datetimerange = '2021-09-07 00:00:00/2021-09-08 23:59:59'
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
    datadir,\
              calibfile=1 )


#code to orgainze files into different directories for pypeit to run through them 
from astropy.io import fits


#Organize Science Files 
for name in glob.glob(datadir + '/lev0/*.fits'):
    f = fits.open(name)
    if f[0].header["SLITNAME"] == 'long_1.0':
        if os.path.exists(datadir + "/" + f[0].header["DATE-OBS"]):
            os.system("cp " + name + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
        else:
            os.system("mkdir " + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
            os.system("cp " + name + " " + datadir + "/" +  f[0].header["DATE-OBS"])
       



date_lt = []

for name in glob.glob(datadir + '/calib/*.fits'):
    f = fits.open(name)
    if f[0].header["SLITNAME"] == 'long_1.0':
        if os.path.exists(datadir + "/" + f[0].header["DATE-OBS"]):
            os.system("cp " + name + " " + datadir + "/" + f[0].header["DATE-OBS"] )
        else:
            os.system("mkdir " + " " + datadir + "/" +  f[0].header["DATE-OBS"] )
            os.system("cp " + name + " " + datadir + "/" + f[0].header["DATE-OBS"])
        date_lt.append(f[0].header["DATE-OBS"])        

date_lt = list(set(date_lt))
 
count_arc = 0 
count_flat = 0 
total = 0


#arcs
for date_dir in date_lt: 

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
			os.system("pypeit_setup -s keck_lris_blue -r  "+ datadir +"/" + date_dir +"/LB -c all")		
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
			os.system('pypeit_sensfunc ' +  file  + ' -o ' + ' /sens_outputs '  + file + '_sens')	        
		print('rm -r ' + datadir + '/*') 
		print('rm -r ' + outdir + '/*')
		print('rm -r ' + 'QA/*')
		print('rm -r ' + 'Callibrations/*')
		print('rm -r ' + 'Science/*')
		print('rm -r ' + 'outputLR/*') 
		print('rm -r ' + 'keck_lris_*/*')   

	except Exception as e: 
		print(e)	
		print("Fluxing did not work")

