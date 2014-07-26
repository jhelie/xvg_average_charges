#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_average_charges_simple', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_average_charges
**************************************************

[ DESCRIPTION ]
 
This script calculate the average of charge densities data contained in several xvg files.

It calculates the avg and (unbiasd) std dev and can deal with NaN.

NB:
the script may give out a warning 'return np.mean(x,axis)/factor', it's ok. it's just
scipy warning us that there were only nans on a row, the result will be a nan as we
expect (see this thread: https://github.com/scipy/scipy/issues/2898).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o	charges_avg	: name of outptut file
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["charges_avg"], help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if len(args.xvgfilenames) == 1:
	print "Error: only 1 data file specified."
	sys.exit(1)
	
for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
		
	global nb_rows
	global nb_cols
	global weights
	global distances
	global data_ions
	global data_peptide
	global data_lipids
	global data_total

	nb_rows = 0
	nb_cols = 0
	weights = np.ones(len(args.xvgfilenames))
		
	for f_index in range(0,len(args.xvgfilenames)):
		#display progress
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.xvgfilenames[f_index]
		with open(filename) as f:
			lines = f.readlines()
		
		#determine legends and nb of lines to skip
		tmp_nb_rows_to_skip = 0
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if "weight" in line:
					if "-> weight = " in line:
						weights[f_index] = float(line.split("-> weight = ")[1])
						if weights[f_index] < 0:
							print "\nError: the weight in file " + str(filename) + " should be a positive number."
							print " -> " + str(line)
							sys.exit(1)
					else:
						print "\nWarning: keyword 'weight' found in the comments of file " + str(filename) + ", but weight not read in as the format '-> weight = ' wasn't found."
		
		#get data
		tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			distances = np.zeros((nb_rows, 1))								#distance from cluster
			data_ions = np.zeros((nb_rows, len(args.xvgfilenames)))			#ion charges density
			data_peptide = np.zeros((nb_rows, len(args.xvgfilenames)))		#peptide charges density
			data_lipids = np.zeros((nb_rows, len(args.xvgfilenames)))		#lipids charges density
			data_total = np.zeros((nb_rows, len(args.xvgfilenames)))		#total charge density
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 0:
			nb_cols = np.shape(tmp_data)[1]
		else:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			distances[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],distances[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		data_ions[:,f_index] = tmp_data[:,1]
		data_peptide[:,f_index] = tmp_data[:,2]
		data_lipids[:,f_index] = tmp_data[:,3]
		data_total[:,f_index] = tmp_data[:,4]

	#replace non sampled values by nan
	data_ions[data_ions == 0] = np.nan
	data_peptide[data_peptide == 0] = np.nan
	data_lipids[data_lipids == 0] = np.nan
	data_total[data_total == 0] = np.nan
	
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_avg():													#DONE

	global avg_charge_ions
	global avg_charge_peptide
	global avg_charge_lipids
	global avg_charge_total
	global std_charge_ions
	global std_charge_peptide
	global std_charge_lipids
	global std_charge_total
				
	avg_charge_ions = np.zeros((nb_rows,1))
	avg_charge_peptide = np.zeros((nb_rows,1))
	avg_charge_lipids = np.zeros((nb_rows,1))
	avg_charge_total = np.zeros((nb_rows,1))
	std_charge_ions = np.zeros((nb_rows,1))
	std_charge_peptide = np.zeros((nb_rows,1))
	std_charge_lipids = np.zeros((nb_rows,1))
	std_charge_total = np.zeros((nb_rows,1))

	#remove nan values of the weights for upper species
	#--------------------------------------------------
	#ions
	weights_ions_nan = np.zeros((nb_rows, 1))	
	weights_ions_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_ions = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_ions[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_ions[r,0] -= 1
	weights_ions_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_ions_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_ions_nan[weights_ions_nan == 0] = 1
	
	#peptide
	weights_peptide_nan = np.zeros((nb_rows, 1))	
	weights_peptide_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_peptide = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_peptide[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_peptide[r,0] -= 1
	weights_peptide_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_peptide_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_peptide_nan[weights_peptide_nan == 0] = 1

	#lipids
	weights_lipids_nan = np.zeros((nb_rows, 1))	
	weights_lipids_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_lipids = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_lipids[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_lipids[r,0] -= 1
	weights_lipids_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_lipids_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_lipids_nan[weights_lipids_nan == 0] = 1

	#total
	weights_total_nan = np.zeros((nb_rows, 1))	
	weights_total_nan_sq = np.zeros((nb_rows, 1))	
	nb_files_total = np.ones((nb_rows, 1)) * len(args.xvgfilenames)
	tmp_weights_nan = np.zeros((nb_rows, len(args.xvgfilenames)))
	for r in range(0, nb_rows):
		tmp_weights_nan[r,:] = weights
		for f_index in range(0, len(args.xvgfilenames)):
			if np.isnan(data_total[r,f_index]):
				tmp_weights_nan[r,f_index] = 0
				nb_files_total[r,0] -= 1
	weights_total_nan[:,0] = np.nansum(tmp_weights_nan, axis = 1)
	weights_total_nan_sq[:,0] = np.nansum(tmp_weights_nan**2, axis = 1)	
	weights_total_nan[weights_total_nan == 0] = 1

	#calculate weighted average taking into account "nan"
	#----------------------------------------------------
	avg_charge_ions[:,0] =  scipy.stats.nanmean(data_ions * weights * nb_files_ions / weights_ions_nan, axis = 1)
	avg_charge_peptide[:,0] =  scipy.stats.nanmean(data_peptide * weights * nb_files_peptide / weights_peptide_nan, axis = 1)
	avg_charge_lipids[:,0] =  scipy.stats.nanmean(data_lipids * weights * nb_files_lipids / weights_lipids_nan, axis = 1)
	avg_charge_total[:,0] =  scipy.stats.nanmean(data_total * weights * nb_files_total / weights_total_nan, axis = 1)
	
	#calculate unbiased weighted std dev taking into account "nan"
	#-------------------------------------------------------------
	tmp_ions = np.zeros((nb_rows, 1))
	tmp_ions[:,0] = np.nansum(weights * (data_ions - avg_charge_ions[:,0:1])**2, axis = 1)			
	tmp_div = np.copy((weights_ions_nan)**2 - weights_ions_nan_sq)
	tmp_div[tmp_div == 0] = 1
	std_charge_ions = np.sqrt(weights_ions_nan / tmp_div * tmp_ions)
	
	tmp_peptide = np.zeros((nb_rows, 1))
	tmp_peptide[:,0] = np.nansum(weights * (data_peptide - avg_charge_peptide[:,0:1])**2, axis = 1)			
	tmp_div = np.copy((weights_peptide_nan)**2 - weights_peptide_nan_sq)
	tmp_div[tmp_div == 0] = 1
	std_charge_peptide = np.sqrt(weights_peptide_nan / tmp_div * tmp_peptide)

	tmp_lipids = np.zeros((nb_rows, 1))
	tmp_lipids[:,0] = np.nansum(weights * (data_lipids - avg_charge_lipids[:,0:1])**2, axis = 1)			
	tmp_div = np.copy((weights_lipids_nan)**2 - weights_lipids_nan_sq)
	tmp_div[tmp_div == 0] = 1
	std_charge_lipids = np.sqrt(weights_lipids_nan / tmp_div * tmp_lipids)

	tmp_total = np.zeros((nb_rows, 1))
	tmp_total[:,0] = np.nansum(weights * (data_total - avg_charge_total[:,0:1])**2, axis = 1)			
	tmp_div = np.copy((weights_total_nan)**2 - weights_total_nan_sq)
	tmp_div[tmp_div == 0] = 1
	std_charge_total = np.sqrt(weights_total_nan / tmp_div * tmp_total)

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_average_charges_simple v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	if np.sum(weights) > len(args.xvgfilenames):
		output_xvg.write("# -> weight = " + str(np.sum(weights)) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"distance from cluster z axis (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"charge densities (e.A-3)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 8\n")
	output_xvg.write("@ s0 legend \"ions (avg)\"\n")
	output_xvg.write("@ s1 legend \"peptide (avg)\"\n")
	output_xvg.write("@ s2 legend \"lipids (avg)\"\n")
	output_xvg.write("@ s3 legend \"total (avg)\"\n")
	output_xvg.write("@ s5 legend \"ions (std)\"\n")
	output_xvg.write("@ s6 legend \"peptide (std)\"\n")
	output_xvg.write("@ s7 legend \"lipids (std)\"\n")
	output_xvg.write("@ s8 legend \"total (avg)\"\n")

	#data
	for r in range(0, nb_rows):
		results = str(distances[r,0])
		results += "	" + "{:.6e}".format(avg_charge_ions[r,0]) + "	" + "{:.6e}".format(avg_charge_peptide[r,0]) + "	" + "{:.6e}".format(avg_charge_lipids[r,0]) + "	" + "{:.6e}".format(avg_charge_total[r,0]) + "	" + "{:.6e}".format(std_charge_ions[r,0]) + "	" + "{:.6e}".format(std_charge_peptide[r,0]) + "	" + "{:.6e}".format(std_charge_lipids[r,0]) + "	" + "{:.6e}".format(std_charge_total[r,0])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\n\nWriting average file..."
calculate_avg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
