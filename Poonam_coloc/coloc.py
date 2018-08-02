#@ File    (label = "Input directory", style = "directory") srcFile
#@ String  (label = "File extension", value=".tif") ext
#@ Float (label = "Channel 1 approx size", value = 4.0) ch1size
#@ Float (label = "Channel 1 relative threshold", value = 3.0) ch1thresh
#@ Float (label = "Channel 2 approx size", value = 4.0) ch2size
#@ Float (label = "Channel 2 relative threshold", value = 3.0) ch2thresh
#@ Float (label = "Maximum colocalisation distance", value = 4.0) coloc

"""

coloc.py 
created by: Erick Martins Ratamero
date: 01/08/18
last updated: 01/08/18

Opens a directory of images and calculates green on red and 
red on green colocalisation. Separates bottom of nucleus calculations
from the rest of the stack.

"""
import os
from java.io import File

from ij import IJ, ImageStack, ImagePlus
from ij.plugin.frame import RoiManager
import math
from ij import WindowManager
from ij.measure import ResultsTable
from loci.plugins import BF
from ij.io import FileSaver 
from ij.process import ImageStatistics as IS  

srcDir = srcFile.getAbsolutePath()

totchannels = 3


def retrieve_channels(image, channels):
	# get stack from current image
	stack = image.getStack()
	#print(image.getStackSize(), stack.getSize())
	final_stack = ImageStack(image.width, image.height)

	for i in range(1, image.getStackSize()+1):
		
		countchannel = i % totchannels
		#print(i, countchannel, countchannel in channels)
		if (countchannel in channels):
			
			myslice = stack.getProcessor(i)
		
			final_stack.addSlice(str(i), myslice)
	return final_stack

def retrieve_dapi(image, channel):
	# get stack from current image
	stack = image.getStack()
	#print(image.getStackSize(), stack.getSize())
	final_stack = ImageStack(image.width, image.height)

	for i in range(1, image.getStackSize()+1):
		
		countchannel = i % totchannels
		#print(i, countchannel, countchannel in channels)
		if (countchannel == channel):
			
			myslice = stack.getProcessor(i)
		
			final_stack.addSlice(str(i), myslice)
	return final_stack

def run_comdet(image):
	IJ.run("Detect Particles", "include two=[Detect in both channels independently] ch1a="+str(ch1size) + " ch1s=" + str(ch1thresh) + " ch2a="+str(ch2size) + " ch2s=" + str(ch2thresh) + " calculate max=" + str(coloc) + " add=Nothing")
	rt = ResultsTable.getResultsTable()
	return rt


def get_red_spots(rt, slices):
	
	spots = [0]*slices
	for count in range(rt.size()):
		channel = int(rt.getValue("Channel",count))
		thisslice = int(rt.getValue("Slice",count))
		if (channel == 1):
			spots[thisslice-1] += 1
	return spots

def get_green_spots(rt, slices):
	spots = [0]*slices
	for count in range(rt.size()):
		channel = int(rt.getValue("Channel",count))
		thisslice = int(rt.getValue("Slice",count))
		if (channel == 2):
			spots[thisslice-1] += 1
	return spots

def get_colocalised(rt, slices):
	spots = [0]*slices
	for count in range(rt.size()):
		coloc = int(rt.getValue("Colocalized",count))
		thisslice = int(rt.getValue("Slice",count))
		if (coloc == 1):
			spots[thisslice-1] += 1
	return spots


# the the list of file names in the input directory
for root, directories, filenames in os.walk(srcDir):
    filenames.sort();

for filename in filenames:
	if not filename.endswith(ext):
		continue
	# generate full file path for opening
	print(os.path.join(srcDir, filename))

	path = os.path.join(srcDir, filename)

	# use the Bioformats importer to open image
	IJ.run("Bio-Formats Importer", "open=" + path + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	image = IJ.getImage()
	directory = srcDir
	channels = [1, 2]
	channel = 3
	channel = channel % totchannels
	twochannel_stack = retrieve_channels(image, channels)
	dapi_stack = retrieve_dapi(image, channel)
	image.close()
	#image = ImagePlus("test", twochannel_stack)
	#fs = FileSaver(image)
	#filepath = directory + "/" + filename + "_twochannel.tif" 
	#fs.saveAsTiff(filepath) 
	image = ImagePlus("test", twochannel_stack).show()
	
	image = IJ.getImage()
	z_slices = twochannel_stack.getSize() / 2
	print("order=xyczt(default) channels=2 slices="+ str(z_slices) + " frames=1 display=Color")
	IJ.run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices="+ str(z_slices) + " frames=1 display=Color")
	
	rt = run_comdet(image)
	rt.save(directory+"/"+filename+"_results.csv" )
	image.close()
	image = IJ.getImage()
	image.close()
	
	red_spots = get_red_spots(rt, z_slices)
	green_spots = get_green_spots(rt, z_slices)
	
	colocalised = get_colocalised(rt, z_slices)
	for i in range(z_slices):
		colocalised[i] = float(colocalised[i] /2)
	print("red: ",red_spots,"green: ", green_spots, "coloc: ", colocalised)

	
	fp = open(directory+"/"+filename+"_summary.csv", "w")
	fp.write("\z-slice,red spots, green spots, colocalised, percentage of red spots colocalising, percentage of green spots colocalising\n")
	for i in range(z_slices):
		if red_spots[i] == 0:
			if green_spots[i] == 0:
				fp.write(str(i)+","+str(red_spots[i])+","+str(green_spots[i])+","+str(colocalised[i])+","+str(0)+","+str(0)+"\n")
			else:
				fp.write(str(i)+","+str(red_spots[i])+","+str(green_spots[i])+","+str(colocalised[i])+","+str(0)+","+str(colocalised[i]/green_spots[i])+"\n")
		else:
			if green_spots[i] == 0:
				fp.write(str(i)+","+str(red_spots[i])+","+str(green_spots[i])+","+str(colocalised[i])+","+str(colocalised[i]/red_spots[i])+","+str(0)+"\n")
			else:
				fp.write(str(i)+","+str(red_spots[i])+","+str(green_spots[i])+","+str(colocalised[i])+","+str(colocalised[i]/red_spots[i])+","+str(colocalised[i]/green_spots[i])+"\n")
	fp.write("\n\n\n")
	fp.write("total red spots, "+str(sum(red_spots))+"\n")
	fp.write("total green spots, "+str(sum(green_spots))+"\n")
	fp.write("total colocalised, "+str(sum(colocalised))+"\n")
	fp.write("percentage of red spots colocalising, "+str(sum(colocalised)/sum(red_spots))+"\n")
	fp.write("percentage of green spots colocalising, "+str(sum(colocalised)/sum(green_spots))+"\n")
	fp.close()
	
	rt.reset()
	w = WindowManager 

	win = w.getWindow("Summary") 
	win.close() 
	win = w.getWindow("Results") 
	win.close() 
	#rt.close()
	