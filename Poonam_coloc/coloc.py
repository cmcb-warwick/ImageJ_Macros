#@ File    (label = "Input directory", style = "directory") srcFile
#@ String  (label = "File extension", value=".tif") ext
#@ Float (label = "Channel 1 approx size", value = 4.0) ch1size
#@ Float (label = "Channel 1 threshold", value = 3.0) ch1thresh
#@ Float (label = "Channel 2 approx size", value = 4.0) ch2size
#@ Float (label = "Channel 2 threshold", value = 3.0) ch2thresh
#@ Float (label = "Maximum colocalisation distance", value = 4.0) coloc
#@ Boolean (label = "Use auto local threshold", value = False) auto_thresh

"""

coloc.py 
created by: Erick Martins Ratamero
date: 01/08/18
last updated: 02/08/18

Opens a directory of images and calculates green on red and 
red on green colocalisation. Separates bottom of nucleus calculations
from the rest of the stack.

"""
import os
from java.io import File

from ij import IJ, ImageStack, ImagePlus, CompositeImage
from ij.plugin.frame import RoiManager
import math
from ij import WindowManager
from ij.measure import ResultsTable
from loci.plugins import BF
from ij.io import FileSaver 
from ij.process import ImageStatistics as IS  
import time
from ij.plugin import HyperStackConverter
from ij.process import ImageConverter

srcDir = srcFile.getAbsolutePath()

totchannels = 3


def safe_div(x,y):
    if y == 0:
        return 0
    return x / y

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
	IJ.run(image,"Detect Particles", " two=[Detect in both channels independently] ch1a="+str(ch1size) + " ch1s=" + str(ch1thresh) + " ch2a="+str(ch2size) + " ch2s=" + str(ch2thresh) + " calculate max=" + str(coloc) + " add=Nothing")
	rt = ResultsTable.getResultsTable()
	return rt


def get_red_spots(rt, slices, image):
	
	spots = [0]*slices
	dapi_spots = [0]*slices
	for count in range(rt.size()):
		channel = int(rt.getValue("Channel",count))
		thisslice = int(rt.getValue("Slice",count))
		X = int(rt.getValue("X_(px)", count))
		Y = int(rt.getValue("Y_(px)", count))
		if (channel == 1):
			spots[thisslice-1] += 1
			if(image.getPixel(X, Y)[0] == 255):
				dapi_spots[thisslice-1] += 1
	return [spots, dapi_spots]

def get_green_spots(rt, slices, image):
	spots = [0]*slices
	dapi_spots = [0]*slices
	for count in range(rt.size()):
		channel = int(rt.getValue("Channel",count))
		thisslice = int(rt.getValue("Slice",count))
		X = int(rt.getValue("X_(px)", count))
		Y = int(rt.getValue("Y_(px)", count))
		if (channel == 2):
			spots[thisslice-1] += 1
			if(image.getPixel(X, Y)[0] == 255):
				dapi_spots[thisslice-1] += 1
	return [spots, dapi_spots]

def get_colocalised(rt, slices, image):
	spots = [0]*slices
	dapi_spots = [0]*slices
	for count in range(rt.size()):
		
		coloc = int(rt.getValue("Colocalized",count))
		thisslice = int(rt.getValue("Slice",count))
		X = int(rt.getValue("X_(px)", count))
		Y = int(rt.getValue("Y_(px)", count))
		#print(count, X, Y, coloc, thisslice, spots, dapi_spots)
		if (coloc == 1):
			spots[thisslice-1] += 1
			if(image.getPixel(X, Y)[0] == 255):
				dapi_spots[thisslice-1] += 1
	return [spots, dapi_spots]


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
	image_dapi = ImagePlus("dapi stack", dapi_stack)
	image_dapi.show()
	image = ImagePlus("two channel stack", twochannel_stack)
	image.show()
	if auto_thresh:
		con = ImageConverter(image)
		con.convertToGray8()
		IJ.run(image, "Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white stack")
		#image = CompositeImage(image_two)
	#image = IJ.getImage()
	z_slices = image.getDimensions()[3] / 2
	
	print("order=xyczt(default) channels=2 slices="+ str(z_slices) + " frames=1 display=Color", image.getDimensions())
	image_two = HyperStackConverter.toHyperStack(image,2,z_slices,1)
	image = CompositeImage(image_two)
	image.show()
	
	rt = run_comdet(image)
	image = IJ.getImage()
	
	rt.save(directory+"/"+filename+"_results.csv" )
	
	image.setDimensions(2, z_slices, 1)
	image.setOpenAsHyperStack(True)
	print(image.isHyperStack(), image.getNChannels(), image.getOverlay())
	#image.flattenStack()
	image.show()
	fs = FileSaver(image)
	filepath = directory + "/" + filename + "_coloc.tiff" 
	
	fs.saveAsTiff(filepath) 
	image.close()

	image = IJ.getImage()
	IJ.run(image_dapi,"Convert to Mask", "method=Default background=Default calculate")

	
	
	#image = IJ.getImage()
	#image.close()
	
	[red_spots, red_spots_dapi] = get_red_spots(rt, z_slices, image_dapi)
	
	[green_spots, green_spots_dapi] = get_green_spots(rt, z_slices, image_dapi)
	
	[colocalised, colocalised_dapi] = get_colocalised(rt, z_slices, image_dapi)
	for i in range(z_slices):
		colocalised[i] = float(colocalised[i] /2)
	for i in range(z_slices):
		colocalised_dapi[i] = float(colocalised_dapi[i] /2)
	print("dapi red: ",red_spots_dapi,"dapi green: ", green_spots_dapi, "dapi coloc: ", colocalised_dapi)
	image.close()
	
	fp = open(directory+"/"+filename+"_summary.csv", "w")
	fp.write("red spots, green spots, colocalised, percentage of red spots colocalising, percentage of green spots colocalising, red spots on DAPI, green spots on DAPI, colocalised on DAPI, percentage of red spots colocalising on DAPI, percentage of green spots colocalising on DAPI\n")
	for i in range(z_slices):
		fp.write(str(i)+","+str(red_spots[i])+","+str(green_spots[i])+","+str(colocalised[i])+","+str(safe_div(colocalised[i],red_spots[i]))+","+str(safe_div(colocalised[i],green_spots[i])))
		fp.write(","+str(red_spots_dapi[i])+","+str(green_spots_dapi[i])+","+str(colocalised_dapi[i])+","+str(safe_div(colocalised_dapi[i],red_spots_dapi[i]))+","+str(safe_div(colocalised_dapi[i],green_spots_dapi[i]))+"\n")
	fp.write("\n\n\n")
	fp.write("total red spots, "+str(sum(red_spots))+"\n")
	fp.write("total green spots, "+str(sum(green_spots))+"\n")
	fp.write("total colocalised, "+str(sum(colocalised))+"\n")
	fp.write("percentage of red spots colocalising, "+str(sum(colocalised)/sum(red_spots))+"\n")
	fp.write("percentage of green spots colocalising, "+str(sum(colocalised)/sum(green_spots))+"\n\n\n")

	fp.write("total red spots on DAPI, "+str(sum(red_spots_dapi))+"\n")
	fp.write("total green spots on DAPI, "+str(sum(green_spots_dapi))+"\n")
	fp.write("total colocalised on DAPI, "+str(sum(colocalised_dapi))+"\n")
	fp.write("percentage of red spots colocalising on DAPI, "+str(safe_div(sum(colocalised_dapi),sum(red_spots_dapi)))+"\n")
	fp.write("percentage of green spots colocalising on DAPI, "+str(safe_div(sum(colocalised_dapi),sum(green_spots_dapi)))+"\n")
	fp.close()

	
	fp = open(directory+"/"+filename+"_bottom.csv", "w")
	fp.write("red spots, green spots, colocalised, percentage of red spots colocalising, percentage of green spots colocalising, red spots on DAPI, green spots on DAPI, colocalised on DAPI, percentage of red spots colocalising on DAPI, percentage of green spots colocalising on DAPI\n")
	fp.write(str(red_spots[-1])+","+str(green_spots[-1])+","+str(colocalised[-1])+","+str(safe_div(colocalised[-1],red_spots[-1]))+","+str(safe_div(colocalised[-1],green_spots[-1])))
	fp.write(","+str(red_spots_dapi[-1])+","+str(green_spots_dapi[-1])+","+str(colocalised_dapi[-1])+","+str(safe_div(colocalised_dapi[-1],red_spots_dapi[-1]))+","+str(safe_div(colocalised_dapi[-1],green_spots_dapi[-1]))+"\n")
	fp.close()
	rt.reset()
	image_dapi.close()
	w = WindowManager 

	win = w.getWindow("Summary") 
	win.close() 
	win = w.getWindow("Results") 
	win.close() 
	#rt.close()
	