"""

distance_spots.py 
created by: Erick Martins Ratamero
date: 29/06/18

From an open image in Fiji with a point ROI,
generate a thresholded stack of channel 3,
detect 3D particles and calculate their
distance to the original point ROI. Finally,
it saves the results as a csv file.

"""




import os
from java.io import File

from ij import IJ, ImageStack, ImagePlus
from ij.plugin.frame import RoiManager
import math
from ij import WindowManager



# returns an instance of ROI Manager (creates one if 
# it doesn't exist)
def get_roi_manager(new=False):
    rm = RoiManager.getInstance()
    if not rm:
        rm = RoiManager()
    if new:
        rm.runCommand("Reset")
    return rm




# gets open image, filename and directory
image = IJ.getImage()
directory = image.getOriginalFileInfo().directory
filename = image.getOriginalFileInfo().fileName
#print(image.getOriginalFileInfo())

# retrieves pixel size for image
cal = image.getCalibration()

size_x = cal.pixelWidth
size_y = cal.pixelHeight
size_z = cal.pixelDepth

# adds point ROI to ROI Manager, get the ROI content
rm = get_roi_manager()
rm.runCommand("Add")
rois = rm.getRoisAsArray()

# just in case, loops over existing ROIs
# and assigns the coordinates of the point ROI
# to variables (taking pixel size into account)

for roi in rois:
	
	z_start = roi.getZPosition()
	z_start = z_start * size_z
	x_start = roi.getXBase()
	x_start = x_start * size_x
	y_start = roi.getYBase()
	y_start = y_start * size_y
	print("centre: ",z_start, x_start, y_start)

	# finally, resets the ROI manager
	rm.runCommand("Reset")


# get stack from current image
stack = image.getStack()
# create empty stack for split channel
red_stack = ImageStack(image.width, image.height)

# now we go through the original image and retrieve each channel 3 slice
for i in range(1, image.getNSlices()+1):
	slice = stack.getProcessor(i*4 -1)

	# then, we assign those slices to the new stack
	red_stack.addSlice(str(i), slice)
# by the end of this for loop, red_stack contains channel 3 only

# we create a new image from that stack and display it
ImagePlus("stack", red_stack).show()	
# then we close the original image
image.close()
		
		
# from now on, "image" refers to the channel 3 only stack
image = IJ.getImage()
# before thresholding, we need to invert it
IJ.run("Invert", "stack")
# then, we create a binary image using MaxEntropy thresholding
IJ.run("Make Binary", "method=MaxEntropy background=Dark");
# now we use the 3D Objects Counter to do 3D particle analysis
IJ.run("3D Objects Counter", "threshold=128 slice=32 min.=10 max.=16777216 exclude_objects_on_edges statistics")
# The results window is called "Statistics for stack". We assign that to rt
rt = WindowManager.getWindow("Statistics for stack").getTextPanel().getOrCreateResultsTable()

# create a new "dist" column in the results table
rt.addValue("dist", 0)
for row in range(rt.size()):

	# for each result, we get XYZ coordinates and then calculate their 
	# "real" value given pixel size
	x = rt.getValue("X", row)
	x = x * size_x
	y = rt.getValue("Y", row)
	y = y * size_y
	z = rt.getValue("Z", row)
	z = z * size_z

	# now, we calculate the "real" distance between this particle
	# and the original point ROI
	dist = math.sqrt( (x - x_start)**2 + (y - y_start)**2 + (z - z_start)**2 )
	# finally, we assign the result to the "dist" column
	rt.setValue("dist", row, dist)

# close the channel 3 stack window
image.close()

# save results table to a csv file
rt.save(directory+filename+".csv" )

# close remaining windows
WindowManager.getWindow("Statistics for stack").close()
WindowManager.getWindow("ROI Manager").close()
 