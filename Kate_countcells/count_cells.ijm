/* Macro for counting cells and green channel pixels
 * v0.0.3
 * 
 * Author: Erick Martins Ratamero
 * Email: E.Martins-Ratamero.1@warwick.ac.uk
 * 
 * www.warwick.ac.uk/camdu
 */
 
//#@ File (label = "Input directory", style = "directory") input
//#@ File (label = "Output directory", style = "directory") output
//#@ String (label = "File suffix", value = ".tif") suffix

 
Dialog.create("Test");

input = getDirectory("Select an input directory");
output = getDirectory("Select an output directory");
suffix = ".tif";
//print(input);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	if(File.exists(output + File.separator + "Summary_bacteria.csv")){
		 File.delete(output + File.separator +"Summary_bacteria.csv");
		 
	}
	bact = File.open(output + File.separator + "Summary_bacteria.csv") ;
	for (i = 0; i < list.length; i++) {
		//if(File.isDirectory(input + File.separator + list[i]))
			//processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i], bact);
	
	}
	File.close(bact);
}
function processFile(input, output, file, bact) {
	print("Processing: " + input + File.separator + file);
	open(input + File.separator + file);

	run("Split Channels");
	selectWindow(file + " (green)");
	run("Auto Threshold", "method=Otsu");
	setOption("BlackBackground", false);
	run("Make Binary");
	run("Create Selection");
	run("Measure");
	print(bact, file + ", " + getResult("Area"));
	close();
	selectWindow(file + " (red)");
	close();
	selectWindow(file + " (blue)");
	run("Median...", "radius=5");
	run("Auto Threshold", "method=Triangle");
	run("Watershed");
	run("Analyze Particles...", "size=2000-Infinity show=Outlines display clear summarize");
	saveAs("Tiff", output + File.separator + File.nameWithoutExtension + "cells.tif");
	saveAs("Results", output + File.separator + File.nameWithoutExtension + "cells.csv");
	close();
	print("Saving to: " + output);
	close();
}
run("Set Measurements...", "area centroid shape redirect=None decimal=3");
processFolder(input);
selectWindow("Summary"); 
saveAs("Text", output + File.separator +"Summary_cells.csv"); 




