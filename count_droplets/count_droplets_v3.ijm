



/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
sameoutput = 0;
if (input == output){
	sameoutput = 1; 
}
suffix = "TxRed"+suffix;
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	print("Directory: "+input);
	
	print("suffix: "+suffix);
	for (i = 0; i < list.length; i++) {
		print("checking file "+list[i]);
		if (sameoutput){
			output = input;
		}
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	
	redfile = file;
	rootfile = split(file,"_");
	rootfile = rootfile[0];
	greenfile = rootfile + "_GFP.tif";
	print("Processing: " + input + File.separator + redfile);
	print("Processing: " + input + File.separator + greenfile);
	
	open(input + File.separator + redfile);
	open(input + File.separator + greenfile);
	print("Saving to: " + output);
	selectWindow(greenfile);
	//calculate total area of GFP - Moments + fill holes
	run("Duplicate...", " ");
	newgreen = getTitle();
	run("Auto Threshold", "method=Moments white");
	run("Fill Holes");
	run("Analyze Particles...", "size=100-Infinity show=Masks summarize");
	cellmask = getTitle();
	print(cellmask);
	selectWindow("Summary"); 
	saveAs("Results", output + "/" + rootfile + "_allareas.csv");
	
   	selectWindow(newgreen);
   	close();
   	selectWindow(rootfile + "_allareas.csv"); 
   	run("Close"); 
	imageCalculator("AND create", cellmask,greenfile);
	//run("Invert");
	cellmask = getTitle();
	print(cellmask);
   	

	selectWindow(redfile);
	run("Auto Threshold", "method=RenyiEntropy white");
	run("Analyze Particles...", "size=500-infinity show=Masks clear summarize");
	selectWindow("Summary"); 
	saveAs("Results", output + "/" + rootfile + "_redareas.csv");
	run("Close"); 

	
	selectWindow(cellmask);
	run("8-bit");
	run("Auto Local Threshold", "method=Bernsen radius=15 parameter_1=0 parameter_2=0 white");
	run("Analyze Particles...", "size=10-100 show=Outlines display clear summarize");
	saveAs("Tiff", output + "/" + rootfile + "_GFP_outlines.tif");
	saveAs("Results", output + "/" + rootfile + "_alldroplets.csv");
	selectWindow("Results"); 
   	run("Close");
	
	imageCalculator("AND create", "Mask of "+redfile,cellmask);
	run("Invert");
	run("Analyze Particles...", "size=10-100 show=Outlines display clear summarize");
	saveAs("Tiff", output + "/" + rootfile + "_TxRed_outlines.tif");
	saveAs("Results", output + "/" + rootfile + "_reddroplets.csv");
	close();
	run("Close All");
	selectWindow("Summary"); 
   	run("Close"); 
   	selectWindow("Results"); 
   	run("Close"); 


	
}
