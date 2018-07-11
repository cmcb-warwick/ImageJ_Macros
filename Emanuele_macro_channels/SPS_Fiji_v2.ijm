/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
output = input
#@ String (label = "File suffix", value = ".tif") suffix


processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	
	for (i = 0; i < list.length; i++) {
		print(list[i]);
		if(File.isDirectory(input + File.separator + list[i])){
			
			processFolder(input + File.separator + list[i]);
		}
		if(endsWith(list[i], suffix)){
			
			output = input;
			
			processFile(input, output, list[i]);
		}
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them -
	// or don't, if you prefer to have an idea of what's happening!
	print("Processing: " + input + File.separator + file);
	run("Bio-Formats Importer", "open=["+ input + File.separator + file+"] color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
	output = input + File.separator + File.nameWithoutExtension;
	print("Making Directory: " + input + File.separator + file + File.separator);
	File.makeDirectory(output);
	while (nImages>0) { 
		  print(nImages);
          selectImage(nImages); 
          print("Selected window " + nImages);
          print("Saving to" + input + File.separator + file + "_c="+(nImages-1) + ".tif");
          saveAs("Tiff", output + File.separator+ file + "_c="+ (nImages-1) + ".tif");
          run("Z Project...", "projection=[Max Intensity] all");
          saveAs("Tiff", output + File.separator+ file + "_c="+ (nImages-2) + "_max.tif");
          close(); 
          selectImage(nImages); 
          close(); 
      } 
	
	
}
