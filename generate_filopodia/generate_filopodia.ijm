filename = getTitle();
dir = getDirectory("image"); 
run("Auto Threshold", "method=Triangle white stack");
run("Fill Holes", "stack");
run("Duplicate...", "duplicate");
rename("original");
selectWindow(filename);
setOption("BlackBackground", true);
run("Erode", "stack");
run("Erode", "stack");
run("Erode", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
imageCalculator("Subtract create stack", "original",filename);

dot = indexOf(filename, "."); 
if (dot >= 0) filename = substring(filename, 0, dot); 
filepath = dir + "/" + filename+"_filopodia.ome.tiff";
run("Bio-Formats Exporter", "save="+filepath);
run("Close All");

