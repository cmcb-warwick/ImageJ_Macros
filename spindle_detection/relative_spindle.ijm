#@ File (label = "Input directory:", style = "directory") input

run("Clear Results");

if (File.exists(input+"/results.csv")){
	File.delete(input+"/results.csv");
}
f = File.open(input+"/results.csv");
print(f, "file,measurement,channel,IntDen,Mean");
processFolder(input);

File.close(f);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], "tif") && lengthOf(list[i])>15)
			processFile(input, list[i]);
	}
}

function processFile(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
	open(input+"/"+file);
	res1 = calculateTotal();
	print(f, file+","+"background,0,"+res1[0]+","+res1[1]);
	print(f, file+","+"background,1,"+res1[2]+","+res1[3]);
	open(input+"/"+file);
	res2 = calculateSpindle(file);
	print(f, file+","+"spindle,0,"+res2[0]+","+res2[1]);
	print(f, file+","+"spindle,1,"+res2[2]+","+res2[3]);

	print(f, file+","+"ratio,0,"+res2[0]/res1[0]+","+res2[1]/res1[1]);
	print(f, file+","+"ratio,1,"+res2[2]/res1[2]+","+res2[3]/res1[3]);
	
	
}













function calculateSpindle(title){
	setTool("multipoint");
	waitForUser("Use the multipoint tool to select 3 points");
	run("Set Measurements...", "area mean standard min fit integrated redirect=None decimal=3");
	
	dir = getDirectory("image"); 
	fitEllipse(title);
	Stack.setPosition(1, 1, 1);
	roiManager("select",0);
	run("Measure");
	results = newArray(4);
	results[0] = getResult("IntDen", 0);
	results[1] = getResult("Mean", 0);
	Stack.setPosition(2, 1, 1);
	//roiManager("select",0);
	run("Measure");
	results[2] = getResult("IntDen", 1);
	results[3] = getResult("Mean", 1);
	if (roiManager("count") > 0){
		roiManager("delete");
	}
	
	saveAs("results", dir+"/"+title+"_spindle.csv");
	run("Close All");
	run("Clear Results");
	return results;
	
}




function fitEllipse(title){
	total = nResults;
	if (roiManager("count") > 0){
		roiManager("delete");
	}
	
	run("Measure");
	x1=getResult("X", 0);
	y1=getResult("Y", 0);
	x2=getResult("X", 1);
	y2=getResult("Y", 1);
	x3=getResult("X", 2);
	y3=getResult("Y", 2);

	//IJ.deleteRows(total - 3, total - 1);
	run("Clear Results");
	run("Select None");
	
	d12 = distance(x1,y1,x2,y2);
	d13 = distance(x1,y1,x3,y3);
	d23 = distance(x2,y2,x3,y3);
	
	getPixelSize(unit,  sizex, sizey);
	//print(d12,d13,d23);
	
	dir = getDirectory("image");
	xcentre = (x1+x2)/2;
	ycentre = (y1+y2)/2;
	
	c = 0.8*d12/2;
	
	a = 0.8*(d13 + d23)/2;
	b = sqrt(a * a + c * c);
	
	xstart = xcentre - a;
	ystart = ycentre - b;
	
	alpha=90+360*angle(x1,x2,y1,y2)/(2*PI);
	
	
	makeOval(xstart/sizex, ystart/sizey, 2*a/sizex, 2*b/sizey);
	//print(alpha);
	run("Rotate...", "  angle="+alpha);
	
	roiManager("Add");
	roiManager("Save", dir+"/"+title+".roi");
	
	run("Clear Results");
	run("Select None");
	

}


function distance(x1,y1,x2,y2){
	xmag = x2-x1;
	ymag = y2-y1;
	prod = xmag * xmag + ymag * ymag;
	return sqrt(prod);
}


function angle(x1,x2,y1,y2){
	diffx = x2 - x1;
	diffy = y2 - y1;
	return atan(diffy/diffx);
}


function calculateTotal(){
	setTool("freehand");
	waitForUser("Use the freehand select tool to select the cell");
	title = getTitle();
	dir = getDirectory("image"); 
	run("Set Measurements...", "area mean standard min fit integrated redirect=None decimal=3");
	run("Create Mask");
	selectWindow("Mask");
	run("Select None");
	
	saveAs("jpeg", dir + "/"+title+"_mask.jpg");
	run("Invert");
	imageCalculator("Transparent-zero create stack", title,"Mask");
	Stack.setPosition(1, 1, 1);
	run("Measure");
	results = newArray(4);
	results[0] = getResult("IntDen", 0);
	results[1] = getResult("Mean", 0);
	
	Stack.setPosition(2, 1, 1);
	run("Measure");
	results[2] = getResult("IntDen", 1);
	results[3] = getResult("Mean", 1);
	saveAs("results", dir+"/"+title+"_background.csv");
	selectWindow("Mask");
	close();
	selectWindow("Result of "+title);
	close();
	run("Select None");
	run("Clear Results");
	return results;
	
}
