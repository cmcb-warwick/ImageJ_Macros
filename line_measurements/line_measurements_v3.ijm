#@ File (label = "Input directory:", style = "directory") input
#@Integer(label="Tubulin Channel:",value=2) tubulin_channel
#@Integer(label="Protein Channel:",value=3) protein_channel

run("Clear Results");


processFolder(input);



function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], "tif")){
			fname = list[i];
			fname = replace(list[i], ".tif", "");
			if (File.exists(input+"/"+fname+".csv")){
				File.delete(input+"/"+fname+".csv");
			}
			f = File.open(input+"/"+fname+".csv");
			print(f, "measurement,line,channel,Length,Mean,StdDev,Min,Max,BX,BY,Width,Height,IntDen,RawIntDen");
			processFile(input, list[i]);
			File.close(f);
		}
	}
}

function processFile(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.

	

	//Draw 3 lines, take measurement on red and green channel, adjust for background

	//Move the 3 lines, take measurement on red and green channel, adjust for background

	open(input+"/"+file);
	bckgrnd = calculateBackground();
	//print(bckgrnd[0], bckgrnd[1]);

	
	first = initLines(f, bckgrnd);

	
	second = moveLines(f, bckgrnd);



	
	run("Close All");

	
	
	
	
}






function calculateBackground(){
	setTool("rectangle");
	Stack.setChannel(protein_channel);
	run("Set Measurements...", "area mean standard min fit integrated redirect=None decimal=3");
	waitForUser("Select a background region:");
	Stack.setChannel(tubulin_channel);
	run("Measure");
	results = newArray(2);
	results[0] = getResult("Mean", 0);
	Stack.setChannel(protein_channel);
	run("Measure");
	results[1] = getResult("Mean", 1);
	run("Clear Results");
	return results;
}



function initLines(f, background){
	setTool("line");
	run("Set Measurements...", "area mean standard min bounding integrated redirect=None decimal=3");
	Stack.setChannel(tubulin_channel);
	waitForUser("Draw line 1:");
	roiManager("Add");
	waitForUser("Draw line 2:");
	roiManager("Add");
	waitForUser("Draw line 3:");
	roiManager("Add");
	roiManager("multi-measure append");
	results = newArray(6);
	for (i=0;i<3;i++){
		
		print(f, "microtubule,"+(i+1)+","+tubulin_channel+","+getResult("Length", i)+","+getResult("Mean", i)-background[0]+","+getResult("StdDev", i)+","+(getResult("Min", i)-background[0])+","+(getResult("Max", i)-background[0])+","+getResult("BX", i)+","+getResult("BY", i)+","+getResult("Width", i)+","+getResult("Height", i)+","+getResult("IntDen", i)+","+getResult("RawIntDen", i));
	}
	
	Stack.setChannel(protein_channel);
	roiManager("multi-measure append");
	for (i=3;i<6;i++){
		print(f, "microtubule,"+(i-2)+","+protein_channel+","+getResult("Length", i)+","+getResult("Mean", i)-background[1]+","+getResult("StdDev", i)+","+(getResult("Min", i)-background[1])+","+(getResult("Max", i)-background[1])+","+getResult("BX", i)+","+getResult("BY", i)+","+getResult("Width", i)+","+getResult("Height", i)+","+getResult("IntDen", i)+","+getResult("RawIntDen", i));
		
	}
	run("Clear Results");
	
	return results;
}


function moveLines(f, background){
	results = newArray(6);
	Stack.setChannel(tubulin_channel);
	total = roiManager("count");
	for (i=0;i<total;i++){
		
		roiManager("Select", i);
		count = i+1;
		Stack.setChannel(tubulin_channel);
		waitForUser("Move line "+count+" :");
		roiManager("Add");
	}
	for (i=0;i<total;i++){
		roiManager("Select", 0);
		roiManager("Delete");
	}
	roiManager("multi-measure append");
	for (i=0;i<3;i++){
		print(f, "background,"+(i+1)+","+tubulin_channel+","+getResult("Length", i)+","+getResult("Mean", i)-background[0]+","+getResult("StdDev", i)+","+(getResult("Min", i)-background[0])+","+(getResult("Max", i)-background[0])+","+getResult("BX", i)+","+getResult("BY", i)+","+getResult("Width", i)+","+getResult("Height", i)+","+getResult("IntDen", i)+","+getResult("RawIntDen", i));
	}
	Stack.setChannel(protein_channel);
	roiManager("multi-measure append");
	for (i=3;i<6;i++){
		print(f, "background,"+(i-2)+","+protein_channel+","+getResult("Length", i)+","+getResult("Mean", i)-background[1]+","+getResult("StdDev", i)+","+(getResult("Min", i)-background[1])+","+(getResult("Max", i)-background[1])+","+getResult("BX", i)+","+getResult("BY", i)+","+getResult("Width", i)+","+getResult("Height", i)+","+getResult("IntDen", i)+","+getResult("RawIntDen", i));
	}
	run("Clear Results");
	roiManager("reset");
	return results;
}


