//Macro allows manual extraction of cell shape and measurement of multiple ROIs in a image or movie. The results files have to be saved manually at each step!


macro "Get cell shape [1]" {

run("Make Binary", "method=MaxEntropy background=Dark calculate");
run("Find Edges", "stack");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
TimeString = "";
if (dayOfMonth<10) {TimeString = TimeString+"0";}
TimeString = TimeString+dayOfMonth+"-"+month+"-"+year+"-";
if (hour<10) {TimeString = TimeString+"0";}
TimeString = TimeString+hour;
if (minute<10) {TimeString = TimeString+"0";}
TimeString = TimeString+minute;
if (second<10) {TimeString = TimeString+"0";}
TimeString = TimeString+second;

//By the end of this, the variable TimeString contains current date and time.

//now, we get the directory where the image is stored...
dir = getDirectory("image"); 

//and then save the Results window as a CSV file at the same dir, with the timestring as part of the file name
setOption("BlackBackground", false);
run("Make Binary", "method=MaxEntropy background=Dark calculate");
run("Find Edges", "stack");
run("Save XY Coordinates...", "background=0 suppress process save="+dir+TimeString+".txt");

setForegroundColor(255, 0, 0);





}

macro "Measure shape features [2]" {

nROI = roiManager("count");

for (i=0;i<nROI;i++) {

roiManager("Select", i);
roiManager("Measure");

setForegroundColor(255, 0, 0);
run("Draw", "slice");

 }
