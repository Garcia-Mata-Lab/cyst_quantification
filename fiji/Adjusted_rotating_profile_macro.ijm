// Macro name: Rotating_Cyst_Profile
// Radially quantifies the distribution of your protein of interest from the cyst centroid outwards
// Author: Gabriel Kreider and Madeline Lovejoy
// Date: 2026-01-21
// Fiji version: 1.54p

// Bio-Formats version: 8.1.1

// Batch processing - choose input and output directory
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
// Code will only open Imaris files
#@ String (label = "File suffix", value = ".ims") suffix

processFolder(input);

// Scans folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
    
// Open images - series 1 
run("Bio-Formats Importer", "open=[" + input + "/" + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

// Renames the image to match the original file name
name = getInfo("image.filename");
rename(name);

// Select myosin channel
Stack.setChannel(2);

// Manual cyst selection
waitForUser("Set Z slider to the center plane of the cyst, draw ROI around one cyst, add to manager, then click OK");

// Goes through each cyst ROI
for (c=0 ; c<roiManager("count"); c++) {
		selectWindow(name);
		roiManager("select", c);
	
run("ROI Manager...", "close");	

// Duplicate only the current slice and channel
currentSlice = getSliceNumber();
// Duplicate for radial intensity measurements
run("Duplicate...", "title=Image channels=2 slices=" + currentSlice);
// Duplicate for calculating cyst diameter
run("Duplicate...", "title=[Image duplicate]");
  
// Cyst diameter calculation
run("Select None");
run("Gaussian Blur...", "sigma=2 stack");
setAutoThreshold("Mean dark no-reset");
run("Threshold...");
waitForUser("Set threshold as heavy as possible, then click OK");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Mean background=Dark black");
run("Fill Holes", "stack");

// Sets measurements to include Feret's diameter value (longest diameter that exists in the cyst)
run("Set Measurements...", "area mean centroid feret's redirect=None decimal=3");
run("Analyze Particles...", "display clear");

// Assigns Feret's diameter and centroid coordinates (µm) as variables
feret_um = getResult("Feret", 0);
x_um = getResult("X", 0);
y_um = getResult("Y", 0);

// Points and lines can only be made with pixel units
// Get pixel size calibration values
getPixelSize(unit, pixelWidth, pixelHeight, voxelDepth);

// Convert from micrometers to pixel units
x_px = x_um / pixelWidth;
y_px = y_um / pixelHeight;
feret_px = feret_um / pixelWidth;  // assumes square pixels

// Select the image used to make the point and line ROIs
selectWindow("Image duplicate");

// Create a point ROI at centroid in pixel coordinates
makePoint(x_px, y_px);
roiManager("Add");

// Create a horizontal line ROI spanning Feret’s diameter
x1 = x_px - (feret_px / 2);
x2 = x_px + (feret_px / 2);
makeLine(x1, y_px, x2, y_px);
roiManager("Add");

// Select the image used to measure the radial intensity
selectWindow("Image");
run("Clear Results");

// Select the line ROI
roiManager("Select", roiManager("count") - 1);

// Take measurements along the line every degree for 360 degrees
 for (c=0; c<360; c++) {
	 profile = getProfile();
	  for (i=0; i<profile.length; i++)
	      setResult(c, i, profile[i]);
	  updateResults;
	  run("Rotate...", "  angle=1");
 }

// Save the intensity measurements
iname = getInfo("image.filename");
dir = getDir("image");
selectWindow("Results");
saveAs("Measurements", dir + name + "_" + iname + "Myo.csv");

// Clears ROI manager
n = roiManager("count");
if (n > 0) {
    roiManager("select", n - 1);
    roiManager("Delete"); 
}
}
}
