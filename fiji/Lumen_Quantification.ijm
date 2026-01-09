// Macro name: Lumen_Quantification.ijm
// Quantifies volume and surface area for individual lumens
// Author: Madeline Lovejoy
// Date: 2026-01-09
// Fiji version: 1.54p

// Bio-Formats version: 8.1.1
// 3D Objects Counter version: 2.0.1

// Batch processing - choose input and output directory
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
// Code will only open Imaris files
#@ String (label = "File suffix", value = ".ims") suffix

// Opens files from input folder
processFolder(input);

// Scans folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
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

// Get image name
name = getInfo("image.filename");
rename(name);

waitForUser("draw ROI around each cyst, add to manager, then click OK");

// Goes through each ROI
for (c=0 ; c<roiManager("count"); c++) {
		selectWindow(name);
		roiManager("select", c);

		run("Duplicate...", "title=Image duplicate");
		run("Split Channels");
		// Pick whichever channel has gp135/apical staining - C2 is an example 
		selectWindow("C2-Image");
		rename("Image4");

// Get the active ROI name
roiManager("select", c);
rName = Roi.getName;

// Lumen Quantification with gp135 channel
selectWindow("Image4");
rename("Image");
run("Gaussian Blur...", "sigma=2 stack");
// Set threshold only heavy enough to make sure that the lumen has a continuous perimeter
setAutoThreshold("Mean dark no-reset");
run("Threshold...");
waitForUser("set threshold");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Mean background=Dark black");
run("Fill Holes", "stack");
run("3D Objects Counter");

// Save the results in output folder as a csv with the image name
selectWindow("Statistics for Image");
saveAs("Results", output + "/"  + name + "_" + rName + "LumenQuant" + ".csv" );

// Close everything
if (isOpen("Objects map of Image")) close("Objects map of Image");
		if (isOpen("Surface map of Image")) close("Surface map of Image");
		if (isOpen("Text")) close("Text");
		close("Image");
		close("Results");
	}

	close("*");
	if (isOpen("ROI Manager")) close("ROI Manager");
}