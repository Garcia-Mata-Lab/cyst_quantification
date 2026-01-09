// Macro name: Nuclei_Cyst_Volume.ijm
// Quantifies overall cyst volume from Hoechst/DAPI-stained cyst z-stack images
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
		// Pick whichever channel has nuclei staining - C3 is an example 
		selectWindow("C3-Image");
		rename("Image");

// Get the active ROI name
roiManager("select", c);
rName = Roi.getName;

// Save image as a tif for Napari nucleus quantification
selectWindow("Image");
run("Select None");
run("Duplicate...", "title=Image duplicate");
saveAs("Tiff", output + "/" + name + "_" + rName + "Nuclei" + ".tif");
close("Image duplicate");

// Cyst volume quantification with nucleus channel
selectWindow("Image");		
run("Gaussian Blur...", "sigma=2 stack");
// Set threshold as heavy as possible without picking up background
setAutoThreshold("Mean dark no-reset");
run("Threshold...");
waitForUser("set threshold");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Mean background=Dark black");
run("Fill Holes", "stack");
run("3D Objects Counter", "threshold=1 slice=100 min.=100 max.=9000000 objects statistics summary");

// Save the results in output folder as a csv with the image name
selectWindow("Statistics for Image");
saveAs("Results", output + "/"  + name + "_" + rName + "CystVol" + ".csv" );
if (isOpen("Objects map of Image")) close("Objects map of Image");
		if (isOpen("Surface map of Image")) close("Surface map of Image");
		if (isOpen("Text")) close("Text");
		close("Image");
		close("Results");
	}

	close("*");
	if (isOpen("ROI Manager")) close("ROI Manager");
}
