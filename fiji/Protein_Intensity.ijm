// Macro name: Protein_Intensity.ijm
// Quantifies intensity of your protein of interest in every z-plane
// Author: Gabriel Kreider and Madeline Lovejoy
// Date: 2026-01-10
// Fiji version: 1.54p

// Bio-Formats version: 8.1.1

// When first running the code, the steps after the batch processing should be executed manually so the thresholding values can be adjusted per experiment

// Batch processing - choose input and output directory
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
// Code will only open Imaris files
#@ String (label = "File suffix", value = ".ims") suffix

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

// Renames the image to match the original file name
name = getInfo("image.filename");
rename(name);

// Manual cyst selection
waitForUser("draw ROI around each cyst, add to manager, then click OK");

// Goes through each cyst ROI
for (c=0 ; c<roiManager("count"); c++) {
		selectWindow(name);
		roiManager("select", c);

		// Duplicates the cyst ROI and splits channels into individual images
		run("Duplicate...", "title=Image duplicate");
		run("Split Channels");
		// The user should change "C3" (channel 3) to whichever channel has staining for the protein of interest 
		selectWindow("C3-Image");
		rename("Image2");

// Get the active ROI name
roiManager("select", c);
rName = Roi.getName;

// Protein intensity quantification
selectWindow("Image2");
// Background subtraction and intensity calculations require the image to be 32-bit
run("32-bit");

// Subtract background only get quantification readings from the protein of interest
// Background ROI will be applied to all z-slices

	waitForUser("draw bkgd ROI then hit enter");
	
      for (i=1; i<=nSlices; i++) {
          setSlice(i);
          // Measures the mean background intensity and subtracts it from the entire slice
          getStatistics(area, mean);
          run("Select None");
          run("Subtract...", "value="+mean);
          run("Restore Selection");
      }
     
run("Select None");

run("Duplicate...", "title=[Binary Mask] duplicate");

// Smooths noise before thresholding
run("Gaussian Blur...", "sigma=2 stack");

setAutoThreshold("Default dark");
// Sets a threshold to determine where signal exists
// These values can stay the same across all experiments and conditions
setThreshold(10, 65535);

setOption("BlackBackground", true);
// Signal = 255, background = 0
run("Convert to Mask", "method=Default background=Dark black");

// Converts mask to signal = 1, background = 0
run("Divide...", "value=255 stack");
// Combines mask and original image
// Keeps the original signal values only where the mask exists
imageCalculator("Multiply create 32-bit stack", "Binary Mask","Image2");


setAutoThreshold("Default dark no-reset");
// These values will need to be adjusted per experiment per protein of interest
// Final threshold removes any residual background or autofluorescence within the mask
setThreshold(103.88, 1e30);

setOption("BlackBackground", true);

//Converts pixels that do not have signal to NaN to exclude from measurements
run("NaN Background", "stack");

// Sets any units to pixels so intensity measurements are consistent across all images
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

// Setting measurements to measure the sum integrated density (intensity) of each slice in the z-stack
run("Set Measurements...", "area mean integrated redirect=None decimal=3");	
run("Measure Stack...");


// Save the results in output folder as a csv with the image name
// The user should replace "EcadInt" with whatever the protein of interest is
saveAs("Results", output + "/"  + name + "_" + rName + "EcadInt" + ".csv" );

// Close everything
if (isOpen("Result of Binary Mask")) close("Result of Binary Mask");
		if (isOpen("Binary Mask")) close("Binary Mask");
		close("Image2");
		close("Results");
	}

	close("*");
	if (isOpen("ROI Manager")) close("ROI Manager");
}