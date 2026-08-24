// Run this macro to automatically change the luts of the dragonfly confocal files to the origenal colors.
// Delete the setMinandMax commands if you do not want to adjust the contrast

run("Make Composite");
//rename("Image");
//run("Z Project...", "projection=[Max Intensity]");
Stack.setChannel(1); run("Magenta"); setMinAndMax(110, 350);
Stack.setChannel(2); run("Blue"); setMinAndMax(104, 1000);
Stack.setChannel(3); run("Yellow"); setMinAndMax(120,650);
Stack.setChannel(4); run("Red");setMinAndMax(110,650);
//run("Split Channels");
//run("Merge Channels...", "c1=C1-MAX_Image c2=C2-MAX_Image c4=C4-MAX_Image create");
//waitForUser("Draw ROI, add to manager, duplicate, copy merge to powerpoint, then click OK");
//run("Split Channels");