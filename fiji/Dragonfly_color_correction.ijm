// Macro name: Dragonfly_color_correction.ijm
// Run this macro to change the LUTs of the Dragonfly confocal files to their orignal colors
// Author: Dr. Gabriel Kreider and Dr. Madeline Lovejoy
// Date: 2026-08-24
// Fiji version: 1.54p

run("Make Composite");
// ***FIRST ADJUST THE CHANNEL VALUES WITH A CTRL CYST IN IMARIS, THEN ADJUST THE COLOR NAMES AND MIN/MAX VALUES IN THIS CODE FOR THE ORDER THEY APPEAR***
Stack.setChannel(1); run("Magenta"); setMinAndMax(110, 350);
Stack.setChannel(2); run("Blue"); setMinAndMax(104, 1000);
Stack.setChannel(3); run("Yellow"); setMinAndMax(120,650);
Stack.setChannel(4); run("Red");setMinAndMax(110,650);
// ***END CHANGE***
