// this macro cuts the wells in the videos acquired with the sony camera
// sizes are adjusted to the proxy smaller videos

image = getImageID();
dir = getDirectory("Choose a Directory to save all images");

// reference selection at first cross
getSelectionCoordinates(xpoints, ypoints);
X0 =  xpoints[0];
Y0 = ypoints[0];

// rectangle size
width = 310;
height = 310;

// condiciones 
condiciones = newArray("WT1", "WT2", "WT3", "WT4", "WT5", "WT6");

// well 1
selectImage(image);
makeRectangle(X0 - width - 12, Y0 - height - 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[0]);
run("8-bit");
saveAs("tiff", dir+condiciones[0]);
close();
run("Collect Garbage");

// well 2
selectImage(image);
makeRectangle(X0 + 12, Y0 - height - 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[1]);
run("8-bit");
saveAs("tiff", dir+condiciones[1]);
close();
run("Collect Garbage");


// well 3
selectImage(image);
makeRectangle(X0 + 346 + 12, Y0 - height - 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[2]);
run("8-bit");
saveAs("tiff", dir+condiciones[2]);
close();
run("Collect Garbage");


// well 4
selectImage(image);
makeRectangle(X0 - width - 12, Y0 + 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[3]);
run("8-bit");
saveAs("tiff", dir+condiciones[3]);
close();
run("Collect Garbage");


// well 5
selectImage(image);
makeRectangle(X0 + 12, Y0 + 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[4]);
run("8-bit");
saveAs("tiff", dir+condiciones[4]);
close();
run("Collect Garbage");

// well 6
selectImage(image);
makeRectangle(X0 + 346 + 12, Y0 + 12, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[5]);
run("8-bit");
saveAs("tiff", dir+condiciones[5]);
close();
run("Collect Garbage");
