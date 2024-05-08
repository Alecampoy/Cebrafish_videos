// this macro cuts the wells in the videos acquired with the sony camera

image = getImageID();
dir = getDirectory("Choose a Directory to save all images");

// reference selection at first cross
getSelectionCoordinates(xpoints, ypoints);
X0 = xpoints[0];
Y0 = ypoints[0];

// rectangle size
proxy = false;
if (proxy == true) {
gap_little = 12;
gap_big = 346;
width = 310;
height = 310;
} else {
gap_little = 36;
gap_big = 523;
width = 451;
height = 451;
}

// condiciones 
//condiciones = newArray("WT1", "WT2", "WT3", "WT4", "WT5", "WT6");
//condiciones = newArray("7", "8", "9", "10", "11", "12");
//condiciones = newArray("13", "14", "15", "16", "17", "18");
condiciones = newArray("19", "20", "21", "22", "23", "24");

// well 1
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 - height - gap_little, width, height);
run("Duplicate...", "duplicate range=200-5800 title="+condiciones[0]);
run("8-bit");
saveAs("tiff", dir+condiciones[0]);
close();
run("Collect Garbage");

// well 2
selectImage(image);
makeRectangle(X0 + gap_little, Y0 - height - gap_little, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[1]);
run("8-bit");
saveAs("tiff", dir+condiciones[1]);
close();
run("Collect Garbage");

// well 3
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 - height - gap_little, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[2]);
run("8-bit");
saveAs("tiff", dir+condiciones[2]);
close();
run("Collect Garbage");

// well 4
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 + gap_little, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[3]);
run("8-bit");
saveAs("tiff", dir+condiciones[3]);
close();
run("Collect Garbage");

// well 5
selectImage(image);
makeRectangle(X0 + gap_little, Y0 + gap_little, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[4]);
run("8-bit");
saveAs("tiff", dir+condiciones[4]);
close();
run("Collect Garbage");

// well 6
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 + gap_little, width, height);
run("Duplicate...", "duplicate range=100-5900 title="+condiciones[5]);
run("8-bit");
saveAs("tiff", dir+condiciones[5]);
close();
run("Collect Garbage");

print("finito");