// this macro cuts the wells in the videos acquired with the sony camera
roiManager("reset");
run("Collect Garbage");
image = getImageID();
dir = getDirectory("Choose a Directory to save all images");

// reference selection at first cross
getSelectionCoordinates(xpoints, ypoints);
X0 = xpoints[0];
Y0 = ypoints[0];

// rectangle size
proxy = false;
if (proxy == true) {
gap_little = 13;
gap_big = 350;
width = 310;
height = 310;
} else {
gap_little = 16;
gap_big = 525;
width = 465;
height = 449;
}


// condiciones 
condiciones = newArray("WT1", "WT2", "WT3", "WT4", "WT5", "WT6");
//condiciones = newArray("WT7", "WT8", "WT9", "WT10", "WT11", "WT12");
//condiciones = newArray("KO44_1", "KO44_2", "KO44_3", "KO44_4", "KO44_5", "KO44_6");
//condiciones = newArray("KO44_7", "KO44_8", "KO44_9", "KO44_10", "KO44_11", "KO44_12");
//condiciones = newArray("KO179_1", "KO179_2", "KO179_3", "KO179_4", "KO179_5", "KO179_6");
//condiciones = newArray("KO179_7", "KO179_8", "KO179_9", "KO179_10", "KO179_11", "KO179_12");
//condiciones = newArray("KO185_1", "KO185_2", "KO185_3", "KO185_4", "KO185_5", "KO185_6");
//condiciones = newArray("KO185_7", "KO185_8", "KO185_9", "KO185_10", "KO185_11", "KO185_12");
//condiciones = newArray("KO44_13", "KO44_14", "KO44_15", "KO44_16", "KO44_17", "KO44_18");
//condiciones = newArray("KO44_19", "KO44_20", "KO44_21", "KO44_22", "KO44_23", "KO44_24");
//condiciones = newArray("KO179_25", "KO179_26", "KO179_27", "KO179_28", "KO179_29", "KO179_30");
//condiciones = newArray("KO179_31", "KO179_32", "KO179_33", "KO179_34", "KO179_35", "KO179_36");

// well 1
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 - height - gap_little, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[0]);
run("8-bit");
saveAs("tiff", dir+condiciones[0]);
close();
run("Collect Garbage");

// well 2
selectImage(image);
makeRectangle(X0 + gap_little+5, Y0 - height - gap_little, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[1]);
run("8-bit");
saveAs("tiff", dir+condiciones[1]);
close();
run("Collect Garbage");

// well 3
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 - height - gap_little -5, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[2]);
run("8-bit");
saveAs("tiff", dir+condiciones[2]);
close();
run("Collect Garbage");

// well 4
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 + gap_little-3, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[3]);
run("8-bit");
saveAs("tiff", dir+condiciones[3]);
close();
run("Collect Garbage");

// well 5
selectImage(image);
makeRectangle(X0 + gap_little, Y0 + gap_little, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[4]);
run("8-bit");
saveAs("tiff", dir+condiciones[4]);
close();
run("Collect Garbage");

// well 6
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 + gap_little+3, width, height);
run("Duplicate...", "duplicate range=500-6100 title="+condiciones[5]);
run("8-bit");
saveAs("tiff", dir+condiciones[5]);
close();
run("Collect Garbage");

print("finito");