// this macro cuts the wells in the videos acquired with the sony camera
roiManager("reset");
image = getImageID();


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
gap_little = 30;
gap_big = 543;
width = 448;
height = 449;
}

// condiciones 
//condiciones = newArray("WT1", "WT2", "WT3", "WT4", "WT5", "WT6");
//condiciones = newArray("7", "8", "9", "10", "11", "12");
//condiciones = newArray("13", "14", "15", "16", "17", "18");
condiciones = newArray("19", "20", "21", "22", "23", "24");

// well 1
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 - height - gap_little, width, height);
roiManager("add");


// well 2
selectImage(image);
makeRectangle(X0 + gap_little, Y0 - height - gap_little, width, height);
roiManager("add");


// well 3
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 - height - gap_little, width, height);
roiManager("add");


// well 4
selectImage(image);
makeRectangle(X0 - width - gap_little, Y0 + gap_little, width, height);
roiManager("add");


// well 5
selectImage(image);
makeRectangle(X0 + gap_little, Y0 + gap_little, width, height);
roiManager("add");


// well 6
selectImage(image);
makeRectangle(X0 + gap_big + gap_little, Y0 + gap_little, width, height);
roiManager("add");
