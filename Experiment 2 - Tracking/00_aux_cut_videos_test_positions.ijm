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

roiManager("show all without labels");
//roiManager("reset");
