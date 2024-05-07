// reference selection at first cross
getSelectionCoordinates(xpoints, ypoints);
X0 =  xpoints[0];
Y0 = ypoints[0];

// rectangle size
width = 310;
height = 310;

// condiciones 
condiciones = newArray("WT", "WT", "WT", "WT", "WT", "WT");

// well 1
makeRectangle(X0 - width - 12, Y0 - height - 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");

// well 2
makeRectangle(X0 + 12, Y0 - height - 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");

// well 3
makeRectangle(X0 + 346 + 12, Y0 - height - 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");

// well 4
makeRectangle(X0 - width - 12, Y0 + 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");

// well 2
makeRectangle(X0 + 12, Y0 + 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");

// well 3
makeRectangle(X0 + 346 + 12, Y0 + 12, width, height);
run("Duplicate...", "title="+condiciones[0]);
run("8-bit");