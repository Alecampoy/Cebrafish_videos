getDimensions(width, height, channels, slices, frames);

run("Z Project...", "projection=Median");
newImage("stack_temp", "32-bit black", width, height, frames);	
imageCalculator("Add stack", "stack_temp","MED_Behave 25hpf 4.lif - 4.tif");
rename("proyection_temp");
imageCalculator("Subtract stack", "proyection_temp","Behave 25hpf 4.lif - 4.tif");
run("16-bit");
run("Gaussian Blur...", "sigma=1 stack");
run("Convert to Mask", "method=RenyiEntropy background=Dark calculate black");