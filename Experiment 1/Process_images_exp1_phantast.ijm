///////////////////////////////////////////////////////////////////////////////////////////////
/* Author: Ale Campoy
 * Microscopy Unit (CABD)
 * Date: 10.05.2022
 * User: Marta Fernandez
 * 	
 * Description: Analysis of moving Cebrafish
 * 
 * Input: Folder with the set of bright field images .tiff coming from Leica magnifier
 *
 * Method: Segmentation of the worm using Phantas.
 * 
 * Output: Segmented images, skeleton and results file
 * 
 *///////////////////////////////////////////////////////////////////////////////////////////////

// to do: 
// debug que la ventana sale negra al completo a partir de linea debub

// Clean previous data in FIJI
run("Close All");
run("Clear Results");
print("\\Clear");
if(roiManager("count") !=0) {roiManager("delete");}

// Set measurements 
run("Options...", "iterations=1 count=1 black");
 // Set black binary bckg
run("Set Measurements...", "area center perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;area;XM;XM;Perim;Circ;Feret;FeretAngle;MinFeret;AR;Round;Solidity;NBranches;AvgBranchLen;MaxBranchLen;BranchLen;EuclideanDist;Time"); // header of the result file

// 0. Select the Folder with the files
dir = getDirectory("Selecciona el directorio con las imagenes en tiff");
list= getFileList (dir);
Results = createFolder(dir, "Results");

Start_time = getTime(); // to inform how long does it take to process the folder
setBatchMode(false);
	// 0.1 Loop to open and process each file
imagen = 0;
for (i=0; i<list.length; i++){
	if (endsWith(list[i], "if")){
		imagen = imagen+1;
	// 0.2 Open and get data
		title=list[i];
		//open(dir+title);
		run("Bio-Formats (Windowless)", "open=["+dir+title+"]");
		run("Select None");
		
// 1. Get dimensions
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pw, ph, pd);
		frame_interval = Stack.getFrameInterval();
		
// 2. Process
		rename("original");
		original = getImageID();
		run("Z Project...", "projection=Median");
		newImage("stack_temp", "32-bit black", width, height, frames);	
		imageCalculator("Add stack", "stack_temp","MED_original");
		close("MED_original");
		rename("proyection_temp");
		imageCalculator("Subtract create stack", "proyection_temp","original");
		run("16-bit");
		run("Gaussian Blur...", "sigma=1 stack");
		debug;
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=RenyiEntropy background=Dark black");
		mask = getImageID();

// 2.3 Segment the worm
	//2.3.1 Loop for every temporal frame
			for (t = 0; t < frames; t++) {
				selectImage(mask);
				Stack.setFrame(t+1);
				run("Duplicate...", "title=binary");
				
	//2.3.2 Get the largtest element
				run("Analyze Particles...", "size=0-Infinity display add");
				rename("binary_temp"); // Output
				binary_temp=getImageID();
				//run("Grays");
				run("Fill Holes");
				selectWindow("Results");
				Area_column = Table.getColumn("Area");
				indices_max = Array.findMaxima(Area_column, 1);
				roiManager("Select", indices_max[0]);
				setBackgroundColor(0, 0, 0);
				run("Clear Outside");
				run("Select None");
				// Clean the results and ROI to later measure again only the largest particle
				run("Clear Results");
				if(roiManager("count") !=0) {roiManager("delete");}
				// binary closing				
				run("Maximum...", "radius=9 stack");
				run("Minimum...", "radius=11 stack");

				
				

	// 2.3.3 Get several features in the frame
				run("Analyze Particles...", "display add");
				selectWindow("Results");
				if (nResults == 1) {
					area = getResultString("Area", 0);
					XM = getResultString("XM", 0);
					YM = getResultString("YM", 0);
					Perim = getResultString("Perim.", 0);
					Circ = getResultString("Circ.", 0);
					Feret = getResultString("Feret", 0);
					FeretAngle = getResultString("FeretAngle", 0);
					MinFeret = getResultString("MinFeret", 0);
					Solidity = getResultString("Solidity", 0);
					AR = getResultString("AR", 0);
					Round = getResultString("Round", 0);
				} else { 
					wait(50);
					roiManager("delete"); // To avoid an error
					area = "NA";
					XM = "NA";
					YM = "NA";
					Perim = "NA";
					Circ = "NA";
					Feret ="NA";
					FeretAngle = "NA";
					MinFeret = "NA";
					Solidity = "NA";
					AR="NA";
					Round="NA";
				
						}

	// 2.3.4 Skeletonize and measure
				selectImage(binary_temp);
				run("Duplicate...", "title=skeleton_temp");
				run("Skeletonize");
				//rename("skeleton_temp");
				skeleton_temp=getImageID();
				if (nResults == 1) { // to avoid erros due to empty images
					run("Clear Results");
					run("Analyze Skeleton (2D/3D)", "prune=none show");
					selectWindow("Results");
					NBranches = getResultString("# Branches", 0);
					AvgBranchLen = getResultString("Average Branch Length", 0);
					MaxBranchLen = getResultString("Maximum Branch Length", 0);
					run("Clear Results");
					selectWindow("Branch information");
					BranchLen = getResultString("Branch length", 0);
					EuclideanDist = getResultString("Euclidean distance", 0);
					V1x = getResultString("V1 x", 0);
					V1y = getResultString("V1 y", 0);
					V2x = getResultString("V2 x", 0);
					V2y = getResultString("V2 y", 0);
					close("Tagged skeleton");
				} else {
					NBranches = "NA";
					AvgBranchLen = "NA";
					MaxBranchLen = "NA";
					BranchLen = "NA";
					EuclideanDist = "NA";
					V1x = 0;
					V1y = 0;
					V2x = 0;
					V2y = 0;
				}
				
	// 2.3.5 Write the results of the frame in the table			
				print((t+1)+";"+area+";"+XM+";"+YM+";"+Perim+";"+Circ+";"+Feret+";"+FeretAngle+";"+MinFeret+";"+AR+";"+Round+";"+Solidity+";"+NBranches+";"+AvgBranchLen+";"+MaxBranchLen+";"+BranchLen+";"+EuclideanDist+";"+frame_interval*t);

// 3. Draw segmentation	on the original at this frame
	// 3.1 Get the frame image to print the results
				selectImage(original);
				Stack.setFrame(t+1);
				run("Duplicate...", "title=result_temp");
				result_temp = getImageID();
				run("RGB Color");
				// If there is one Roi, then draw 
				if(roiManager("count") >= 1) {
	//3.2 Draw the segmented worm on the original image
					selectImage(result_temp);
					roiManager("Select", 0);
					setForegroundColor(255, 255, 0); // draw in yellow
					run("Draw", "slice");
	// 3.3 Draw Feret
					//run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1.0000");
					List.setMeasurements;
					x1 = List.getValue("FeretX")*pw;
					y1 = List.getValue("FeretY")*pw;
					length = List.getValue("Feret");
					degrees = List.getValue("FeretAngle");
					if (degrees>90){degrees -= 180;}
					angle = degrees*PI/180;
					x2 = x1 + cos(angle)*length;
					y2 = y1 - sin(angle)*length;
					setForegroundColor(255, 0, 0);  // draw in red
					drawLine(x1/pw, y1/pw, x2/pw, y2/pw); // functions needs arguments in pixels
					// 3.3 Draw the Skeleton on the original image
					selectImage(skeleton_temp);
					run("Create Selection");
					selectImage(result_temp);
					run("Restore Selection");
					setForegroundColor(0, 255, 255);  // draw in skyblue
					run("Fill", "slice");
					run("Select None");
					// 3.4 Draw the Euclidean Distance
					setForegroundColor(0, 255, 0);  // draw in green
					drawLine(V1x, V1y, V2x, V2y); // functions needs arguments in pixels
					roiManager("Delete");
				} 
				close("skeleton_temp");
				close("binary_temp");
				run("Clear Results");
				// Create the result as stack
				if (t==0) {
					selectImage(result_temp);
					rename("Stack_Result");
				} else {
					run("Concatenate...", "  title=Stack_Result open image1=Stack_Result image2=result_temp image3=[-- None --]");
				}
			} // End of loop for every frame
			
// 4. Save the results
	// 4.1 Save the result stack image
		selectWindow("Stack_Result");
		run("Scale...", "x=0.5 y=0.5 z=1.0 interpolation=Bilinear fill process create");
		rename(title+"_result");
		saveAs("Tiff", Results+title+"_segment.tif");

	// 4.2 Save results and clean for the next image
		selectWindow("Log");
		saveAs("Text", Results+title+"_Results.csv");
		print("\\Clear");
		print("Frame;area;XM;XM;Perim;Circ;Feret;FeretAngle;MinFeret;AR;Round;Solidity;NBranches;AvgBranchLen;MaxBranchLen;BranchLen;EuclideanDist;Time");
		run("Close All");
		run("Clear Results");
	}
}
setBatchMode(false);

// Macro is finished						
print("\\Clear");
print("Terminado");
Finish_time = getTime();
Time_used = Finish_time - Start_time;
print("It took =", Time_used/1000, "second to finish the proccess");
//Functions

function createFolder(dir, name) {
	mydir = dir+name+File.separator;
	File.makeDirectory(mydir);
	if(!File.exists(mydir)){
		print("Unable to create the folder");
	}
	return mydir;
}
