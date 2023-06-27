///////////////////////////////////////////////////////////////////////////////////////////////
/* Author: Ale Campoy
 * Microscopy Unit (CABD)
 * Date: 27/06/2023
 * User: Marta Fernandez	
 * 	
 * Description: Tracks the position of a Cebra fish inside a well and measures the time it stays in contact the edges of the well
 * 
 * Input: Folder with the set of bright field images .mp4 at 720p & 5fps from an Iphone 
 *
 * Method: see segmentation details
 * 
 * Output: Segmented images in folder and result file
 * 
 *///////////////////////////////////////////////////////////////////////////////////////////////


// 0.0 Clean previous data in FIJI
run("Close All");
run("Clear Results");
print("\\Clear");
if(roiManager("count") !=0) {roiManager("delete");}

// 0.1 Set measurements
run("Options...", "iterations=1 count=1 black"); // Set black binary bckg
run("Set Measurements...", "area center perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;X;Y;Frame;Mean-Distance;Time;"); // header of the result file in the Log window

// 1 Select the Folder with the files
dir = getDirectory("Select the folder with the .mp4 movies");
list= getFileList (dir);
Results = createFolder(dir, "Results");

Start_time = getTime(); // to inform how long does it take to process the folder
setBatchMode(false);

//  Loop to open and process each file
for (i=0; i<list.length; i++){
	if (endsWith(list[i], ".mp4")){
	
	// 1.2 Open and get data
		title=list[i];
		run("Movie (FFMPEG)...", "choose=["+dir+title+"] first_frame=1001 last_frame=5500");
		rename("original");
		original = getImageID();
				
	// 1.3 Get dimensions
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pw, ph, pd);
		frame_interval = Stack.getFrameInterval();
		
// 2. Process
	// 2.1 generate the distance map
		selectImage(original);
		Stack.setFrame(frames/2);
		run("Duplicate...", "title=well_edge");
		run("Gaussian Blur...", "sigma=2");
		// run("Enhance Contrast...", "saturated=0.001 normalize process_all"); // optional if wand does not work properly
		doWand(width/2, height/2, 12.0, "Legacy");
		run("Create Mask");
		run("Morphological Filters", "operation=Opening element=Disk radius=80");
		run("Distance Map");
		distance_map=getImageID();
		rename("distance_map");
		
	// 2.2 process of the original image
		selectImage(original);
		run("8-bit");
		run("Gaussian Blur...", "sigma=1 stack"); // opcional
		run("Z Project...", "projection=Median");
		rename("median_proyection");
		run("Invert");
		imageCalculator("Add stack", "original","median_proyection");
		run("Invert", "stack");
		
	//2.1 Loop for every temporal frame to detect the points
		for (t = 0; t < frames; t++) {
				
	// 2.3.3 Get features in the frame
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
					roiManager("delete"); // To avoid an error if ROI Manager has several ROI
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

				close("binary_temp");
				run("Clear Results");
				// Create the result as stack concatenation
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
// Macro is finished. Print time						
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
