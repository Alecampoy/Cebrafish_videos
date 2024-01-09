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
run("Options...", "iterations=1 count=1 black");
 // Set black binary bckg
setBackgroundColor(0, 0, 0);
run("Set Measurements...", "mean perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;X;Y;Mean-Distance;Time;"); // header of the result file in the Log window


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
		run("Movie (FFMPEG)...", "choose=["+dir+title+"] first_frame=500 last_frame=4701"); // 14 minutes video
		rename("original");
		run("8-bit");
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
		run("Gamma...", "value=0.41");
		run("Gaussian Blur...", "sigma=2");
		// run("Enhance Contrast...", "saturated=0.001 normalize process_all"); // optional if wand does not work properly
		doWand(width/2, height/2, 4.0, "Legacy");
		run("Fit Circle");
		roiManager("Add");
		run("Create Mask");
		run("Distance Map");
		distance_map=getImageID();
		rename("distance_map");
		
	// 2.2 process of the original image
		selectImage(original);
		run("Gaussian Blur...", "sigma=1 stack"); 
		run("Z Project...", "projection=Median");
		rename("median_proyection");
		run("Invert");
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=-8");
		run("Gaussian Blur...", "sigma=8"); // elimina cualquier rastro del pez en caso de que este tanto tiempo quieto que aparezca en la proyeccion mediana
		imageCalculator("Add stack", "original","median_proyection");
		selectImage(original);
		run("Invert", "stack");
		// limpio fuera del pocillo para evitar que se detecte debris que ocurre en el video
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=8");
		run("Clear Outside", "stack");
		run("Select None");
		
	//2.3 Loop for every temporal frame to detect the points
		 	for (t = 0; t < slices; t++) {
			selectImage(original);
			run("Select None");
			Stack.setSlice(t+1);
			run("Find Maxima...", "prominence=83 output=[Point Selection]");
			selectImage(distance_map);
			run("Restore Selection");
			run("Measure");
			run("Select None");
							
	// 2.3.3 Get features in the frame
			selectWindow("Results");
				if (nResults == 1) {
					if (t>0) {X_0=X; Y_0=Y;} //to draw a line later
					frame = t+1;
					X = getResultString("X", 0);
					Y = getResultString("Y", 0);
					Distance_edge = getResultString("Mean", 0); // the distance to the edge is the mean value of the distance map
					time = frame/5; // 5 fps
					
				} else { 
					wait(50);
					frame = t+1;
					X = "NA";
					Y = "NA";
					Distance_edge = "NA";
					time = frame/5; // 5 fps
				
}

	// 2.3.5 Write the results of the frame in the table			
				print(frame+";"+X+";"+Y+";"+Distance_edge+";"+time);

// 3. Draw result on the original at this frame
	// 3.1 Get the frame image to print the results
				selectImage(original);
				Stack.setSlice(t+1);
				run("Duplicate...", "title=result_temp");
				result_temp = getImageID();
				run("RGB Color");
				// If there is one Roi, then draw 
				if (nResults == 1 && t>0){
					selectImage(result_temp);
					setForegroundColor(255, 0, 0);  // draw in red
					// pinto una linea del movimiento entre frame y frame
					drawLine(X_0/pw, Y_0/pw, X/pw, Y/pw); // functions needs arguments in pixels
				}
				run("Scale...", "x=0.6 y=0.6 z=1.0 interpolation=Bilinear fill process create"); 
				result_temp_2 = getImageID();
				close("result_temp");
				selectImage(result_temp_2);
				rename("result_temp");
			
				// Create the result as stack concatenation
				if (t==0) {
					selectImage(result_temp_2);
					rename("Stack_Result");
				} else {
					run("Concatenate...", "  title=Stack_Result open image1=Stack_Result image2=result_temp image3=[-- None --]");
				}

			run("Clear Results");
			
			} // End of loop for every frame

// 4. Save the results
	// 4.1 Save the result stack image
		selectWindow("Stack_Result");
		rename(title+"_result");
		saveAs("Tiff", Results+title+"_tracking.tif");
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", Results+title+"_tracking_projection.tif");
		
		
	// 4.2 Save results and clean for the next image
		selectWindow("Log");
		saveAs("Text", Results+title+"_Results.csv");
		print("\\Clear");
		print("Frame;X;Y;Mean-Distance;Time;");  // header of the result file  in the Log window
		run("Close All");
		roiManager("delete"); 
		run("Clear Results");
	}
}		

setBatchMode(false);
// Macro is finished. Print time						
print("\\Clear");
print("Terminado");
Finish_time = getTime();
Time_used = Finish_time - Start_time;
print("It took =", Time_used/60000, "minutes to finish the proccess");


//Functions
function createFolder(dir, name) {
	mydir = dir+name+File.separator;
	File.makeDirectory(mydir);
	if(!File.exists(mydir)){
		print("Unable to create the folder");
	}
	return mydir;
}
