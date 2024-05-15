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

STRICT = true; // argumento de find maxima
if (STRICT == true) {strict_value = 10;} // quizas poner el argumento  strict para que no aparezca ningun punto si el pez no se detecta

// 0.0 Clean previous data in FIJI
run("Close All");
run("Clear Results");
print("\\Clear");
roiManager("reset");
Start_time = getTime(); // to inform how long does it take to process the folder
setBatchMode(false);

// 0.1 Set measurements
run("Options...", "iterations=1 count=1 black");
 // Set black binary bckg
setBackgroundColor(0, 0, 0);
run("Set Measurements...", "area mean perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;X;Y;Mean-Distance;Time"); // header of the result file in the Log window

// Parent Folder to process all batches

dir_parent = getDirectory("Select the folder where each folder is a batch");
list_parent =  getFileList (dir_parent); // lista de las carpetas en dir_parent

for (j = 0; j<list_parent.length; j++) { // loop en las carpetas de los batches, llega hasta abajo
	dir = dir_parent+list_parent[j]; // carpeta con las peliculas

	// Folder with the files
	list= getFileList(dir); // lista con los archivos que se van a procesar
	if (STRICT == true) {Results = createFolder(dir, "Results_strict_"+strict_value);}
	else {Results = createFolder(dir, "Results");}

// 1. Loop to open and process each file
	for (i=0; i<list.length; i++){
		if (endsWith(list[i], ".tif")){
	
		// 1.2 Open and get data
		title=list[i];
		open(dir+title);
		run("Select None");
		rename("original");
		original = getImageID();		
				
		// 1.3 Get dimensions
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pw, ph, pd);
		frame_interval = 1/6; // number of frames of the acquisition videos per second
				
// 2. Process files
	// 2.1 generate the distance map
		selectImage(original);
		Stack.setSlice(slices/2);
		run("Duplicate...", "title=well_edge");
		run("Gaussian Blur...", "sigma=1");
		run("Enhance Contrast...", "saturated=0.40 normalize"); // ojo: tienen que estar las imagenes limpias por fuera del pocillo. el cartel perturba esta ejecucion
		run("Gamma...", "value=1.42");
		run("Gaussian Blur...", "sigma=1");
		run("Subtract Background...", "rolling=50 light");
		wand=24;
		doWand(width/2, height/2, wand, "4-connected");
		roiManager("Add");
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
		//roiManager("Select", 1);
		//run("Enlarge...", "enlarge=-8");
		//run("Gaussian Blur...", "sigma=8"); // elimina cualquier rastro del pez en caso de que este tanto tiempo quieto que aparezca en la proyeccion mediana
		imageCalculator("Add stack", "original","median_proyection");
		selectImage(original);
		run("Invert", "stack");
		// limpio fuera del pocillo para evitar que se detecte debris que ocurre en el video
		roiManager("Select", 1);
		run("Enlarge...", "enlarge=18");
		mean_bck = getValue("Mean raw");
		setBackgroundColor(mean_bck, mean_bck, mean_bck);
		run("Clear Outside", "stack"); 
		run("Select None");
		
	//2.3 Loop for every temporal frame to detect the points
		 	for (t = 0; t < slices; t++) {
			selectImage(original);
			run("Select None");
			Stack.setSlice(t+1);
			wait(32);
			if (STRICT == true) {run("Find Maxima...", "prominence="+strict_value+" strict output=[Point Selection]");
				// start selecting rois
				// keep only the littles roi, probably the fish
				run("Measure");
				if (nResults >=2 && t !=0) {
					roiManager("reset"); // clean the roi for further filtering of the big particles
					run("Clear Results"); // cleaned for analyzed particles
					selectImage(original);
					run("Select None");
					Stack.setSlice(t+1);
					run("Find Maxima...", "prominence="+strict_value+" strict output=[Maxima Within Tolerance]");
					temp_2_results = getImageID();
					run("Analyze Particles...", "size=0-Infinity display add");
					selectWindow("Results");
					Area_column = Table.getColumn("Area");
					indices_min = Array.findMinima(Area_column, 0);
					if (Area_column[indices_min[0]] < 108) {
						roiManager("Select", indices_min[0]);
						setBackgroundColor(0, 0, 0);
						run("Clear Outside");
						run("Select None");
						run("Find Maxima...", "prominence=100 output=[Point Selection]");
						}
					selectImage(temp_2_results);
					close();
					}
				// common measurement of distance map
					run("Clear Results");
					selectImage(distance_map);
					wait(15);			
					run("Restore Selection");
					run("Measure");
					run("Select None");

			} // Finish strict
			
			else { // NO STRICT
				run("Find Maxima...", "prominence=100 output=[Point Selection]");	// no strict, considers the fish the brightest point
				selectImage(distance_map);
				wait(32);
				run("Restore Selection");
				run("Measure");
				run("Select None");
				}


	// 2.3.3 Get features in the frame
			if (t==0) {X="NA"; X_0="NA";} // to avoid an error in following 'if'. It behaves different as in python where only the first argument is evaluated if false
			selectWindow("Results");
				if (nResults == 1) {
					wait(32);	
					if (t>0 && X != "NA") {X_0=X; Y_0=Y;} //to draw a line later, avoiding the NA to draw a line not appropiate.
					frame = t+1;
					X = getResultString("X", 0);
					Y = getResultString("Y", 0);
					Distance_edge = getResultString("Mean", 0); // the distance to the edge is the mean value of the distance map
					time = frame/5; // 5 fps
					
				} else { 
					wait(32);
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
				run("Select None");
				run("Duplicate...", "title=result_temp");
				result_temp = getImageID();
				run("RGB Color");
				// If there is one Roi, then draw 
				if (nResults == 1 && t>0 && X_0 != "NA" && X != "NA"){
					selectImage(result_temp);
					setForegroundColor(255, 0, 0);  // draw in red
					// pinto una linea del movimiento entre frame y frame
					drawLine(X_0/pw, Y_0/pw, X/pw, Y/pw); // functions needs arguments in pixels
					//run("Scale...", "x=0.6 y=0.6 z=1.0 interpolation=Bilinear fill process create"); 
					//result_temp_2 = getImageID();
					//close("result_temp");
					//selectImage(result_temp_2);
					//rename("result_temp");
				}
				// Create the result as stack concatenation
				if (t==0) {
					setForegroundColor(255, 255, 0); // Draw in yellow
					selectImage(result_temp);
					roiManager("Select", 0);
					run("Draw", "slice");
					run("Select None");
					//run("Scale...", "x=0.6 y=0.6 z=1.0 interpolation=Bilinear fill process create"); 
					rename("Stack_Result");
					roiManager("reset"); // clean the roi for further filtering of the big particles
				} else {
					imageCalculator("Max", "Stack_Result","result_temp");
					close("result_temp");
					//run("Concatenate...", "  title=Stack_Result open image1=Stack_Result image2=result_temp image3=[-- None --]");
				}

			run("Clear Results");
			
			} // End of loop for every frame

// 4. Save the results
	// 4.1 Save the result stack image
		selectWindow("Stack_Result");
		rename(title+"_result");
		//saveAs("Tiff", Results+title+"_tracking.tif");
		//run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", Results+title+"_tracking_projection.tif");
		
		
	// 4.2 Save results and clean for the next image
		selectWindow("Log");
		saveAs("Text", Results+title+"_Results.csv");
		print("\\Clear");
		print("Frame;X;Y;Mean-Distance;Time");  // header of the result file  in the Log window
		close("*");
		roiManager("reset"); 
		run("Clear Results");
		run("Collect Garbage");
	}
}		

} // cierre loop de carpeta_parent

setBatchMode(false);
// Macro is finished. Print time						
print("\\Clear");
print("Terminado");
Finish_time = getTime();
Time_used = Finish_time - Start_time;
print("It took =", Time_used/(60000 * 60), "h to finish the proccess");


//Functions
function createFolder(dir, name) {
	mydir = dir+name+File.separator;
	File.makeDirectory(mydir);
	if(!File.exists(mydir)){
		print("Unable to create the folder");
	}
	return mydir;
}
