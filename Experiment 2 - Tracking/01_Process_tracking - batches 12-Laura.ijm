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
if (STRICT == true) {strict_value = 11;} // quizas poner el argumento  strict para que no aparezca ningun punto si el pez no se detecta

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
run("Set Measurements...", "area mean centroid center perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;X;Y;Mean-Distance;Time;Pocillo_XM;Pocillo_YM;Pocillo_diam"); // header of the result file in the Log window

// Parent Folder to process all batches

dir_parent = getDirectory("Select the folder where each folder is a batch");
list_parent =  getFileList (dir_parent); // lista de las carpetas en dir_parent

for (j = 0; j<list_parent.length; j++) { // loop en las carpetas de los batches, llega hasta abajo
	dir = dir_parent+list_parent[j]; // carpeta con las peliculas
	// Folder with the files
	list= getFileList(dir); // lista con los archivos que se van a procesar
	if (STRICT == true) {Results = createFolder(dir, substring(list_parent[j],0,lengthOf(list_parent[j])-10)+"_Results_strict_"+strict_value);}
	else {Results = createFolder(dir, list_parent[j]+"Results");}

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
		List.setMeasurements;
  		pocillo_XM = List.getValue("XM");
		pocillo_YM = List.getValue("YM");
		pocillo_diam = List.getValue("Feret");
		roiManager("Add");
		run("Fit Circle");
		roiManager("Add");
		roiManager("Select", roiManager("count")-1); // selecciona la ultima para que se ejecute la acción
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
		run("Enlarge...", "enlarge=12");
		mean_bck = getValue("Mean raw");
		setBackgroundColor(mean_bck, mean_bck, mean_bck);
		run("Clear Outside", "stack"); 
		run("Select None");
		
		// 2.2.1 proyección para pintar el resultado
		run("Z Project...", "projection=[Max Intensity]");
		run("RGB Color");
		result_proyection = getImageID();
		rename("result_proyection");
		
		//2.3 Loop for every temporal frame to detect the points
		for (t = 0; t < slices; t++) {
			selectImage(original);
			run("Select None");
			Stack.setSlice(t+1);
			wait(5);
			if (STRICT == true) {
				run("Find Maxima...", "prominence="+strict_value+" strict output=[Point Selection]");
				wait(5);
				// compruebo el numero de puntos-maxima que ha detectado. Quiero garantizar que hay uno 1
				getSelectionCoordinates(x, y);
				n_results = x.length;
				if (n_results >=2 && t !=0) { // creo que este if es para quedarme solo con un elemento en caso que aparezcan varios, sospecho que no se usa        
					roiManager("reset"); 
					run("Clear Results"); 
					selectImage(original);
					run("Select None");
					Stack.setSlice(t+1);
					run("Find Maxima...", "prominence="+strict_value+" strict output=[Maxima Within Tolerance]");
					temp_2_results = getImageID();
					run("Analyze Particles...", "size=0-Infinity display add");
					selectWindow("Results");
					wait(5);
					Area_column = Table.getColumn("Area");
					if (Area_column.length == 1) {
					indices_min = newArray(0); // Manually create array with index 0
					indices_min[0]=0;
					} else {
						  indices_min = Array.findMinima(Area_column, 0, 0);
						}
						if (Area_column[indices_min[0]] < 108) {
							roiManager("Select", indices_min[0]);
							setBackgroundColor(0, 0, 0);
							run("Clear Outside");
							run("Select None");
							run("Find Maxima...", "prominence=100 output=[Point Selection]");
							}
						selectImage(temp_2_results);
						close();
				} // fin del if sospechoso
			} // Finish strict = true
			
			else { // NO STRICT considers the brightest point
				run("Find Maxima...", "prominence=100 output=[Point Selection]");
				getSelectionCoordinates(x, y);
				n_results = x.length;
			}

			// measures the points
			selectImage(distance_map);
			wait(5);			
			run("Restore Selection");
			List.setMeasurements; // equivalent to run("Measure");
			run("Select None");

		// 2.4 Get the measurements avoiding errors
			if (t==0) {X="NA"; X_0="NA";} // to avoid an error in following 'if' clauses to paint the result
			// selectWindow("Results");
			if (n_results == 1) { // medida correcta cuando tenemos solo un Punto
				wait(5);	
				if (t>0 && X != "NA") {X_0=X; Y_0=Y;} //to draw a line later, avoiding the NA to draw a line
				X = toString(List.getValue("X"));
				Y = toString(List.getValue("Y"));
				Distance_edge = toString(List.getValue("Mean")); // mean represent the distance to the edge is the mean value of the distance map
			} else { // more than 1 selection
				wait(5);
				X = "NA";
				Y = "NA";
				Distance_edge = "NA";
			}
			frame = t+1;
			time = frame*frame_interval;

			// 2.5 Write the results of the frame in the table			
			print(frame+";"+X+";"+Y+";"+Distance_edge+";"+time+";"+pocillo_XM+";"+pocillo_YM+";"+pocillo_diam);

// 3. Draw result on the original at this frame
			// 3.1 Get the result projection and print results
			selectImage(result_proyection);
			// If there is one Roi, then pinto una linea del movimiento entre frame y frame
			if (t==0) {
				setForegroundColor(255, 255, 0); // Draw in yellow
				roiManager("Select", 0);
				run("Draw", "slice");
				run("Select None");
			} else {setForegroundColor(255, 0, 0);}  // draw in red
			// pinto solo si la medida ha sido correcta
			if (n_results == 1 && t>0 && X_0 != "NA" && X != "NA"){
				setForegroundColor(255, 0, 0);  // draw in red
				drawLine(X_0/pw, Y_0/pw, X/pw, Y/pw); // functions needs arguments in pixels
			}
		} // End of loop for every frame

// 4. Save the results
		// 4.1 Save the result stack image
		selectWindow("result_proyection");
		rename(title+"_result");
		saveAs("Tiff", Results+title+"_tracking_projection.tif");
		
		// 4.2 Save results and clean for the next image
		selectWindow("Log");
		saveAs("Text", Results+title+"_Results.csv");
		print("\\Clear");
		print("Frame;X;Y;Mean-Distance;Time;Pocillo_XM;Pocillo_YM;Pocillo_diam"); // header of the result file in the Log window
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
