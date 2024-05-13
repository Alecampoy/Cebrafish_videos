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
run("Set Measurements...", "area mean perimeter fit shape feret's area_fraction stack redirect=None decimal=2");
print("Frame;X;Y;Mean-Distance;Time"); // header of the result file in the Log window

// Parent Folder to process all batches

dir_parent = getDirectory("Select the folder where each folder is a batch");
list_parent =  getFileList (dir_parent); // lista de las carpetas en dir_parent

for (j = 0; j<list_parent.length; j++) { // loop en las carpetas de los batches, llega hasta abajo
	dir = dir_parent+list_parent[j]; // carpeta con las peliculas


	// Folder with the files
	list= getFileList(dir);
 // lista con los archivos que se van a procesar
	if (STRICT == true) {Results = createFolder(dir, "Results_strict_"+strict_value);
}
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
		

// 4. Save the results
	// 4.1 Save the result stack image
		selectWindow("well_edge");
		roiManager("Select", 0);
		run("Flatten");
		rename(title+"_result");
		//saveAs("Tiff", Results+title+"_tracking.tif");
		//run("Z Project...", "projection=[Max Intensity]");
		saveAs("Jpeg", Results+title+"_tracking_projection.tif");
		close("*");
		roiManager("reset"); 
		run("Clear Results");
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
