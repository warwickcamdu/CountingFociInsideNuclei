////FIXED PARAMETERS:


//Channels:

//Channel number with stain/fluorophore for nuclei
nucleiChannel = 3;

//Channel number with stain/fluorophore for foci
fociChannel = 1;


//Project preferences:

//Preferred Z project method for each channel
nucleiProjectPreference = "Max Intensity"
maximaProjectPreference = "Max Intensity"


//Background subtraction:

// Size of rolling ball used for background subtraction
backgroundSubtractionBallPixelRadius = 50;


//Thresholding:

//Preferred thresholding method for the nuclei
nucleiThresholdMethod = "Triangle";


//Batch mode
setBatchMode(true);


//MAIN METHOD:

//Input directory where images are stored
#@ File (label = "Input directory", style = "directory") input

//Output directory for ROIs and CSVs to be saved
#@ File (label = "Output directory", style = "directory") output

//File suffix for images
#@ String (label = "File suffix", value = ".tif") suffix

//Nucleus size parameters
#@ Integer (label="Minimum nucleus size (micron^2)", style="slider", min=0, max=2000, value=500) nucleiMin
#@ Integer (label="Maximum nucleus size (micron^2)", style="slider", min=0, max=2000, value=500) nucleiMax

//Tick box to use standard deviation to choose most in-focus slice for nuclei detection
//(default is Z projection of the whole stack)
#@ Boolean (label="Select most in-focus nuclei slice") select_focus

//Foci prominence threshold
#@ Integer (label="Prominence of foci", style="slider",min=0, max=300, value=50) foci_prominence



//Count of images processed
imageCount = 0;

//Process the input folder
processFolder(input);

//Completion dialogue
Dialog.create("Process Finished");
Dialog.addMessage("Foci counting completed for " + imageCount + " images.\nYou can find results as CSV files and image files in" +output);
Dialog.show();



//Function for cycling through folder, finding images and processing them
function processFolder(input)
{
	//Get and arrange list of files/directories in the input folder
	list = getFileList(input);
	list = Array.sort(list);
	
	
	//Cycle through the list
	for (i = 0; i < list.length; i++) {
		//If the item is a directory, process this directory
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		//If the item is an image with the given suffix, then process the image
		if(endsWith(list[i], suffix))
		{
			processFile(input, output, list[i]);
			imageCount++;
		}
	}
}


//Find nuclei and foci, and create output with foci counts per nucleus.
//Please note, folder and file names cannot contain square brackets.
function processFile(input, output, file)
{
	//Use Bioformats to import file
	importString = "open=[" + input + File.separator + file + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack";
	run("Bio-Formats Importer", importString);
	
	//Name the image currentImage
	currentImage = getTitle();
	
	//Get variables for the dimensions of the image
	getDimensions(width, height, channels, slices, frames);
	
	
	//Get nuclei as ROIs
	identifyNuclei(currentImage, nucleiChannel);
	
	//Create resutls table for nuclei
	Table.create("Nuclei Results");
	
	//Numbers will be linked to nucleus number in the ROIs and other results
	Table.showRowNumbers(true);
	
	
	//Measure characteristics of all nuclei in each channel
	for (i = 1; i <= channels; i ++)
	{
		selectWindow(currentImage);
		Stack.setChannel(i);
		
		//Ensures all ROIs are measured
		roiManager("Deselect");
		
		roiManager("Measure");
		
		intDenResults = getResultsAsArray("RawIntDen");
		
		selectWindow("Nuclei Results");
		Table.setColumn("Channel " + i + " Raw Integrated Density", intDenResults);
		intDenResults = getResultsAsArray("IntDen");
		
		selectWindow("Nuclei Results");
		Table.setColumn("Channel " + i + " Integrated Density", intDenResults);
		intDenResults = newArray(1);
		
		meanResults = getResultsAsArray("Mean");
		
		selectWindow("Nuclei Results");
		Table.setColumn("Channel " + i + " Mean", meanResults);
		meanResults = newArray(1);
		
		medianResults = getResultsAsArray("Median");
		
		selectWindow("Nuclei Results");
		Table.setColumn("Channel " + i + " Median", medianResults);
		medianResults = newArray(1);
		
		//Wipe the results table
		run("Clear Results");
		updateResults();
	}
	
	//Now get the number of nuclei from number of rows in
	//the table we created from the ROIs
	numNuclei = Table.size;
	
	//Make array that allows all nucleus ROIs to be selected
	nucleusROIs = newArray(numNuclei);
	for (i = 0; i < numNuclei; i++)
	{
		nucleusROIs[i] = i;
	}
	
	//Empty array to populate column for number of foci inside nucleus
	fociColumn = newArray(numNuclei);
	
	//Get maxima as point selection ROI
	identifyMaxima(currentImage, fociChannel);
	
	//Index used to select the ROI containing the foci
	maximaROIIndex = numNuclei;
	
	//Select all nucleus ROIs and combine them into a single ROI
	roiManager("select",nucleusROIs);
	roiManager("Combine");
	roiManager("Add");
	
	//Index used to select ROI containing all nuclei 
	nucleiROIIndex = numNuclei+1;
	
	//Select all nuclei and maxima
	maximaAndNuclei = newArray(maximaROIIndex,nucleiROIIndex);
	
	roiManager("select",maximaAndNuclei);

	
	//AND filter creates new ROI with only maxima that are inside nuclei
	roiManager("AND");
	
	
	//Check if there are any foci inside nuclei. If so, add to ROI manager
	type = selectionType();
	
	if (type != -1)
	{
		roiManager("Add");
	}
	
	else
	{
		print("Error - no foci inside nuclei");
	}
	
	
	//Get list of the co-ordinates of all foci inside nuclei
	getSelectionCoordinates(maximaX, maximaY);
	
	//Count number of such foci
	numMaxima = maximaX.length;
	
	//Create table for results regarding the foci
	Table.create("Foci Results");
	
	//Empty arrays to populate columns in foci table
	nucleiColumn = newArray(numMaxima);
	
	xCoordsColumn = newArray(numMaxima);
	
	yCoordsColumn = newArray(numMaxima);
	
	intensityColumns = newArray(numMaxima*channels);
	
	
	//Show row numbers in table
	Table.showRowNumbers(true);
	
	//Variable for reading and storing intensities
	intensity = 0;
	
	//Variable for stepping through foci in each nucleus
	fociCount = 0;
	
	
	//Index for selecting ROI of all foci in a particular nucleus
	subsetIndex = numNuclei + 3;

	//Counter for which row of foci table to write in
	fociRow = 0;
	
	//Step through nuclei and their foci and record all data in the two tables
	for (i = 0; i < numNuclei; i++)
	{
		//Create ROI of foci in only the ith nucleus
		roiManager("Select", newArray(i,numNuclei));
		roiManager("AND");
		
		//If there are foci in the nucleus, record how many and the data for each focus
		type = selectionType();
		if (type != -1)
		{
			//Add the ROI to the manager
			roiManager("Add");
			
			roiManager("Select", subsetIndex);
			
			//Get co-ordinates of points
			getSelectionCoordinates(xPoints, yPoints);
			
			//Count number of foci and add to the relevant table column
			fociCount = xPoints.length;
			fociColumn[i] = fociCount;
			
			//Step through foci
			for (j = 0; j < fociCount; j++)
			{
				//Record which nucleus the focus belongs to
				nucleiColumn[fociRow] = i+1;
				
				//Record the co-ordinates
				xCoordsColumn[fociRow] = xPoints[j];
				yCoordsColumn[fociRow] = yPoints[j];
				
				//Step through channels
				for (k = 0; k < channels; k++)
				{
					//Select the given channel in the image
					selectWindow(currentImage);
					Stack.setChannel(k+1);
					
					//Get the intensity at the focus
					intensity = getPixel(xPoints[j], yPoints[j]);
					
					//Which position in the array to save the intensity in
					index = (k*numMaxima) + fociRow;
					
					//Save the intensity
					intensityColumns[index] = intensity;
				}
				fociRow ++;
			}
			//Delete the ROI for this nucleus's foci
			roiManager("delete");
		}
		
		//If no foci, simply record the zero foci in the table
		else
		{
			fociColumn[i] = 0;
		}
	}
	
	//Place the foci column in the table regarding nuclei
	selectWindow("Nuclei Results");
	Table.setColumn("Number of Foci", fociColumn);
	fociColumn = newArray(1);
	
	//Place the various columns in the table regarding foci
	selectWindow("Foci Results");
	Table.setColumn("Nucleus", nucleiColumn);
	nucleiColumn = newArray(1);
	Table.setColumn("X Co-ordinate", xCoordsColumn);
	xCoordsColumn = newArray(1);
	Table.setColumn("Y Co-ordinate", yCoordsColumn);
	yCoordsColumn = newArray(1);
	
	//Use loop to create intensity columns - this system means
	//the script works for images with different numbers of channels
	for (i = 0; i < channels; i ++)
	{
		column = Array.slice(intensityColumns, i*numMaxima, ((i+1)*numMaxima));
		Table.setColumn("Channel " + i+1 + " Intensity", column);
	}
	column = newArray(1);
	
	
	//Save results

	saveResults(currentImage, "Nuclei Results", "nuclei", "csv");

	saveResults(currentImage, "Foci Results", "foci", "csv");

	
	//Close the tables
	close("Nuclei Results");
	close("Foci Results");
	
	//Save the ROIs - including nuclei, foci and their combination
	saveROIs(currentImage);
	
	//Tidy up
	run("Clear Results");
	roiManager("Deselect");
	roiManager("Delete");
	close(currentImage);
}


//Identify nuclei in nuclei channel and save as ROIs
function identifyNuclei(imageTitle, channel)
{
	//Get a single image to work with (no slices or channels)
	if (select_focus)
	{
		duplicateChannelAndFocus(imageTitle, channel);
	}
	else
	{
		duplicateChannelAndProjectStack(imageTitle, channel, nucleiProjectPreference);
	}
	
	//Get title of the new single frame image
	nucleiDuplicate = getTitle();
	
	//Subtract the background
	backgroundSubtractionString = "rolling=" + backgroundSubtractionBallPixelRadius;
	run("Subtract Background...", backgroundSubtractionString);
	
	//Segment the image to identify nuclei
	thresholdString = nucleiThresholdMethod + " dark no-reset";
	setAutoThreshold(thresholdString);
	run("Convert to Mask");
	
	//Fill in holes in the mask
	run("Fill Holes");
	
	//Run watershed to split overlapping or touching nuclei
	run("Watershed");
	
	//Find all particles fitting the user-defined size bounds,
	//and add to the ROI manager
	anaPartString = "size=[" + nucleiMin + " - " + nucleiMax + "] exclude add";
	run("Analyze Particles...", anaPartString);
	
	//Close the single frame image
	close(nucleiDuplicate);
}

//Identify foci as maxima (points of high intensity)
function identifyMaxima(imageTitle, channel)
{
	//Get single image to work with using a Z project
	duplicateChannelAndProjectStack(imageTitle, channel, maximaProjectPreference);
	
	//Get title of the new single frame image
	maximaDuplicate = getTitle();
	
	//Subtract the background
	backgroundSubtractionString = "rolling=" + backgroundSubtractionBallPixelRadius;
	run("Subtract Background...", backgroundSubtractionString);
	
	//Find maxima using the prominence value given by the user
	findMaximaString = "prominence=" + foci_prominence + " output=[Point Selection]";
	run("Find Maxima...", findMaximaString);
	
	//Add maxima point selection to the ROI manager
	roiManager("Add");
	
	//Close the single frame image
	close(maximaDuplicate);
}


//Extract a single channel from a multi-channel image and then do a Z project
function duplicateChannelAndProjectStack(imageTitle, channel, projectPreference)
{
	//Duplicate specific channel from image to extract it as a new image
	duplicateChannel(imageTitle, channel);
	
	//Do Z project to get single frame if the image is a Z-stack
	getDimensions(width, height, channels, slices, frames);
	if (slices > 1)
	{
		projectString = "projection=[" + projectPreference + "]";
		stackTitle = getTitle();
		run("Z Project...", projectString);
		zProjectTitle = getTitle();
		close(stackTitle);
		
		//Select the new image before finishing
		selectWindow(zProjectTitle);
	}
}

//Extract a single channel from a multi-channel image and then select most in focus Z-slice
function duplicateChannelAndFocus(imageTitle, channel)
{
	//Duplicate specific channel from image to extract it as a new image
	duplicateChannel(imageTitle, channel);
	getDimensions(width, height, channels, slices, frames);
	
	//Measure standard deviations to select most in focus Z-slice and duplicate this
	//as new image
	if (slices > 1)
	{
		stackTitle = getTitle();
		
		//Highest standard deviation observed so far
		maxStDev = 0;
		
		//Initial slice assumed to be most in focus before beginning loop
		focussedSlice = 1;
		
		//Cycle through slices
		for (i = 1; i <= slices; i++)
		{
			setSlice(i);
			stDev = getValue("StdDev");
			
			//If this slice's standard deviation is higher than the current
			//maximum, then update the maximum and the most focussed slice
			if (stDev > maxStDev)
			{
				maxStDev = stDev;
				focussedSlice = i;
			}
		}
		
		//Duplicate the most focussed slice
		duplicateSlice(stackTitle, focussedSlice);
		focussedTitle = getTitle();
		
		close(stackTitle);
		
		selectWindow(focussedTitle);
	}
}

//Duplicate single channel from an image to make new image
function duplicateChannel(imageTitle, channel)
{
	selectWindow(imageTitle);
	duplicationString = "duplicate channels=" + channel;
	run("Duplicate...", duplicationString);
}

//Duplicate single Z-slice from an image to make new image
function duplicateSlice(imageTitle, slice)
{
	selectWindow(imageTitle);
	duplicationString = "slices=" + slice;
	run("Make Substack...", duplicationString);

}

//Save results as CSV, with image name, appended with "_foci_per_nucleus"
function saveResults(imageTitle, table, fileName, extension)
{
	outCSVName = imageTitle + "_" + fileName + "." + extension;
	fullPath = output + "/" + outCSVName;
	selectWindow(table);
	saveAs(table, fullPath);
	close(outCSVName);
}

//Save ROIs in a zip file
function saveROIs(imageTitle)
{
	outROIName = output + "/" + imageTitle + "_ROIs.zip";
	run("Select None");
	roiManager("Save", outROIName);
}

//Turn a column from the results window into an ImageJ array
function getResultsAsArray(columnTitle)
{
	//Count number of rows in results table
	numResults = nResults;
	
	//Create empty array of correct length
	resultsArray = newArray(numResults);
	
	//Loop through results
	for (i = 0; i < numResults; i++)
	{
		//Write result in row i to element i of the array
		resultsArray[i] = getResult(columnTitle, i);
	}
	
	return resultsArray;
}