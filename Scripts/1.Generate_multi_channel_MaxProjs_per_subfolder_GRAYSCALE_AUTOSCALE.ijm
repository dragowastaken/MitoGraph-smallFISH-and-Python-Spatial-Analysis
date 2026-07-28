// ============================================================================
// Generate one grayscale, autoscaled multi-channel MaxProjs.tif per folder
//
// Choose a parent directory. The macro recursively visits that folder and all
// subfolders. For every folder that directly contains input .tif/.tiff files,
// it creates that folder's own MaxProjs.tif from only those TIFFs.
//
// Output: <each folder containing TIFFs>/MaxProjs.tif
//
// Added behavior:
//   1. Output is displayed with grayscale LUTs.
//   2. Output display range is reset/autoscaled before saving.
//      - For 8-bit images: display min=0, max=255.
//      - For 16-bit images: display min/max are set from the smallest and
//        largest pixel values in the output channel's histogram.
//   3. The saved TIFF stores the grayscale LUT/display settings, so opening the
//      output in Fiji should not require manually changing LUT/autoscale.
// ============================================================================

rootFolder = getDirectory("Choose the parent directory containing TIFF folders");

print("Root folder: " + rootFolder);
foldersDone = processFolderTree(rootFolder, 0);

showMessage("Done", "Created grayscale autoscaled MaxProjs.tif in " + foldersDone + " folder(s).\n\nParent folder:\n" + rootFolder);

// ============================================================================
// Recursively visit folders. Each folder is processed independently.
// ============================================================================
function processFolderTree(dir, foldersDone) {
	// First process this folder's own TIFFs, if any.
	nDirect = countDirectTiffs(dir);
	if (nDirect > 0) {
		processOneFolder(dir, nDirect);
		foldersDone++;
	}

	// Then recurse into subfolders.
	list = getFileList(dir);
	for (i = 0; i < list.length; i++) {
		path = dir + list[i];
		if (File.isDirectory(path)) {
			foldersDone = processFolderTree(path, foldersDone);
		}
	}
	return foldersDone;
}

// ============================================================================
// Process only the TIFFs directly inside one folder.
// Does NOT combine TIFFs from child folders.
// ============================================================================
function processOneFolder(dir, nTiffs) {
	firstTiff = findFirstDirectTiff(dir);

	if (firstTiff == "")
		return;

	open(firstTiff);
	firstTitle = getTitle();
	getDimensions(w, h, cCount, zCount, tCount);
	close();

	if (w < 1 || h < 1 || cCount < 1) {
		showMessage("Could not read image dimensions",
			"ImageJ could not read valid dimensions from this TIFF:\n" +
			firstTiff + "\n\nWidth=" + w + ", height=" + h + ", channels=" + cCount +
			"\n\nTry opening this file manually in Fiji/ImageJ to confirm it imports correctly.");
		exit();
	}

	print("------------------------------------------------------------");
	print("Processing folder: " + dir);
	print("TIFFs in this folder: " + nTiffs);
	print("First TIFF: " + firstTiff);
	print("Width x height: " + w + " x " + h);
	print("Channels: " + cCount);

	// If an older MaxProjs window is still open for any reason, close it.
	if (isOpen("MaxProjs")) {
		selectWindow("MaxProjs");
		close();
	}

	// Create one output hyperstack for this folder:
	// channels = source channels; slices = one max projection per TIFF.
	newImage("MaxProjs", "16-bit black", w, h, cCount, nTiffs, 1);

	processed = 0;
	list = getFileList(dir);
	for (i = 0; i < list.length; i++) {
		name = list[i];
		path = dir + name;

		if (!File.isDirectory(path) && isInputTiff(name)) {
			processed++;

			open(path);
			origTitle = getTitle();
			getDimensions(currW, currH, currC, currZ, currT);

			print("  " + processed + "/" + nTiffs + ": " + name);

			if (currW != w || currH != h || currC != cCount) {
				close();
				selectWindow("MaxProjs");
				close();
				showMessage("Dimension mismatch",
					"This file does not match the first TIFF in the folder and cannot be added:\n" +
					path + "\n\nExpected: " + w + " x " + h + ", channels=" + cCount +
					"\nFound: " + currW + " x " + currH + ", channels=" + currC);
				exit();
			}

			if (currZ < 1)
				currZ = 1;

			for (c = 1; c <= cCount; c++) {
				selectWindow(origTitle);

				// Duplicate this channel across all z slices.
				run("Duplicate...", "title=Temp_Channel duplicate channels=" + c + " slices=1-" + currZ + " frames=1");

				// Max project this duplicated channel.
				run("Z Project...", "start=1 stop=" + currZ + " projection=[Max Intensity]");

				// Copy projection.
				selectWindow("MAX_Temp_Channel");
				run("Copy");
				close();

				// Close duplicated channel stack.
				selectWindow("Temp_Channel");
				close();

				// Paste projection into the correct channel and output slice.
				selectWindow("MaxProjs");
				Stack.setChannel(c);
				Stack.setSlice(processed);
				run("Paste");
				setMetadata("Label", name);
			}

			selectWindow(origTitle);
			close();
		}
	}

	selectWindow("MaxProjs");

	// Make the output open as grayscale and autoscaled, equivalent to using
	// Image > Lookup Tables > Grays and then Reset in Brightness/Contrast.
	applyGrayscaleAndResetAutoscale("MaxProjs", cCount, nTiffs);

	outFile = dir + "MaxProjs.tif";
	print("Saving grayscale autoscaled folder output to: " + outFile);
	saveAs("Tiff", outFile);
	close();
}

// ============================================================================
// Apply grayscale LUT and reset-style display autoscaling.
// ============================================================================
function applyGrayscaleAndResetAutoscale(title, cCount, nSlices) {
	selectWindow(title);


	for (c = 1; c <= cCount; c++) {
		selectWindow(title);
		Stack.setChannel(c);
		Stack.setSlice(1);

		// Force grayscale LUT for this channel.
		run("Grays");

		// Reset-style autoscale:
		// 8-bit: min=0, max=255.
		// 16-bit: min/max equal smallest/largest pixel values in this channel
		// across all output slices.
		if (bitDepth() == 8) {
			setMinAndMax(0, 255);
		} else {
			channelMin = 1/0;
			channelMax = -1/0;

			for (s = 1; s <= nSlices; s++) {
				Stack.setChannel(c);
				Stack.setSlice(s);
				getStatistics(area, mean, min, max, std);

				if (min < channelMin)
					channelMin = min;
				if (max > channelMax)
					channelMax = max;
			}

			Stack.setChannel(c);
			Stack.setSlice(1);

			if (channelMax > channelMin) {
				setMinAndMax(channelMin, channelMax);
			} else {
				// Fallback for empty/constant images.
				setMinAndMax(0, 65535);
			}

			print("  Channel " + c + " grayscale autoscale min=" + channelMin + ", max=" + channelMax);
		}
	}

	selectWindow(title);
	Stack.setChannel(1);
	Stack.setSlice(1);
}

// ============================================================================
// Direct-folder helpers
// ============================================================================
function findFirstDirectTiff(dir) {
	list = getFileList(dir);
	for (i = 0; i < list.length; i++) {
		path = dir + list[i];
		if (!File.isDirectory(path) && isInputTiff(list[i]))
			return path;
	}
	return "";
}

function countDirectTiffs(dir) {
	count = 0;
	list = getFileList(dir);
	for (i = 0; i < list.length; i++) {
		path = dir + list[i];
		if (!File.isDirectory(path) && isInputTiff(list[i]))
			count++;
	}
	return count;
}

function isInputTiff(name) {
	lower = toLowerCase(name);

	// Skip existing output files if the macro is rerun.
	if (lower == "maxprojs.tif" || lower == "maxprojs.tiff")
		return false;

	return endsWith(lower, ".tif") || endsWith(lower, ".tiff");
}
