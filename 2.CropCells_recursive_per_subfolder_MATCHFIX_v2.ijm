// ==========================================================================
// CropCells_recursive_per_subfolder.ijm
// Adapted from Matheus Viana's CropCells.ijm.
//
// Select a parent directory. The macro recursively checks that directory and
// all subfolders. Each folder that contains both RoiSet.zip and MaxProjs.tif
// is processed independently, preserving the original CropCells behavior:
//   - crops cells from the original TIFF z-stacks using ROIs from RoiSet.zip
//   - saves cropped cell stacks into a "cells" subfolder inside that same folder
// ===========================================================================

// Defining the size in pixels of the single cell z-stacks
_xy = getNumber("Single cell image size in pixels", 200);

_RootFolder = getDirectory("Choose a parent directory");

setBatchMode(true);

processedFolders = 0;
processedCells = 0;
skippedFolders = 0;

run("ROI Manager...");
processDirectoryRecursive(_RootFolder);

setBatchMode(false);

print("Cell cropping is complete.");
print("Folders processed: " + processedFolders);
print("Cells cropped: " + processedCells);
print("Folders skipped because they lacked RoiSet.zip or MaxProjs.tif: " + skippedFolders);


function processDirectoryRecursive(dir) {
    // Process this directory if it contains the expected files.
    processOneFolder(dir);

    // Recurse into subdirectories, except the output folder.
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;
        if (File.isDirectory(path)) {
            cleanName = stripTrailingSlash(name);
            if (cleanName != "cells") {
                processDirectoryRecursive(ensureTrailingSlash(path));
            }
        }
    }
}


function processOneFolder(dir) {
    roiPath = dir + "RoiSet.zip";
    maxPath = dir + "MaxProjs.tif";

    if (!File.exists(roiPath) || !File.exists(maxPath)) {
        skippedFolders++;
        return;
    }

    print("Processing folder: " + dir);

    File.makeDirectory(dir + "cells");

    roiManager("Reset");
    roiManager("Open", roiPath);

    if (roiManager("count") == 0) {
        print("  Skipping: RoiSet.zip contains no ROIs.");
        return;
    }

    open(maxPath);
    MAXP = getImageID;

    // For each ROI (cell)
    for (roi = 0; roi < roiManager("count"); roi++) {
        roiManager("Select", roi);
        _FileName = getInfo("slice.label");
        _FileName = cleanTiffBaseName(_FileName);

        // Match the ROI/MaxProjs slice label back to the original TIFF.
        // Some ND2-derived TIFFs have labels like:
        //   "sample.nd2 - sample.nd2"
        // while the actual TIFF file is:
        //   "sample.nd2 - sample.nd2 (series 03).tif"
        // Therefore, use exact matching first, then prefix matching, then
        // fall back to the only TIFF in the folder if there is exactly one.
        originalPath = findOriginalTiffForLabel(dir, _FileName);
        if (originalPath == "") {
            print("  Skipping ROI " + roi + ": could not find original TIFF for label '" + _FileName + "'.");
            continue;
        }

        open(originalPath);
        ORIGINAL = getImageID;

        run("Restore Selection");

        newImage("CELL", "16-bit Black", _xy, _xy, nSlices);
        CELL = getImageID;

        // Estimating the noise distribution around the ROI
        max_ai = 0;
        slice_max_ai = 1;
        for (s = 1; s <= nSlices; s++) {
            selectImage(MAXP);

            selectImage(ORIGINAL);
            setSlice(s);
            run("Restore Selection");
            run("Make Band...", "band=5");
            getStatistics(area, mean, min, max, std);
            run("Restore Selection");
            run("Copy");

            selectImage(CELL);
            setSlice(s);
            run("Select None");
            run("Add...", "value=" + mean + " slice");
            run("Add Specified Noise...", "slice standard=" + 0.5 * std);
            run("Paste");

            getStatistics(area, mean, min, max, std);
            if (mean > max_ai) {
                max_ai = mean;
                slice_max_ai = s;
            }
        }

        run("Select None");
        resetMinAndMax();

        save(dir + "cells/" + _FileName + "_" + IJ.pad(roi, 3) + ".tif");
        processedCells++;

        selectImage(CELL); close();
        selectImage(ORIGINAL); close();
    }

    selectImage(MAXP); close();
    processedFolders++;
}



function findOriginalTiffForLabel(dir, label) {
    labelBase = trimString(cleanTiffBaseName(label));
    labelLower = toLowerCase(labelBase);

    exactTif = dir + labelBase + ".tif";
    if (File.exists(exactTif))
        return exactTif;

    exactTiff = dir + labelBase + ".tiff";
    if (File.exists(exactTiff))
        return exactTiff;

    list = getFileList(dir);
    nTiffs = 0;
    onlyTiff = "";

    // Prefix/contains matching is needed for Bio-Formats series exports.
    for (j = 0; j < list.length; j++) {
        name = list[j];
        path = dir + name;

        if (!File.isDirectory(path) && isInputTiff(name)) {
            nTiffs++;
            onlyTiff = path;

            fileBase = trimString(cleanTiffBaseName(name));
            fileLower = toLowerCase(fileBase);

            if (fileLower == labelLower)
                return path;

            // Handles: label="...nd2 - ...nd2" and file="...nd2 - ...nd2 (series 03)"
            if (startsWith(fileLower, labelLower))
                return path;

            // Handles the reverse case if labels include text not present in filenames.
            if (startsWith(labelLower, fileLower))
                return path;

            // Last-resort substring match, still restricted to TIFFs in this exact folder.
            if (indexOf(fileLower, labelLower) >= 0 || indexOf(labelLower, fileLower) >= 0)
                return path;
        }
    }

    // If this Series folder contains only one raw TIFF, use it.
    if (nTiffs == 1)
        return onlyTiff;

    return "";
}

function isInputTiff(name) {
    lower = toLowerCase(name);
    if (lower == "maxprojs.tif" || lower == "maxprojs.tiff")
        return false;
    if (endsWith(lower, ".tif")) return true;
    if (endsWith(lower, ".tiff")) return true;
    return false;
}

function trimString(s) {
    while (lengthOf(s) > 0 && substring(s, 0, 1) == " ")
        s = substring(s, 1, lengthOf(s));

    while (lengthOf(s) > 0 && substring(s, lengthOf(s) - 1, lengthOf(s)) == " ")
        s = substring(s, 0, lengthOf(s) - 1);

    return s;
}

function cleanTiffBaseName(label) {
    // File.getName safely removes folder components without regex parsing.
    label = File.getName(label);

    lower = toLowerCase(label);
    if (endsWith(lower, ".tiff")) {
        label = substring(label, 0, lengthOf(label) - 5);
    } else if (endsWith(lower, ".tif")) {
        label = substring(label, 0, lengthOf(label) - 4);
    }

    trimmed = trimString(label);
    return trimmed;
}


function ensureTrailingSlash(path) {
    if (endsWith(path, "/") || endsWith(path, "\\\\")) {
        return path;
    }
    return path + File.separator;
}


function stripTrailingSlash(path) {
    while (endsWith(path, "/") || endsWith(path, "\\\\")) {
        path = substring(path, 0, lengthOf(path) - 1);
    }
    return path;
}
