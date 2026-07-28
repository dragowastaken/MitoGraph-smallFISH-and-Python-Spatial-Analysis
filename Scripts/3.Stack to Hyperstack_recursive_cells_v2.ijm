// Batch convert stacks to grayscale hyperstacks with user-specified dimensions
// Recursive version: choose a parent folder, then process every subfolder named "cells".
// Each cells folder gets its own Hyperstacks_Grayscale output folder.

parentDir = getDirectory("Choose parent folder containing subfolders with cells folders");

// Ask for hyperstack dimensions once, then apply to every cells folder
Dialog.create("Hyperstack settings");
Dialog.addNumber("Channels:", 5);
Dialog.addNumber("Slices (Z):", 21);
Dialog.addNumber("Frames (T):", 1);
Dialog.addChoice("Stack order:", newArray("xyczt(default)", "xyzct", "xyztc"), "xyczt(default)");
Dialog.show();

channels = Dialog.getNumber();
slices = Dialog.getNumber();
frames = Dialog.getNumber();
order = Dialog.getChoice();
expected = channels * slices * frames;

processedFolders = 0;
processedFiles = 0;
skippedFiles = 0;

setBatchMode(true);
processTree(parentDir);
setBatchMode(false);

print("Done.");
print("Cells folders processed: " + processedFolders);
print("Files converted: " + processedFiles);
print("Files skipped: " + skippedFiles);

function processTree(dir) {
    list = getFileList(dir);

    // If this directory itself is named cells, process it and do not recurse into its output folder.
    if (isCellsFolder(dir)) {
        processCellsFolder(dir);
    }

    // Continue searching subfolders, but skip output folders.
    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;

        if (File.isDirectory(path)) {
            cleanName = stripTrailingSlash(name);
            lowerName = toLowerCase(cleanName);

            if (lowerName == "hyperstacks_grayscale")
                continue;

            processTree(path);
        }
    }
}

function processCellsFolder(inputDir) {
    outputDir = inputDir + "Hyperstacks_Grayscale/";
    File.makeDirectory(outputDir);

    print("Processing cells folder: " + inputDir);
    processedFolders++;

    list = getFileList(inputDir);

    for (j = 0; j < list.length; j++) {
        filename = list[j];
        path = inputDir + filename;

        if (File.isDirectory(path))
            continue;

        lowerFilename = toLowerCase(filename);

        if (endsWith(lowerFilename, ".tif") || endsWith(lowerFilename, ".tiff")) {
            // Do not re-process anything already in an output folder; this loop only sees direct children,
            // but keep this check in case a file is named like an output artifact.
            if (startsWith(lowerFilename, "hyperstack")) {
                // This is only a filename check, not required for normal output folders.
            }

            open(path);
            title = getTitle();
            n = nSlices;

            if (n != expected) {
                print("  Skipping " + title + ": stack has " + n + " slices, expected " + expected);
                close();
                skippedFiles++;
                continue;
            }

            run("Stack to Hyperstack...",
                "order=" + order +
                " channels=" + channels +
                " slices=" + slices +
                " frames=" + frames +
                " display=Grayscale");

            // Ensure it is not saved as a composite color image
            run("Grays");

            saveAs("Tiff", outputDir + title);
            close();
            processedFiles++;
        }
    }
}

function isCellsFolder(dir) {
    d = dir;
    // Return numeric 1/0 for use inside if().
    if (endsWith(d, "/") || endsWith(d, "\\\\"))
        d = substring(d, 0, lengthOf(d) - 1);

    slash1 = lastIndexOf(d, "/");
    slash2 = lastIndexOf(d, "\\\\");
    slash = slash1;
    if (slash2 > slash)
        slash = slash2;

    if (slash >= 0)
        base = substring(d, slash + 1);
    else
        base = d;

    if (toLowerCase(base) == "cells")
        return 1;
    else
        return 0;
}

function stripTrailingSlash(name) {
    out = name;
    while (endsWith(out, "/") || endsWith(out, "\\\\"))
        out = substring(out, 0, lengthOf(out) - 1);
    return out;
}
