// Recursive Subtract Background - OVERWRITE ORIGINALS
// WARNING: This macro overwrites the original TIFF files.
// Make a backup before running.

macro "Recursive Subtract Background Overwrite Originals" {

    parentDir = getDirectory("Choose parent folder");
    if (parentDir == "") exit("No folder selected.");

    Dialog.create("Subtract Background Settings");
    Dialog.addNumber("Rolling ball radius (pixels):", 30);
    Dialog.addCheckbox("Light background", false);
    Dialog.addCheckbox("Separate colors", false);
    Dialog.addCheckbox("Create background only", false);
    Dialog.addCheckbox("Sliding paraboloid", false);
    Dialog.addCheckbox("Disable smoothing", false);
    Dialog.addCheckbox("Process stack slices", true);
    Dialog.addMessage("WARNING: This version overwrites the original TIFF files.");
    Dialog.show();

    radius = Dialog.getNumber();
    lightBackground = Dialog.getCheckbox();
    separateColors = Dialog.getCheckbox();
    createBackground = Dialog.getCheckbox();
    slidingParaboloid = Dialog.getCheckbox();
    disableSmoothing = Dialog.getCheckbox();
    processStack = Dialog.getCheckbox();

    processed = 0;
    skipped = 0;

    processFolder(parentDir);

    print("Recursive background subtraction complete.");
    print("Files overwritten: " + processed);
    print("Files skipped: " + skipped);
}


function processFolder(dir) {
    list = getFileList(dir);

    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;

        if (File.isDirectory(path)) {
            cleanName = replace(name, "/", "");
            cleanName = replace(cleanName, "\\", "");

            // Skip output folders from older versions if present.
            if (cleanName == "Background_Subtracted") {
                continue;
            }

            processFolder(path);
        }
        else {
            lower = toLowerCase(name);

            if (endsWith(lower, ".tif") || endsWith(lower, ".tiff")) {

                // Optional safety skip: avoid reprocessing files already marked as bgsub.
                // Delete this block if you truly want to reprocess every TIFF.
                if (indexOf(lower, "_bgsub") >= 0) {
                    skipped++;
                    continue;
                }

                processImageOverwrite(path);
            }
        }
    }
}


function processImageOverwrite(path) {
    open(path);
    title = getTitle();

    command = "rolling=" + radius;

    if (lightBackground) command = command + " light";
    if (separateColors) command = command + " separate";
    if (createBackground) command = command + " create";
    if (slidingParaboloid) command = command + " sliding";
    if (disableSmoothing) command = command + " disable";
    if (processStack) command = command + " stack";

    run("Subtract Background...", command);

    // OVERWRITE original file
    saveAs("Tiff", path);

    close();
    processed++;
}
