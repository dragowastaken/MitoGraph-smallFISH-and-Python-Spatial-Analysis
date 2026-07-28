// Recursive batch split/extract selected channels from hyperstacks
// Select a parent directory. The macro recursively finds folders named
// "Hyperstacks_Grayscale" and extracts selected channels from TIFF hyperstacks
// directly inside each matching folder.
// Saves each selected channel into named folders:
//   Hyperstacks_Grayscale/Extracted_Channels/<your channel name>/
// Uses Split Channels instead of Duplicate to avoid display/border artifacts.

parentDir = getDirectory("Choose parent folder to search recursively");

Dialog.create("Extract channels");
Dialog.addString("Channels to extract, separated by commas:", "1");
Dialog.addMessage("Examples: 1 or 1,3,5");
Dialog.show();

channelText = Dialog.getString();
channelText = replace(channelText, " ", "");
channelsToExtract = split(channelText, ",");

// Ask for a descriptive output-folder name for each selected channel.
// These names are used instead of generic folders like "Channel 1".
Dialog.create("Name extracted channel folders");
for (nameIndex = 0; nameIndex < channelsToExtract.length; nameIndex++) {
    cText = channelsToExtract[nameIndex];
    if (cText == "")
        continue;
    Dialog.addString("Folder name for channel " + cText + ":", "Channel " + cText);
}
Dialog.addMessage("Avoid using these filename characters: \\ / : * ? \" < > |");
Dialog.show();

channelNames = newArray(channelsToExtract.length);
for (nameIndex = 0; nameIndex < channelsToExtract.length; nameIndex++) {
    cText = channelsToExtract[nameIndex];
    if (cText == "") {
        channelNames[nameIndex] = "";
        continue;
    }
    rawName = Dialog.getString();
    rawName = trimStringSafe(rawName);
    if (rawName == "")
        rawName = "Channel " + cText;
    cleanName = sanitizeFolderName(rawName);
    if (cleanName == "")
        cleanName = "Channel " + cText;
    channelNames[nameIndex] = cleanName;
}

processedFolders = 0;
processedFiles = 0;
skippedFiles = 0;

setBatchMode(true);
processDirectory(parentDir);
setBatchMode(false);

print("Done. Hyperstacks_Grayscale folders processed: " + processedFolders);
print("Files processed: " + processedFiles);
print("Files skipped: " + skippedFiles);

function processDirectory(dir) {
    if (isHyperstackFolder(dir) == 1) {
        processHyperstackFolder(dir);
        return;
    }

    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;

        if (File.isDirectory(path)) {
            cleanName = stripTrailingSlash(name);
            lowerName = toLowerCaseSafe(cleanName);

            // Do not recurse through previously generated output folders.
            if (lowerName == "extracted_channels")
                continue;
            if (startsWith(lowerName, "channel "))
                continue;

            processDirectory(path);
        }
    }
}

function processHyperstackFolder(inputDir) {
    print("Processing Hyperstacks_Grayscale folder: " + inputDir);

    outputRoot = inputDir + "Extracted_Channels/";
    File.makeDirectory(outputRoot);

    localProcessed = 0;
    localSkipped = 0;

    list = getFileList(inputDir);

    for (i = 0; i < list.length; i++) {
        filename = list[i];
        path = inputDir + filename;

        if (File.isDirectory(path))
            continue;

        lowerFilename = toLowerCaseSafe(filename);

        if (!(endsWith(lowerFilename, ".tif") || endsWith(lowerFilename, ".tiff")))
            continue;

        open(path);

        originalTitle = getTitle();
        getDimensions(width, height, channels, slices, frames);

        if (width < 1 || height < 1 || channels < 1) {
            print("  Skipping " + filename + ": invalid dimensions or channel count.");
            if (isOpen(originalTitle)) {
                selectWindow(originalTitle);
                close();
            }
            skippedFiles++;
            localSkipped++;
            continue;
        }

        // Force grayscale display before splitting.
        run("Grays");

        // Split all channels into separate grayscale images.
        run("Split Channels");

        // Save only requested channels.
        for (cIndex = 0; cIndex < channelsToExtract.length; cIndex++) {

            if (channelsToExtract[cIndex] == "")
                continue;

            c = parseInt(channelsToExtract[cIndex]);

            if (c < 1 || c > channels) {
                print("  Skipping channel " + c + " for " + originalTitle +
                      ": image has only " + channels + " channel(s).");
                continue;
            }

            channelDir = outputRoot + channelNames[cIndex] + "/";
            File.makeDirectory(channelDir);

            // Split Channels names windows like C1-filename.tif, C2-filename.tif, etc.
            splitTitle = "C" + c + "-" + originalTitle;

            if (isOpen(splitTitle)) {
                selectWindow(splitTitle);
                run("Grays");
                saveAs("Tiff", channelDir + originalTitle);
            } else {
                print("  Could not find split channel window: " + splitTitle);
            }
        }

        // Close all split channel windows from this image.
        for (c = 1; c <= channels; c++) {
            splitTitle = "C" + c + "-" + originalTitle;
            if (isOpen(splitTitle)) {
                selectWindow(splitTitle);
                close();
            }
        }

        // Close original if still open.
        if (isOpen(originalTitle)) {
            selectWindow(originalTitle);
            close();
        }

        processedFiles++;
        localProcessed++;
    }

    if (localProcessed > 0) {
        processedFolders++;
        print("  Extracted channels from " + localProcessed + " file(s). Output: " + outputRoot);
    } else {
        print("  No TIFF files processed in this folder.");
    }

    if (localSkipped > 0)
        print("  Skipped " + localSkipped + " file(s) in this folder.");
}

function isHyperstackFolder(dir) {
    folderName = getLastFolderName(dir);
    folderName = toLowerCaseSafe(folderName);
    if (folderName == "hyperstacks_grayscale")
        return 1;
    return 0;
}

function getLastFolderName(path) {
    p = stripTrailingSlash(path);
    lastSlash = -1;
    for (j = 0; j < lengthOf(p); j++) {
        ch = substring(p, j, j + 1);
        if (ch == "/" || ch == "\\")
            lastSlash = j;
    }
    if (lastSlash >= 0) {
        out = substring(p, lastSlash + 1, lengthOf(p));
        return out;
    }
    return p;
}

function stripTrailingSlash(s) {
    out = s;
    while (lengthOf(out) > 0) {
        last = substring(out, lengthOf(out) - 1, lengthOf(out));
        if (last == "/" || last == "\\")
            out = substring(out, 0, lengthOf(out) - 1);
        else
            break;
    }
    return out;
}

function toLowerCaseSafe(s) {
    // ImageJ macro supports toLowerCase() in current Fiji/ImageJ versions.
    out = toLowerCase(s);
    return out;
}

function trimStringSafe(s) {
    out = s;
    while (lengthOf(out) > 0) {
        first = substring(out, 0, 1);
        if (first == " " || first == "\t")
            out = substring(out, 1, lengthOf(out));
        else
            break;
    }
    while (lengthOf(out) > 0) {
        last = substring(out, lengthOf(out) - 1, lengthOf(out));
        if (last == " " || last == "\t")
            out = substring(out, 0, lengthOf(out) - 1);
        else
            break;
    }
    return out;
}

function sanitizeFolderName(s) {
    out = "";
    for (k = 0; k < lengthOf(s); k++) {
        ch = substring(s, k, k + 1);
        if (ch == "/" || ch == "\\" || ch == ":" || ch == "*" || ch == "?" || ch == "\"" || ch == "<" || ch == ">" || ch == "|")
            out = out + "_";
        else
            out = out + ch;
    }
    out = trimStringSafe(out);
    return out;
}
