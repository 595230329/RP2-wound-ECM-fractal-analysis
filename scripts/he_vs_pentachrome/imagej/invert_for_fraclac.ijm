// FracLac (2015Sep090313a9330) treats black as foreground regardless of the
// ImageJ "Black background" option (log shows "changed inverted LUT on ...").
// This writes inverted copies of the masks (collagen = black) for FracLac only.
// The original masks (collagen = white) are kept for compare_masks.py.
input = getDirectory("Folder with original masks");
output = getDirectory("Empty folder for inverted copies");
list = getFileList(input);
setBatchMode(true);
for (i = 0; i < list.length; i++) {
    name = toLowerCase(list[i]);
    if (endsWith(name, ".tif") || endsWith(name, ".tiff")) {
        open(input + list[i]);
        run("Invert");
        saveAs("Tiff", output + list[i]);
        close();
    }
}
setBatchMode(false);
print("Inverted " + list.length + " files.");
