/**
 * export_masks.groovy
 *
 * Exports binary collagen masks for axis-aligned (unrotated) square ROIs.
 * Run with Automate > Run for project on both H&E and pentachrome images.
 *
 * ROIs:   annotations whose name contains '9M' (e.g. '9M06 NORMAL SKIN 1').
 * Output: <project>/masks/H&E/<ROI name>_mask.tif  or  <project>/masks/Penta/<ROI name>_mask.tif
 *         8-bit, collagen = 255, background = 0, at OUT_PIXEL_SIZE µm/px.
 * Detections left inside the ROIs from earlier runs are removed before classifying.
 */
import qupath.lib.images.servers.LabeledImageServer
import qupath.lib.regions.RegionRequest

// ===================== Settings =====================
def HE_KEYWORD        = 'H&E'                      // text found only in H&E image names
def HE_CLASSIFIER     = 'HE_Collagen_v1'
def PENTA_CLASSIFIER  = 'Pentachrome_Collagen_v1'  // yellow only = collagen
def COLLAGEN_CLASS    = 'collagen'                 // class name in both classifiers
double OUT_PIXEL_SIZE = 0.2431                     // µm/px, both stains
// ====================================================

def imageData = getCurrentImageData()
def server = imageData.getServer()
def imageName = getProjectEntry().getImageName()
def isHE = imageName.toUpperCase().contains(HE_KEYWORD.toUpperCase())
def classifierName = isHE ? HE_CLASSIFIER : PENTA_CLASSIFIER
println "${imageName} -> ${isHE ? 'H&E' : 'Pentachrome'} (${classifierName})"

// Named ROIs only
def rois = getAnnotationObjects().findAll { it.getName()?.toUpperCase()?.contains('9M') }
if (rois.isEmpty()) {
    println "  No named ROIs found - skipped."
    return
}

// Remove detections from earlier runs inside these ROIs
def oldDetections = getDetectionObjects().findAll { rois.contains(it.getParent()) }
if (!oldDetections.isEmpty()) {
    removeObjects(oldDetections, true)
    println "  Removed ${oldDetections.size()} old detections."
}

// Classify inside the ROIs
selectObjects(rois)
createDetectionsFromPixelClassifier(classifierName, 0, 0)
resetSelection()

// Build and save the masks
def outDir = buildFilePath(PROJECT_BASE_DIR, 'masks', isHE ? 'H&E' : 'Penta')
mkdirs(outDir)
double downsample = OUT_PIXEL_SIZE / server.getPixelCalibration().getAveragedPixelSize()

def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.BLACK)
    .downsample(downsample)
    .useDetections()
    .addLabel(COLLAGEN_CLASS, 255, ColorTools.WHITE)
    .multichannelOutput(false)
    .build()

rois.each { roi ->
    def request = RegionRequest.createInstance(labelServer.getPath(), downsample, roi.getROI())
    def path = buildFilePath(outDir, roi.getName() + '_mask.tif')
    writeImageRegion(labelServer, request, path)
    println "  Saved: ${roi.getName()}_mask.tif"
}
