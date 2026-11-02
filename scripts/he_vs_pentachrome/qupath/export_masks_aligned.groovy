/**
 * export_masks_aligned.groovy
 *
 * Exports H&E collagen masks for rotated square ROIs, resampled into the
 * pentachrome ROI frame (pentachrome = reference). Used for ROIs drawn on
 * pentachrome and propagated to the registered H&E section, where they appear rotated.
 * Run with Automate > Run for project on both stains:
 *   - pentachrome images (9M05-9M09): only save reference check images (<ROI>_penta.png)
 *   - H&E images (9M10-9M14): classify, then save the mask and an H&E check image (<ROI>_HE.png)
 *
 * Output: <project>/masks/H&E_aligned/<ROI name>_mask.tif   (8-bit, collagen = 255)
 *         <project>/masks/check_aligned/<ROI name>_HE.png / _penta.png
 *
 * START_BY_IMAGE sets which corner of the rotated ROI becomes the top-left corner,
 * so the H&E mask has the same orientation as the pentachrome mask. It was chosen
 * per H&E image by comparing the check images.
 * Sampling: nearest neighbour at output-pixel centres.
 */
import qupath.lib.images.servers.LabeledImageServer
import qupath.lib.regions.RegionRequest
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

// ===================== Settings =====================
def HE_CLASSIFIER     = 'HE_Collagen_v1'
def COLLAGEN_CLASS    = 'collagen'
double OUT_PIXEL_SIZE = 0.2431     // µm/px, same as pentachrome masks
double ROI_SIZE_UM    = 500        // side of the pentachrome ROI
def START_BY_IMAGE    = ['9M10': 1, '9M11': 0, '9M13': 1, '9M14': 0]   // H&E code -> start corner
// ====================================================

// Nearest-neighbour resampling of a rotated square into an n x n grid.
// (ox, oy): start corner; (exx, exy), (eyx, eyy): one output pixel along each side.
@groovy.transform.CompileStatic
static int[] resample(int[] src, int sw, int sh, int n,
                      double ox, double oy, double exx, double exy,
                      double eyx, double eyy, int fill) {
    int[] out = new int[n * n]
    for (int v = 0; v < n; v++) {
        for (int u = 0; u < n; u++) {
            double fu = u + 0.5d
            double fv = v + 0.5d
            int x = (int) Math.floor(ox + fu * exx + fv * eyx)
            int y = (int) Math.floor(oy + fu * exy + fv * eyy)
            out[v * n + u] = (x >= 0 && y >= 0 && x < sw && y < sh) ? src[y * sw + x] : fill
        }
    }
    return out
}

def imageData = getCurrentImageData()
def server = imageData.getServer()
def imageName = getProjectEntry().getImageName()
def isHE = (imageName.toUpperCase() =~ /9M1[0-4]/).find()
int START_CORNER = 0
def codeMatch = (imageName.toUpperCase() =~ /9M1[0-4]/)
if (codeMatch.find()) START_CORNER = START_BY_IMAGE.get(codeMatch.group(), 0)
println "${imageName} -> ${isHE ? 'H&E' : 'Pentachrome (check images only)'}, START_CORNER = ${START_CORNER}"

def rois = getAnnotationObjects().findAll { it.getName()?.toUpperCase()?.contains('9M') }
if (rois.isEmpty()) { println "  No named ROIs found - skipped."; return }

double fullPx = server.getPixelCalibration().getAveragedPixelSize()
double downsample = OUT_PIXEL_SIZE / fullPx
int N = (int) Math.round(ROI_SIZE_UM / OUT_PIXEL_SIZE)

def checkDir = buildFilePath(PROJECT_BASE_DIR, 'masks', 'check_aligned')
mkdirs(checkDir)

// ---------- Pentachrome: reference check images only ----------
if (!isHE) {
    rois.each { r ->
        def req = RegionRequest.createInstance(server.getPath(), downsample, r.getROI())
        ImageIO.write(server.readRegion(req), 'png',
            new File(buildFilePath(checkDir, r.getName() + '_penta.png')))
        println "  Saved check: ${r.getName()}_penta.png"
    }
    return
}

// ---------- H&E: classify inside the rotated ROIs ----------
def oldDet = getDetectionObjects().findAll { rois.contains(it.getParent()) }
if (!oldDet.isEmpty()) removeObjects(oldDet, true)
selectObjects(rois)
createDetectionsFromPixelClassifier(HE_CLASSIFIER, 0, 0)
resetSelection()

def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.BLACK)
    .downsample(downsample)
    .useDetections()
    .addLabel(COLLAGEN_CLASS, 255, ColorTools.WHITE)
    .multichannelOutput(false)
    .build()

def outDir = buildFilePath(PROJECT_BASE_DIR, 'masks', 'H&E_aligned')
mkdirs(outDir)

for (r in rois) {
    def roi = r.getROI()

    // Four distinct corners of the rotated square
    def c = []
    for (p in roi.getAllPoints()) {
        if (!c.any { Math.abs(it.getX() - p.getX()) < 1e-6 && Math.abs(it.getY() - p.getY()) < 1e-6 })
            c << p
    }
    if (c.size() != 4) { println "  ${r.getName()}: ${c.size()} corners (expected 4) - skipped"; continue }

    def p0 = c[START_CORNER % 4]
    def p1 = c[(START_CORNER + 1) % 4]
    def p3 = c[(START_CORNER + 3) % 4]
    double exx = (p1.getX() - p0.getX()) / N, exy = (p1.getY() - p0.getY()) / N
    double eyx = (p3.getX() - p0.getX()) / N, eyy = (p3.getY() - p0.getY()) / N
    if (exx * eyy - exy * eyx < 0) {   // corner order reversed: swap axes to avoid a mirror image
        double tx = exx, ty = exy
        exx = eyx; exy = eyy; eyx = tx; eyy = ty
    }
    double angle = Math.toDegrees(Math.atan2(exy, exx))
    double sideUm = Math.hypot(exx, exy) * N * fullPx

    // Region around the rotated square
    double pad = 4 * downsample
    int x0 = (int) Math.max(0, Math.floor(roi.getBoundsX() - pad))
    int y0 = (int) Math.max(0, Math.floor(roi.getBoundsY() - pad))
    int x1 = (int) Math.min(server.getWidth(),  Math.ceil(roi.getBoundsX() + roi.getBoundsWidth()  + pad))
    int y1 = (int) Math.min(server.getHeight(), Math.ceil(roi.getBoundsY() + roi.getBoundsHeight() + pad))
    double ox = (p0.getX() - x0) / downsample, oy = (p0.getY() - y0) / downsample
    double ax = exx / downsample, ay = exy / downsample
    double bx = eyx / downsample, by = eyy / downsample

    // Mask
    def lab = labelServer.readRegion(RegionRequest.createInstance(labelServer.getPath(), downsample, x0, y0, x1 - x0, y1 - y0))
    int sw = lab.getWidth(), sh = lab.getHeight()
    int[] src = lab.getRaster().getSamples(0, 0, sw, sh, 0, (int[]) null)
    int[] out = resample(src, sw, sh, N, ox, oy, ax, ay, bx, by, 0)
    def mask = new BufferedImage(N, N, BufferedImage.TYPE_BYTE_GRAY)
    mask.getRaster().setSamples(0, 0, N, N, 0, out)
    ImageIO.write(mask, 'tif', new File(buildFilePath(outDir, r.getName() + '_mask.tif')))

    // RGB check image in pentachrome orientation
    def rgb = server.readRegion(RegionRequest.createInstance(server.getPath(), downsample, x0, y0, x1 - x0, y1 - y0))
    int rw = rgb.getWidth(), rh = rgb.getHeight()
    int[] rsrc = rgb.getRGB(0, 0, rw, rh, null, 0, rw)
    int[] rout = resample(rsrc, rw, rh, N, ox, oy, ax, ay, bx, by, -1)
    def chk = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB)
    chk.setRGB(0, 0, N, N, rout, 0, N)
    ImageIO.write(chk, 'png', new File(buildFilePath(checkDir, r.getName() + '_HE.png')))

    println String.format("  Saved: %s_mask.tif | angle = %.1f deg | side = %.0f um", r.getName(), angle, sideUm)
}