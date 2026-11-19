/**
 * export_masks_aligned_penta.groovy
 *
 * The reverse of export_masks_aligned.groovy: H&E is the reference. Used for
 * wound ROIs that were redrawn on H&E with anatomical criteria and propagated to
 * pentachrome, where they appear rotated. The pentachrome mask is resampled into
 * the H&E ROI frame. NORMAL SKIN ROIs are excluded so their masks are not overwritten.
 * Run with Automate > Run for project on both stains:
 *   - H&E images (9M10-9M14): only save reference check images (<ROI>_HE.png)
 *   - pentachrome images: classify, then save the mask and a pentachrome check image (<ROI>_penta.png)
 *
 * Output: <project>/masks/Penta_aligned/<ROI name>_mask.tif   (8-bit, collagen = 255)
 *         <project>/masks/check_aligned_penta/<ROI name>_HE.png / _penta.png
 *
 * Start corner: for each ROI, the corner whose first edge points most nearly along +x
 * becomes the top-left corner, so the pentachrome mask keeps the orientation of the H&E
 * ROI. START_OVERRIDE fixes the corner for an image instead (e.g. a section mounted
 * rotated by 90 deg). assemble_masks.py checks the orientation of every image.
 * Sampling: nearest neighbour at output-pixel centres.
 *
 * Each saved ROI is also logged to <project>/masks/roi_angles.csv (appended; the
 * last row per ROI and stain is used by assemble_masks.py).
 */
import qupath.lib.images.servers.LabeledImageServer
import qupath.lib.regions.RegionRequest
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

// ===================== Settings =====================
def PENTA_CLASSIFIER  = 'Pentachrome_Collagen_v1'  // yellow only = collagen
def COLLAGEN_CLASS    = 'collagen'
double OUT_PIXEL_SIZE = 0.2431                     // µm/px, same as other masks
double ROI_SIZE_UM    = 500
def START_OVERRIDE    = [:]   // pentachrome code -> fixed start corner, e.g. ['9M09': 0]; empty = automatic
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
def codeMatch = (imageName.toUpperCase() =~ /9M0[5-9]/)
def overrideCorner = codeMatch.find() ? START_OVERRIDE.get(codeMatch.group()) : null
println "${imageName} -> ${isHE ? 'H&E (check images only)' : 'Pentachrome'}, start corner ${overrideCorner == null ? 'automatic' : overrideCorner}"

// Wound ROIs only
def rois = getAnnotationObjects().findAll {
    def n = it.getName()?.toUpperCase()
    n?.contains('9M') && !n.contains('NORMAL')
}
if (rois.isEmpty()) { println "  No wound ROIs found - skipped."; return }

double fullPx = server.getPixelCalibration().getAveragedPixelSize()
double downsample = OUT_PIXEL_SIZE / fullPx
int N = (int) Math.round(ROI_SIZE_UM / OUT_PIXEL_SIZE)

def checkDir = buildFilePath(PROJECT_BASE_DIR, 'masks', 'check_aligned_penta')
mkdirs(checkDir)

// ---------- H&E: reference check images only ----------
if (isHE) {
    rois.each { r ->
        def req = RegionRequest.createInstance(server.getPath(), downsample, r.getROI())
        ImageIO.write(server.readRegion(req), 'png',
            new File(buildFilePath(checkDir, r.getName() + '_HE.png')))
        println "  Saved check: ${r.getName()}_HE.png"
    }
    return
}

// ---------- Pentachrome: classify inside the rotated ROIs ----------
def oldDet = getDetectionObjects().findAll { rois.contains(it.getParent()) }
if (!oldDet.isEmpty()) removeObjects(oldDet, true)
selectObjects(rois)
createDetectionsFromPixelClassifier(PENTA_CLASSIFIER, 0, 0)
resetSelection()

def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.BLACK)
    .downsample(downsample)
    .useDetections()
    .addLabel(COLLAGEN_CLASS, 255, ColorTools.WHITE)
    .multichannelOutput(false)
    .build()

def outDir = buildFilePath(PROJECT_BASE_DIR, 'masks', 'Penta_aligned')
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

    // Start corner: the corner whose first edge points most nearly along +x, so the output
    // keeps the orientation of the reference ROI. START_OVERRIDE fixes it for an image instead.
    int startCorner = 0
    double exx = 0, exy = 0, eyx = 0, eyy = 0, bestAbs = Double.MAX_VALUE
    for (int s = 0; s < 4; s++) {
        if (overrideCorner != null && s != overrideCorner) continue
        def q0 = c[s], q1 = c[(s + 1) % 4], q3 = c[(s + 3) % 4]
        double ax_ = (q1.getX() - q0.getX()) / N, ay_ = (q1.getY() - q0.getY()) / N
        double bx_ = (q3.getX() - q0.getX()) / N, by_ = (q3.getY() - q0.getY()) / N
        if (ax_ * by_ - ay_ * bx_ < 0) {   // corner order reversed: swap axes to avoid a mirror image
            double tx = ax_, ty = ay_
            ax_ = bx_; ay_ = by_; bx_ = tx; by_ = ty
        }
        double a = Math.abs(Math.toDegrees(Math.atan2(ay_, ax_)))
        if (a < bestAbs) { bestAbs = a; startCorner = s; exx = ax_; exy = ay_; eyx = bx_; eyy = by_ }
    }
    def p0 = c[startCorner]
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

    // RGB check image in H&E orientation
    def rgb = server.readRegion(RegionRequest.createInstance(server.getPath(), downsample, x0, y0, x1 - x0, y1 - y0))
    int rw = rgb.getWidth(), rh = rgb.getHeight()
    int[] rsrc = rgb.getRGB(0, 0, rw, rh, null, 0, rw)
    int[] rout = resample(rsrc, rw, rh, N, ox, oy, ax, ay, bx, by, -1)
    def chk = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB)
    chk.setRGB(0, 0, N, N, rout, 0, N)
    ImageIO.write(chk, 'png', new File(buildFilePath(checkDir, r.getName() + '_penta.png')))

    println String.format("  Saved: %s_mask.tif | angle = %.1f deg | side = %.0f um | start corner %d", r.getName(), angle, sideUm, startCorner)

    // Angle log read by assemble_masks.py to decide which stain is the reference for each ROI
    def angleLog = new File(buildFilePath(PROJECT_BASE_DIR, 'masks', 'roi_angles.csv'))
    if (!angleLog.exists()) angleLog.text = 'Stain,Image,ROI,Angle_deg,Side_um\n'
    angleLog << String.format('%s,"%s","%s",%.3f,%.1f%n', 'Penta', imageName, r.getName(), angle, sideUm)
}
