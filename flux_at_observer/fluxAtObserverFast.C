#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TThread.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TGeoVolume.h"
#include "TGLViewer.h"
#include "TPad.h"

#include "AFocalSurface.h"
#include "AGeoUtil.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ABorderSurfaceCondition.h"
#include "ARayShooter.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
#include <queue>
#include <condition_variable>

const double cm = AOpticsManager::cm();
const double nm = AOpticsManager::nm();
const double THETA_MAX = 170.;     // governs size of exit port
const int MAX_REFLECTIONS = 50000;

const double INNER_RADIUS = 100.1*cm;
const double OUTER_RADIUS = 101*cm;
const double REFLECTANCE = 0.99;
const double ROUGHNESS = 0.01;

struct Detector {
    double x, y, z;           // position
    double nx, ny, nz;        // normal vector
    double width, height;     // dimensions
    int hitCount;             // number of rays that hit the detector
    
    Detector(double w = 10*cm, double h = 10*cm) 
        : width(w), height(h), hitCount(0) {}
    
    /* Sets the detector position and orientation using spherical coordinates
     * @param theta: polar angle from negative z-axis (in degrees)
     *              0° = on axis with exit port
     *              90° = perpendicular to exit port axis
     * @param phi: azimuthal angle in x-y plane (in degrees)
     *            0° = +x axis
     *            90° = +y axis
     * @param radius: distance from exit port center
     */
    void setPosition(double theta, double phi, double radius) {
        // Convert spherical to cartesian coordinates
        double theta_rad = theta * M_PI / 180.0;
        double phi_rad = phi * M_PI / 180.0;
        
        // Calculate position relative to world origin
        x = radius * sin(theta_rad) * cos(phi_rad);
        y = radius * sin(theta_rad) * sin(phi_rad);
        z = -100*cm - radius * cos(theta_rad);  // Offset by exit port position
        
        // Calculate normal vector pointing to (0,0,-100)
        double dx = x - 0;
        double dy = y - 0;
        double dz = z - (-100*cm);
        double mag = sqrt(dx*dx + dy*dy + dz*dz);
        
        nx = -dy/mag;  // Negative because we want to point towards the exit port
        ny = dx/mag;
        nz = dz/mag;
    }
    
    bool checkIntersection(const ARay* ray) {
        Double_t lastPoint[3];
        Double_t direction[3];
        ray->GetLastPoint(lastPoint);
        ray->GetDirection(direction);
        
        // Calculate dot product of ray direction with detector normal
        double dot = direction[0]*nx + direction[1]*ny + direction[2]*nz;
        if (fabs(dot) < 1e-10) return false;  // Ray is parallel to detector plane
        
        // Calculate intersection with detector plane
        double dx = lastPoint[0] - x;  // Changed sign: point - center
        double dy = lastPoint[1] - y;
        double dz = lastPoint[2] - z;
        
        // Distance along ray to intersection point
        double t = -(dx*nx + dy*ny + dz*nz) / dot;  // Added negative sign
        
        // Calculate intersection point
        double ix = lastPoint[0] + direction[0] * t;
        double iy = lastPoint[1] + direction[1] * t;
        double iz = lastPoint[2] + direction[2] * t;
        
        // Vector from detector center to intersection point
        double rx = ix - x;
        double ry = iy - y;
        double rz = iz - z;
        
        // Compute distance in detector plane
        // Using cross product with normal to get perpendicular components
        double ux = ny*rz - nz*ry;  // First perpendicular direction
        double uy = nz*rx - nx*rz;
        double uz = nx*ry - ny*rx;
        
        double r2 = ux*ux + uy*uy + uz*uz;  // Square of distance in detector plane
        
        return r2 <= (width/2)*(width/2);
    }

    TGeoVolume* CreateGeometry() const {
        // Create a disk instead of a box
        TGeoTube* detDisk = new TGeoTube("detDisk", 
            0,           // inner radius
            width/2,     // outer radius (using width as diameter)
            1*cm         // half-length in z (thickness)
        );
        TGeoVolume* detVol = new TGeoVolume("detector", detDisk);
        detVol->SetLineColor(kBlue);
        detVol->SetLineWidth(2);
        detVol->SetFillColor(kCyan);
        detVol->SetTransparency(20);
        return detVol;
    }
    
    void AddToGeometry(TGeoVolume* top) const {
        TGeoVolume* detVol = CreateGeometry();
        
        // Calculate rotation matrix to point detector normal to exit port
        // We need to rotate from the default orientation (normal along z-axis)
        // to point towards (0,0,-100)
        
        // First convert normal vector to spherical coordinates
        double theta_rot = acos(-nz) * 180/M_PI;  // Note the negative nz
        double phi_rot = atan2(-ny, -nx) * 180/M_PI;  // Note the negatives
        
        // Create rotation that first rotates around Y axis (theta)
        // then around Z axis (phi)
        TGeoRotation* rot = new TGeoRotation();
        rot->SetAngles(phi_rot, theta_rot, 0);
        
        TGeoTranslation* trans = new TGeoTranslation(x, y, z);
        TGeoCombiTrans* transform = new TGeoCombiTrans(*trans, *rot);
        
        top->AddNode(detVol, 1, transform);
    }
};

class RayPool {
private:
    std::vector<ARay*> rays;
    size_t currentIndex;
    
public:
    RayPool(size_t size = 1000) : currentIndex(0) {
        rays.reserve(size);
        for (size_t i = 0; i < size; i++) {
            rays.push_back(new ARay(i, 400*nm, 0, 0, 0, 0, 0, 0, 0));
        }
    }
    
    ~RayPool() {
        for (auto ray : rays) {
            delete ray;
        }
    }
    
    ARay* getRay() {
        if (currentIndex >= rays.size()) {
            size_t i = rays.size();
            rays.push_back(new ARay(i, 400*nm, 0, 0, 0, 0, 0, 0, 0));
        }
        
        return rays[currentIndex++];
    }
    
    void reset() {
        currentIndex = 0;
    }
};

void setupOpticsManager(AOpticsManager* manager, int maxReflections = MAX_REFLECTIONS, 
                         double roughness = ROUGHNESS, double reflectance = REFLECTANCE, double thetaMax = THETA_MAX,
                         bool drawRays = false) {
    // Set ray trace limit
    manager->SetLimit(maxReflections);
    
    // Create slightly larger world box to ensure ray origins are inside
    TGeoBBox* box = new TGeoBBox("box", 300*cm, 300*cm, 300*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    
    // Create spherical mirror with exit port
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", INNER_RADIUS, OUTER_RADIUS, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    mirror->SetReflectance(reflectance);
    
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(true);
    condition->SetGaussianRoughness(roughness);
    
    world->AddNode(mirror, 1);
    
    // Use fewer segments when not visualizing
    manager->SetNsegments(drawRays ? 60: 20);
    
    // Close geometry
    manager->CloseGeometry();
    
    // Manually voxelize if not drawing rays
    if (!drawRays) {
        world->Voxelize("");
    }
    
    // Verify geometry is properly closed
    if (!manager->IsClosed()) {
        std::cerr << "Warning: Geometry not properly closed!" << std::endl;
        manager->CloseGeometry();
    }
}

bool isRayPassingThroughExitPort(ARay* ray, double exitPortZ) {
    Double_t lastPoint[3];
    ray->GetLastPoint(lastPoint);
    return lastPoint[2] < exitPortZ;
}

// Fix traceRays function to update detector's hit count
int traceRays(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false, int maxPoints = MAX_REFLECTIONS) {
    // Use the hardcoded maxPoints parameter that matches what you set with SetLimit()
    // instead of trying to retrieve it from the manager
    
    std::vector<TPolyLine3D*> rays;
    if (drawRays) {
        rays.reserve(n);
    }
    int hitCount = 0;
    
    for (int i = 0; i < n; ++i) {
        // Ray parameters
        double x = -60*cm;
        ARay* ray = new ARay(i, 660*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        
        manager->TraceNonSequential(*ray);
        
        // Process ray results
        if (drawRays) {
            // Visualization code - unchanged
        } else {
            Double_t lastPoint[3];
            ray->GetLastPoint(lastPoint);
            if (lastPoint[2] < exitPortZ) {
                if (detector.checkIntersection(ray)) {
                    detector.hitCount++;
                    hitCount++;
                }
            }
        }

        // Fix: Use maxPoints parameter instead of fLimit
        if (ray->IsRunning() && ray->GetNpoints() >= maxPoints) {
            ray->Suspend();
        }

        delete ray;
    }
    return hitCount;
}

// Parallelized version of traceRays
int traceRaysParallel(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false, double x = -60*cm, double y = 0*cm, 
    double z = -80*cm, double dirX = 5, double dirY = 2, double dirZ = 0) {
    int hitCount = 0;
    
    // Create an ARayArray for batch processing
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        ARay* ray = new ARay(i, 660*nm, x, y, z, 0, dirX, dirY, dirZ);
        rayArray->Add(ray);
    }
    
    // Let ROBAST handle multi-threading internally
    manager->TraceNonSequential(rayArray);
    
    // Process results - this doesn't need to be parallel since tracing is done
    TObjArray* stopped = rayArray->GetStopped();
    TObjArray* exited = rayArray->GetExited();
    
    // Count hits from stopped rays
    for (int i = 0; i <= stopped->GetLast(); i++) {
        ARay* ray = (ARay*)stopped->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ && detector.checkIntersection(ray)) {
            hitCount++;
            detector.hitCount++;
        }
    }
    
    // Count hits from exited rays
    for (int i = 0; i <= exited->GetLast(); i++) {
        ARay* ray = (ARay*)exited->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ && detector.checkIntersection(ray)) {
            hitCount++;
            detector.hitCount++;
        }
    }
    
    // Clean up
    delete rayArray;
    
    return hitCount;
}

// Modified function to trace rays and check intersection with two detectors
int traceRaysParallelTwofold(AOpticsManager* manager, int n, double exitPortZ, 
                            Detector& detector1, Detector& detector2, 
                            bool drawRays = false, double x = -60*cm, double y = 0*cm, 
                            double z = -80*cm, double dirX = 5, double dirY = 2, double dirZ = 0) {
    int hitCount = 0;
    
    // Create an ARayArray for batch processing
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        ARay* ray = new ARay(i, 660*nm, x, y, z, 0, dirX, dirY, dirZ);
        rayArray->Add(ray);
    }
    
    // Let ROBAST handle multi-threading internally
    manager->TraceNonSequential(rayArray);
    
    // Process results - check both detectors for each ray
    TObjArray* stopped = rayArray->GetStopped();
    TObjArray* exited = rayArray->GetExited();
    
    // Process stopped rays
    for (int i = 0; i <= stopped->GetLast(); i++) {
        ARay* ray = (ARay*)stopped->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ) {
            // Check first detector
            if (detector1.checkIntersection(ray)) {
                detector1.hitCount++;
                hitCount++;
            }
            
            // Check second detector
            if (detector2.checkIntersection(ray)) {
                detector2.hitCount++;
                hitCount++;
            }
        }
    }
    
    // Process exited rays
    for (int i = 0; i <= exited->GetLast(); i++) {
        ARay* ray = (ARay*)exited->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ) {
            // Check first detector
            if (detector1.checkIntersection(ray)) {
                detector1.hitCount++;
                hitCount++;
            }
            
            // Check second detector
            if (detector2.checkIntersection(ray)) {
                detector2.hitCount++;
                hitCount++;
            }
        }
    }
    
    // Clean up
    delete rayArray;
    
    return hitCount;
}

// Function to get a unique filename if the target file already exists
std::string getUniqueFilename(const std::string& basePath) {
    // First check if file exists
    FILE* file = fopen(basePath.c_str(), "r");
    if (!file) {
        // File doesn't exist, so we can use this name
        return basePath;
    }
    fclose(file);
    
    // If file exists, find a unique name
    std::string directory;
    std::string filename;
    
    // Parse the path
    size_t lastSlash = basePath.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        directory = basePath.substr(0, lastSlash + 1);
        filename = basePath.substr(lastSlash + 1);
    } else {
        directory = "";
        filename = basePath;
    }
    
    // Split filename into stem and extension
    size_t lastDot = filename.find_last_of('.');
    std::string stem;
    std::string extension;
    if (lastDot != std::string::npos) {
        stem = filename.substr(0, lastDot);
        extension = filename.substr(lastDot);
    } else {
        stem = filename;
        extension = "";
    }
    
    // Find unique filename by incrementing counter
    int counter = 1;
    std::string newPath;
    do {
        newPath = directory + stem + "_" + std::to_string(counter) + extension;
        counter++;
        file = fopen(newPath.c_str(), "r");
        if (file) {
            fclose(file);
        } else {
            // Found a name that doesn't exist
            break;
        }
    } while (true);
    
    return newPath;
}

// Thread-safe queue for processing positions
class ThreadSafeQueue {
private:
    std::queue<std::pair<int, int>> data;
    std::mutex mutex;
    std::condition_variable cv;
    bool done = false;

public:
    bool pop(std::pair<int, int>& item) {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this] { return !data.empty() || done; });
        
        if (data.empty() && done) return false;
        
        item = data.front();
        data.pop();
        return true;
    }
    
    void push(std::pair<int, int> item) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            data.push(item);
        }
        cv.notify_one();
    }
    
    void finish() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            done = true;
        }
        cv.notify_all();
    }
};

// Move this struct outside the function to fix linkage issues
// Add this after the ThreadSafeQueue class definition and before sweepDetectorOptimized
struct WorkChunk {
    int startIdx;
    int count;
};

// Helper function to get current time formatted as string
std::string getCurrentTimeString() {
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%H:%M:%S", timeinfo);
    return std::string(timeBuffer);
}

// Add timing points to the sweepDetectorTwofold function
void sweepDetectorTwofold(bool notify = true, const char* saveFolder = "results", int threads = -1, 
    double srcX = -60*cm, double srcY = 0*cm, double srcZ = -80*cm, double dirX = 5,
    double dirY = 2, double dirZ = 0, double thetaMax = THETA_MAX) {
    
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Starting sweep setup" << std::endl;
    TStopwatch setupTimer;
    setupTimer.Start();
    
    // Initialize ROOT
    TThread::Initialize();
    
    // Track start and finish times
    time_t startTime = time(nullptr);
    time_t finishTime;
    
    // Create the geometry manager
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Set up geometry
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Starting geometry setup" << std::endl;
    TStopwatch geometryTimer;
    geometryTimer.Start();
    
    setupOpticsManager(manager, MAX_REFLECTIONS, ROUGHNESS, REFLECTANCE, thetaMax, false);
    
    geometryTimer.Stop();
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Geometry setup completed in " 
              << geometryTimer.RealTime() << " seconds" << std::endl;
    
    // Configure internal threading for ray tracing
    if (std::thread::hardware_concurrency() > 1) {
        int maxRayThreads = std::max(1, std::min(4, (int)std::thread::hardware_concurrency()));
        manager->SetMaxThreads(maxRayThreads);
        std::cout << "Using " << maxRayThreads << " threads for ray tracing" << std::endl;
    }
    
    int n = 50000; // rays per position
    double exitPortZ = -100*cm;
    
    // Configure angular binning
    const int nThetaBins = 180;
    const int nPhiBins = 90;
    
    // Create directory for results
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    int result = system(mkdirCmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Could not create directory: " << saveFolder << std::endl;
    } else {
        std::cout << "Using directory: " << saveFolder << std::endl;
    }
    
    // Generate filename based on parameters with twofold indicator
    std::string filename = "fluxmap_twofold_" + std::to_string(n) + "rays_" + 
                          std::to_string(nThetaBins) + "x" + std::to_string(nPhiBins) + 
                          "_src" + std::to_string(int(srcX/cm)) + "_" + 
                          std::to_string(int(srcY/cm)) + "_" + 
                          std::to_string(int(srcZ/cm)) + ".csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    
    // Get unique filename if already exists
    fullPath = getUniqueFilename(fullPath);
    
    // Open CSV file
    std::ofstream csvFile(fullPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
        return;
    }
    
    // Create 2D histogram for the angular sweep
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    // Create two detectors - use a larger size to catch more rays
    Detector detector1(40*cm, 40*cm);
    Detector detector2(40*cm, 40*cm);
    int completedPositions = 0;
    
    // Write CSV metadata
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Writing file metadata" << std::endl;
    
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    // Write metadata (header content unchanged)
    csvFile << "# Flux Map Data (Twofold Method) - Generated: " << timeBuffer << std::endl;
    csvFile << "# Number of rays per position: " << n << std::endl;
    csvFile << "# Detector dimensions: " << detector1.width/cm << "cm x " << detector1.height/cm << "cm" << std::endl;
    csvFile << "# Sphere inner radius: " << INNER_RADIUS/cm << "cm" << std::endl;
    csvFile << "# Sphere outer radius: " << OUTER_RADIUS/cm << "cm" << std::endl;
    csvFile << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvFile << "# Theta bins: " << nThetaBins << std::endl;
    csvFile << "# Phi bins: " << nPhiBins << std::endl;
    csvFile << "# Mirror reflectance: " << REFLECTANCE << std::endl;
    csvFile << "# Gaussian roughness: " << ROUGHNESS << std::endl;
    csvFile << "# Lambertian scattering: enabled" << std::endl;
    csvFile << "# Source position (x,y,z): " << srcX/cm << "cm, " << srcY/cm << "cm, " << srcZ/cm << "cm" << std::endl;
    csvFile << "# Source direction (x,y,z): " << dirX << ", " << dirY << ", " << dirZ << std::endl;
    csvFile << "# Max reflections: " << MAX_REFLECTIONS << std::endl;
    csvFile << "# Method: Twofold (two detectors 180° apart)" << std::endl;
    csvFile << "theta,phi,fraction" << std::endl;  // Header
    
    setupTimer.Stop();
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Sweep setup completed in " 
              << setupTimer.RealTime() << " seconds" << std::endl;
    
    int totalPositions = nThetaBins * nPhiBins;
    // We only need to run half the phi angles
    int actualRuns = nThetaBins * (nPhiBins/2);
    std::cout << "\nStarting twofold detector sweep with " << n << " rays per position "
              << "(" << actualRuns << " runs for " << totalPositions << " total positions)..." << std::endl;
    
    // Start timer for overall progress
    TStopwatch timer;
    timer.Start();
    
    // Initialize timers for measuring angle completion times
    TStopwatch thetaTimer;
    std::vector<double> thetaTimes;
    thetaTimes.reserve(nThetaBins);

    // For tracking recent point times for better time estimates
    std::deque<double> recentPointTimes;
    const int maxRecentPoints = 20;

    // For debug - track total ray hits
    int totalHitRays = 0;

    // Process each position - using half the phi angles
    for(int i = 0; i < nThetaBins; i++) {
        thetaTimer.Start();
        double theta = (i + 0.5) * 90.0/nThetaBins;
        std::cout << "\n[DEBUG TIME " << getCurrentTimeString() << "] Starting theta angle: " << theta << "°" << std::endl;
        
        // Only need to process half the phi range because we're using two detectors
        for(int j = 0; j < nPhiBins/2; j++) {
            double phi1 = (j + 0.5) * 360.0/nPhiBins;
            double phi2 = phi1 + 180.0;  // Second detector 180° opposite
            
            // Ensure phi2 is in range [0, 360)
            if (phi2 >= 360.0) phi2 -= 360.0;
            
            // Display positions being processed
            std::cout << "\n[DEBUG TIME " << getCurrentTimeString() << "] Processing positions (θ=" << theta 
                      << "°, φ1=" << phi1 << "°, φ2=" << phi2 << "°)" << std::endl;
            
            // Reset detectors and position them
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Setting up detectors" << std::endl;
            TStopwatch detectorSetupTimer;
            detectorSetupTimer.Start();
            
            detector1.hitCount = 0;
            detector2.hitCount = 0;
            detector1.setPosition(theta, phi1, 100*cm);
            detector2.setPosition(theta, phi2, 100*cm);
            
            detectorSetupTimer.Stop();
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Detector setup completed in " 
                      << detectorSetupTimer.RealTime() << " seconds" << std::endl;
            
            // Debug output for detector positions
            std::cout << "  Detector 1 position (x,y,z): (" 
                      << detector1.x/cm << ", " << detector1.y/cm << ", " << detector1.z/cm 
                      << ") cm" << std::endl;
            std::cout << "  Detector 2 position (x,y,z): (" 
                      << detector2.x/cm << ", " << detector2.y/cm << ", " << detector2.z/cm 
                      << ") cm" << std::endl;
            
            // Start timing this individual point
            TStopwatch pointTimer;
            pointTimer.Start();

            // Process rays with twofold detector checking
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Starting ray tracing for this position" << std::endl;
            
            int hitCount = traceRaysParallelTwofold(manager, n, exitPortZ, detector1, detector2, 
                                                   false, srcX, srcY, srcZ, dirX, dirY, dirZ);
            totalHitRays += hitCount;
            
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Calculating fractions and updating results" << std::endl;
            TStopwatch resultUpdateTimer;
            resultUpdateTimer.Start();
            
            // Calculate fractions for both detectors
            double fraction1 = double(detector1.hitCount)/double(n);
            double fraction2 = double(detector2.hitCount)/double(n);
            
            // Find index for detector 2 position in histogram
            int j2 = j + nPhiBins/2;
            if (j2 >= nPhiBins) j2 -= nPhiBins;
            
            // Update histogram and write to CSV for both positions
            fluxMap->SetBinContent(i+1, j+1, fraction1);
            fluxMap->SetBinContent(i+1, j2+1, fraction2);
            
            csvFile << std::fixed << std::setprecision(6) 
                   << theta << "," << phi1 << "," << fraction1 << std::endl;
            csvFile << std::fixed << std::setprecision(6) 
                   << theta << "," << phi2 << "," << fraction2 << std::endl;
            
            // Flush to ensure data is written even if program crashes
            csvFile.flush();
            
            resultUpdateTimer.Stop();
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Results update completed in " 
                      << resultUpdateTimer.RealTime() << " seconds" << std::endl;
            
            // After the positions are completed:
            pointTimer.Stop();
            double pointTime = pointTimer.RealTime();
            recentPointTimes.push_back(pointTime);
            if (recentPointTimes.size() > maxRecentPoints) {
                recentPointTimes.pop_front();
            }

            // Update position counter - we've completed 2 positions
            completedPositions += 2;
            
            // Print progress (unchanged)
            double percentComplete = (100.0 * completedPositions) / totalPositions;
            std::cout << "  Progress: " << std::fixed << std::setprecision(1) 
                      << percentComplete << "% (" << completedPositions << "/" << totalPositions 
                      << " positions)" << std::endl;
            std::cout << "  Detector 1 (φ=" << phi1 << "°): " << detector1.hitCount << "/" << n 
                      << " = " << std::setprecision(8) << fraction1 << std::endl;
            std::cout << "  Detector 2 (φ=" << phi2 << "°): " << detector2.hitCount << "/" << n 
                      << " = " << std::setprecision(8) << fraction2 << std::endl;

            // Calculate estimated time based on average of recent points
            if (recentPointTimes.size() > 5) {
                double avgTimePerPoint = 0;
                for (double time : recentPointTimes) {
                    avgTimePerPoint += time;
                }
                avgTimePerPoint /= recentPointTimes.size();
                
                // We're doing two positions per run, so divide remaining positions by 2
                int remainingRuns = (totalPositions - completedPositions) / 2;
                double remainingSeconds = avgTimePerPoint * remainingRuns;
                
                // Calculate estimated end time
                time_t currentTime = time(nullptr);
                time_t estimatedEndTime = currentTime + static_cast<time_t>(remainingSeconds);
                
                // Format times for display
                int hours = static_cast<int>(remainingSeconds / 3600);
                int minutes = static_cast<int>((remainingSeconds - hours * 3600) / 60);
                int seconds = static_cast<int>(remainingSeconds - hours * 3600 - minutes * 60);
                
                char endTimeBuffer[80];
                struct tm* timeinfo = localtime(&estimatedEndTime);
                strftime(endTimeBuffer, sizeof(endTimeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
                
                std::cout << "  Estimated remaining time: ";
                if (hours > 0) std::cout << hours << "h ";
                if (hours > 0 || minutes > 0) std::cout << minutes << "m ";
                std::cout << seconds << "s (ETA: " << endTimeBuffer << ")" << std::endl;
            }
        }
        
        // Record theta angle completion time
        thetaTimer.Stop();
        double thetaTime = thetaTimer.RealTime();
        thetaTimes.push_back(thetaTime);
        
        std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Completed theta=" << theta 
                  << "° in " << thetaTime << " seconds" << std::endl;
        
        // Calculate estimated remaining time (unchanged)
        if (i < nThetaBins - 1) {
            double avgThetaTime = 0;
            for (double t : thetaTimes) {
                avgThetaTime += t;
            }
            avgThetaTime /= thetaTimes.size();
            
            double remainingTime = avgThetaTime * (nThetaBins - i - 1);
            int hours = static_cast<int>(remainingTime / 3600);
            int minutes = static_cast<int>((remainingTime - hours * 3600) / 60);
            int seconds = static_cast<int>(remainingTime - hours * 3600 - minutes * 60);
            
            std::cout << "Estimated remaining time for sweep: ";
            if (hours > 0) std::cout << hours << "h ";
            if (hours > 0 || minutes > 0) std::cout << minutes << "m ";
            std::cout << seconds << std::endl;
        }
    }

    // After sweep completes, add final metadata to file
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Finalizing sweep and saving results" << std::endl;
    TStopwatch finalizingTimer;
    finalizingTimer.Start();
    
    timer.Stop();
    double realTime = timer.RealTime();
    double cpuTime = timer.CpuTime();
    
    // Finalize CSV file and canvas creation (unchanged)
    finishTime = time(nullptr);
    char finishTimeBuffer[80];
    struct tm* finishTimeInfo = localtime(&finishTime);
    strftime(finishTimeBuffer, sizeof(finishTimeBuffer), "%Y-%m-%d %H:%M:%S", finishTimeInfo);
    
    csvFile << "# Sweep completed at: " << finishTimeBuffer << std::endl;
    csvFile << "# Total execution time: " << realTime << " seconds" << std::endl;
    csvFile << "# Total ray hits: " << totalHitRays << " out of " << (n * totalPositions) << std::endl;
    csvFile.close();
    
    std::cout << "\nFlux map data saved to '" << fullPath << "'" << std::endl;
    std::cout << "Sweep completed in " << realTime << " seconds (wall clock)" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;
    std::cout << "Efficiency gain: ~2x (processed " << totalPositions << " positions with " 
              << actualRuns << " simulation runs)" << std::endl;
    
    // Create canvas and draw the histogram
    TCanvas* c = new TCanvas("c", "Detector Flux Map (Twofold Method)", 1000, 800);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    
    fluxMap->GetZaxis()->SetTitle("Fraction of rays detected");
    fluxMap->Draw("COLZ");
    
    // Add color scale
    gPad->Update();
    TPaletteAxis* palette = (TPaletteAxis*)fluxMap->GetListOfFunctions()->FindObject("palette");
    if(palette) {
        palette->SetX1NDC(0.9);
        palette->SetX2NDC(0.92);
    }

    finalizingTimer.Stop();
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Finalization completed in " 
              << finalizingTimer.RealTime() << " seconds" << std::endl;
    
    if (notify) {
        // NOTIFICATION SYSTEM
        std::cout << "\n***** TWOFOLD SWEEP COMPLETE *****\n" << std::endl;
        std::cout << '\a' << std::endl;
    }
    
    // Clean up manager
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Cleaning up manager" << std::endl;
    delete manager;
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Sweep completed" << std::endl;
}

void visualizeDetector(double detTheta = 45.0, double detPhi = 0.0) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    double angle = 150.0;
    // Setup geometry using our standard function
    setupOpticsManager(manager, MAX_REFLECTIONS, ROUGHNESS, REFLECTANCE, angle, true);
    
    // Create and position detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(detTheta, detPhi, 100*cm);
    
    // Add detector to geometry - Fix type mismatch
    TGeoVolume* vol = manager->GetTopVolume();
    AOpticalComponent* world = dynamic_cast<AOpticalComponent*>(vol);
    if (!world) {
        std::cerr << "Error: Top volume is not an AOpticalComponent" << std::endl;
        return;
    }
    detector.AddToGeometry(world);
    
    // Close the geometry after all components are added
    manager->CloseGeometry();
    
    // Source position
    double srcX = -60*cm;
    double srcY = 0*cm;
    double srcZ = -75*cm;
    
    // Direction
    double dirX = 5;
    double dirY = 0;
    double dirZ = 0;
    
    // Create source marker
    TGeoSphere* sourceSphere = new TGeoSphere("sourceSphere", 0, 3*cm);
    TGeoVolume* sourceVol = new TGeoVolume("source", sourceSphere);
    sourceVol->SetLineColor(kRed);
    sourceVol->SetFillColor(kRed);
    world->AddNode(sourceVol, 1, new TGeoTranslation(srcX, srcY, srcZ));
    
    // Set up the 3D viewer
    gStyle->SetCanvasPreferGL(true);
    world->Draw("ogl");
    
    // Get viewer and set initial properties
    TGLViewer* glv = (TGLViewer*)gPad->GetViewer3D();
    if (glv) {
        glv->SetStyle(TGLRnrCtx::kWireFrame);
        glv->UpdateScene();
    }
    
    c->Update();
    gSystem->ProcessEvents();
    gSystem->Sleep(500); // Give time for the viewer to initialize
    
    // Create a transparent pad for rays
    TPad* rayPad = new TPad("rayPad", "Ray Visualization", 0, 0, 1, 1);
    rayPad->SetFillStyle(4000); // Transparent
    rayPad->Draw();
    rayPad->cd();
    
    // Draw rays with our parallel tracing function
    int n = 1000; // Start with fewer rays for visualization
    double exitPortZ = -100*cm;
    detector.hitCount = 0;
    
    // We need to trace the rays first, then manually create polylines
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        ARay* ray = new ARay(i, 660*nm, srcX, srcY, srcZ, 0, dirX, dirY, dirZ);
        rayArray->Add(ray);
    }
    
    // Trace rays
    manager->TraceNonSequential(rayArray);
    
    // Process results and draw rays
    TObjArray* stopped = rayArray->GetStopped();
    TObjArray* exited = rayArray->GetExited();
    
    int exitCount = 0;
    int hitCount = 0;
    int suspendedCount = 0;
    int absorbedCount = 0;
    std::vector<TPolyLine3D*> rayLines; // Store to prevent garbage collection
    
    // Process stopped rays
    for (int i = 0; i <= stopped->GetLast(); i++) {
        ARay* ray = (ARay*)stopped->At(i);
        if (!ray) continue;
        
        TPolyLine3D* pol = ray->MakePolyLine3D();
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (ray->IsSuspended()) {
            pol->SetLineColor(kMagenta);
            suspendedCount++;
        } else if (ray->IsAbsorbed()) {
            pol->SetLineColor(kBlack);
            absorbedCount++;
        } else if (lastPoint[2] < exitPortZ) {
            exitCount++;
            if (detector.checkIntersection(ray)) {
                pol->SetLineColor(kGreen);
                detector.hitCount++;
                hitCount++;
            } else {
                pol->SetLineColor(kYellow);
            }
        } else {
            pol->SetLineColor(kRed);
            if (lastPoint[2] >= exitPortZ) {
                // This is a "reflected back" ray
                pol->SetLineColor(kRed);
                pol->SetLineWidth(5);  // Make thicker
                
                // Add a marker at the endpoint
                TMarker* marker = new TMarker(lastPoint[0], lastPoint[1], 20);
                marker->SetMarkerColor(kRed);
                marker->SetMarkerStyle(20);  // Solid circle
                marker->SetMarkerSize(2);
                marker->Draw();
                
                std::cout << "Red ray endpoint: (" << lastPoint[0]/cm << ", " 
                          << lastPoint[1]/cm << ", " << lastPoint[2]/cm << ") cm" << std::endl;
            }
        }
        
        pol->SetLineWidth(3);
        pol->Draw();
        rayLines.push_back(pol);
    }
    
    // Process exited rays
    for (int i = 0; i <= exited->GetLast(); i++) {
        ARay* ray = (ARay*)exited->At(i);
        if (!ray) continue;
        
        TPolyLine3D* pol = ray->MakePolyLine3D();
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ) {
            exitCount++;
            if (detector.checkIntersection(ray)) {
                pol->SetLineColor(kGreen);
                detector.hitCount++;
                hitCount++;
            } else {
                pol->SetLineColor(kYellow);
            }
        } else {
            pol->SetLineColor(kRed);
        }
        
        pol->SetLineWidth(3);
        pol->Draw();
        rayLines.push_back(pol);
    }
    
    // Clean up ray array
    delete rayArray;
    
    rayPad->Modified();
    rayPad->Update();
    c->Update();
    
    // Print statistics
    std::cout << "\nDetector Information:" << std::endl;
    std::cout << "Total rays: " << n << std::endl;
    std::cout << "Rays exiting port: " << exitCount << std::endl;
    std::cout << "Rays hitting detector: " << hitCount << "/" << exitCount << " of exiting rays" << std::endl;
    std::cout << "Rays suspended: " << suspendedCount << std::endl;
    std::cout << "Rays absorbed: " << absorbedCount << std::endl;
    std::cout << "Rays reflected back: " << (n - exitCount - suspendedCount - absorbedCount) << std::endl;
    std::cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << std::endl;
    
    // Add legend
    TPaveText* legend = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    legend->SetFillColor(kWhite);
    legend->SetTextAlign(12);
    legend->AddText("Green: Hit detector");
    legend->AddText("Yellow: Exit port, miss detector");
    legend->AddText("Red: Did not exit");
    legend->AddText("Magenta: Suspended (max reflections)");
    legend->AddText("Black: Absorbed");
    legend->Draw();
    
    c->Update();
}

// New optimized function using the trace-once approach
void sweepDetectorTraceOnce(bool notify = true, const char* saveFolder = "results", int threads = -1, 
    double srcX = -60*cm, double srcY = 0*cm, double srcZ = -80*cm, double dirX = 5,
    double dirY = 2, double dirZ = 0, double thetaMax = THETA_MAX) {
    
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Starting sweep setup" << std::endl;
    TStopwatch setupTimer;
    setupTimer.Start();
    
    // Initialize ROOT and geometry as before
    TThread::Initialize();
    time_t startTime = time(nullptr);
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    setupOpticsManager(manager, MAX_REFLECTIONS, ROUGHNESS, REFLECTANCE, thetaMax, false);
    
    // Configure threading
    if (std::thread::hardware_concurrency() > 1) {
        int maxRayThreads = std::max(1, std::min(4, (int)std::thread::hardware_concurrency()));
        manager->SetMaxThreads(maxRayThreads);
        std::cout << "Using " << maxRayThreads << " threads for ray tracing" << std::endl;
    }
    
    // Configure parameters
    int n = 100000; // rays total
    double exitPortZ = -100*cm;
    const int nThetaBins = 180;
    const int nPhiBins = 90;
    int totalPositions = nThetaBins * nPhiBins;
    
    // Create output directory and file
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    system(mkdirCmd.c_str());
    
    std::string filename = "fluxmap_traceonce_" + std::to_string(n) + "rays_" + 
                          std::to_string(nThetaBins) + "x" + std::to_string(nPhiBins) + 
                          "_src" + std::to_string(int(srcX/cm)) + "_" + 
                          std::to_string(int(srcY/cm)) + "_" + 
                          std::to_string(int(srcZ/cm)) + ".csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    fullPath = getUniqueFilename(fullPath);
    
    // Buffer for CSV data
    std::stringstream csvBuffer;
    
    // Write metadata
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    csvBuffer << "# Flux Map Data (Trace-Once Method) - Generated: " << timeBuffer << std::endl;
    csvBuffer << "# Number of rays: " << n << std::endl;
    csvBuffer << "# Detector dimensions: 40cm x 40cm" << std::endl;
    csvBuffer << "# Sphere inner radius: " << INNER_RADIUS/cm << "cm" << std::endl;
    csvBuffer << "# Sphere outer radius: " << OUTER_RADIUS/cm << "cm" << std::endl;
    csvBuffer << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvBuffer << "# Theta bins: " << nThetaBins << std::endl;
    csvBuffer << "# Phi bins: " << nPhiBins << std::endl;
    csvBuffer << "# Mirror reflectance: " << REFLECTANCE << std::endl;
    csvBuffer << "# Gaussian roughness: " << ROUGHNESS << std::endl;
    csvBuffer << "# Lambertian scattering: enabled" << std::endl;
    csvBuffer << "# Source position (x,y,z): " << srcX/cm << "cm, " << srcY/cm << "cm, " << srcZ/cm << "cm" << std::endl;
    csvBuffer << "# Source direction (x,y,z): " << dirX << ", " << dirY << ", " << dirZ << std::endl;
    csvBuffer << "# Max reflections: " << MAX_REFLECTIONS << std::endl;
    csvBuffer << "# Method: Trace-Once (single trace, multiple detector positions)" << std::endl;
    csvBuffer << "theta,phi,fraction" << std::endl;
    
    // Create histogram for results
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map (Trace-Once Method);#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    // ----- TRACE RAYS ONCE -----
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Tracing all rays once" << std::endl;
    TStopwatch rayTraceTimer;
    rayTraceTimer.Start();
    
    // Create an ARayArray for batch processing
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        ARay* ray = new ARay(i, 660*nm, srcX, srcY, srcZ, 0, dirX, dirY, dirZ);
        rayArray->Add(ray);
    }
    
    // Trace all rays at once
    manager->TraceNonSequential(rayArray);
    
    // Structure to store ray endpoints for intersection testing
    struct RayEndpoint {
        double start[3];  // Starting point of final segment
        double end[3];    // End point (ray exit position)
        double dir[3];    // Direction of final segment
    };
    std::vector<RayEndpoint> exitedRays;
    
    // Process stopped rays
    TObjArray* stopped = rayArray->GetStopped();
    for (int i = 0; i <= stopped->GetLast(); i++) {
        ARay* ray = (ARay*)stopped->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ) {
            // We only need the final ray segment for detector intersection
            RayEndpoint endpoint;
            
            // Get last two points to determine the final segment
            int nPoints = ray->GetNpoints();
            if (nPoints >= 2) {
                // Get second-to-last point as segment start
                Double_t secondLastPoint[3];
                ray->GetPoint(nPoints - 2, secondLastPoint);
                
                // Copy points to our storage
                endpoint.start[0] = secondLastPoint[0];
                endpoint.start[1] = secondLastPoint[1];
                endpoint.start[2] = secondLastPoint[2];
                
                endpoint.end[0] = lastPoint[0];
                endpoint.end[1] = lastPoint[1];
                endpoint.end[2] = lastPoint[2];
                
                // Calculate direction from the two points
                double dx = lastPoint[0] - secondLastPoint[0];
                double dy = lastPoint[1] - secondLastPoint[1];
                double dz = lastPoint[2] - secondLastPoint[2];
                
                // Normalize direction
                double mag = sqrt(dx*dx + dy*dy + dz*dz);
                endpoint.dir[0] = dx/mag;
                endpoint.dir[1] = dy/mag;
                endpoint.dir[2] = dz/mag;
                
                // Store this endpoint
                exitedRays.push_back(endpoint);
            }
        }
    }
    
    // Process exited rays (similar to stopped rays)
    TObjArray* exited = rayArray->GetExited();
    for (int i = 0; i <= exited->GetLast(); i++) {
        ARay* ray = (ARay*)exited->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (lastPoint[2] < exitPortZ) {
            // Same process as for stopped rays
            RayEndpoint endpoint;
            
            int nPoints = ray->GetNpoints();
            if (nPoints >= 2) {
                Double_t secondLastPoint[3];
                ray->GetPoint(nPoints - 2, secondLastPoint);
                
                endpoint.start[0] = secondLastPoint[0];
                endpoint.start[1] = secondLastPoint[1];
                endpoint.start[2] = secondLastPoint[2];
                
                endpoint.end[0] = lastPoint[0];
                endpoint.end[1] = lastPoint[1];
                endpoint.end[2] = lastPoint[2];
                
                double dx = lastPoint[0] - secondLastPoint[0];
                double dy = lastPoint[1] - secondLastPoint[1];
                double dz = lastPoint[2] - secondLastPoint[2];
                
                double mag = sqrt(dx*dx + dy*dy + dz*dz);
                endpoint.dir[0] = dx/mag;
                endpoint.dir[1] = dy/mag;
                endpoint.dir[2] = dz/mag;
                
                exitedRays.push_back(endpoint);
            }
        }
    }
    
    // Clean up ray array as we no longer need the original ray objects
    delete rayArray;
    
    rayTraceTimer.Stop();
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Ray tracing completed in " 
              << rayTraceTimer.RealTime() << " seconds" << std::endl;
    std::cout << "Total rays exiting port: " << exitedRays.size() << " out of " << n << std::endl;
    
    // ----- TEST AGAINST ALL DETECTOR POSITIONS -----
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Starting detector position sweep" << std::endl;
    TStopwatch sweepTimer;
    sweepTimer.Start();
    
    // Use vector to store results
    std::vector<std::tuple<double, double, double>> results;
    results.reserve(1000); // Reserve space for batched results
    
    int completedPositions = 0;
    
    // Process each position
    for(int i = 0; i < nThetaBins; i++) {
        double theta = (i + 0.5) * 90.0/nThetaBins;
        
        for(int j = 0; j < nPhiBins; j++) {
            double phi = (j + 0.5) * 360.0/nPhiBins;
            
            // Create detector at this position
            Detector detector(40*cm, 40*cm);
            detector.setPosition(theta, phi, 100*cm);
            
            // Count hits for this detector position
            int hitCount = 0;
            
            // Optimized intersection testing
            for (const auto& endpoint : exitedRays) {
                // Create temporary ray for testing
                ARay tempRay(0, 660*nm, 
                            endpoint.start[0], endpoint.start[1], endpoint.start[2], 
                            0, 
                            endpoint.dir[0], endpoint.dir[1], endpoint.dir[2]);
                
                // Check intersection
                if (detector.checkIntersection(&tempRay)) {
                    hitCount++;
                }
            }
            
            // Calculate fraction
            double fraction = double(hitCount) / double(n);
            
            // Store result
            results.push_back(std::make_tuple(theta, phi, fraction));
            
            // Update histogram
            fluxMap->SetBinContent(i+1, j+1, fraction);
            
            // Update progress counter
            completedPositions++;
            
            // Print progress every 100 positions
            if (completedPositions % 100 == 0 || completedPositions == totalPositions) {
                double percentComplete = (100.0 * completedPositions) / totalPositions;
                std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Progress: " 
                          << std::fixed << std::setprecision(1) << percentComplete << "% (" 
                          << completedPositions << "/" << totalPositions << " positions)" << std::endl;
            }
        }
        
        // Write accumulated results every 10 theta angles
        if (i % 10 == 9 || i == nThetaBins - 1) {
            std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Writing batch of results to CSV" << std::endl;
            
            for (const auto& result : results) {
                csvBuffer << std::fixed << std::setprecision(6)
                         << std::get<0>(result) << ","
                         << std::get<1>(result) << ","
                         << std::get<2>(result) << std::endl;
            }
            
            // Write to file
            std::ofstream csvFile(fullPath, std::ios::trunc);
            if (csvFile.is_open()) {
                csvFile << csvBuffer.str();
                csvFile.close();
            } else {
                std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
            }
            
            // Clear results vector to free memory
            results.clear();
            results.reserve(1000);
        }
    }
    
    sweepTimer.Stop();
    double sweepTime = sweepTimer.RealTime();
    double rayTime = rayTraceTimer.RealTime();
    double totalTime = setupTimer.RealTime() + rayTime + sweepTime;
    
    // After sweep completes, create visualization and print stats
    std::cout << "[DEBUG TIME " << getCurrentTimeString() << "] Finalizing results" << std::endl;
    
    // Create canvas and draw the histogram
    TCanvas* c = new TCanvas("c", "Detector Flux Map (Trace-Once Method)", 1000, 800);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    
    fluxMap->GetZaxis()->SetTitle("Fraction of rays detected");
    fluxMap->Draw("COLZ");
    
    // Add color scale
    gPad->Update();
    TPaletteAxis* palette = (TPaletteAxis*)fluxMap->GetListOfFunctions()->FindObject("palette");
    if(palette) {
        palette->SetX1NDC(0.9);
        palette->SetX2NDC(0.92);
    }
    
    // Add final metadata to CSV
    time_t finishTime = time(nullptr);
    char finishTimeBuffer[80];
    struct tm* finishTimeInfo = localtime(&finishTime);
    strftime(finishTimeBuffer, sizeof(finishTimeBuffer), "%Y-%m-%d %H:%M:%S", finishTimeInfo);
    
    // Append to CSV file
    std::ofstream csvFileAppend(fullPath, std::ios::app);
    if (csvFileAppend.is_open()) {
        csvFileAppend << "# Sweep completed at: " << finishTimeBuffer << std::endl;
        csvFileAppend << "# Total execution time: " << totalTime << " seconds" << std::endl;
        csvFileAppend << "# Ray tracing time: " << rayTime << " seconds" << std::endl;
        csvFileAppend << "# Detector sweep time: " << sweepTime << " seconds" << std::endl;
        csvFileAppend << "# Total rays exiting port: " << exitedRays.size() << " out of " << n << std::endl;
        csvFileAppend.close();
    }
    
    std::cout << "\nFlux map data saved to '" << fullPath << "'" << std::endl;
    std::cout << "Ray tracing completed in " << rayTime << " seconds" << std::endl;
    std::cout << "Detector sweep completed in " << sweepTime << " seconds" << std::endl;
    std::cout << "Total execution time: " << totalTime << " seconds" << std::endl;
    
    if (notify) {
        // NOTIFICATION SYSTEM
        std::cout << "\n***** TRACE-ONCE SWEEP COMPLETE *****\n" << std::endl;
        std::cout << '\a' << std::endl;
    }
    
    // Clean up manager
    delete manager;
}

// Updated visualization function with option to show only red rays
void visualizeDetector(double detTheta = 45.0, double detPhi = 0.0, bool onlyShowRedRays = false) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    double angle = 150.0;
    // Setup geometry using our standard function
    setupOpticsManager(manager, MAX_REFLECTIONS, ROUGHNESS, REFLECTANCE, angle, true);
    
    // Create and position detector
    Detector detector(20*cm, 20*cm);
    detector.setPosition(detTheta, detPhi, 100*cm);
    
    // Add detector to geometry
    TGeoVolume* vol = manager->GetTopVolume();
    AOpticalComponent* world = dynamic_cast<AOpticalComponent*>(vol);
    if (!world) {
        std::cerr << "Error: Top volume is not an AOpticalComponent" << std::endl;
        return;
    }
    detector.AddToGeometry(world);
    
    // Close the geometry after all components are added
    manager->CloseGeometry();
    
    // Source position
    double srcX = -60*cm;
    double srcY = 0*cm;
    double srcZ = -75*cm;
    
    // Direction
    double dirX = 5;
    double dirY = 0;
    double dirZ = 0;
    
    // Create source marker
    TGeoSphere* sourceSphere = new TGeoSphere("sourceSphere", 0, 3*cm);
    TGeoVolume* sourceVol = new TGeoVolume("source", sourceSphere);
    sourceVol->SetLineColor(kRed);
    sourceVol->SetFillColor(kRed);
    world->AddNode(sourceVol, 1, new TGeoTranslation(srcX, srcY, srcZ));
    
    // Set up the 3D viewer
    gStyle->SetCanvasPreferGL(true);
    world->Draw("ogl");
    
    // Get viewer and set initial properties
    TGLViewer* glv = (TGLViewer*)gPad->GetViewer3D();
    if (glv) {
        glv->SetStyle(TGLRnrCtx::kWireFrame);
        glv->UpdateScene();
    }
    
    c->Update();
    gSystem->ProcessEvents();
    gSystem->Sleep(500); // Give time for the viewer to initialize
    
    // Create a transparent pad for rays
    TPad* rayPad = new TPad("rayPad", "Ray Visualization", 0, 0, 1, 1);
    rayPad->SetFillStyle(4000); // Transparent
    rayPad->Draw();
    rayPad->cd();
    
    // Draw rays with our parallel tracing function
    int n = onlyShowRedRays ? 5000 : 1000; // More rays when showing only red rays
    double exitPortZ = -100*cm;
    detector.hitCount = 0;
    
    // We need to trace the rays first, then manually create polylines
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        ARay* ray = new ARay(i, 660*nm, srcX, srcY, srcZ, 0, dirX, dirY, dirZ);
        rayArray->Add(ray);
    }
    
    // Trace rays
    manager->TraceNonSequential(rayArray);
    
    // Process results and draw rays
    TObjArray* stopped = rayArray->GetStopped();
    TObjArray* exited = rayArray->GetExited();
    
    int exitCount = 0;
    int hitCount = 0;
    int suspendedCount = 0;
    int absorbedCount = 0;
    int reflectedBackCount = 0;
    std::vector<TPolyLine3D*> rayLines; // Store to prevent garbage collection
    
    std::cout << "\nProcessing stopped rays..." << std::endl;
    // Process stopped rays
    for (int i = 0; i <= stopped->GetLast(); i++) {
        ARay* ray = (ARay*)stopped->At(i);
        if (!ray) continue;
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        bool isRed = false;
        
        if (ray->IsSuspended()) {
            suspendedCount++;
            if (onlyShowRedRays) continue; // Skip if only showing red rays
        } else if (ray->IsAbsorbed()) {
            absorbedCount++;
            if (onlyShowRedRays) continue; // Skip if only showing red rays
        } else if (lastPoint[2] < exitPortZ) {
            exitCount++;
            if (detector.checkIntersection(ray)) {
                detector.hitCount++;
                hitCount++;
            }
            if (onlyShowRedRays) continue; // Skip if only showing red rays
        } else {
            // This is a red ray (reflected back)
            reflectedBackCount++;
            isRed = true;
        }
        
        // If we're only showing red rays and this isn't red, skip
        if (onlyShowRedRays && !isRed) continue;
        
        TPolyLine3D* pol = ray->MakePolyLine3D();
        
        if (isRed) {
            // Enhanced styling for red rays
            pol->SetLineColor(kRed);
            pol->SetLineWidth(onlyShowRedRays ? 4 : 3);
            
            // If only showing red rays, add endpoint markers
            if (onlyShowRedRays) {
                // Add a marker at the last point
                TMarker* marker = new TMarker(lastPoint[0]/cm, lastPoint[1]/cm, 20);
                marker->SetMarkerColor(kRed);
                marker->SetMarkerStyle(20);  // Solid circle
                marker->SetMarkerSize(1.5);
                marker->Draw();
                
                std::cout << "Red ray endpoint: (" << lastPoint[0]/cm << ", " 
                          << lastPoint[1]/cm << ", " << lastPoint[2]/cm << ") cm" << std::endl;
            }
        } else {
            // Standard coloring for other rays
            if (ray->IsSuspended()) {
                pol->SetLineColor(kMagenta);
            } else if (ray->IsAbsorbed()) {
                pol->SetLineColor(kBlack);
            } else if (detector.checkIntersection(ray)) {
                pol->SetLineColor(kGreen);
            } else {
                pol->SetLineColor(kYellow);
            }
        }
        
        pol->SetLineWidth(2);
        pol->Draw();
        rayLines.push_back(pol);
    }
    
    // Process exited rays (only if not in red-only mode)
    if (!onlyShowRedRays) {
        for (int i = 0; i <= exited->GetLast(); i++) {
            ARay* ray = (ARay*)exited->At(i);
            if (!ray) continue;
            
            TPolyLine3D* pol = ray->MakePolyLine3D();
            
            Double_t lastPoint[3];
            ray->GetLastPoint(lastPoint);
            
            if (lastPoint[2] < exitPortZ) {
                exitCount++;
                if (detector.checkIntersection(ray)) {
                    pol->SetLineColor(kGreen);
                    detector.hitCount++;
                    hitCount++;
                } else {
                    pol->SetLineColor(kYellow);
                }
            } else {
                pol->SetLineColor(kRed);
                reflectedBackCount++;
            }
            
            pol->SetLineWidth(2);
            pol->Draw();
            rayLines.push_back(pol);
        }
    }
    
    // Clean up ray array
    delete rayArray;
    
    rayPad->Modified();
    rayPad->Update();
    c->Update();
    
    // Print statistics
    std::string mode = onlyShowRedRays ? "RED RAYS ONLY MODE" : "STANDARD VISUALIZATION";
    std::cout << "\n=== " << mode << " ===" << std::endl;
    std::cout << "Detector Information:" << std::endl;
    std::cout << "Total rays: " << n << std::endl;
    std::cout << "Rays exiting port: " << exitCount << std::endl;
    std::cout << "Rays hitting detector: " << hitCount << "/" << exitCount << " of exiting rays" << std::endl;
    std::cout << "Rays suspended: " << suspendedCount << std::endl;
    std::cout << "Rays absorbed: " << absorbedCount << std::endl;
    std::cout << "Rays reflected back: " << reflectedBackCount << std::endl;
    std::cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << std::endl;
    
    // Add legend appropriate to the mode
    TPaveText* legend = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    legend->SetFillColor(kWhite);
    legend->SetTextAlign(12);
    
    if (onlyShowRedRays) {
        legend->AddText("RED RAYS ONLY MODE");
        legend->AddText("Red: Rays that did not exit port");
        legend->AddText(Form("Count: %d / %d rays", reflectedBackCount, n));
        legend->AddText(Form("%.1f%% of total rays", 100.0 * reflectedBackCount / n));
    } else {
        legend->AddText("Green: Hit detector");
        legend->AddText("Yellow: Exit port, miss detector");
        legend->AddText("Red: Did not exit port");
        legend->AddText("Magenta: Suspended (max reflections)");
        legend->AddText("Black: Absorbed");
    }
    
    legend->Draw();
    
    c->Update();
}

// Helper function specifically for examining red rays in detail
void showRedRaysOnly() {
    visualizeDetector(45.0, 0.0, true);
}

void sweepSeries(){
    // Run sweeps for different source angles - using parameters instead of hardcoded values
    const double srcX = -60*cm;
    const double srcY = 0*cm;
    const double srcZ = -75*cm;
    const double dirXBase = 5;
    const double portAngle = 164.0; // degrees
    
    // Create a results directory with the base parameters in the name
    std::string baseFolder = "portAngleSweep_04_03_" + 
                            std::to_string(int(srcX/cm)) + "_" +
                            std::to_string(int(srcY/cm)) + "_" +
                            std::to_string(int(srcZ/cm)) + "_" +
                            std::to_string(int(portAngle));
    
    // Series different source angles
    // sweepDetector(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 0, 0);
    // sweepDetector(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 2, 0);
    // sweepDetector(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 4, 0);
    // sweepDetector(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 6, 0);
    // sweepDetector(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 8, 0);

    // Series different exit port sizes
    int n = 5;
    for (size_t i = 0; i < n; i++)
    {
        sweepDetectorTraceOnce(false, baseFolder.c_str(), 1, srcX, srcY, srcZ, dirXBase, 0, 0, portAngle);
    }
    
    // Notify user when all sweeps are complete
    std::cout << "\n***** ALL SWEEP SERIES COMPLETE *****\n" << std::endl;
    std::cout << '\a' << std::endl;  // Sound alert
}
