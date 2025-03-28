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
const double thetaMax = 170.;     // governs size of exit port
const int MAX_REFLECTIONS = 50000;

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
        
        // FIX: Calculate normal vector that actually points TO the exit port
        double dx = 0 - x;            // Vector FROM detector TO exit port
        double dy = 0 - y;
        double dz = -100*cm - z;
        double mag = sqrt(dx*dx + dy*dy + dz*dz);
        
        nx = dx/mag;  // Normal now correctly points toward exit port
        ny = dy/mag;
        nz = dz/mag;
    }
    
    bool checkIntersection(const ARay* ray) {
        static const double radiusSquared = (width/2) * (width/2);
        static const double maxDistance = 10*cm; // Maximum distance to consider
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        // Calculate distance from last point to detector plane
        double signedDistance = (lastPoint[0] - x)*nx + 
                                (lastPoint[1] - y)*ny + 
                                (lastPoint[2] - z)*nz;
        
        // Only consider points very close to or behind the detector plane
        if (signedDistance > maxDistance) return false;
        
        // Project last point onto detector plane
        double projected_x = lastPoint[0] - signedDistance*nx;
        double projected_y = lastPoint[1] - signedDistance*ny;
        double projected_z = lastPoint[2] - signedDistance*nz;
        
        // Calculate distance from projection to detector center
        double dx = projected_x - x;
        double dy = projected_y - y;
        double dz = projected_z - z;
        
        // Remove normal component to get projection on detector plane
        double dot = dx*nx + dy*ny + dz*nz;
        double px = dx - dot*nx;
        double py = dy - dot*ny;
        double pz = dz - dot*nz;
        
        // Debug output for first few rays to see where they're projecting
        // static int debugCount = 0;
        // if (debugCount < 5) {
            // std::cout << "Ray projected to (" << projected_x/cm << ", " 
                    //   << projected_y/cm << ", " << projected_z/cm << ") cm" << std::endl;
            // std::cout << "Distance from center: " << sqrt(px*px + py*py + pz*pz)/cm << " cm" 
                    //   << " (detector radius: " << width/2/cm << " cm)" << std::endl;
            // debugCount++;
        // }
        
        return (px*px + py*py + pz*pz) <= radiusSquared;
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
                         double roughness = 0.75, double reflectance = 0.99, 
                         bool drawRays = false) {
    // Set ray trace limit
    manager->SetLimit(maxReflections);
    
    // Create slightly larger world box to ensure ray origins are inside
    TGeoBBox* box = new TGeoBBox("box", 300*cm, 300*cm, 300*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    
    // Create spherical mirror with exit port
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    mirror->SetReflectance(reflectance);
    
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    // condition->EnableLambertian(false);
    condition->EnableLambertian(true);
    // condition->SetGaussianRoughness(roughness);
    condition->SetGaussianRoughness(0.01);
    
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
int traceRaysParallel(AOpticsManager* manager, int n, double exitPortZ, Detector& detector, bool drawRays = false) {
    int hitCount = 0;
    
    // Create an ARayArray for batch processing
    ARayArray* rayArray = new ARayArray();
    
    // Fill the array with rays
    for (int i = 0; i < n; ++i) {
        // Ray parameters
        double x = -60*cm;
        ARay* ray = new ARay(i, 660*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
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

int traceRaysWithStats(AOpticsManager* manager, int n, double exitPortZ, Detector& detector) {
    int hitCount = 0;
    int exitCount = 0;
    int suspendedCount = 0;
    int absorbedCount = 0;
    
    for (int i = 0; i < n; ++i) {
        // Use EXACT same ray parameters as the working code
        // Ray parameters
        double x = -60*cm;
        ARay* ray = new ARay(i, 660*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        
        manager->TraceNonSequential(*ray);
        
        if (ray->IsSuspended()) {
            suspendedCount++;
        } else if (ray->IsAbsorbed()) {
            absorbedCount++;
        }
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        // First check if ray exits port, count separately
        if (lastPoint[2] < exitPortZ) {
            exitCount++;
            
            // Then check if it also hits detector
            if (detector.checkIntersection(ray)) {
                hitCount++;
            }
        }
        
        delete ray;
    }
    
    std::cout << "Total rays: " << n << std::endl;
    std::cout << "Rays exiting port: " << exitCount << std::endl;
    std::cout << "Rays hitting detector: " << hitCount << std::endl;
    std::cout << "Rays suspended: " << suspendedCount << std::endl;
    std::cout << "Rays absorbed: " << absorbedCount << std::endl;
    std::cout << "Rays reflected back: " << (n - exitCount - suspendedCount - absorbedCount) << std::endl;
    
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

// Updated sweepDetector function with parallel position processing
void sweepDetector(bool notify = true, const char* saveFolder = "results", int threads = -1) {
    // Initialize ROOT
    TThread::Initialize();
    
    // Inform user we're using single-threaded mode
    std::cout << "Using single-threaded processing mode for maximum stability" << std::endl;
    
    // Create the geometry manager
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Set up geometry
    setupOpticsManager(manager, MAX_REFLECTIONS, 0.75, 0.99, false);
    
    // Configure internal threading for ray tracing
    if (std::thread::hardware_concurrency() > 1) {
        int maxRayThreads = std::max(1, std::min(4, (int)std::thread::hardware_concurrency()));
        manager->SetMaxThreads(maxRayThreads);
        std::cout << "Using " << maxRayThreads << " threads for ray tracing" << std::endl;
    }
    
    // Ray count - configurable
    int n = 50000;
    double exitPortZ = -100*cm;
    
    // Configure angular binning
    const int nThetaBins = 45;
    const int nPhiBins = 90;
    
    // Create directory and set up CSV file
    std::string filename = "detector_data_" + std::to_string(n) + "rays_" + std::to_string(nThetaBins * nPhiBins) + "points.csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    int result = system(mkdirCmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Could not create directory: " << saveFolder << std::endl;
    } else {
        std::cout << "Using directory: " << saveFolder << std::endl;
    }
    
    fullPath = getUniqueFilename(fullPath);
    
    // Open CSV file
    std::ofstream csvFile(fullPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
        return;
    }
    
    // Write metadata as comments in CSV
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    csvFile << "# Detector Data - Generated: " << timeBuffer << std::endl;
    csvFile << "# Number of rays: " << n << std::endl;
    csvFile << "# Detector dimensions: 20cm x 20cm" << std::endl;
    csvFile << "# Sphere inner radius: 100.1cm" << std::endl;
    csvFile << "# Sphere outer radius: 101cm" << std::endl;
    csvFile << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvFile << "# Theta bins: " << nThetaBins << std::endl;
    csvFile << "# Phi bins: " << nPhiBins << std::endl;
    csvFile << "# Mirror reflectance: 0.99" << std::endl;
    csvFile << "# Gaussian roughness: 0.75" << std::endl;
    csvFile << "# Lambertian scattering: enabled" << std::endl;
    csvFile << "theta,phi,fraction" << std::endl;
    
    // Create 2D histogram for the angular sweep
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    int totalPositions = nThetaBins * nPhiBins;
    std::cout << "\nStarting detector sweep with " << n << " rays per position "
              << "(" << totalPositions << " positions total)..." << std::endl;
    
    // Initialize detector and counter
    Detector detector(20*cm, 20*cm);
    int completedPositions = 0;
    
    // Start timer for overall progress
    TStopwatch timer;
    timer.Start();
    
    // Process each position sequentially
    for(int i = 0; i < nThetaBins; i++) {
        for(int j = 0; j < nPhiBins; j++) {
            double theta = (i + 0.5) * 90.0/nThetaBins;
            double phi = (j + 0.5) * 360.0/nPhiBins;
            
            // Display position being processed
            std::cout << "Processing position (θ=" << theta << "°, φ=" << phi << "°)..." << std::endl;
            
            // Reset detector and position it
            detector.hitCount = 0;
            detector.setPosition(theta, phi, 100*cm);
            
            // Process rays
            TStopwatch posTimer;
            posTimer.Start();
            
            traceRaysParallel(manager, n, exitPortZ, detector, false);
            
            posTimer.Stop();
            double fraction = double(detector.hitCount)/double(n);
            
            // Update histogram and write to CSV
            fluxMap->SetBinContent(i+1, j+1, fraction);
            csvFile << std::fixed << std::setprecision(6) 
                   << theta << "," << phi << "," << fraction << std::endl;
            
            // Update position counter
            completedPositions++;
            
            // Display result and progress
            std::cout << "  Result: " << detector.hitCount << "/" << n << " rays hit detector ("
                     << std::fixed << std::setprecision(6) << fraction * 100 << "%)" << std::endl;
            std::cout << "  Position time: " << posTimer.RealTime() << " seconds" << std::endl;
            
            // Progress reporting
            double percentComplete = (100.0 * completedPositions) / totalPositions;
            double elapsedTime = timer.RealTime();
            
            // Calculate time per position for estimation
            double avgTimePerPos = elapsedTime / completedPositions;
            double remainingTime = avgTimePerPos * (totalPositions - completedPositions);
            
            // Apply conservative factor based on progress
            if (percentComplete < 1.0) {
                remainingTime *= 1.2;
            } else if (percentComplete < 5.0) {
                remainingTime *= 1.1;
            }
            
            // Format elapsed and remaining times
            int elapsed_hours = static_cast<int>(elapsedTime / 3600);
            int elapsed_minutes = static_cast<int>((elapsedTime - elapsed_hours * 3600) / 60);
            int elapsed_seconds = static_cast<int>(elapsedTime - elapsed_hours * 3600 - elapsed_minutes * 60);
            
            int hours = static_cast<int>(remainingTime / 3600);
            int minutes = static_cast<int>((remainingTime - hours * 3600) / 60);
            int seconds = static_cast<int>(remainingTime - hours * 3600 - minutes * 60);
            
            // Display progress information
            std::cout << "  Progress: " << std::fixed << std::setprecision(1) 
                      << percentComplete << "% (" << completedPositions << "/" << totalPositions 
                      << " positions)" << std::endl;
            
            std::cout << "  Elapsed: ";
            if (elapsed_hours > 0) std::cout << elapsed_hours << "h ";
            if (elapsed_hours > 0 || elapsed_minutes > 0) std::cout << elapsed_minutes << "m ";
            std::cout << elapsed_seconds << "s" << std::endl;
            
            std::cout << "  Remaining: ";
            if (hours > 0) std::cout << hours << "h ";
            if (hours > 0 || minutes > 0) std::cout << minutes << "m ";
            std::cout << seconds << "s" << std::endl;
            
            // Show estimated completion time
            time_t finishTime = now + static_cast<time_t>(remainingTime);
            struct tm* finishTm = localtime(&finishTime);
            strftime(timeBuffer, sizeof(timeBuffer), "%H:%M:%S (%Y-%m-%d)", finishTm);
            std::cout << "  Est. completion at: " << timeBuffer << std::endl;
            
            std::cout << std::endl;
        }
    }
    
    // After sweep completes, add final metadata to file
    timer.Stop();
    double realTime = timer.RealTime();
    double cpuTime = timer.CpuTime();
    
    csvFile << "# Total execution time: " << realTime << " seconds" << std::endl;
    csvFile.close();
    
    std::cout << "Data saved to " << fullPath << std::endl;
    std::cout << "Sweep completed in " << realTime << " seconds (wall clock)" << std::endl;
    std::cout << "CPU time used: " << cpuTime << " seconds" << std::endl;
    
    // Create canvas and draw the histogram
    TCanvas* c = new TCanvas("c", "Detector Flux Map", 1000, 800);
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

    if (notify) {
        // NOTIFICATION SYSTEM - Play a sound when sweep is complete
        std::cout << "\n***** SWEEP COMPLETE *****\n" << std::endl;
        
        // Terminal bell (works in most terminals including WSL)
        std::cout << '\a' << std::endl;
    }
    
    // Clean up manager
    delete manager;
}

// Add this function to provide a more efficient sweepDetector

void sweepDetectorOptimized(bool notify = true, const char* saveFolder = "results", int threads = -1) {
    // Initialize ROOT threading
    TThread::Initialize();
    
    // Determine thread count - use fewer threads if they're causing slowdowns
    int maxPositionThreads = (threads > 0) ? threads : 1;
    if (maxPositionThreads > 4) maxPositionThreads = 4;  // Cap at 4 threads max
    
    std::cout << "Using " << maxPositionThreads << " threads for position processing" << std::endl;
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Set up geometry BEFORE enabling threading
    setupOpticsManager(manager, MAX_REFLECTIONS, 0.75, 0.99, false);
    
    // Only set max threads AFTER geometry is closed
    if (std::thread::hardware_concurrency() > 1) {
        int maxRayThreads = std::max(1, (int)std::thread::hardware_concurrency()/2);
        manager->SetMaxThreads(maxRayThreads);
        std::cout << "Using " << maxRayThreads << " threads for each ray tracing task" << std::endl;
    }
    
    // Ray count - configurable
    int n = 50000;
    double exitPortZ = -100*cm;
    
    // Configure angular binning
    const int nThetaBins = 180;
    const int nPhiBins = 90;
    
    // File setup (unchanged)
    std::string filename = "detector_data_" + std::to_string(n) + "rays.csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    
    // Create directory
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    int result = system(mkdirCmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Could not create directory: " << saveFolder << std::endl;
    } else {
        std::cout << "Using directory: " << saveFolder << std::endl;
    }
    
    fullPath = getUniqueFilename(fullPath);
    
    // Open CSV file
    std::mutex csvMutex;
    std::ofstream csvFile(fullPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
        return;
    }
    
    // Write metadata as comments
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    csvFile << "# Detector Data - Generated: " << timeBuffer << std::endl;
    csvFile << "# Number of rays: " << n << std::endl;
    csvFile << "# Detector dimensions: 20cm x 20cm" << std::endl;
    csvFile << "# Sphere inner radius: 100.1cm" << std::endl;
    csvFile << "# Sphere outer radius: 101cm" << std::endl;
    csvFile << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvFile << "# Theta bins: " << nThetaBins << std::endl;
    csvFile << "# Phi bins: " << nPhiBins << std::endl;
    csvFile << "# Mirror reflectance: 0.99" << std::endl;
    csvFile << "# Gaussian roughness: 0.75" << std::endl;
    csvFile << "# Lambertian scattering: enabled" << std::endl;
    csvFile << "theta,phi,fraction" << std::endl;
    
    // Create histogram
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    int totalPositions = nThetaBins * nPhiBins;
    std::cout << "\nStarting detector sweep with " << n << " rays per position "
              << "(" << totalPositions << " positions total)..." << std::endl;
    
    // Use chunked processing - divide the space into chunks that threads can work on
    int chunksPerThread = 20;  // Each thread gets 20 chunks to process
    int totalChunks = maxPositionThreads * chunksPerThread;
    int positionsPerChunk = (totalPositions + totalChunks - 1) / totalChunks;
    
    std::cout << "Dividing work into " << totalChunks << " chunks (" 
              << positionsPerChunk << " positions per chunk)" << std::endl;
    
    // Prepare worker queue with chunks instead of individual positions
    struct WorkChunk {
        int startIdx;
        int count;
    };
    
    std::vector<WorkChunk> chunks;
    int remaining = totalPositions;
    int startIdx = 0;
    
    while (remaining > 0) {
        int count = std::min(positionsPerChunk, remaining);
        chunks.push_back({startIdx, count});
        remaining -= count;
        startIdx += count;
    }
    
    std::atomic<int> chunkIndex(0);
    std::atomic<int> completedPositions(0);
    std::mutex consoleMutex;
    
    // Shared timer
    TStopwatch timer;
    timer.Start();
    
    // Worker function for chunk processing
    auto processChunk = [&](int threadId) {
        // Each thread processes chunks of positions
        Detector detector(20*cm, 20*cm);
        
        while (true) {
            // Get next chunk
            int currentChunkIdx = chunkIndex++;
            if (currentChunkIdx >= chunks.size()) break;
            
            WorkChunk chunk = chunks[currentChunkIdx];
            int startPos = chunk.startIdx;
            int endPos = startPos + chunk.count - 1;
            
            {
                std::lock_guard<std::mutex> lock(consoleMutex);
                std::cout << "Thread " << threadId << " processing chunk " << currentChunkIdx
                         << " (positions " << startPos << " to " << endPos << ")..." << std::endl;
            }
            
            for (int pos = startPos; pos <= endPos; pos++) {
                int i = pos / nPhiBins;
                int j = pos % nPhiBins;
                
                double theta = (i + 0.5) * 90.0/nThetaBins;
                double phi = (j + 0.5) * 360.0/nPhiBins;
                
                // Reset detector and position it
                detector.hitCount = 0;
                detector.setPosition(theta, phi, 100*cm);
                
                // Process rays
                TStopwatch posTimer;
                posTimer.Start();
                
                // No need for mutex here - each thread has its own detector
                traceRaysParallel(manager, n, exitPortZ, detector, false);
                
                posTimer.Stop();
                double fraction = double(detector.hitCount)/double(n);
                
                // Update histogram and CSV (thread-safe)
                {
                    std::lock_guard<std::mutex> lock(csvMutex);
                    fluxMap->SetBinContent(i+1, j+1, fraction);
                    csvFile << std::fixed << std::setprecision(6) 
                           << theta << "," << phi << "," << fraction << std::endl;
                }
                
                // Update progress
                int current = ++completedPositions;
                double percentComplete = (100.0 * current) / totalPositions;
                
                // Occasional progress update (not for every position to reduce output)
                if (current % 10 == 0 || current == totalPositions || current <= 10) {
                    std::lock_guard<std::mutex> lock(consoleMutex);
                    std::cout << "Position (θ=" << theta << "°, φ=" << phi << "°): " 
                             << detector.hitCount << "/" << n << " rays (" 
                             << std::fixed << std::setprecision(2) << fraction * 100 << "%)" 
                             << " in " << posTimer.RealTime() << "s" << std::endl;
                    
                    // Full progress report
                    if (current % 50 == 0 || current == totalPositions) {
                        // Calculate estimated time
                        double elapsedTime = timer.RealTime();
                        double avgTimePerPos = elapsedTime / current;
                        double remainingTime = avgTimePerPos * (totalPositions - current);
                        
                        // Apply safety factor based on progress
                        if (percentComplete < 5.0) remainingTime *= 1.5;
                        else if (percentComplete < 20.0) remainingTime *= 1.2;
                        
                        // Format times
                        int hours = static_cast<int>(remainingTime / 3600);
                        int minutes = static_cast<int>((remainingTime - hours * 3600) / 60);
                        int seconds = static_cast<int>(remainingTime - hours * 3600 - minutes * 60);
                        
                        std::cout << "\nProgress: " << std::fixed << std::setprecision(1) 
                                 << percentComplete << "% (" << current << "/" << totalPositions 
                                 << " positions)" << std::endl;
                        
                        if (hours > 0) {
                            std::cout << "Remaining time: " << hours << "h " 
                                     << minutes << "m " << seconds << "s" << std::endl;
                        } else if (minutes > 0) {
                            std::cout << "Remaining time: " << minutes << "m " 
                                     << seconds << "s" << std::endl;
                        } else {
                            std::cout << "Remaining time: " << seconds << "s" << std::endl;
                        }
                        
                        // Show estimated completion time
                        time_t finishTime = now + static_cast<time_t>(remainingTime);
                        struct tm* finishTm = localtime(&finishTime);
                        strftime(timeBuffer, sizeof(timeBuffer), "%H:%M:%S (%Y-%m-%d)", finishTm);
                        std::cout << "Est. completion at: " << timeBuffer << "\n" << std::endl;
                    }
                }
            }
        }
    };
    
    // Start worker threads
    std::vector<std::thread> workers;
    for (int t = 0; t < maxPositionThreads; t++) {
        workers.emplace_back(processChunk, t);
    }
    
    // Wait for all threads to finish
    for (auto& thread : workers) {
        thread.join();
    }
    
    // After sweep completes
    timer.Stop();
    double realTime = timer.RealTime();
    
    csvFile << "# Total execution time: " << realTime << " seconds" << std::endl;
    csvFile.close();
    
    std::cout << "Data saved to " << fullPath << std::endl;
    std::cout << "Sweep completed in " << realTime << " seconds (wall clock)" << std::endl;
    
    // Create canvas and draw
    TCanvas* c = new TCanvas("c", "Detector Flux Map", 1000, 800);
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

    if (notify) {
        // Play a sound when sweep is complete
        std::cout << "\n***** SWEEP COMPLETE *****\n" << std::endl;
        std::cout << '\a' << std::endl;
    }
    
    delete manager;
}

// Add a new function with improved monitoring capabilities

void sweepDetectorWithMonitor(bool notify = true, const char* saveFolder = "results", int threads = -1) {
    // Initialize ROOT threading
    TThread::Initialize();
    
    // Determine thread count - use fewer threads if they're causing slowdowns
    int maxPositionThreads = (threads > 0) ? threads : 1;
    if (maxPositionThreads > 4) maxPositionThreads = 4;  // Cap at 4 threads max
    
    std::cout << "Using " << maxPositionThreads << " threads for position processing" << std::endl;
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Set up geometry BEFORE enabling threading
    setupOpticsManager(manager, MAX_REFLECTIONS, 0.75, 0.99, false);
    
    // Only set max threads AFTER geometry is closed
    if (std::thread::hardware_concurrency() > 1) {
        int maxRayThreads = std::max(1, (int)std::thread::hardware_concurrency()/2);
        manager->SetMaxThreads(maxRayThreads);
        std::cout << "Using " << maxRayThreads << " threads for each ray tracing task" << std::endl;
    }
    
    // Ray count and binning configuration
    int n = 50000;
    double exitPortZ = -100*cm;
    const int nThetaBins = 180;
    const int nPhiBins = 90;
    
    // File setup
    std::string filename = "detector_data_" + std::to_string(n) + "rays.csv";
    std::string fullPath = std::string(saveFolder) + "/" + filename;
    std::string mkdirCmd = "mkdir -p \"" + std::string(saveFolder) + "\"";
    int result = system(mkdirCmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Could not create directory: " << saveFolder << std::endl;
    } else {
        std::cout << "Using directory: " << saveFolder << std::endl;
    }
    
    fullPath = getUniqueFilename(fullPath);
    
    // Open CSV file and write metadata
    std::mutex csvMutex;
    std::ofstream csvFile(fullPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open file " << fullPath << " for writing." << std::endl;
        return;
    }
    
    time_t now = time(nullptr);
    struct tm* timeinfo = localtime(&now);
    char timeBuffer[80];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    csvFile << "# Detector Data - Generated: " << timeBuffer << std::endl;
    csvFile << "# Number of rays: " << n << std::endl;
    csvFile << "# Detector dimensions: 20cm x 20cm" << std::endl;
    csvFile << "# Sphere inner radius: 100.1cm" << std::endl;
    csvFile << "# Sphere outer radius: 101cm" << std::endl;
    csvFile << "# Exit port angle: " << thetaMax << " degrees" << std::endl;
    csvFile << "# Theta bins: " << nThetaBins << std::endl;
    csvFile << "# Phi bins: " << nPhiBins << std::endl;
    csvFile << "# Mirror reflectance: 0.99" << std::endl;
    csvFile << "# Gaussian roughness: 0.75" << std::endl;
    csvFile << "# Lambertian scattering: enabled" << std::endl;
    csvFile << "theta,phi,fraction" << std::endl;
    
    // Create histogram
    TH2D* fluxMap = new TH2D("fluxMap", "Detector Flux Map;#theta (deg);#phi (deg)", 
                            nThetaBins, 0, 90, nPhiBins, 0, 360);
    
    int totalPositions = nThetaBins * nPhiBins;
    std::cout << "\nStarting detector sweep with " << n << " rays per position "
              << "(" << totalPositions << " positions total)..." << std::endl;
    
    // Set up chunking for work distribution
    int chunksPerThread = 20;  // Each thread gets 20 chunks to process
    int totalChunks = maxPositionThreads * chunksPerThread;
    int positionsPerChunk = (totalPositions + totalChunks - 1) / totalChunks;
    
    std::cout << "Dividing work into " << totalChunks << " chunks (" 
              << positionsPerChunk << " positions per chunk)" << std::endl;
    
    // Prepare worker queue with chunks
    struct WorkChunk {
        int startIdx;
        int count;
    };
    
    std::vector<WorkChunk> chunks;
    int remaining = totalPositions;
    int startIdx = 0;
    
    while (remaining > 0) {
        int count = std::min(positionsPerChunk, remaining);
        chunks.push_back({startIdx, count});
        remaining -= count;
        startIdx += count;
    }
    
    // Shared variables and synchronization primitives
    std::atomic<int> chunkIndex(0);
    std::atomic<int> completedPositions(0);
    std::atomic<bool> processingComplete(false);
    std::mutex consoleMutex;
    
    // Performance tracking variables
    struct ThreadStats {
        int threadId;
        int completedChunks;
        int completedPositions;
        double totalProcessingTime;
        double lastChunkStartTime;
        bool isProcessing;
    };
    
    std::vector<ThreadStats> threadStats(maxPositionThreads);
    std::mutex statsLock;
    
    // Initialize thread stats
    for (int i = 0; i < maxPositionThreads; i++) {
        threadStats[i] = {i, 0, 0, 0.0, 0.0, false};
    }
    
    // Overall timer
    TStopwatch timer;
    timer.Start();
    
    // Worker function for chunk processing
    auto processChunk = [&](int threadId) {
        // Each thread processes chunks of positions
        Detector detector(20*cm, 20*cm);
        
        while (true) {
            // Get next chunk
            int currentChunkIdx = chunkIndex++;
            if (currentChunkIdx >= chunks.size()) break;
            
            // Record chunk start time
            {
                std::lock_guard<std::mutex> lock(statsLock);
                threadStats[threadId].isProcessing = true;
                threadStats[threadId].lastChunkStartTime = timer.RealTime();
            }
            
            WorkChunk chunk = chunks[currentChunkIdx];
            int startPos = chunk.startIdx;
            int endPos = startPos + chunk.count - 1;
            
            {
                std::lock_guard<std::mutex> lock(consoleMutex);
                std::cout << "Thread " << threadId << " processing chunk " << currentChunkIdx
                          << " (positions " << startPos << " to " << endPos << ")..." << std::endl;
            }
            
            double chunkStartTime = timer.RealTime();
            int positionsInChunk = 0;
            
            for (int pos = startPos; pos <= endPos; pos++) {
                int i = pos / nPhiBins;
                int j = pos % nPhiBins;
                
                double theta = (i + 0.5) * 90.0/nThetaBins;
                double phi = (j + 0.5) * 360.0/nPhiBins;
                
                // Reset detector and position it
                detector.hitCount = 0;
                detector.setPosition(theta, phi, 100*cm);
                
                // Process rays
                TStopwatch posTimer;
                posTimer.Start();
                
                // No need for mutex - each thread has its own detector
                traceRaysParallel(manager, n, exitPortZ, detector, false);
                
                posTimer.Stop();
                double fraction = double(detector.hitCount)/double(n);
                
                // Update histogram and CSV (thread-safe)
                {
                    std::lock_guard<std::mutex> lock(csvMutex);
                    fluxMap->SetBinContent(i+1, j+1, fraction);
                    csvFile << std::fixed << std::setprecision(6) 
                           << theta << "," << phi << "," << fraction << std::endl;
                }
                
                // Update progress
                int current = ++completedPositions;
                positionsInChunk++;
                
                // Occasional progress update (not for every position)
                if (current % 10 == 0 || current == totalPositions || current <= 5) {
                    std::lock_guard<std::mutex> lock(consoleMutex);
                    std::cout << "Thread " << threadId << " - Position (θ=" << theta << "°, φ=" << phi << "°): " 
                             << detector.hitCount << "/" << n << " rays (" 
                             << std::fixed << std::setprecision(2) << fraction * 100 << "%)" 
                             << " in " << posTimer.RealTime() << "s" << std::endl;
                }
            }
            
            double chunkTime = timer.RealTime() - chunkStartTime;
            
            // Update thread statistics after chunk completion
            {
                std::lock_guard<std::mutex> lock(statsLock);
                threadStats[threadId].completedChunks++;
                threadStats[threadId].completedPositions += positionsInChunk;
                threadStats[threadId].totalProcessingTime += chunkTime;
                threadStats[threadId].isProcessing = false;
            }
        }
    };
    
    // Monitoring thread function
    auto monitorProgress = [&]() {
        int lastReportedPositions = 0;
        
        // Wait a little before starting to monitor
        std::this_thread::sleep_for(std::chrono::seconds(5));
        
        // Continue monitoring until processing is complete
        while (!processingComplete) {
            // Sleep between updates to avoid excessive reporting
            std::this_thread::sleep_for(std::chrono::seconds(10));
            
            int current;
            std::vector<ThreadStats> currentStats;
            
            // Make thread-safe copy of the stats
            {
                std::lock_guard<std::mutex> lock(statsLock);
                currentStats = threadStats;
                current = completedPositions;
            }
            
            // Only report if we've made progress
            if (current == lastReportedPositions) {
                continue;
            }
            
            lastReportedPositions = current;
            double percentComplete = (100.0 * current) / totalPositions;
            double elapsedTime = timer.RealTime();
            
            // Calculate average position processing time across active threads
            double totalPositionsProcessed = 0;
            double totalThreadTime = 0;
            double activeThreadCount = 0;
            
            for (const auto& stats : currentStats) {
                if (stats.completedPositions > 0) {
                    totalPositionsProcessed += stats.completedPositions;
                    totalThreadTime += stats.totalProcessingTime;
                    activeThreadCount++;
                }
            }
            
            double avgTimePerPosition = totalPositionsProcessed > 0 ? 
                totalThreadTime / totalPositionsProcessed : 0;
            
            // Calculate completion estimate based on actual thread performance
            double effectiveThreads = std::min(maxPositionThreads, (int)activeThreadCount);
            
            // If we don't have enough data yet, use a conservative estimate
            if (effectiveThreads < 1 || avgTimePerPosition == 0) {
                effectiveThreads = 1;
                avgTimePerPosition = elapsedTime / (current > 0 ? current : 1);
            }
            
            // Calculate estimated remaining time
            double remainingPositions = totalPositions - current;
            double positionsPerThread = remainingPositions / effectiveThreads;
            double rawRemainingTime = positionsPerThread * avgTimePerPosition;
            
            // Apply safety factor based on progress
            double remainingTime;
            if (percentComplete < 5.0) {
                remainingTime = rawRemainingTime * 1.5; // Very early - be very conservative
            } else if (percentComplete < 20.0) {
                remainingTime = rawRemainingTime * 1.2; // Early - be somewhat conservative
            } else {
                remainingTime = rawRemainingTime;       // Later - use the raw estimate
            }
            
            // Format for display
            int hours = static_cast<int>(remainingTime / 3600);
            int minutes = static_cast<int>((remainingTime - hours * 3600) / 60);
            int seconds = static_cast<int>(remainingTime - hours * 3600 - minutes * 60);
            
            int elapsed_hours = static_cast<int>(elapsedTime / 3600);
            int elapsed_minutes = static_cast<int>((elapsedTime - elapsed_hours * 3600) / 60);
            int elapsed_seconds = static_cast<int>(elapsedTime - elapsed_hours * 3600 - elapsed_minutes * 60);
            
            // Display detailed progress report
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::cout << "\n=== MONITOR: PROGRESS UPDATE ===" << std::endl;
            std::cout << "Progress: " << std::fixed << std::setprecision(1) 
                      << percentComplete << "% (" << current << "/" << totalPositions 
                      << " positions)" << std::endl;
            
            std::cout << "Elapsed time: ";
            if (elapsed_hours > 0) std::cout << elapsed_hours << "h ";
            if (elapsed_hours > 0 || elapsed_minutes > 0) std::cout << elapsed_minutes << "m ";
            std::cout << elapsed_seconds << "s" << std::endl;
            
            std::cout << "Average position time: " << avgTimePerPosition << "s" << std::endl;
            std::cout << "Effective threads: " << effectiveThreads << std::endl;
            
            std::cout << "Remaining time: ";
            if (hours > 0) std::cout << hours << "h ";
            if (hours > 0 || minutes > 0) std::cout << minutes << "m ";
            std::cout << seconds << "s" << std::endl;
            
            time_t finishTime = time(nullptr) + static_cast<time_t>(remainingTime);
            struct tm* finishTm = localtime(&finishTime);
            strftime(timeBuffer, sizeof(timeBuffer), "%H:%M:%S (%Y-%m-%d)", finishTm);
            std::cout << "Estimated completion at: " << timeBuffer << std::endl;
            
            // Show per-thread statistics
            std::cout << "\nPer-thread statistics:" << std::endl;
            for (const auto& stats : currentStats) {
                std::cout << "Thread " << stats.threadId 
                          << ": " << stats.completedPositions << " positions in " 
                          << stats.totalProcessingTime << "s";
                if (stats.completedPositions > 0) {
                    std::cout << " (avg: " << stats.totalProcessingTime / stats.completedPositions 
                              << "s/position)";
                }
                std::cout << std::endl;
            }
            std::cout << "==========================\n" << std::endl;
        }
    };
    
    // Start worker threads
    std::vector<std::thread> workers;
    for (int t = 0; t < maxPositionThreads; t++) {
        workers.emplace_back(processChunk, t);
    }
    
    // Start monitoring thread
    std::thread monitorThread(monitorProgress);
    
    // Wait for all worker threads to finish
    for (auto& thread : workers) {
        thread.join();
    }
    
    // Signal that processing is complete and wait for monitor thread
    processingComplete = true;
    monitorThread.join();
    
    // After sweep completes
    timer.Stop();
    double realTime = timer.RealTime();
    
    csvFile << "# Total execution time: " << realTime << " seconds" << std::endl;
    csvFile.close();
    
    std::cout << "Data saved to " << fullPath << std::endl;
    std::cout << "Sweep completed in " << realTime << " seconds (wall clock)" << std::endl;
    
    // Create canvas and draw
    TCanvas* c = new TCanvas("c", "Detector Flux Map", 1000, 800);
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

    if (notify) {
        // Play a sound when sweep is complete
        std::cout << "\n***** SWEEP COMPLETE *****\n" << std::endl;
        std::cout << '\a' << std::endl;
    }
    
    delete manager;
}


void visualizeDetector(double theta = 45.0, double phi = 0.0) {
    // Create canvas with square proportions
    TCanvas* c = new TCanvas("cvis", "Detector Visualization", 800, 800);
    c->cd();
    
    AOpticsManager* manager = new AOpticsManager("manager", "spherical shell");
    
    // Create standard geometry
    TGeoBBox* box = new TGeoBBox("box", 200*cm, 200*cm, 200*cm);
    AOpticalComponent* world = new AOpticalComponent("world", box);
    manager->SetTopVolume(world);
    
    TGeoSphere* sphere = new TGeoSphere("sphereWithExitPort", 100.1*cm, 101*cm, 0., thetaMax);
    AMirror* mirror = new AMirror("mirror", sphere);
    
    // Add the border surface condition for mirror reflectivity
    ABorderSurfaceCondition* condition = new ABorderSurfaceCondition(world, mirror);
    condition->EnableLambertian(false);
    // condition->SetGaussianRoughness(0.00000001);  
    
    world->AddNode(mirror, 1);
    
    // Add detector to geometry
    Detector detector(20*cm, 20*cm);
    detector.setPosition(theta, phi, 100*cm);
    detector.AddToGeometry(world);
    
    // Configure and close geometry
    manager->SetNsegments(60);
    manager->CloseGeometry();
    
    // Draw with OpenGL
    gStyle->SetCanvasPreferGL(true);
    world->Draw("ogl");
    
    // Update viewer safely
    TGLViewer* glv = (TGLViewer*)gPad->GetViewer3D();
    if (glv) {
        glv->SetStyle(TGLRnrCtx::kWireFrame);
        glv->UpdateScene();
        c->Update();
    }
    
    // Let the viewer initialize completely before ray tracing
    gSystem->ProcessEvents();
    
    // Now trace rays with statistics
    int n = 1000;
    double exitPortZ = -100*cm;
    
    // First run traceRaysWithStats to get statistics
    std::cout << "\nTracing " << n << " rays with detailed statistics:" << std::endl;
    int hitCount = traceRaysWithStats(manager, n, exitPortZ, detector);
    
    // Reset counters for visualization
    detector.hitCount = 0;
    int exitCount = 0;
    int suspendedCount = 0;
    int absorbedCount = 0;
    
    // Use individual rays for better stability
    for (int i = 0; i < n; ++i) {
        // Match parameters with working code
        // Ray parameters
        double x = -60*cm;
        ARay* ray = new ARay(i, 660*nm, x, 0*cm, -80*cm, 0, 5, 0, 0);
        
        manager->TraceNonSequential(*ray);
        
        TPolyLine3D* pol = ray->MakePolyLine3D();
        
        Double_t lastPoint[3];
        ray->GetLastPoint(lastPoint);
        
        if (ray->IsSuspended()) {
            pol->SetLineColor(kMagenta);  // Suspended rays in magenta
            suspendedCount++;
        } else if (ray->IsAbsorbed()) {
            pol->SetLineColor(kBlack);    // Absorbed rays in black
            absorbedCount++;
        } else if (lastPoint[2] < exitPortZ) {
            exitCount++;
            if (detector.checkIntersection(ray)) {
                pol->SetLineColor(kGreen);
                detector.hitCount++;
            } else {
                pol->SetLineColor(kYellow);
            }
        } else {
            pol->SetLineColor(kRed);    // Didn't exit
        }
        
        pol->SetLineWidth(2);
        pol->Draw();
        
        // Add debugging to see actual internal calculations
        // std::cout << "Detector at: (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << std::endl;
        // std::cout << "Ray exit at: (" << lastPoint[0]/cm << ", " << lastPoint[1]/cm << ", " << lastPoint[2]/cm << ") cm" << std::endl;

        delete ray;
        
        // Update periodically
        if (i % 10 == 0) {
            c->Update();
            gSystem->ProcessEvents();
        }
    }
    
    c->Update();
    
    // Print complete detector information
    std::cout << "\nDetector Information:" << std::endl;
    std::cout << "Total rays: " << n << std::endl;
    std::cout << "Rays exiting port: " << exitCount << std::endl;
    std::cout << "Rays hitting detector: " << detector.hitCount << "/" << exitCount << " of exiting rays" << std::endl;
    std::cout << "Rays suspended: " << suspendedCount << std::endl;
    std::cout << "Rays absorbed: " << absorbedCount << std::endl;
    std::cout << "Rays reflected back: " << (n - exitCount - suspendedCount - absorbedCount) << std::endl;
    std::cout << "Angular position: theta = " << theta << "°, phi = " << phi << "°" << std::endl;
    std::cout << "Position (x,y,z): (" << detector.x/cm << ", " << detector.y/cm << ", " << detector.z/cm << ") cm" << std::endl;
    
    // Legend for ray colors
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