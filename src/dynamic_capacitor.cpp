#include <embree4/rtcore.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

// Physical constants for capacitor calculation
const double EPSILON_0 = 8.854e-12; // F/m (permittivity of free space)
const double GLYCERIN_RELATIVE_PERMITTIVITY = 42.28; // Relative permittivity of glycerin

struct Triangle {
    float v0[3], v1[3], v2[3];
};

struct Mesh {
    std::vector<Triangle> triangles;
};

struct Vec3 {
    float x, y, z;
    Vec3(float x_ = 0, float y_ = 0, float z_ = 0) : x(x_), y(y_), z(z_) {}
};

struct Matrix4x4 {
    float m[4][4];
    
    Matrix4x4() {
        // Initialize as identity matrix
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }
    
    Vec3 transformPoint(const Vec3& point) const {
        float x = m[0][0] * point.x + m[0][1] * point.y + m[0][2] * point.z + m[0][3];
        float y = m[1][0] * point.x + m[1][1] * point.y + m[1][2] * point.z + m[1][3];
        float z = m[2][0] * point.x + m[2][1] * point.y + m[2][2] * point.z + m[2][3];
        return Vec3(x, y, z);
    }
    
    void print() const {
        std::cout << "Matrix:" << std::endl;
        for (int i = 0; i < 4; i++) {
            std::cout << "  [";
            for (int j = 0; j < 4; j++) {
                std::cout << std::setw(10) << std::fixed << std::setprecision(6) << m[i][j];
                if (j < 3) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
    }
};

struct CapacitorResult {
    int step;
    double capacitance_pF;
    double minDistance_mm;
    double maxDistance_mm;
    double totalArea_mm2;
    int hits;
    int misses;
    Vec3 translation;
    
    void print() const {
        std::cout << "Step " << step << ": " 
                  << capacitance_pF << " pF, "
                  << "gap: " << minDistance_mm << " mm, "
                  << "area: " << totalArea_mm2 << " mm^2, "
                  << "hits: " << hits << "/" << (hits + misses) << std::endl;
    }
};

// Load transformation matrices from file
std::vector<Matrix4x4> loadTransformations(const std::string& filename) {
    std::vector<Matrix4x4> matrices;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return matrices;
    }
    
    std::string line;
    int numMatrices = 0;
    
    // Read header
    while (std::getline(file, line)) {
        if (line.find("MATRICES") == 0) {
            std::istringstream iss(line);
            std::string keyword;
            iss >> keyword >> numMatrices;
            break;
        }
    }
    
    std::cout << "Loading " << numMatrices << " matrices from " << filename << std::endl;
    
    // Read matrices
    for (int m = 0; m < numMatrices; m++) {
        Matrix4x4 matrix;
        
        // Skip "MATRIX x" line
        std::getline(file, line);
        
        // Read 4 rows of matrix
        for (int i = 0; i < 4; i++) {
            if (!std::getline(file, line)) {
                std::cerr << "Error reading matrix " << m << " row " << i << std::endl;
                return matrices;
            }
            
            std::istringstream iss(line);
            for (int j = 0; j < 4; j++) {
                if (!(iss >> matrix.m[i][j])) {
                    std::cerr << "Error parsing matrix " << m << " element [" << i << "][" << j << "]" << std::endl;
                    return matrices;
                }
            }
        }
        
        matrices.push_back(matrix);
        
        // Skip empty line
        std::getline(file, line);
    }
    
    file.close();
    std::cout << "Successfully loaded " << matrices.size() << " transformation matrices" << std::endl;
    return matrices;
}

// Apply transformation to mesh
Mesh transformMesh(const Mesh& originalMesh, const Matrix4x4& transform) {
    Mesh transformedMesh;
    transformedMesh.triangles.reserve(originalMesh.triangles.size());
    
    for (const auto& tri : originalMesh.triangles) {
        Triangle newTri;
        
        Vec3 v0(tri.v0[0], tri.v0[1], tri.v0[2]);
        Vec3 v1(tri.v1[0], tri.v1[1], tri.v1[2]);
        Vec3 v2(tri.v2[0], tri.v2[1], tri.v2[2]);
        
        Vec3 tv0 = transform.transformPoint(v0);
        Vec3 tv1 = transform.transformPoint(v1);
        Vec3 tv2 = transform.transformPoint(v2);
        
        newTri.v0[0] = tv0.x; newTri.v0[1] = tv0.y; newTri.v0[2] = tv0.z;
        newTri.v1[0] = tv1.x; newTri.v1[1] = tv1.y; newTri.v1[2] = tv1.z;
        newTri.v2[0] = tv2.x; newTri.v2[1] = tv2.y; newTri.v2[2] = tv2.z;
        
        transformedMesh.triangles.push_back(newTri);
    }
    
    return transformedMesh;
}

// Calculate triangle area and center
void getTriangleInfo(const Triangle& tri, float center[3], float normal[3], float& area) {
    // Calculate center (average of vertices)
    for (int i = 0; i < 3; i++) {
        center[i] = (tri.v0[i] + tri.v1[i] + tri.v2[i]) / 3.0f;
    }
    
    // Calculate normal and area using cross product
    float edge1[3], edge2[3];
    for (int i = 0; i < 3; i++) {
        edge1[i] = tri.v1[i] - tri.v0[i];
        edge2[i] = tri.v2[i] - tri.v0[i];
    }
    
    // Cross product for normal
    normal[0] = edge1[1] * edge2[2] - edge1[2] * edge2[1];
    normal[1] = edge1[2] * edge2[0] - edge1[0] * edge2[2];
    normal[2] = edge1[0] * edge2[1] - edge1[1] * edge2[0];
    
    // Area is half the magnitude of cross product
    float length = sqrtf(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    area = length / 2.0f;
    
    // Normalize the normal
    if (length > 0) {
        normal[0] /= length;
        normal[1] /= length;
        normal[2] /= length;
    }
}

// Load mesh into Embree scene
RTCScene createScene(RTCDevice device, const Mesh& mesh) {
    RTCScene scene = rtcNewScene(device);
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
                                                       RTC_FORMAT_FLOAT3, 3 * sizeof(float),
                                                       mesh.triangles.size() * 3);
    
    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0,
                                                           RTC_FORMAT_UINT3, 3 * sizeof(unsigned),
                                                           mesh.triangles.size());
    
    for (size_t i = 0; i < mesh.triangles.size(); i++) {
        vertices[i*9 + 0] = mesh.triangles[i].v0[0];
        vertices[i*9 + 1] = mesh.triangles[i].v0[1];
        vertices[i*9 + 2] = mesh.triangles[i].v0[2];
        
        vertices[i*9 + 3] = mesh.triangles[i].v1[0];
        vertices[i*9 + 4] = mesh.triangles[i].v1[1];
        vertices[i*9 + 5] = mesh.triangles[i].v1[2];
        
        vertices[i*9 + 6] = mesh.triangles[i].v2[0];
        vertices[i*9 + 7] = mesh.triangles[i].v2[1];
        vertices[i*9 + 8] = mesh.triangles[i].v2[2];
        
        indices[i*3 + 0] = i*3 + 0;
        indices[i*3 + 1] = i*3 + 1;
        indices[i*3 + 2] = i*3 + 2;
    }
    
    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);
    
    return scene;
}

// Calculate capacitor for a single configuration
CapacitorResult calculateCapacitorStep(const Mesh& positiveMesh, const Mesh& negativeMesh, int step) {
    RTCDevice device = rtcNewDevice(nullptr);
    RTCScene negativeScene = createScene(device, negativeMesh);
    
    const float maxDistance = 10.0f; // 10mm max search distance
    
    CapacitorResult result;
    result.step = step;
    result.capacitance_pF = 0.0;
    result.minDistance_mm = 1e9;
    result.maxDistance_mm = 0;
    result.totalArea_mm2 = 0;
    result.hits = 0;
    result.misses = 0;
    result.translation = Vec3(0, 0, 0);
    
    // Process each triangle in positive model
    for (size_t i = 0; i < positiveMesh.triangles.size(); i++) {
        float center[3], normal[3], area_mm2;
        getTriangleInfo(positiveMesh.triangles[i], center, normal, area_mm2);
        
        // Setup ray
        RTCRayHit rayhit;
        rayhit.ray.org_x = center[0];
        rayhit.ray.org_y = center[1];
        rayhit.ray.org_z = center[2];
        rayhit.ray.dir_x = normal[0];
        rayhit.ray.dir_y = normal[1];
        rayhit.ray.dir_z = normal[2];
        rayhit.ray.tnear = 0.0f;
        rayhit.ray.tfar = maxDistance;
        rayhit.ray.mask = -1;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        
        rtcIntersect1(negativeScene, &rayhit);
        
        float distance_mm = (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) ? 
                            rayhit.ray.tfar : maxDistance;
        
        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
            result.hits++;
            result.minDistance_mm = std::min(result.minDistance_mm, (double)distance_mm);
            result.maxDistance_mm = std::max(result.maxDistance_mm, (double)distance_mm);
            result.totalArea_mm2 += area_mm2;
        } else {
            result.misses++;
        }
    }
    
    // Calculate capacitance for parallel plates with glycerin dielectric
    if (result.hits > 0 && result.minDistance_mm < 10.0) {
        double area_m2 = result.totalArea_mm2 * 1e-6;  // Convert mm^2 to m^2
        double distance_m = result.minDistance_mm * 1e-3;  // Convert mm to m
        result.capacitance_pF = (EPSILON_0 * GLYCERIN_RELATIVE_PERMITTIVITY * area_m2 / distance_m) * 1e12;  // Convert to pF
    }
    
    rtcReleaseScene(negativeScene);
    rtcReleaseDevice(device);
    
    return result;
}

// Simple OBJ loader with meter to millimeter conversion
Mesh loadOBJ(const std::string& filename) {
    Mesh mesh;
    
    FILE* file = fopen(filename.c_str(), "r");
    if (!file) {
        std::cerr << "Failed to open: " << filename << std::endl;
        return mesh;
    }
    
    char line[256];
    std::vector<Vec3> verts;
    
    std::cout << "Loading " << filename << " and converting from meters to millimeters..." << std::endl;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == 'v' && line[1] == ' ') {
            float x, y, z;
            sscanf(line, "v %f %f %f", &x, &y, &z);
            // CONVERT METERS TO MILLIMETERS
            verts.push_back(Vec3(x * 1000.0f, y * 1000.0f, z * 1000.0f));
        }
        else if (line[0] == 'f' && line[1] == ' ') {
            int i1, i2, i3;
            if (sscanf(line, "f %d %d %d", &i1, &i2, &i3) == 3) {
                Triangle tri;
                tri.v0[0] = verts[i1-1].x; tri.v0[1] = verts[i1-1].y; tri.v0[2] = verts[i1-1].z;
                tri.v1[0] = verts[i2-1].x; tri.v1[1] = verts[i2-1].y; tri.v1[2] = verts[i2-1].z;
                tri.v2[0] = verts[i3-1].x; tri.v2[1] = verts[i3-1].y; tri.v2[2] = verts[i3-1].z;
                mesh.triangles.push_back(tri);
            }
            else if (sscanf(line, "f %d/%*d/%*d %d/%*d/%*d %d/%*d/%*d", &i1, &i2, &i3) == 3 ||
                     sscanf(line, "f %d//%*d %d//%*d %d//%*d", &i1, &i2, &i3) == 3) {
                Triangle tri;
                tri.v0[0] = verts[i1-1].x; tri.v0[1] = verts[i1-1].y; tri.v0[2] = verts[i1-1].z;
                tri.v1[0] = verts[i2-1].x; tri.v1[1] = verts[i2-1].y; tri.v1[2] = verts[i2-1].z;
                tri.v2[0] = verts[i3-1].x; tri.v2[1] = verts[i3-1].y; tri.v2[2] = verts[i3-1].z;
                mesh.triangles.push_back(tri);
            }
        }
    }
    
    fclose(file);
    std::cout << "Loaded " << mesh.triangles.size() << " triangles from " << filename << std::endl;
    return mesh;
}

// Process all time steps for a sensor
void processSensor(const std::string& sensorName, const Mesh& positiveMesh, const Mesh& negativeMesh) {
    std::string transformFile = "transformations_" + sensorName + ".txt";
    std::vector<Matrix4x4> transformations = loadTransformations(transformFile);
    
    if (transformations.empty()) {
        std::cout << "No transformations found for sensor " << sensorName << std::endl;
        return;
    }
    
    std::cout << "\n=== Processing Sensor " << sensorName << " (Glycerin Dielectric) ===" << std::endl;
    std::cout << "Time steps: " << transformations.size() << std::endl;
    
    std::vector<CapacitorResult> results;
    std::string outputFile = "capacitance_" + sensorName + ".csv";
    std::ofstream csvFile(outputFile);
    
    // CSV header
    csvFile << "Step,Capacitance_pF,MinDistance_mm,MaxDistance_mm,Area_mm2,Hits,Misses,Translation_X,Translation_Y,Translation_Z" << std::endl;
    
    for (size_t step = 0; step < transformations.size(); step++) {
        // Transform the negative mesh (A2 relative to A1, so A1 stays fixed)
        Mesh transformedNegative = transformMesh(negativeMesh, transformations[step]);
        
        // Calculate capacitance
        CapacitorResult result = calculateCapacitorStep(positiveMesh, transformedNegative, step);
        
        // Extract translation from matrix
        result.translation.x = transformations[step].m[0][3];
        result.translation.y = transformations[step].m[1][3];
        result.translation.z = transformations[step].m[2][3];
        
        results.push_back(result);
        
        // Write to CSV
        csvFile << step << "," << result.capacitance_pF << "," << result.minDistance_mm << ","
                << result.maxDistance_mm << "," << result.totalArea_mm2 << "," << result.hits << ","
                << result.misses << "," << result.translation.x << "," << result.translation.y << ","
                << result.translation.z << std::endl;
        
        // Print progress for first few and every 10th step
        if (step < 3 || step % 10 == 0) {
            result.print();
        }
    }
    
    csvFile.close();
    
    // Summary statistics
    double avgCapacitance = 0;
    double minCapacitance = 1e9;
    double maxCapacitance = 0;
    int validResults = 0;
    
    for (const auto& r : results) {
        if (r.capacitance_pF > 0) {
            avgCapacitance += r.capacitance_pF;
            minCapacitance = std::min(minCapacitance, r.capacitance_pF);
            maxCapacitance = std::max(maxCapacitance, r.capacitance_pF);
            validResults++;
        }
    }
    
    if (validResults > 0) {
        avgCapacitance /= validResults;
        
        std::cout << "\n=== Sensor " << sensorName << " Summary (Glycerin εᵣ=" << GLYCERIN_RELATIVE_PERMITTIVITY << ") ===" << std::endl;
        std::cout << "Average capacitance: " << avgCapacitance << " pF" << std::endl;
        std::cout << "Min capacitance: " << minCapacitance << " pF" << std::endl;
        std::cout << "Max capacitance: " << maxCapacitance << " pF" << std::endl;
        std::cout << "Capacitance variation: " << (maxCapacitance - minCapacitance) << " pF" << std::endl;
        std::cout << "Valid results: " << validResults << "/" << results.size() << std::endl;
    }
    
    std::cout << "Results saved to: " << outputFile << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "=== Dynamic Capacitor Calculator (Glycerin Dielectric) ===" << std::endl;
    std::cout << "Using glycerin relative permittivity: " << GLYCERIN_RELATIVE_PERMITTIVITY << std::endl;
    
    // Load base meshes
    std::cout << "\nLoading base meshes..." << std::endl;
    Mesh positiveMesh = loadOBJ("positive_model.obj");
    Mesh negativeMesh = loadOBJ("negative_model.obj");
    
    if (positiveMesh.triangles.empty() || negativeMesh.triangles.empty()) {
        std::cerr << "Error: Could not load base models!" << std::endl;
        return 1;
    }
    
    // Process each sensor
    std::vector<std::string> sensors = {"A", "B", "C"};
    
    for (const std::string& sensor : sensors) {
        processSensor(sensor, positiveMesh, negativeMesh);
    }
    
    std::cout << "\n=== All sensors processed ===" << std::endl;
    std::cout << "Output files generated:" << std::endl;
    for (const std::string& sensor : sensors) {
        std::cout << "  - capacitance_" << sensor << ".csv" << std::endl;
    }
    
    return 0;
}