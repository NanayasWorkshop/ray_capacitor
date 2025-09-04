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

struct SensorConfig {
    std::string sensorName;
    std::string model1Name;
    std::string model2Name;
    std::string transformFile;
    Mesh model1Mesh;
    Mesh model2Mesh;
    std::vector<Matrix4x4> transformations;
};

struct RayData {
    float origin[3];
    float direction[3];
    float hit_point[3];
    bool hit;
    float distance;
};

// Export ray data for animation visualization
void exportRayDataForStep(const std::vector<RayData>& rays, const Mesh& posMesh, const Mesh& negMesh, 
                         const std::string& sensorName, int step) {
    std::string filename = "ray_data_" + sensorName + "_step_" + 
                          (step < 10 ? "00" : (step < 100 ? "0" : "")) + std::to_string(step) + ".txt";
    
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "  -> ERROR: Could not create file " << filename << std::endl;
        return;
    }
    
    // Export positive mesh
    file << "POSITIVE_MESH\n";
    for (const auto& tri : posMesh.triangles) {
        file << tri.v0[0] << " " << tri.v0[1] << " " << tri.v0[2] << " ";
        file << tri.v1[0] << " " << tri.v1[1] << " " << tri.v1[2] << " ";
        file << tri.v2[0] << " " << tri.v2[1] << " " << tri.v2[2] << "\n";
    }
    
    // Export negative mesh
    file << "NEGATIVE_MESH\n";
    for (const auto& tri : negMesh.triangles) {
        file << tri.v0[0] << " " << tri.v0[1] << " " << tri.v0[2] << " ";
        file << tri.v1[0] << " " << tri.v1[1] << " " << tri.v1[2] << " ";
        file << tri.v2[0] << " " << tri.v2[1] << " " << tri.v2[2] << "\n";
    }
    
    // Export rays
    file << "RAYS\n";
    for (const auto& ray : rays) {
        file << ray.origin[0] << " " << ray.origin[1] << " " << ray.origin[2] << " ";
        file << ray.hit_point[0] << " " << ray.hit_point[1] << " " << ray.hit_point[2] << " ";
        file << (ray.hit ? 1 : 0) << " " << ray.distance << "\n";
    }
    file.close();
}

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

// Calculate capacitor for a single configuration with ray data collection
CapacitorResult calculateCapacitorStep(const Mesh& positiveMesh, const Mesh& negativeMesh, int step, 
                                     bool exportRays, const std::string& sensorName) {
    RTCDevice device = rtcNewDevice(nullptr);
    RTCScene negativeScene = createScene(device, negativeMesh);
    
    const float maxDistance = 5.0f; // 5mm max search distance (REDUCED FROM 10mm)
    
    CapacitorResult result;
    result.step = step;
    result.capacitance_pF = 0.0;
    result.minDistance_mm = 1e9;
    result.maxDistance_mm = 0;
    result.totalArea_mm2 = 0;
    result.hits = 0;
    result.misses = 0;
    result.translation = Vec3(0, 0, 0);
    
    // Collect ray data for visualization if requested (ONLY HITS)
    std::vector<RayData> rayDataCollection;
    
    // Process each triangle in positive model
    for (size_t i = 0; i < positiveMesh.triangles.size(); i++) {
        float center[3], normal[3], area_mm2;
        getTriangleInfo(positiveMesh.triangles[i], center, normal, area_mm2);
        
        // TRY BOTH DIRECTIONS - shoot rays in both normal directions
        bool hit_found = false;
        float best_distance = maxDistance;
        float best_hit_point[3];
        float best_direction[3];
        
        // Direction 1: Original flipped normal
        {
            RTCRayHit rayhit;
            rayhit.ray.org_x = center[0];
            rayhit.ray.org_y = center[1];
            rayhit.ray.org_z = center[2];
            rayhit.ray.dir_x = -normal[0];
            rayhit.ray.dir_y = -normal[1];
            rayhit.ray.dir_z = -normal[2];
            rayhit.ray.tnear = 0.0f;
            rayhit.ray.tfar = maxDistance;
            rayhit.ray.mask = -1;
            rayhit.ray.flags = 0;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            
            rtcIntersect1(negativeScene, &rayhit);
            
            if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID && rayhit.ray.tfar < best_distance) {
                hit_found = true;
                best_distance = rayhit.ray.tfar;
                best_hit_point[0] = center[0] + (-normal[0]) * rayhit.ray.tfar;
                best_hit_point[1] = center[1] + (-normal[1]) * rayhit.ray.tfar;
                best_hit_point[2] = center[2] + (-normal[2]) * rayhit.ray.tfar;
                best_direction[0] = -normal[0];
                best_direction[1] = -normal[1];
                best_direction[2] = -normal[2];
            }
        }
        
        // Direction 2: Positive normal direction
        {
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
            
            if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID && rayhit.ray.tfar < best_distance) {
                hit_found = true;
                best_distance = rayhit.ray.tfar;
                best_hit_point[0] = center[0] + normal[0] * rayhit.ray.tfar;
                best_hit_point[1] = center[1] + normal[1] * rayhit.ray.tfar;
                best_hit_point[2] = center[2] + normal[2] * rayhit.ray.tfar;
                best_direction[0] = normal[0];
                best_direction[1] = normal[1];
                best_direction[2] = normal[2];
            }
        }
        
        // Update results based on best hit (if any)
        if (hit_found) {
            result.hits++;
            result.minDistance_mm = std::min(result.minDistance_mm, (double)best_distance);
            result.maxDistance_mm = std::max(result.maxDistance_mm, (double)best_distance);
            result.totalArea_mm2 += area_mm2;
            
            // ONLY EXPORT RAYS THAT HIT
            if (exportRays) {
                RayData rd;
                rd.origin[0] = center[0];
                rd.origin[1] = center[1];
                rd.origin[2] = center[2];
                rd.direction[0] = best_direction[0];
                rd.direction[1] = best_direction[1];
                rd.direction[2] = best_direction[2];
                rd.hit = true;
                rd.distance = best_distance;
                rd.hit_point[0] = best_hit_point[0];
                rd.hit_point[1] = best_hit_point[1];
                rd.hit_point[2] = best_hit_point[2];
                
                rayDataCollection.push_back(rd);
            }
        } else {
            result.misses++;
            // DO NOT EXPORT MISS RAYS
        }
    }
    
    // Calculate capacitance for parallel plates with glycerin dielectric
    if (result.hits > 0 && result.minDistance_mm < 5.0) {  // Updated max distance check
        double area_m2 = result.totalArea_mm2 * 1e-6;  // Convert mm^2 to m^2
        double distance_m = result.minDistance_mm * 1e-3;  // Convert mm to m
        result.capacitance_pF = (EPSILON_0 * GLYCERIN_RELATIVE_PERMITTIVITY * area_m2 / distance_m) * 1e12;  // Convert to pF
    }
    
    // Export ray data for animation if requested
    if (exportRays) {
        exportRayDataForStep(rayDataCollection, positiveMesh, negativeMesh, sensorName, step);
    }
    
    rtcReleaseScene(negativeScene);
    rtcReleaseDevice(device);
    
    return result;
}

// Enhanced OBJ loader with KeyShot compatibility
Mesh loadOBJ(const std::string& filename) {
    Mesh mesh;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open: " << filename << std::endl;
        return mesh;
    }
    
    std::string line;
    std::vector<Vec3> verts;
    int lineNumber = 0;
    
    std::cout << "Loading " << filename << " and converting from meters to millimeters..." << std::endl;
    
    while (std::getline(file, line)) {
        lineNumber++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == 'm' || line[0] == 'g' || line[0] == 'u') {
            continue;
        }
        
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        
        if (token == "v") {
            // Parse vertex with robust floating point parsing
            double x, y, z;
            if (iss >> x >> y >> z) {
                // MANUAL UNIT CONTROL:
                const bool KEYSHOT_UNITS_MM = true;
                
                if (KEYSHOT_UNITS_MM) {
                    // KeyShot format: already in mm, use as-is
                    verts.push_back(Vec3(x, y, z));
                    if (lineNumber == 7) {
                        std::cout << "Units: Using KeyShot mm format, vertex: (" << x << ", " << y << ", " << z << ")" << std::endl;
                    }
                } else {
                    // SolidEdge format: in meters, convert to mm
                    verts.push_back(Vec3(x * 1000.0f, y * 1000.0f, z * 1000.0f));
                    if (lineNumber == 7) {
                        std::cout << "Units: Converting meters to mm, vertex: (" << x*1000 << ", " << y*1000 << ", " << z*1000 << ")" << std::endl;
                    }
                }
            } else {
                std::cout << "Warning: Could not parse vertex on line " << lineNumber << ": " << line << std::endl;
            }
        }
        else if (token == "f") {
            // Parse face indices
            std::vector<int> indices;
            std::string faceToken;
            
            while (iss >> faceToken) {
                // Handle various face formats: "1", "1/2", "1/2/3", "1//3"
                size_t firstSlash = faceToken.find('/');
                std::string vertexIndex = (firstSlash != std::string::npos) ? 
                                         faceToken.substr(0, firstSlash) : faceToken;
                
                int index = std::stoi(vertexIndex);
                indices.push_back(index);
            }
            
            // Only process triangular faces
            if (indices.size() >= 3) {
                // Handle both positive and negative indices (KeyShot compatibility)
                int i1 = (indices[0] > 0) ? indices[0] - 1 : (int)verts.size() + indices[0];
                int i2 = (indices[1] > 0) ? indices[1] - 1 : (int)verts.size() + indices[1];
                int i3 = (indices[2] > 0) ? indices[2] - 1 : (int)verts.size() + indices[2];
                
                // Bounds checking
                if (i1 >= 0 && i1 < (int)verts.size() && 
                    i2 >= 0 && i2 < (int)verts.size() && 
                    i3 >= 0 && i3 < (int)verts.size()) {
                    
                    Triangle tri;
                    tri.v0[0] = verts[i1].x; tri.v0[1] = verts[i1].y; tri.v0[2] = verts[i1].z;
                    tri.v1[0] = verts[i2].x; tri.v1[1] = verts[i2].y; tri.v1[2] = verts[i2].z;
                    tri.v2[0] = verts[i3].x; tri.v2[1] = verts[i3].y; tri.v2[2] = verts[i3].z;
                    mesh.triangles.push_back(tri);
                } else {
                    std::cout << "Warning: Face indices out of bounds on line " << lineNumber 
                              << " (vertices loaded: " << verts.size() << ", indices: " 
                              << i1+1 << " " << i2+1 << " " << i3+1 << ")" << std::endl;
                }
                
                // Handle quads by creating a second triangle
                if (indices.size() >= 4) {
                    int i4 = (indices[3] > 0) ? indices[3] - 1 : (int)verts.size() + indices[3];
                    if (i4 >= 0 && i4 < (int)verts.size()) {
                        Triangle tri2;
                        tri2.v0[0] = verts[i1].x; tri2.v0[1] = verts[i1].y; tri2.v0[2] = verts[i1].z;
                        tri2.v1[0] = verts[i3].x; tri2.v1[1] = verts[i3].y; tri2.v1[2] = verts[i3].z;
                        tri2.v2[0] = verts[i4].x; tri2.v2[1] = verts[i4].y; tri2.v2[2] = verts[i4].z;
                        mesh.triangles.push_back(tri2);
                    }
                }
            } else {
                std::cout << "Warning: Face with insufficient vertices on line " << lineNumber << std::endl;
            }
        }
    }
    
    file.close();
    std::cout << "Loaded " << verts.size() << " vertices and " << mesh.triangles.size() 
              << " triangles from " << filename << std::endl;
    return mesh;
}

// Process a single model against the stationary negative plate
void processModel(const std::string& modelName, const Mesh& modelMesh, const Mesh& stationaryNegative,
                  const std::vector<Matrix4x4>& transformations, const std::string& sensorGroup, bool exportAnimation) {
    
    std::cout << "\n=== Processing Model " << modelName << " (Sensor Group " << sensorGroup << ") ===" << std::endl;
    std::cout << "Time steps: " << transformations.size() << std::endl;
    if (exportAnimation) {
        std::cout << "Animation export: ENABLED" << std::endl;
    }
    
    std::vector<CapacitorResult> results;
    std::string outputFile = "capacitance_" + modelName + ".csv";
    std::ofstream csvFile(outputFile);
    
    // CSV header
    csvFile << "Step,Capacitance_pF,MinDistance_mm,MaxDistance_mm,Area_mm2,Hits,Misses,Translation_X,Translation_Y,Translation_Z" << std::endl;
    
    for (size_t step = 0; step < transformations.size(); step++) {
        // Transform the moving model
        Mesh transformedModel = transformMesh(modelMesh, transformations[step]);
        
        // Calculate capacitance against stationary negative plate
        CapacitorResult result = calculateCapacitorStep(transformedModel, stationaryNegative, step, 
                                                       exportAnimation, modelName);
        
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
        
        std::cout << "\n=== Model " << modelName << " Summary (Glycerin εᵣ=" << GLYCERIN_RELATIVE_PERMITTIVITY << ") ===" << std::endl;
        std::cout << "Average capacitance: " << avgCapacitance << " pF" << std::endl;
        std::cout << "Min capacitance: " << minCapacitance << " pF" << std::endl;
        std::cout << "Max capacitance: " << maxCapacitance << " pF" << std::endl;
        std::cout << "Capacitance variation: " << (maxCapacitance - minCapacitance) << " pF" << std::endl;
        std::cout << "Valid results: " << validResults << "/" << results.size() << std::endl;
    }
    
    std::cout << "Results saved to: " << outputFile << std::endl;
    
    if (exportAnimation) {
        std::cout << "Animation data exported: ray_data_" << modelName << "_step_XXX.txt" << std::endl;
    }
}

// Initialize sensor configurations
std::vector<SensorConfig> initializeSensors() {
    std::vector<SensorConfig> sensors;
    
    std::vector<std::string> sensorGroups = {"A", "B", "C"};
    
    for (const std::string& group : sensorGroups) {
        SensorConfig config;
        config.sensorName = group;
        config.model1Name = group + "1";
        config.model2Name = group + "2";
        config.transformFile = "transformations_" + group + ".txt";
        
        // Load transformation matrices
        config.transformations = loadTransformations(config.transformFile);
        
        if (config.transformations.empty()) {
            std::cout << "Warning: No transformations found for sensor " << group << std::endl;
            continue;
        }
        
        // Load model meshes
        std::string model1File = config.model1Name + "_model.obj";
        std::string model2File = config.model2Name + "_model.obj";
        
        config.model1Mesh = loadOBJ(model1File);
        config.model2Mesh = loadOBJ(model2File);
        
        if (config.model1Mesh.triangles.empty()) {
            std::cout << "Warning: Could not load " << model1File << std::endl;
            continue;
        }
        
        if (config.model2Mesh.triangles.empty()) {
            std::cout << "Warning: Could not load " << model2File << std::endl;
            continue;
        }
        
        sensors.push_back(config);
    }
    
    return sensors;
}

int main(int argc, char** argv) {
    bool exportAnimation = (argc > 1 && std::string(argv[1]) == "-anim");
    
    std::cout << "=== Multi-Sensor Dynamic Capacitor Calculator (Glycerin Dielectric) ===" << std::endl;
    std::cout << "Using glycerin relative permittivity: " << GLYCERIN_RELATIVE_PERMITTIVITY << std::endl;
    
    if (exportAnimation) {
        std::cout << "\n*** ANIMATION MODE ENABLED ***" << std::endl;
        std::cout << "Will export ray data for each time step" << std::endl;
        std::cout << "Use: python multi_sensor_animated_viewer.py animate" << std::endl;
    }
    
    std::cout << "\nExpected model files:" << std::endl;
    std::cout << "  - A1_model.obj, A2_model.obj" << std::endl;
    std::cout << "  - B1_model.obj, B2_model.obj" << std::endl;
    std::cout << "  - C1_model.obj, C2_model.obj" << std::endl;
    std::cout << "  - stationary_negative.obj" << std::endl;
    std::cout << "\nExpected transformation files:" << std::endl;
    std::cout << "  - transformations_A.txt" << std::endl;
    std::cout << "  - transformations_B.txt" << std::endl;
    std::cout << "  - transformations_C.txt" << std::endl;
    
    // Load stationary negative plate (same for all measurements)
    std::cout << "\nLoading stationary negative plate..." << std::endl;
    Mesh stationaryNegative = loadOBJ("stationary_negative.obj");
    
    if (stationaryNegative.triangles.empty()) {
        std::cerr << "Error: Could not load stationary_negative.obj!" << std::endl;
        return 1;
    }
    
    // Initialize sensor configurations
    std::vector<SensorConfig> sensors = initializeSensors();
    
    if (sensors.empty()) {
        std::cerr << "Error: No valid sensor configurations found!" << std::endl;
        return 1;
    }
    
    // Process each model in each sensor
    for (const auto& sensor : sensors) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "Processing Sensor Group " << sensor.sensorName << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
        // Process Model 1 (e.g., A1)
        processModel(sensor.model1Name, sensor.model1Mesh, stationaryNegative, 
                    sensor.transformations, sensor.sensorName, exportAnimation);
        
        // Process Model 2 (e.g., A2)  
        processModel(sensor.model2Name, sensor.model2Mesh, stationaryNegative,
                    sensor.transformations, sensor.sensorName, exportAnimation);
    }
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "=== All sensors processed ===" << std::endl;
    std::cout << "Output files generated:" << std::endl;
    for (const auto& sensor : sensors) {
        std::cout << "  - capacitance_" << sensor.model1Name << ".csv" << std::endl;
        std::cout << "  - capacitance_" << sensor.model2Name << ".csv" << std::endl;
    }
    
    if (exportAnimation) {
        std::cout << "\nAnimation files generated:" << std::endl;
        std::cout << "  - ray_data_[sensor]_step_[XXX].txt files" << std::endl;
        std::cout << "\nTo view animation:" << std::endl;
        std::cout << "  python multi_sensor_animated_viewer.py animate" << std::endl;
        std::cout << "To view single step:" << std::endl;
        std::cout << "  python multi_sensor_animated_viewer.py [step_number]" << std::endl;
    }
    
    std::cout << std::string(60, '=') << std::endl;
    
    return 0;
}