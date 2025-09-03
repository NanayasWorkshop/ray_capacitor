#include <embree4/rtcore.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <string>

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

struct RayData {
    float origin[3];
    float direction[3];
    float hit_point[3];
    bool hit;
    float distance;
};

// Physical constants for capacitor calculation
const double EPSILON_0 = 8.854e-12; // F/m (permittivity of free space)

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

// Export ray data for visualization
void exportRayData(const std::vector<RayData>& rays, const Mesh& posMesh, const Mesh& negMesh) {
    std::ofstream file("ray_data.txt");
    
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

// Main calculation function
void calculateCapacitor(const Mesh& positiveMesh, const Mesh& negativeMesh, bool visualize = false) {
    RTCDevice device = rtcNewDevice(nullptr);
    RTCScene negativeScene = createScene(device, negativeMesh);
    
    const float maxDistance = 10.0f; // 10mm max search distance
    
    double totalSum = 0.0;
    double minDistance = 1e9;
    double maxDistance_found = 0;
    double totalArea = 0;
    int hits = 0;
    int misses = 0;
    
    std::vector<RayData> rayDataCollection;
    
    std::cout << "\n=== Ray Tracing Analysis ===" << std::endl;
    
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
        
        // Collect ray data for visualization
        if (visualize) {
            RayData rd;
            rd.origin[0] = center[0];
            rd.origin[1] = center[1];
            rd.origin[2] = center[2];
            rd.direction[0] = normal[0];
            rd.direction[1] = normal[1];
            rd.direction[2] = normal[2];
            
            if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                rd.hit = true;
                rd.distance = rayhit.ray.tfar;
                rd.hit_point[0] = center[0] + normal[0] * rayhit.ray.tfar;
                rd.hit_point[1] = center[1] + normal[1] * rayhit.ray.tfar;
                rd.hit_point[2] = center[2] + normal[2] * rayhit.ray.tfar;
            } else {
                rd.hit = false;
                rd.distance = maxDistance;
                rd.hit_point[0] = center[0] + normal[0] * 2.0f;
                rd.hit_point[1] = center[1] + normal[1] * 2.0f;
                rd.hit_point[2] = center[2] + normal[2] * 2.0f;
            }
            
            rayDataCollection.push_back(rd);
        }
        
        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
            hits++;
            minDistance = std::min(minDistance, (double)distance_mm);
            maxDistance_found = std::max(maxDistance_found, (double)distance_mm);
            totalArea += area_mm2;
        } else {
            misses++;
        }
        
        // Detailed output for first few triangles
        if (i < 3) {
            std::cout << "Triangle " << i << ":" << std::endl;
            std::cout << "  Center: (" << center[0] << ", " << center[1] << ", " << center[2] << ") mm" << std::endl;
            std::cout << "  Normal: (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << std::endl;
            std::cout << "  Area: " << area_mm2 << " mm^2" << std::endl;
            std::cout << "  Distance: " << distance_mm << " mm" << std::endl;
            std::cout << "  Hit: " << (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID ? "Yes" : "No") << std::endl;
        }
        
        // Convert area from mm^2 to m^2 for the sum
        double area_m2 = area_mm2 * 1e-6;
        totalSum += (area_m2 * distance_mm);
    }
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total triangles processed: " << positiveMesh.triangles.size() << std::endl;
    std::cout << "Hits: " << hits << ", Misses: " << misses << std::endl;
    std::cout << "Min distance: " << minDistance << " mm" << std::endl;
    std::cout << "Max distance: " << maxDistance_found << " mm" << std::endl;
    std::cout << "Total area with hits: " << totalArea << " mm^2" << std::endl;
    std::cout << "Sum (area x distance): " << totalSum << " m^2*mm" << std::endl;
    
    // Calculate capacitance for parallel plates
    if (hits > 0 && minDistance < 10.0) {
        double area_m2 = totalArea * 1e-6;  // Convert mm^2 to m^2
        double distance_m = minDistance * 1e-3;  // Convert mm to m
        double capacitance = EPSILON_0 * area_m2 / distance_m;
        
        std::cout << "\n=== Capacitor Calculation ===" << std::endl;
        std::cout << "Effective plate area: " << totalArea << " mm^2 (" << area_m2 << " m^2)" << std::endl;
        std::cout << "Gap distance: " << minDistance << " mm (" << distance_m << " m)" << std::endl;
        std::cout << "Capacitance: " << capacitance * 1e12 << " pF" << std::endl;
        std::cout << "Expected for 64 mm^2 at 0.2mm: 2.833 pF" << std::endl;
    }
    
    if (visualize) {
        exportRayData(rayDataCollection, positiveMesh, negativeMesh);
        std::cout << "\nVisualization data exported to ray_data.txt" << std::endl;
        std::cout << "To visualize, run:" << std::endl;
        std::cout << "  cd capacitor_viz" << std::endl;
        std::cout << "  .\\venv\\Scripts\\activate" << std::endl;
        std::cout << "  python capacitor_viewer.py" << std::endl;
    }
    
    rtcReleaseScene(negativeScene);
    rtcReleaseDevice(device);
}

// Simple OBJ loader - NOW WITH METER TO MILLIMETER CONVERSION
Mesh loadOBJ(const std::string& filename) {
    Mesh mesh;
    
    FILE* file = fopen(filename.c_str(), "r");
    if (!file) {
        std::cerr << "Failed to open: " << filename << std::endl;
        return mesh;
    }
    
    char line[256];
    std::vector<Vec3> verts;
    
    std::cout << "Loading and converting from meters to millimeters..." << std::endl;
    
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

int main(int argc, char** argv) {
    bool visualize = (argc > 1 && std::string(argv[1]) == "-v");
    
    std::cout << "=== Parallel Plate Capacitor Ray Tracer ===" << std::endl;
    std::cout << "Converting OBJ files from meters to millimeters" << std::endl;
    std::cout << "Expected for 8x8mm plates with 0.2mm gap: ~2.833 pF" << std::endl;
    
    if (visualize) {
        std::cout << "Visualization mode enabled" << std::endl;
    }
    std::cout << std::endl;
    
    // Load models
    std::cout << "Loading positive model..." << std::endl;
    Mesh positiveMesh = loadOBJ("positive_model.obj");
    
    std::cout << "Loading negative model..." << std::endl;
    Mesh negativeMesh = loadOBJ("negative_model.obj");
    
    if (positiveMesh.triangles.empty() || negativeMesh.triangles.empty()) {
        std::cerr << "Error: Could not load models!" << std::endl;
        return 1;
    }
    
    calculateCapacitor(positiveMesh, negativeMesh, visualize);
    
    return 0;
}