#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

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
    
    std::cout << "Loading " << filename << "..." << std::endl;
    
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
                const bool KEYSHOT_UNITS_MM = true;  // Set to false if models are in meters
                
                if (KEYSHOT_UNITS_MM) {
                    // KeyShot format: already in mm, use as-is
                    verts.push_back(Vec3(x, y, z));
                } else {
                    // SolidEdge format: in meters, convert to mm
                    verts.push_back(Vec3(x * 1000.0f, y * 1000.0f, z * 1000.0f));
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
            }
        }
    }
    
    file.close();
    std::cout << "Loaded " << verts.size() << " vertices and " << mesh.triangles.size() 
              << " triangles from " << filename << std::endl;
    return mesh;
}

void printModelDimensions(const Mesh& mesh, const std::string& modelName) {
    if (mesh.triangles.empty()) {
        std::cout << "ERROR: " << modelName << " has no triangles!" << std::endl;
        return;
    }
    
    float minX = 1e9, maxX = -1e9;
    float minY = 1e9, maxY = -1e9; 
    float minZ = 1e9, maxZ = -1e9;
    
    // Check all vertices
    for (const auto& tri : mesh.triangles) {
        // Check all 3 vertices of each triangle
        float verts[9] = {tri.v0[0], tri.v0[1], tri.v0[2],
                         tri.v1[0], tri.v1[1], tri.v1[2], 
                         tri.v2[0], tri.v2[1], tri.v2[2]};
        
        for (int i = 0; i < 9; i += 3) {
            minX = std::min(minX, verts[i]);     maxX = std::max(maxX, verts[i]);
            minY = std::min(minY, verts[i+1]);   maxY = std::max(maxY, verts[i+1]);
            minZ = std::min(minZ, verts[i+2]);   maxZ = std::max(maxZ, verts[i+2]);
        }
    }
    
    float sizeX = maxX - minX;
    float sizeY = maxY - minY;
    float sizeZ = maxZ - minZ;
    float centerX = (minX + maxX) / 2.0f;
    float centerY = (minY + maxY) / 2.0f;
    float centerZ = (minZ + maxZ) / 2.0f;
    
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << "MODEL: " << modelName << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Dimensions: " << sizeX << " x " << sizeY << " x " << sizeZ << " mm" << std::endl;
    std::cout << "Center: (" << centerX << ", " << centerY << ", " << centerZ << ") mm" << std::endl;
    std::cout << "Bounds: X[" << minX << " to " << maxX << "]" << std::endl;
    std::cout << "        Y[" << minY << " to " << maxY << "]" << std::endl;
    std::cout << "        Z[" << minZ << " to " << maxZ << "]" << std::endl;
    std::cout << "Triangles: " << mesh.triangles.size() << std::endl;
    
    // Expected size check
    bool sizeOK = (sizeX >= 15.0f && sizeX <= 25.0f) && 
                  (sizeY >= 15.0f && sizeY <= 25.0f) && 
                  (sizeZ >= 5.0f && sizeZ <= 12.0f);
    
    std::cout << "\nExpected: ~20 x 20 x 8 mm" << std::endl;
    std::cout << "Status: " << (sizeOK ? "✓ DIMENSIONS OK" : "✗ UNEXPECTED SIZE") << std::endl;
    
    if (!sizeOK) {
        std::cout << "\nPOSSIBLE ISSUES:" << std::endl;
        if (sizeX < 1.0f || sizeY < 1.0f || sizeZ < 1.0f) {
            std::cout << "- Model too small: Units might be in meters, need conversion" << std::endl;
        }
        if (sizeX > 100.0f || sizeY > 100.0f || sizeZ > 100.0f) {
            std::cout << "- Model too big: Units already converted or wrong scale" << std::endl;
        }
    }
}

int main() {
    std::cout << "=== Model Dimension Checker ===" << std::endl;
    std::cout << "Verifying model sizes and unit conversions" << std::endl;
    
    // List of models to check
    std::vector<std::string> modelFiles = {
        "A1_model.obj", "A2_model.obj",
        "B1_model.obj", "B2_model.obj", 
        "C1_model.obj", "C2_model.obj",
        "stationary_negative.obj"
    };
    
    int modelsFound = 0;
    for (const auto& filename : modelFiles) {
        Mesh mesh = loadOBJ(filename);
        if (!mesh.triangles.empty()) {
            printModelDimensions(mesh, filename);
            modelsFound++;
        } else {
            std::cout << "\nWARNING: Could not load " << filename << std::endl;
        }
    }
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "SUMMARY: Checked " << modelsFound << "/" << modelFiles.size() << " model files" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    if (modelsFound == 0) {
        std::cout << "No models found! Make sure you're running from the bin directory." << std::endl;
    }
    
    return 0;
}