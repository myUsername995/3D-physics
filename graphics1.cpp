#include <glad/glad.h>
#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3/SDL_opengl.h>
#include <SDL3_ttf/SDL_ttf.h>

std::mutex printMutex;

// SDL_opengl defines these goofy things that match my variable names
#undef near
#undef far

using namespace std;

#define WINDOW_HEIGHT 800
#define WINDOW_WIDTH 800
#define PI 3.14159

struct mat4x4 {
    float m[4][4] = { 0 };
};

struct vec2d {
    double x = 0;
    double y = 0;
};

struct vec3d {
    double x = 0;
    double y = 0;
    double z = 0;
    float w = 1;

    vec3d& operator/=(float k) {
        x /= k;
        y /= k;
        z /= k;
        return *this;
    }
};

vec3d operator-(const vec3d& a) {
    return {-a.x, -a.y, -a.z};
}

vec3d operator-(const vec3d& a, const vec3d& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

vec3d operator+(const vec3d& a, const vec3d& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

vec3d operator*(const vec3d& a, const float& k) {
    return {a.x * k, a.y * k, a.z * k};
}

vec3d operator*(const vec3d& a, const vec3d& k) {
    return {a.x * k.x, a.y * k.y, a.z * k.z};
}

struct color {
    float r, g, b, a;
};

color operator+(const color& a, const color& b){
    float red = a.r + b.r;
    float green = a.g + b.g;
    float blue = a.b + b.b;

    return {red, green, blue, 255};
}

color operator-(const color& a, const color& b){
    float red = a.r - b.r;
    float green = a.g - b.g;
    float blue = a.b - b.b;

    return {red, green, blue, 255};
}

color operator*(const color& a, const float& k){
    float red = a.r * k;
    float green = a.g * k;
    float blue = a.b * k;

    return {red, green, blue, 255};
}

struct material {

    std::string name;
    vec3d ambient = {1.0, 1.0, 1.0};
    vec3d diffuse;
    vec3d specular;
    vec3d emitting;
    float shininess;

};

struct triangle {
    int p[3];                // Indexes of the vertices
    vec3d normals[3];
    color vertCols[3];
    material mat;
    bool cullBackface = true;      // Some triangles shouldn't be culled!!!
};

// Used for temporary triangles
struct triangleP {
    vec3d p[3];              // Actual points
    vec3d normals[3];        // Store vertex normal information
    color vertCols[3];
    material mat;
};

struct quad {
    vec3d p[4];
};

struct Link {
    int startIndex = 0;
    int endIndex = 0;
    int index[2];
    double distance;
    bool visual;
};

struct Spring {
    int startIndex = 0;
    int endIndex = 0;
    int index[2];             // The two objects
    float k;                 // The stiffness
    float restingDist;       // The distance where spring is at a resting position
};

// An object is just a collection of triangles with some properties and physics are applied to it
struct object {

    int startIndex = 0;
    int endIndex = 0;
    int startTri = 0; 
    int endTri = 0;

    int velocityVec = -1;               // Index into the obj array, used to display the velocity vector
    vec3d midP, prevMid;
    vec3d acceleration = {0, 0, 0};
    vec3d color;
    float radius;

    bool drawSphere = true;
    bool anchored = false;
    bool selfCollision = true;          // Collide within it's own object?
    bool doCollide = true;              // Collide with any objects?

    int objIndex;               // Reserve 0 for normal balls

};

struct physicsFace {

    int index[3];       // Indexes into the objects array
    int normal;         // Index into the objects array to the point an unit normal away
    color color;
    bool cullBackface = true;
    int startIndex, endIndex;

};

struct mesh {
    vector<vec2d> textures;
    vector<vec3d> vertices;
    vector<triangle> tris;
};

struct rectangle {
    vector<quad> faces;     // Store geometry data
    unordered_set<int> trisIndex;  // Stores the triangles which are inside the box

    vec3d minVert, maxVert;        // Stores min and max vertex
};

struct Plane {
    vec3d normal;
    float d;
};

struct Edge{
    int y_min, y_max;
    double x, inv_slope;
};

struct DynamicMesh {
    GLuint VAO;
    GLuint VBO;
    GLsizei maxVertexCount;
};

struct Vertex {
    vec3d pos;
    vec3d normal;
    color col; // if you have vertex colors
};

// Create a global meshObj
mesh meshObj;
unordered_map<std::string, material> materials;

// Shader functions
std::string LoadShaderSource(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open shader file: " << filepath << std::endl;
        return "";
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string source = buffer.str();

    if (source.empty()) {
        std::cerr << "Shader file is empty: " << filepath << std::endl;
    }

    return source;
}

DynamicMesh CreateDynamicMesh(size_t maxVertices) {
    DynamicMesh mesh;

    glGenVertexArrays(1, &mesh.VAO);
    glBindVertexArray(mesh.VAO);

    glGenBuffers(1, &mesh.VBO);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.VBO);
    glBufferData(GL_ARRAY_BUFFER, maxVertices * 6 * sizeof(float), nullptr, GL_DYNAMIC_DRAW); // allocate empty buffer

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    mesh.maxVertexCount = static_cast<GLsizei>(maxVertices);
    return mesh;
}

void UpdateDynamicMesh(DynamicMesh& mesh, const std::vector<float>& vertexData) {
    //glBindVertexArray(mesh.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, vertexData.size() * sizeof(float), vertexData.data());
}

void DrawDynamicMeshTriangle(const DynamicMesh& mesh, GLuint shaderProgram, GLsizei vertexCount){
    glUseProgram(shaderProgram);

    glBindVertexArray(mesh.VAO);
    glDrawArrays(GL_TRIANGLES, 0, vertexCount);
}

void DrawDynamicMeshLine(const DynamicMesh& mesh, GLuint shaderProgram, GLsizei vertexCount) {
    glUseProgram(shaderProgram);
    glBindVertexArray(mesh.VAO);
    glDrawArrays(GL_LINES, 0, vertexCount);
}

vec3d meshFindIndex(int index){
    return meshObj.vertices.at(index);
}

// Returns the index it was added to
int meshAddToVertex(vec3d vertex){
    meshObj.vertices.push_back(vertex);

    return meshObj.vertices.size() - 1;
}

void meshAddToTriangle(triangle tri){
    meshObj.tris.push_back(tri);
}

int meshFindVertex(vec3d vertex){
    for (int i = 0; i < meshObj.vertices.size(); i++){
        vec3d vert = meshObj.vertices.at(i);
        // Check if its the same vertex
        if (vertex.x == vert.x && vertex.y == vert.y && vertex.z == vert.z){
            return i;
        }
    }

    return 0;
}

triangleP meshConvert(const triangle& t) {
    triangleP tp;
    for (int i = 0; i < 3; ++i) {
        tp.p[i] = meshObj.vertices[t.p[i]];       // Convert index to actual point
        tp.vertCols[i] = t.vertCols[i];
    }
    return tp;
}

// Load an MTL file into a materials array
void loadMTL(const string& filepath){
    const string folderFilename = "Textures//" + filepath;
    ifstream file(folderFilename);

    if (!file){
        cout << "Coulnd't open MTL file\n";
        exit(1);
    }

    string line;

    material cur;
    while (getline(file, line)){
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "newmtl"){
            if (!cur.name.empty()){
                materials[cur.name] = cur;
            }

            iss >> cur.name;
        }
        else if (prefix == "Ns"){
            float nsValue;
            iss >> nsValue;

            cur.shininess = nsValue;
        }
        else if (prefix == "Ka"){
            float r, g, b;
            iss >> r >> g >> b;
            cur.ambient = {r, g, b};
        }
        else if (prefix == "Ks"){
            float r, g, b;
            iss >> r >> g >> b;
            cur.specular = {r, g, b};
        }
        else if (prefix == "Kd"){
            float r, g, b;
            iss >> r >> g >> b;
            cur.diffuse = {r, g, b};
        }
        else if (prefix == "Ke"){
            float r, g, b;
            iss >> r >> g >> b;
            cur.emitting = {r, g, b};
        }
    }

    if (!cur.name.empty()) {
        materials[cur.name] = cur;
    }
}

void parseObj(const std::string& fileName, vec3d& rotation){
    const string folderFilename = "objFiles//" + fileName;
    ifstream inFile(folderFilename, ios::in);

    // Make it so you can parse multiple obj files and still keep them right
    size_t startIndex = meshObj.vertices.size();

    if (!inFile) {
        cerr << "Error: Could not open file " << folderFilename << endl;
        exit(1);
    }

    string line;
    material cur;
    line.reserve(128);  // Reserve space to reduce reallocations

    while (getline(inFile, line)){
        if (line.empty()) continue;
        istringstream iss(line);
        string prefix;

        iss >> prefix;

        if (prefix == "v"){
            float x, y, z;

            iss >> x >> y >> z;
            meshObj.vertices.push_back({x, -y, z});
        }
        // Defines how much you should rotate for the map to be horizontal
        else if (prefix == "rot"){
            iss >> rotation.x >> rotation.y >> rotation.z;
        }
        else if (prefix == "mtllib"){
            string fileName;
            iss >> fileName;

            loadMTL(fileName);
        }
        else if (prefix == "usemtl"){
            // Get the name of the material
            string name;
            iss >> name;

            // Look up in the hashtable
            cur = materials.at(name);
        }
        else if (prefix == "f"){
            array<string, 3> faces;
            iss >> faces[0] >> faces[1] >> faces[2];

            triangle newTri;
            for (int i = 0; i < 3; i++){
                string v, t;
                int firstSlash = faces[i].find('/');
                int secondSlash; int secondIndex;

                // Check if there is a second slash "//"
                if (faces[i].at(firstSlash + 1) == '/'){
                    secondSlash = faces[i].find('/', firstSlash + 2);
                    secondIndex = firstSlash + 2;
                }
                else {
                    secondSlash = faces[i].find('/', firstSlash + 1);
                    secondIndex = firstSlash + 1;
                }

                // Case where there is no slash in the file
                if (firstSlash == 0){
                    // Just a vertex index and no texture index
                    v = faces[i];

                    newTri.p[i] = startIndex + stoi(v) - 1;
                }
                else {
                    v = faces[i].substr(0, firstSlash);
                    t = faces[i].substr(secondIndex, secondSlash - secondIndex);

                    newTri.p[i] = startIndex + stoi(v) - 1;
                }
            }

            newTri.mat = cur;
            meshAddToTriangle(newTri);
        }
    }
}

void fillTriangle(SDL_Renderer* renderer, vec3d v0, vec3d v1, vec3d v2) {
    // Sort by y
    if (v1.y < v0.y) std::swap(v0, v1);
    if (v2.y < v1.y) std::swap(v1, v2);
    if (v1.y < v0.y) std::swap(v0, v1);

    // Inverse slopes
    auto invSlope = [&](const vec3d& A, const vec3d& B) {
        return (B.x - A.x) / (B.y - A.y);
    };
    double inv01 = (v1.y != v0.y) ? invSlope(v0, v1) : 0.0;
    double inv02 = (v2.y != v0.y) ? invSlope(v0, v2) : 0.0;
    double inv12 = (v2.y != v1.y) ? invSlope(v1, v2) : 0.0;

    int yStart = (int)std::ceil(v0.y);
    int yMid   = (int)std::ceil(v1.y);
    int yEnd   = (int)std::ceil(v2.y);

    double xLeft  = v0.x;
    double xRight = v0.x;

    // Top half
    for (int y = yStart; y < yMid; ++y) {
        SDL_RenderLine(renderer, xLeft, y, xRight, y);
        xLeft  += inv01;
        xRight += inv02;
    }

    // Bottom half
    xLeft = v1.x;  // restart at the middle vertex
    for (int y = yMid; y < yEnd; ++y) {
        SDL_RenderLine(renderer, xLeft, y, xRight, y);
        xLeft  += inv12;
        xRight += inv02;
    }
}

vec3d rectangleFindMax(const rectangle& boundingB){
    // Find the minimum x y and maximum x y cordinates of the meshcube
    float maxX = -INFINITY; float maxY = -INFINITY; float maxZ = -INFINITY;

    // Loop through all faces
    for (quad face : boundingB.faces){
        // Loop through all vertices
        for (int i = 0; i < 4; i++){
            vec3d vert = face.p[i];

            // Update min / max variables
            maxX = fmax(maxX, vert.x);
            maxY = fmax(maxY, vert.y);
            maxZ = fmax(maxZ, vert.z);
        }
    }

    return {maxX, maxY, maxZ};
}

vec3d rectangleFindMin(const rectangle& boundingB) {
    // Start with positive infinity so any vertex will be smaller
    float minX = INFINITY;
    float minY = INFINITY;
    float minZ = INFINITY;

    // Loop through all faces
    for (const quad& face : boundingB.faces) {
        // Loop through all vertices
        for (int i = 0; i < 4; i++) {
            const vec3d& vert = face.p[i];

            // Update min values
            minX = fmin(minX, vert.x);
            minY = fmin(minY, vert.y);
            minZ = fmin(minZ, vert.z);
        }
    }

    return {minX, minY, minZ};
}

bool isTriangleInsideFrustum(const triangle& tri, const vector<Plane>& frustumPlanes) {
    for (const Plane& plane : frustumPlanes) {
        int outsideCount = 0;

        for (const vec3d& v : {meshFindIndex(tri.p[0]), meshFindIndex(tri.p[1]), meshFindIndex(tri.p[2])}) {
            float distance = plane.normal.x * v.x +
                             plane.normal.y * v.y +
                             plane.normal.z * v.z + plane.d;
                             
            if (distance < 0) {
                outsideCount++;
            }
        }

        // All vertices are outside this plane → triangle is outside
        if (outsideCount == 3) {
            return false;
        }
    }

    return true; // Triangle is at least partially inside
}

rectangle shiftRectangle(const rectangle& rect, float distanceX, float distanceY, float distanceZ) {
    rectangle shifted = rect;  // Make a copy

    for (quad& face : shifted.faces) {
        for (int i = 0; i < 4; i++) {
            face.p[i].x += distanceX;
            face.p[i].y += distanceY;
            face.p[i].z += distanceZ;
        }
    }

    return shifted;
}

// Declare this for the createEnclosingRectangle function
vec3d matrixMultiplyVector(const mat4x4& m, const vec3d& i);

void createEnclosingRectangle(const mesh& meshCube, rectangle& boundingB){
    // Find the minimum x y and maximum x y cordinates of the meshcube
    float minX = INFINITY; float maxX = -INFINITY;
    float minY = INFINITY; float maxY = -INFINITY;
    float minZ = INFINITY; float maxZ = -INFINITY;

    // Loop through all triangles
    for (triangle tri : meshCube.tris){
        // Loop through all vertices
        for (int i = 0; i < 3; i++){
            vec3d vert = meshFindIndex(tri.p[i]);

            // Update min / max variables
            minX = fmin(minX, vert.x);
            maxX = fmax(maxX, vert.x);
            minY = fmin(minY, vert.y);
            maxY = fmax(maxY, vert.y);
            minZ = fmin(minZ, vert.z);
            maxZ = fmax(maxZ, vert.z);
        }
    }

    // Create the 8 corner vertices of the bounding box
    vec3d v000 = {minX, minY, minZ}; // left-bottom-back
    vec3d v001 = {minX, minY, maxZ}; // left-bottom-front
    vec3d v010 = {minX, maxY, minZ}; // left-top-back
    vec3d v011 = {minX, maxY, maxZ}; // left-top-front
    vec3d v100 = {maxX, minY, minZ}; // right-bottom-back
    vec3d v101 = {maxX, minY, maxZ}; // right-bottom-front
    vec3d v110 = {maxX, maxY, minZ}; // right-top-back
    vec3d v111 = {maxX, maxY, maxZ}; // right-top-front

    // Create the 6 faces (quads)
    quad front  = {v011, v111, v101, v001}; // +Z
    quad back   = {v100, v110, v010, v000}; // -Z
    quad left   = {v010, v011, v001, v000}; // -X
    quad right  = {v101, v111, v110, v100}; // +X
    quad bottom = {v001, v101, v100, v000}; // -Y
    quad top    = {v110, v111, v011, v010}; // +Y

    // Store in boundingB
    boundingB.faces.clear();
    boundingB.faces.push_back(front);
    boundingB.faces.push_back(back);
    boundingB.faces.push_back(left);
    boundingB.faces.push_back(right);
    boundingB.faces.push_back(bottom);
    boundingB.faces.push_back(top);

    boundingB.minVert = {minX, minY, minZ};
    boundingB.maxVert = {maxX, maxY, maxZ};
}

// subdivide only along X (width) and Z (depth)
vector<rectangle> makeGrid(const rectangle& original, int detail){
    double totalW = original.maxVert.x - original.minVert.x;
    double totalD = original.maxVert.z - original.minVert.z;

    double cellW = totalW / detail;
    double cellD = totalD / detail;

    vector<rectangle> out;
    out.reserve(detail * detail);

    for (int ix = 0; ix < detail; ++ix)
    {
        for (int iz = 0; iz < detail; ++iz)
        {
            rectangle cell;

            // lower‐corner
            cell.minVert.x = original.minVert.x + ix * cellW;
            cell.minVert.y = original.minVert.y;                      // keep full Y
            cell.minVert.z = original.minVert.z + iz * cellD;

            // upper‐corner
            cell.maxVert.x = original.minVert.x + (ix + 1) * cellW;
            cell.maxVert.y = original.maxVert.y;                      // keep full Y
            cell.maxVert.z = original.minVert.z + (iz + 1) * cellD;

            // Ensure we have exactly 6 quads
            cell.faces.resize(6);

            // Alias for readability
            const auto& mn = cell.minVert;
            const auto& mx = cell.maxVert;

            // Build the eight corners of the AABB
            vec3d v000{ mn.x, mn.y, mn.z };  // left-bottom-back
            vec3d v001{ mn.x, mn.y, mx.z };  // left-bottom-front
            vec3d v010{ mn.x, mx.y, mn.z };  // left-top-back
            vec3d v011{ mn.x, mx.y, mx.z };  // left-top-front
            vec3d v100{ mx.x, mn.y, mn.z };  // right-bottom-back
            vec3d v101{ mx.x, mn.y, mx.z };  // right-bottom-front
            vec3d v110{ mx.x, mx.y, mn.z };  // right-top-back
            vec3d v111{ mx.x, mx.y, mx.z };  // right-top-front

            // 0) Front  (+Z)
            cell.faces[0].p[0] = v011;
            cell.faces[0].p[1] = v111;
            cell.faces[0].p[2] = v101;
            cell.faces[0].p[3] = v001;

            // 1) Back   (–Z)
            cell.faces[1].p[0] = v100;
            cell.faces[1].p[1] = v110;
            cell.faces[1].p[2] = v010;
            cell.faces[1].p[3] = v000;

            // 2) Left   (–X)
            cell.faces[2].p[0] = v010;
            cell.faces[2].p[1] = v011;
            cell.faces[2].p[2] = v001;
            cell.faces[2].p[3] = v000;

            // 3) Right  (+X)
            cell.faces[3].p[0] = v101;
            cell.faces[3].p[1] = v111;
            cell.faces[3].p[2] = v110;
            cell.faces[3].p[3] = v100;

            // 4) Bottom (–Y)
            cell.faces[4].p[0] = v000;
            cell.faces[4].p[1] = v100;
            cell.faces[4].p[2] = v101;
            cell.faces[4].p[3] = v001;

            // 5) Top    (+Y)
            cell.faces[5].p[0] = v010;
            cell.faces[5].p[1] = v110;
            cell.faces[5].p[2] = v111;
            cell.faces[5].p[3] = v011;

            out.push_back(cell);
        }
    }

    return out;
}

void addBoundingBoxToMesh(mesh& meshCube, const vector<rectangle> boundingBs){
    // Go through all the bounding boxes
    for (const rectangle& rect : boundingBs){
        // Go through all the faces of the boundingBs
        for (const quad& face : rect.faces){
            // Split the faces into 2 triangles
            vec3d v1 = face.p[0];
            vec3d v2 = face.p[1];
            vec3d v3 = face.p[2];
            vec3d v4 = face.p[3];

            // Add these new vertices to the vertex array
            int i1 = meshAddToVertex(v1);
            int i2 = meshAddToVertex(v2);
            int i3 = meshAddToVertex(v3);
            int i4 = meshAddToVertex(v4);
            
            // Triangles
            triangle tri1 = {i1, i2, i3};
            triangle tri2 = {i1, i3, i4};

            meshAddToTriangle(tri1);
            meshAddToTriangle(tri2);
        }
    }
}

vec3d vectorAdd(const vec3d& v1, const vec3d& v2){
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

vec3d vectorSub(const vec3d& v1, const vec3d& v2){
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

vec3d vectorMul(const vec3d& v1, const float k){
    return {v1.x * k, v1.y * k, v1.z * k};
}

vec3d vectorDiv(const vec3d& v1, const float k){
    return {v1.x / k, v1.y / k, v1.z / k};
}


float vectorDot(const vec3d& v1, const vec3d& v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float vectorLength(const vec3d& v){
    return sqrtf(vectorDot(v, v));
}

vec3d vectorNormalise(const vec3d& v){
    float l = vectorLength(v);
    return {v.x / l, v.y /l, v.z / l};
}

vec3d vectorCross(const vec3d& v1, const vec3d& v2) {
    vec3d v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}


vec3d matrixMultiplyVector(const mat4x4& m, const vec3d& i) {
    vec3d v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

    return v;
}

inline void matrixMultiplyVectorFast(vec3d& __restrict out, const mat4x4& __restrict m, const vec3d& __restrict i) {
    vec3d v;
    out.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    out.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    out.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    out.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
}

mat4x4 matrixMakeIdentity(){
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;

    return matrix;
}

mat4x4 matrixMakeRotationX(const float angleRad){
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = cosf(angleRad);
    matrix.m[1][2] = sinf(angleRad);
    matrix.m[2][1] = -sinf(angleRad);
    matrix.m[2][2] = cosf(angleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}


mat4x4 matrixMakeRotationY(const float angleRad){
    mat4x4 matrix;
    matrix.m[0][0] = cosf(angleRad);
    matrix.m[0][2] = sinf(angleRad);
    matrix.m[2][0] = -sinf(angleRad);
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = cosf(angleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixMakeRotationZ(const float angleRad){
    mat4x4 matrix;
    matrix.m[0][0] = cosf(angleRad);
    matrix.m[0][1] = sinf(angleRad);
    matrix.m[1][0] = -sinf(angleRad);
    matrix.m[1][1] = cosf(angleRad);
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}


mat4x4 matrixMakeTranslation(float x, float y, float z){
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    matrix.m[3][0] = x;
    matrix.m[3][1] = y;
    matrix.m[3][2] = z;

    return matrix;
}

mat4x4 matrixMakeProjection(float fovDegrees, float aspectRatio, float near, float far){
    float fovRad = 1.0f / tanf(fovDegrees * 0.5f / 180.0f * 3.14159f);

    mat4x4 matProj;
    matProj.m[0][0] = aspectRatio * fovRad;
    matProj.m[1][1] = fovRad;
    matProj.m[2][2] = far / (far - near);
    matProj.m[3][2] = (-far * near) / (far - near);
    matProj.m[2][3] = 1.0f;
    matProj.m[3][3] = 0.0f;

    return matProj;
}

mat4x4 matrixTranspose(const mat4x4& m) {
    mat4x4 result;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            result.m[i][j] = m.m[j][i];
    return result;
}

mat4x4 matrixPointAt(const vec3d& pos, const vec3d& target, const vec3d& up){
    vec3d newForward = vectorSub(target, pos);
    newForward = vectorNormalise(newForward);

    vec3d a = vectorMul(newForward, vectorDot(up, newForward));
    vec3d newUp = vectorSub(up, a);
    newUp = vectorNormalise(newUp);

    vec3d newRight = vectorCross(newUp, newForward);

    mat4x4 matrix;
    matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixQuickInverse(const mat4x4& m){
    mat4x4 matrix;
    matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixMultiplyMatrix(const mat4x4& m1, const mat4x4& m2) {
    mat4x4 matrix;

    for (int c = 0; c < 4; c++) {
        for (int r = 0; r < 4; r++) {
            matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
        }
    }

    return matrix;
}

std::array<Plane, 6> extractFrustumPlanes(const mat4x4& m, bool normalize = true) {
    std::array<Plane, 6> planes;

    // Left
    planes[0].normal.x = m.m[0][3] + m.m[0][0];
    planes[0].normal.y = m.m[1][3] + m.m[1][0];
    planes[0].normal.z = m.m[2][3] + m.m[2][0];
    planes[0].d        = m.m[3][3] + m.m[3][0];

    // Right
    planes[1].normal.x = m.m[0][3] - m.m[0][0];
    planes[1].normal.y = m.m[1][3] - m.m[1][0];
    planes[1].normal.z = m.m[2][3] - m.m[2][0];
    planes[1].d        = m.m[3][3] - m.m[3][0];

    // Bottom
    planes[2].normal.x = m.m[0][3] + m.m[0][1];
    planes[2].normal.y = m.m[1][3] + m.m[1][1];
    planes[2].normal.z = m.m[2][3] + m.m[2][1];
    planes[2].d        = m.m[3][3] + m.m[3][1];

    // Top
    planes[3].normal.x = m.m[0][3] - m.m[0][1];
    planes[3].normal.y = m.m[1][3] - m.m[1][1];
    planes[3].normal.z = m.m[2][3] - m.m[2][1];
    planes[3].d        = m.m[3][3] - m.m[3][1];

    // Near
    planes[4].normal.x = m.m[0][3] + m.m[0][2];
    planes[4].normal.y = m.m[1][3] + m.m[1][2];
    planes[4].normal.z = m.m[2][3] + m.m[2][2];
    planes[4].d        = m.m[3][3] + m.m[3][2];

    // Far
    planes[5].normal.x = m.m[0][3] - m.m[0][2];
    planes[5].normal.y = m.m[1][3] - m.m[1][2];
    planes[5].normal.z = m.m[2][3] - m.m[2][2];
    planes[5].d        = m.m[3][3] - m.m[3][2];

    // Normalize the planes
    if (normalize) {
        for (auto& plane : planes) {
            float length = sqrtf(plane.normal.x * plane.normal.x +
                                 plane.normal.y * plane.normal.y +
                                 plane.normal.z * plane.normal.z);
            plane.normal = vectorDiv(plane.normal, length);
            plane.d /= length;
        }
    }

    return planes;
}

bool frustumBoundingBoxIntersection(const std::array<Plane, 6>& frustumPlanes, const vec3d& boxMin, const vec3d& boxMax) {
    for (const auto& plane : frustumPlanes) {
        // Compute the positive vertex (furthest in direction of plane normal)
        vec3d positiveVertex = {
            plane.normal.x >= 0 ? boxMax.x : boxMin.x,
            plane.normal.y >= 0 ? boxMax.y : boxMin.y,
            plane.normal.z >= 0 ? boxMax.z : boxMin.z
        };

        // Plane equation: Ax + By + Cz + D
        float distance = plane.normal.x * positiveVertex.x +
                         plane.normal.y * positiveVertex.y +
                         plane.normal.z * positiveVertex.z + plane.d;

        if (distance < 0) {
            // Box is completely outside this plane
            return false;
        }
    }
    return true; // Box intersects or is inside the frustum
}

void boundForce(object& object){
    vec3d vel = object.midP - object.prevMid;

    // Floor
    if (object.midP.y + object.radius > 0) {
        object.midP.y = -object.radius;
        vel.y *= -0.9f; 
        object.prevMid = object.midP - vel;
    }

    // Sides
    if (object.midP.x - object.radius < 0){
        object.midP.x = object.radius;
        vel.x *= -0.9f;
        object.prevMid = object.midP - vel;
    }
    if (object.midP.x + object.radius > 380){
        object.midP.x = 380 - object.radius;
        vel.x *= -0.9f;
        object.prevMid = object.midP - vel;
    }
    if (object.midP.z - object.radius < 0){
        object.midP.z = object.radius;
        vel.z *= -0.9f;
        object.prevMid = object.midP - vel;
    }
    if (object.midP.z + object.radius > 380){
        object.midP.z = 380 - object.radius;
        vel.z *= -0.9f;
        object.prevMid = object.midP - vel;
    }
}

void springForce(vector<object>& objects, Spring spring){
    object& circle1 = objects.at(spring.index[0]);
    object& circle2 = objects.at(spring.index[1]);

    vec3d delta = circle2.midP - circle1.midP;
    float dist = vectorLength(delta);
    if (dist < 0.0001f) return;

    float diff = (dist - spring.restingDist) / dist;
    vec3d correction = delta * spring.k * 0.5f * diff;

    if (!circle1.anchored) circle1.midP = vectorAdd(circle1.midP, correction);
    if (!circle2.anchored) circle2.midP = vectorSub(circle2.midP, correction);
}

void linkForce(vector<object>& objects, Link link){
    object& circle1 = objects.at(link.index[0]);
    object& circle2 = objects.at(link.index[1]);

    vec3d linkAxis = circle2.midP - circle1.midP;
    double currentDist = vectorLength(linkAxis);

    if (currentDist < 1e-6) return;

    vec3d correction = vectorMul(linkAxis, 1 / currentDist) * (currentDist - link.distance);

    if (circle1.anchored){
        circle2.midP = vectorSub(circle2.midP, correction);
    }
    else if (circle2.anchored){
        circle1.midP = vectorAdd(circle1.midP, correction);
    }
    else {
        circle1.midP = vectorAdd(circle1.midP, correction * 0.5f);
        circle2.midP = vectorSub(circle2.midP, correction * 0.5f);
    }
}

void collisionForce(object& circle1, object& circle2) {
    if (circle1.midP.x == circle2.midP.x &&
        circle1.midP.y == circle2.midP.y &&
        circle1.midP.z == circle2.midP.z) return;

    vec3d delta = circle1.midP - circle2.midP;
    float distSq = vectorDot(delta, delta);
    float minDist = circle1.radius + circle2.radius;
    float minDistSq = minDist * minDist;

    if (distSq > minDistSq) return; // no collision

    float dist = sqrt(distSq);
    float overlap = minDist - dist;

    // velocities
    vec3d velocity1 = circle1.midP - circle1.prevMid;
    vec3d velocity2 = circle2.midP - circle2.prevMid;

    // normalized collision axis
    vec3d collisionAxis = vectorMul(delta, 1 / dist);;

    // normal velocity components
    float v1n = vectorDot(velocity1, collisionAxis);
    float v2n = vectorDot(velocity2, collisionAxis);

    // swap normal components
    vec3d newV1 = velocity1 + collisionAxis * (v2n - v1n);
    vec3d newV2 = velocity2 + collisionAxis * (v1n - v2n);

    newV1 = vectorMul(newV1, 0.9f);
    newV2 = vectorMul(newV2, 0.9f);

    // offsets for overlap resolution
    float offset1, offset2;
    if (circle1.anchored) {
        offset1 = 0; offset2 = overlap; newV1 = {0,0};
    } else if (circle2.anchored) {
        offset1 = overlap; offset2 = 0; newV2 = {0,0};
    } else {
        offset1 = overlap * 0.5f;
        offset2 = overlap * 0.5f;
    }

    // separate circles
    circle1.prevMid = circle1.midP + collisionAxis * offset1;
    circle2.prevMid = circle2.midP - collisionAxis * offset2;

    // apply new positions
    circle1.midP = circle1.prevMid + newV1;
    circle2.midP = circle2.prevMid + newV2;
}

// Return closest point on triangle (v0,v1,v2) to point p
vec3d closestPointOnTriangle(const vec3d& p, const vec3d& v0, const vec3d& v1, const vec3d& v2) {
    // Edges
    vec3d ab = v1 - v0;
    vec3d ac = v2 - v0;
    vec3d ap = p  - v0;

    // Barycentric coordinates
    float d1 = vectorDot(ab, ap);
    float d2 = vectorDot(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f) return v0; // barycentric (1,0,0)

    vec3d bp = p - v1;
    float d3 = vectorDot(ab, bp);
    float d4 = vectorDot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) return v1; // barycentric (0,1,0)

    float vc = d1*d4 - d3*d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        float v = d1 / (d1 - d3);
        return v0 + ab * v; // between v0 and v1
    }

    vec3d cp = p - v2;
    float d5 = vectorDot(ab, cp);
    float d6 = vectorDot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) return v2; // barycentric (0,0,1)

    float vb = d5*d2 - d1*d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        float w = d2 / (d2 - d6);
        return v0 + ac * w; // between v0 and v2
    }

    float va = d3*d6 - d5*d4;
    if (va <= 0.0f && (d4-d3) >= 0.0f && (d5-d6) >= 0.0f) {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return v1 + (v2 - v1) * w; // between v1 and v2
    }

    // Inside face region: use barycentrics
    float denom = 1.0f / (va + vb + vc);
    float u = va * denom;
    float v = vb * denom;
    float w = 1.0f - u - v;
    return v0*u + v1*v + v2*w;
}

void collisionForceFace(vector<object>& objects, object& circle1, vector<physicsFace>& faces) {
    double smallestDist = INFINITY;
    int smallestIndex = -1;

    // See which face is the closest to the ball
    for (int i = 0; i < faces.size(); i++){
        auto& face = faces[i];

        if (circle1.objIndex == objects[face.index[0]].objIndex) continue;

        object v0 = objects[face.index[0]];
        object v1 = objects[face.index[1]];
        object v2 = objects[face.index[2]];

        vec3d curClosest = closestPointOnTriangle(circle1.midP, v0.midP, v1.midP, v2.midP);
        vec3d prevClosest = closestPointOnTriangle(circle1.prevMid, v0.prevMid, v1.prevMid, v2.prevMid);

        // The direction to the closest point, where the ball will collide with the plane
        vec3d curVector = curClosest - circle1.midP;
        vec3d prevVector = prevClosest - circle1.prevMid;

        double curDist = vectorLength(curVector);
        double prevDist = vectorLength(prevVector);

        // See if the ball is close enough to the closest point to touch, if not then its not colliding
        if (curDist > circle1.radius && prevDist > circle1.radius) continue;

        // Only apply the closest point of collision
        if (curDist < smallestDist){
            smallestDist = curDist;
            smallestIndex = i;
        }
    }

    // No collisions
    if (smallestIndex == -1) return;

    object& v0 = objects[faces[smallestIndex].index[0]];
    object& v1 = objects[faces[smallestIndex].index[1]];
    object& v2 = objects[faces[smallestIndex].index[2]];

    // Calculate resolution vector, we do this by taking the planes and balls velocity,
    // and applying a random equation I got from chatGPT to calculate the resolution vector accuretaly

    vec3d planeVelocity = vectorDiv((v0.midP - v0.prevMid + v1.midP - v1.prevMid + v2.midP - v2.prevMid), 3.0f);
    vec3d planeNormal = vectorNormalise(vectorCross(v0.midP - v1.midP, v0.midP - v2.midP));

    // See if the ball actually collides with the plane
    // What to check: if the ball is on one side of the plane in the prevMid phase, and then on the other side
    // at the midP phase, then the ball has collided with the plane.

    vec3d curClosest = closestPointOnTriangle(circle1.midP, v0.midP, v1.midP, v2.midP);
    vec3d prevClosest = closestPointOnTriangle(circle1.prevMid, v0.prevMid, v1.prevMid, v2.prevMid);

    // The direction to the closest point, where the ball will collide with the plane
    vec3d curVector = curClosest - circle1.midP;
    vec3d prevVector = prevClosest - circle1.prevMid;

     // The closest point on the SPHERE to the plane, use THIS to detect collisions with the actual ball
    vec3d curEdgePoint = circle1.midP + vectorMul(vectorNormalise(curVector), circle1.radius);
    vec3d prevEdgePoint = circle1.prevMid + vectorMul(vectorNormalise(prevVector), circle1.radius);

    vec3d curEdgeVector = vectorNormalise(curClosest - curEdgePoint);
    vec3d prevEdgeVector = vectorNormalise(prevClosest - prevEdgePoint);

    bool curDirection = vectorDot(curEdgeVector, planeNormal) < 0;
    bool prevDirection = vectorDot(prevEdgeVector, planeNormal) < 0;

    if (curDirection == prevDirection) return;

    // Distance from current sphere center to plane
    double distToPlane = vectorDot(circle1.midP - v0.midP, -planeNormal);

    vec3d ballVelocity = circle1.midP - circle1.prevMid;

    // If inside the plane, push it out along the normal
    double penetrationDepth = circle1.radius - fabs(distToPlane);
    vec3d correctionDir = (distToPlane < 0) ? planeNormal : -planeNormal;
    circle1.midP = vectorAdd(circle1.midP, correctionDir * penetrationDepth);

    // Assume a mass of 1 for the ball and 10 for the plane
    float ballMass = 1.0f;
    float planeMass = 1.0f;

    // Get the normal and tangent velocity of each body
    vec3d normalV1 = -planeNormal * vectorDot(ballVelocity, -planeNormal);
    vec3d tangentV1 = ballVelocity - normalV1;

    vec3d normalV2 = -planeNormal * vectorDot(planeVelocity, -planeNormal);
    vec3d tangentV2 = planeVelocity - normalV2;

    // The velocities only change along the normal direction, so calculate that
    vec3d newNormalV1 = vectorDiv((normalV1 * (ballMass - planeMass)) + (normalV2 * 2 * planeMass), ballMass + planeMass);
    vec3d newNormalV2 = vectorDiv((normalV2 * (planeMass - ballMass)) + (normalV1 * 2 * ballMass), ballMass + planeMass);

    vec3d newVel1 = tangentV1 + newNormalV1;
    vec3d newVel2 = tangentV2 + newNormalV2;

    // Update the ball
    circle1.prevMid = circle1.midP - newVel1;

    // Update the planes
    if (!v0.anchored) v0.midP = v0.prevMid + newVel2;
    if (!v1.anchored) v1.midP = v1.prevMid + newVel2;
    if (!v2.anchored) v2.midP = v2.prevMid + newVel2;
}

int nearestSphere(const vector<object>& objects, vec3d forward, vec3d camera){
    // The highest dot product (the most similar vector) with the normal vector is the object the user clicked
    float maxDot = -2;
    int maxObj = 0;
    forward = vectorNormalise(forward);
    for (int i = 0; i < objects.size(); i++){
        object cur = objects[i];

        // Vector from camera to the current sphere
        vec3d objVec = vectorNormalise(cur.midP - camera);

        float dot = vectorDot(objVec, forward);

        if (dot > maxDot){
            maxDot = dot;
            maxObj = i;
        }
    }

    return maxObj;
}

void initDrawSprings(const vector<object>& objects, vector<Spring>& springs, int startIndex){
    float springRadius = 0.5;

    for (int i = startIndex; i < springs.size(); i++){
        Spring& spring = springs[i];

        object object1 = objects[spring.index[0]];
        object object2 = objects[spring.index[1]];

        spring.startIndex = meshObj.vertices.size();

        vec3d normal = vectorNormalise(object2.midP - object1.midP);
        normal = vectorMul(normal, springRadius);

        // Pick an arbitrary vector not parallel to normal
        vec3d arbitrary = (fabs(normal.x) > 0.9) ? vec3d{0,1,0} : vec3d{1,0,0};

        // Cross products to build an orthogonal basis
        vec3d u = vectorNormalise(vectorCross(normal, arbitrary));
        vec3d v = vectorCross(normal, u);

        // First object
        vec3d center = object1.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, springRadius * cos(theta)) 
                        + vectorMul(v, springRadius * sin(theta));

            // Push back the points
            meshObj.vertices.push_back(point);
        }

        // Second object
        center = object2.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, springRadius * cos(theta)) 
                        + vectorMul(v, springRadius * sin(theta));

            // Push back these points too
            meshObj.vertices.push_back(point);
        }

        // The average of each of the nodes
        vec3d springCol = vectorDiv((object1.color + object2.color), 2.0f);

        // Triangulate the links
        for (int i = 0; i < 4; i++) {
            int iNext = (i + 1) % 4;

            triangle tri1;
            tri1.p[0] = spring.startIndex + i;
            tri1.p[1] = spring.startIndex + 4 + i;
            tri1.p[2] = spring.startIndex + 4 + iNext;
            tri1.mat.diffuse = springCol;
            tri1.mat.shininess = 50;
            tri1.mat.specular = {0.5, 0.5, 0.5};
            tri1.cullBackface = false;

            triangle tri2;
            tri2.p[0] = spring.startIndex + i;
            tri2.p[1] = spring.startIndex + 4 + iNext;
            tri2.p[2] = spring.startIndex + iNext;
            tri2.mat.diffuse = springCol;
            tri2.mat.shininess = 50;
            tri2.mat.specular = {0.5, 0.5, 0.5};
            tri2.cullBackface = false;

            meshObj.tris.push_back(tri1);
            meshObj.tris.push_back(tri2);
        }
    }
}

void drawSprings(const vector<object>& objects, vector<Spring>& springs){
    float springRadius = 0.5;

    for (auto& spring : springs){
        object object1 = objects[spring.index[0]];
        object object2 = objects[spring.index[1]];

        vec3d normal = vectorNormalise(object2.midP - object1.midP);
        normal = vectorMul(normal, springRadius);

        // Pick an arbitrary vector not parallel to normal
        vec3d arbitrary = (fabs(normal.x) > 0.9) ? vec3d{0,1,0} : vec3d{1,0,0};

        // Cross products to build an orthogonal basis
        vec3d u = vectorNormalise(vectorCross(normal, arbitrary));
        vec3d v = vectorCross(normal, u);

        // First object
        vec3d center = object1.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, springRadius * cos(theta)) 
                        + vectorMul(v, springRadius * sin(theta));

            // Update the points (assume theyre already initialised)
            meshObj.vertices.at(spring.startIndex + i) = point;
        }

        // Second object
        center = object2.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, springRadius * cos(theta)) 
                        + vectorMul(v, springRadius * sin(theta));

            // Update these points too, but different index
            meshObj.vertices.at(spring.startIndex + i + 4) = point;
        }
    }
}

void initDrawLinks(const vector<object>& objects, vector<Link>& links, int startIndex){
    float linkRadius = 0.5;

    for (int i = startIndex; i < links.size(); i++){
        Link& link = links[i];

        object object1 = objects[link.index[0]];
        object object2 = objects[link.index[1]];

        link.startIndex = meshObj.vertices.size();

        vec3d normal = vectorNormalise(object2.midP - object1.midP);
        normal = vectorMul(normal, linkRadius);

        // Pick an arbitrary vector not parallel to normal
        vec3d arbitrary = (fabs(normal.x) > 0.9) ? vec3d{0,1,0} : vec3d{1,0,0};

        // Cross products to build an orthogonal basis
        vec3d u = vectorNormalise(vectorCross(normal, arbitrary));
        vec3d v = vectorCross(normal, u);

        // First object
        vec3d center = object1.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, linkRadius * cos(theta)) 
                        + vectorMul(v, linkRadius * sin(theta));

            // Push back the points
            meshObj.vertices.push_back(point);
        }

        // Second object
        center = object2.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, linkRadius * cos(theta)) 
                        + vectorMul(v, linkRadius * sin(theta));

            // Push back these points too
            meshObj.vertices.push_back(point);
        }

        vec3d linkCol = vectorDiv((object1.color + object2.color), 2.0f);

        // Triangulate the links
        for (int i = 0; i < 4; i++) {
            int iNext = (i + 1) % 4;

            triangle tri1;
            tri1.p[0] = link.startIndex + i;
            tri1.p[1] = link.startIndex + 4 + i;
            tri1.p[2] = link.startIndex + 4 + iNext;
            tri1.mat.diffuse = linkCol;
            tri1.mat.shininess = 50;
            tri1.mat.specular = {0.5, 0.5, 0.5};
            tri1.cullBackface = false;

            triangle tri2;
            tri2.p[0] = link.startIndex + i;
            tri2.p[1] = link.startIndex + 4 + iNext;
            tri2.p[2] = link.startIndex + iNext;
            tri2.mat.diffuse = linkCol;
            tri2.mat.shininess = 50;
            tri2.mat.specular = {0.5, 0.5, 0.5};
            tri2.cullBackface = false;

            meshObj.tris.push_back(tri1);
            meshObj.tris.push_back(tri2);
        }
    }
}

void drawLinks(const vector<object>& objects, vector<Link>& links){
    float linkRadius = 0.5;

    for (auto& link : links){
        object object1 = objects[link.index[0]];
        object object2 = objects[link.index[1]];

        vec3d normal = vectorNormalise(object2.midP - object1.midP);
        normal = vectorMul(normal, linkRadius);

        // Pick an arbitrary vector not parallel to normal
        vec3d arbitrary = (fabs(normal.x) > 0.9) ? vec3d{0,1,0} : vec3d{1,0,0};

        // Cross products to build an orthogonal basis
        vec3d u = vectorNormalise(vectorCross(normal, arbitrary));
        vec3d v = vectorCross(normal, u);

        // First object
        vec3d center = object1.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, linkRadius * cos(theta)) 
                        + vectorMul(v, linkRadius * sin(theta));

            // Update the points (assume theyre already initialised)
            meshObj.vertices.at(link.startIndex + i) = point;
        }

        // Second object
        center = object2.midP;

        for (int i = 0; i < 4; i++) {
            float theta = 2.0f * PI * i / 4;
            vec3d point = center 
                        + vectorMul(u, linkRadius * cos(theta)) 
                        + vectorMul(v, linkRadius * sin(theta));

            // Update these points too, but different index
            meshObj.vertices.at(link.startIndex + i + 4) = point;
        }
    }
}

// Lets you pick any position and radius for a sphere
void transformUnitSphere(object& sphere){
    // Change the object based on the radius
    for (int i = sphere.startIndex; i < sphere.endIndex; i++){
        meshObj.vertices.at(i) = vectorMul(meshObj.vertices.at(i), sphere.radius);

        //sphere.midP.y = -sphere.midP.y;
        meshObj.vertices.at(i) = vectorAdd(meshObj.vertices.at(i), sphere.midP);
    }
}

vec3d calculateNormal(vec3d side1, vec3d side2){
    return vectorNormalise(vectorCross(side1, side2));
}

void createObj(vector<object>& objects, vec3d midP, vec3d color, float radius, int objIndex, bool drawSphere, 
                 bool anchored, bool selfCollision, bool doCollide){

    object newObj;
    newObj.drawSphere = drawSphere;
    vec3d rotation = {0, 0, 0};
    
    if (newObj.drawSphere){
        newObj.startIndex = meshObj.vertices.size();
        newObj.startTri = meshObj.tris.size();
        parseObj("unit_sphere.obj", rotation);
        newObj.endTri = meshObj.tris.size();
        newObj.endIndex = meshObj.vertices.size();
    }

    newObj.midP = midP;
    newObj.prevMid = newObj.midP;
    newObj.acceleration = {0, 0, 0};
    newObj.anchored = anchored;
    newObj.radius = radius;
    newObj.color = color;
    newObj.objIndex = objIndex;
    newObj.selfCollision = selfCollision;
    newObj.doCollide = doCollide;
    transformUnitSphere(newObj);

    objects.push_back(newObj);
}

void createSpring(vector<object>& objects, vector<Spring>& springs, int index1, int index2, float k){
    Spring newSpring;
    newSpring.index[0] = index1;
    newSpring.index[1] = index2;
    newSpring.k = k;
    newSpring.restingDist = vectorLength(objects[index1].midP - objects[index2].midP);
    springs.push_back(newSpring);
}

void createLink(vector<object>& objects, vector<Link>& links, int index1, int index2, bool visual){
    Link newLink;
    newLink.index[0] = index1;
    newLink.index[1] = index2;
    newLink.distance = vectorLength(objects[index1].midP - objects[index2].midP);
    newLink.visual = visual;
    links.push_back(newLink);
}

void createFace(vector<physicsFace>& faces, int index1, int index2, int index3, color color, bool cullBackface){
    physicsFace newFace;
    newFace.index[0] = index1;
    newFace.index[1] = index2;
    newFace.index[2] = index3;
    newFace.color = color;
    newFace.cullBackface = cullBackface;

    faces.push_back(newFace);
}

Vertex vertexLerp(const Vertex& a, const Vertex& b, float t) {
    Vertex out;
    out.pos = a.pos + (b.pos - a.pos) * t;
    out.normal = vectorNormalise(a.normal + (b.normal - a.normal) * t);
    out.col = a.col + (b.col - a.col) * t;
    return out;
}

Vertex vertexIntersectPlane(const vec3d& plane_p, const vec3d& plane_n,
                            const Vertex& a, const Vertex& b, float& t){
    
    vec3d n = vectorNormalise(plane_n);
    float d = -vectorDot(n, plane_p);
    float ad = vectorDot(a.pos, n);
    float bd = vectorDot(b.pos, n);
    float denom = (bd - ad);
    const float EPS = 1e-6f;
    if (fabsf(denom) < EPS) {
        t = 0.0f; // or 0.5f; choose a convention
        return a; // segment parallel (or coincident with plane); pick an endpoint or handle specially
    }
    t = (-d - ad) / denom;
    return vertexLerp(a, b, t);
}

int triangleClipAgainstPLane(vec3d plane_p, vec3d plane_n, triangleP& in_tri, triangleP& out_tri1, triangleP& out_tri2){
    
    plane_n = vectorNormalise(plane_n);

    // Return if a point is inside or outside of a plane
    auto dist = [&](vec3d &p)
    {
        vec3d n = vectorNormalise(p);
        return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - vectorDot(plane_n, plane_p));
    };

    Vertex insidePoints[3]; int insidePointCount = 0;
    Vertex outsidePoints[3]; int outsidePointCount = 0;

    // Copy info from the input triangle
    Vertex p1, p2, p3;
    p1.pos = in_tri.p[0];
    p1.col = in_tri.vertCols[0];
    p1.normal = in_tri.normals[0];

    p2.pos = in_tri.p[1];
    p2.col = in_tri.vertCols[1];
    p2.normal = in_tri.normals[1];

    p3.pos = in_tri.p[2];
    p3.col = in_tri.vertCols[2];
    p3.normal = in_tri.normals[2];

    // Calculate whether a point is outside or inside
    float d0 = dist(p1.pos);
    float d1 = dist(p2.pos);
    float d2 = dist(p3.pos);

    if (d0 >= 0){
        insidePoints[insidePointCount++] = p1;
    }
    else {
        outsidePoints[outsidePointCount++] = p1;
    }

    if (d1 >= 0){
        insidePoints[insidePointCount++] = p2;
    }
    else {
        outsidePoints[outsidePointCount++] = p2;
    }

    if (d2 >= 0){
        insidePoints[insidePointCount++] = p3;
    }
    else {
        outsidePoints[outsidePointCount++] = p3;
    }

    // Classify the triangles
    if (insidePointCount == 0){
        return 0;
    }

    if (insidePointCount == 3){
        out_tri1 = in_tri;

        return 1;
    }

    // Forming a new triangle
    if (insidePointCount == 1 && outsidePointCount == 2){
        // Keep the inside point, because its valid
        out_tri1.p[0] = insidePoints[0].pos;
        out_tri1.normals[0] = insidePoints[0].normal;
        out_tri1.vertCols[0] = insidePoints[0].col;

        float t1, t2;
        Vertex intersect1 = vertexIntersectPlane(plane_p, plane_n, insidePoints[0], outsidePoints[0], t1);
        Vertex intersect2 = vertexIntersectPlane(plane_p, plane_n, insidePoints[0], outsidePoints[1], t2);

        // Update the other two points
        out_tri1.p[1] = intersect1.pos;
        out_tri1.normals[1] = intersect1.normal;
        out_tri1.vertCols[1] = intersect1.col;

        out_tri1.p[2] = intersect2.pos;
        out_tri1.normals[2] = intersect2.normal;
        out_tri1.vertCols[2] = intersect2.col;

        return 1;
    }

    // Forming a quad
    if (insidePointCount == 2 && outsidePointCount == 1){
        // Add the two new triangles
        out_tri1.p[0] = insidePoints[0].pos;
        out_tri1.normals[0] = insidePoints[0].normal;
        out_tri1.vertCols[0] = insidePoints[0].col;

        out_tri1.p[1] = insidePoints[1].pos;
        out_tri1.normals[1] = insidePoints[1].normal;
        out_tri1.vertCols[1] = insidePoints[1].col;

        float t1, t2;
        Vertex intersect1 = vertexIntersectPlane(plane_p, plane_n, insidePoints[0], outsidePoints[0], t1);
        Vertex intersect2 = vertexIntersectPlane(plane_p, plane_n, insidePoints[1], outsidePoints[0], t2);

        out_tri1.p[2] = intersect1.pos;
        out_tri1.normals[2] = intersect1.normal;
        out_tri1.vertCols[2] = intersect1.col;

        // Add the other triangle
        out_tri2.p[0] = insidePoints[1].pos;
        out_tri2.normals[0] = insidePoints[1].normal;
        out_tri2.vertCols[0] = insidePoints[1].col;

        out_tri2.p[1] = out_tri1.p[2];
        out_tri2.normals[1] = out_tri1.normals[2];
        out_tri2.vertCols[1] = out_tri1.vertCols[2];

        out_tri2.p[2] = intersect2.pos;
        out_tri2.normals[2] = intersect2.normal;
        out_tri2.vertCols[2] = intersect2.col;

        return 2;
    }

    return 0;
}

int findContainingBox(const vec3d& pt, const vec3d& minB, const vec3d& maxB, int nx, int nz) {
    // Compute per-axis cell sizes (only x and z)
    vec3d fullSize {
        maxB.x - minB.x,
        minB.y,
        maxB.z - minB.z
    };
    vec3d cellSize {
        fullSize.x / nx,
        maxB.y,
        fullSize.z / nz
    };

    // Map pt → (i,k)
    int i = int((pt.x - minB.x) / cellSize.x);
    int k = int((pt.z - minB.z) / cellSize.z);

    // Clamp to valid range
    i = std::clamp(i, 0, nx-1);
    k = std::clamp(k, 0, nz-1);

    // Flatten into 1D index (no y indexing)
    return k + i * nx;
}

void putTriangleInBoundingBox(const triangle& tri, int triIndex, vector<rectangle>& boundingBs,
                              const vec3d& globalMin, const vec3d& globalMax, int detail){

    array<vec3d, 3> vertices;
    vertices[0] = meshFindIndex(tri.p[0]);
    vertices[1] = meshFindIndex(tri.p[1]);
    vertices[2] = meshFindIndex(tri.p[2]);

    for (int i = 0; i < 3; i++){
        int boxIdx = findContainingBox(vertices[i], globalMin, globalMax, detail, detail);
        boundingBs.at(boxIdx).trisIndex.insert(triIndex);
    }
}

vec3d reflectVector(vec3d in, vec3d normal){
    return vectorSub(-in, vectorMul(normal, 2 * vectorDot(normal, -in)));
}

void render3D(vector<triangleP>& tris, const vector<rectangle>& boundingBs, const mesh& meshCube, const mat4x4& matProj, 
              const mat4x4& matWorld, const mat4x4& matView, const vec3d& camera, const mat4x4& matCamera){

    std::array<triangleP, 64> clippedTris;
    int clipped_track = 0;

    std::array<triangleP, 64> newTris;
    int new_track = 0;

    int maxSize = 0;
    for (rectangle box : boundingBs){
        maxSize += box.trisIndex.size();
    }

    tris.reserve(maxSize);

    // Extract the frustum planes from the projection matrix
    std::array<Plane, 6> frustumPlanesView = extractFrustumPlanes(matProj);
    std::array<Plane, 6> frustumPlanesWorld;

    mat4x4 matCombined = matrixMultiplyMatrix(matView, matProj);

    // Transform frustum planes to world space
    mat4x4 transposeView = matrixTranspose(matView);
    
    for (int i = 0; i < 6; i++){
        Plane plane = frustumPlanesView[i];
        vec3d planeVec = {plane.normal.x, plane.normal.y, plane.normal.z, plane.d};

        vec3d worldPlaneVec = matrixMultiplyVector(transposeView, planeVec);

        Plane worldPlane;
        worldPlane.normal = { worldPlaneVec.x, worldPlaneVec.y, worldPlaneVec.z };
        worldPlane.d = worldPlaneVec.w;

        frustumPlanesWorld[i] = worldPlane;
    }

    array<vec3d, 8> corners;

    // Loop over the bounding boxes and perform frustum clipping, and only draw triangles inside those bounding boxes
    for (const auto& box : boundingBs){
        // Rotate the boxes into world space
        vec3d minVert = matrixMultiplyVector(matWorld, box.minVert);
        vec3d maxVert = matrixMultiplyVector(matWorld, box.maxVert);

        // Test box against frustum
        if (!frustumBoundingBoxIntersection(frustumPlanesWorld, minVert, maxVert)) continue;

        for (int index : box.trisIndex){
            triangle tri = meshCube.tris[index];

            triangleP triTransformed;
            triangle triProjected;

            matrixMultiplyVectorFast(triTransformed.p[0], matWorld, meshCube.vertices[tri.p[0]]);
            matrixMultiplyVectorFast(triTransformed.p[1], matWorld, meshCube.vertices[tri.p[1]]);
            matrixMultiplyVectorFast(triTransformed.p[2], matWorld, meshCube.vertices[tri.p[2]]);

            // Copy the normal data
            matrixMultiplyVectorFast(triTransformed.normals[0], matWorld, tri.normals[0]);
            matrixMultiplyVectorFast(triTransformed.normals[1], matWorld, tri.normals[1]);
            matrixMultiplyVectorFast(triTransformed.normals[2], matWorld, tri.normals[2]);

            // Use the cross product to get surface normal
            vec3d normal, line1, line2;
            line1 = vectorSub(triTransformed.p[1], triTransformed.p[0]);
            line2 = vectorSub(triTransformed.p[2], triTransformed.p[0]);

            normal = vectorCross(line1, line2);
            normal = vectorNormalise(normal);

            vec3d cameraRay = vectorSub(triTransformed.p[0], camera);

            if (!tri.cullBackface || vectorDot(normal, cameraRay) > 0){
                // Clip against frustum planes
                clippedTris[0] = triTransformed;
                clipped_track = 1;

                bool switchNormal = !tri.cullBackface && vectorDot(normal, cameraRay) < 0;

                for (const auto& plane : frustumPlanesWorld) {
                    new_track = 0;
                    for (int i = 0; i < clipped_track; i++) {
                        triangleP t = clippedTris[i];
                        triangleP out1, out2;
                        int n = triangleClipAgainstPLane(plane.normal * -plane.d, plane.normal, t, out1, out2);
                        if (n == 1) newTris[new_track++] = out1;
                        if (n == 2) {
                            newTris[new_track++] = out1;
                            newTris[new_track++] = out2;
                        }
                    }
                    copy(newTris.begin(), newTris.begin() + new_track, clippedTris.begin());

                    clipped_track = new_track;
                    if (clipped_track == 0){
                        break;
                    }
                }

                SDL_Color finalColors[3];
                // Project and place into an array
                for (int c_track = 0; c_track < clipped_track; c_track++) {
                    triangleP clipped = clippedTris[c_track];
                    vec3d lightPos = {190, 1000, 190};

                    float ambient = 0.3f;     // baseline ambient (reduce to avoid oversaturation)
                    float emission = 0.2f;    // mild emission

                    for (int i = 0; i < 3; i++){
                        vec3d normal = switchNormal ? -clipped.normals[i] : clipped.normals[i];

                        vec3d lightDir = vectorNormalise(vectorSub(lightPos, clipped.p[i]));
                        vec3d viewDir  = vectorNormalise(vectorSub(camera, clipped.p[i]));
                        vec3d reflectDir = reflectVector(-lightDir, normal);

                        // Ambient + Emission using material color
                        vec3d vAmbient = vectorMul(tri.mat.diffuse, ambient);
                        vec3d vEmissive = vectorMul(tri.mat.diffuse, emission);

                        // Diffuse
                        float diffuse = fmax(vectorDot(normal, lightDir), 0.0f);
                        vec3d vDiffuse = vectorMul(tri.mat.diffuse, diffuse);

                        // Specular
                        float specular = pow(fmax(vectorDot(viewDir, reflectDir), 0.0f), tri.mat.shininess);
                        vec3d vSpecular = vectorMul(tri.mat.specular, specular);

                        // Combine all
                        vec3d final = vectorAdd(vectorAdd(vAmbient, vEmissive), vectorAdd(vDiffuse, vSpecular));

                        // Clamp final after scaling to 0-255
                        finalColors[i].r = fmin(final.x * 255.0f, 255.0f);
                        finalColors[i].g = fmin(final.y * 255.0f, 255.0f);
                        finalColors[i].b = fmin(final.z * 255.0f, 255.0f);
                    }


                    triangleP triProjected;

                    matrixMultiplyVectorFast(triProjected.p[0], matCombined, clipped.p[0]);
                    matrixMultiplyVectorFast(triProjected.p[1], matCombined, clipped.p[1]);
                    matrixMultiplyVectorFast(triProjected.p[2], matCombined, clipped.p[2]);

                    triProjected.p[0] = vectorDiv(triProjected.p[0], triProjected.p[0].w);
                    triProjected.p[1] = vectorDiv(triProjected.p[1], triProjected.p[1].w);
                    triProjected.p[2] = vectorDiv(triProjected.p[2], triProjected.p[2].w);

                    vec3d offsetView = {1, 1, 0};
                    triProjected.p[0] = vectorAdd(triProjected.p[0], offsetView);
                    triProjected.p[1] = vectorAdd(triProjected.p[1], offsetView);
                    triProjected.p[2] = vectorAdd(triProjected.p[2], offsetView);

                    triProjected.p[0].x *= 0.5f * WINDOW_WIDTH;
                    triProjected.p[0].y *= 0.5f * WINDOW_HEIGHT;
                    triProjected.p[1].x *= 0.5f * WINDOW_WIDTH;
                    triProjected.p[1].y *= 0.5f * WINDOW_HEIGHT;
                    triProjected.p[2].x *= 0.5f * WINDOW_WIDTH;
                    triProjected.p[2].y *= 0.5f * WINDOW_HEIGHT;

                    // Pass color information per vertex
                    for (int i = 0; i < 3; i++){
                        triProjected.vertCols[i].r = finalColors[i].r;
                        triProjected.vertCols[i].g = finalColors[i].g;
                        triProjected.vertCols[i].b = finalColors[i].b;
                        triProjected.vertCols[i].a = 255;
                    }

                    tris.push_back(triProjected);
                }
            }
        }
    }
}

int SDL_main(int argc, char* argv[]){
    SDL_Window* window;

    SDL_Init(SDL_INIT_VIDEO);

    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    window = SDL_CreateWindow("3D graphics", WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_OPENGL);

    // Integrate OpenGL
    SDL_GLContext context = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, context);

    // Load OpenGL functions using GLAD
    int out = gladLoadGLLoader((GLADloadproc)SDL_GL_GetProcAddress);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // Default: pass if incoming depth < stored depth

    // Compile vertex shader
    std::string vertexShaderCode = LoadShaderSource("shaderCode\\shader.vert");
    const char* vertexShaderSource = vertexShaderCode.c_str();

    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    // Compile fragment shader
    std::string fragShaderCode = LoadShaderSource("shaderCode\\shader.frag");
    const char* fragmentShaderSource = fragShaderCode.c_str();

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    // Link shaders into a program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    // Cleanup shaders (they're linked now)
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    bool run = true;
    SDL_Event event;

    vector<object> objects;
    vector<Link> links;
    vector<Spring> springs;
    vector<physicsFace> faces;
    
    vector<rectangle> boundingBs;

    // Rotate the map so that it looks horizontal
    meshObj.vertices.resize(1);
    vec3d rotation = {0, 0, 0};
    parseObj("grid.obj", rotation);

    // CUBE
    // int numX = 4;
    // int numY = 4;
    // int numZ = 4;
    // int radius = 1;

    // for (int x = 0; x < numX; x++){
    //     for (int y = 0; y < numY; y++){
    //         for (int z = 0; z < numZ; z++){
    //             object newObj;
    //             newObj.drawSphere = false;
                
    //             if (newObj.drawSphere){
    //                 newObj.startIndex = meshObj.vertices.size();
    //                 newObj.startTri = meshObj.tris.size();
    //                 parseObj("unit_sphere.obj", rotation);
    //                 newObj.endTri = meshObj.tris.size();
    //                 newObj.endIndex = meshObj.vertices.size();
    //             }

    //             newObj.midP = {60.0 + x * (radius * 2 + 1), 
    //                            -100.0 - y * (radius * 2 + 1),
    //                            60.0 + z * (radius * 2 + 1)};

    //             newObj.prevMid = newObj.midP;
    //             newObj.acceleration = {0, 0, 0};
    //             newObj.anchored = false;
    //             newObj.radius = radius;
    //             newObj.color = {0, 127, 255};
    //             newObj.objIndex = 1;
    //             newObj.selfCollision = false;
    //             transformUnitSphere(newObj);

    //             objects.push_back(newObj);
    //         }
    //     }
    // }

    // // Connect the indices
    // auto getFlatIndex = [numX, numY](int x, int y, int z) {
    //     return x + y * numX + z * (numX * numY);
    // };

    // for (int x = 0; x < numX; x++){
    //     for (int y = 0; y < numY; y++){
    //         for (int z = 0; z < numZ; z++){
    //             // Current index
    //             int i = getFlatIndex(x, y, z);

    //             // Right neighbour
    //             if (x + 1 < numX){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x+1, y, z);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Bottom neighbour
    //             if (y + 1 < numY){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x, y+1, z);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Depth neighbour
    //             if (z + 1 < numZ){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x, y, z+1);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Diagonal neighbour
    //             if (x + 1 < numX && y + 1 < numY){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x+1, y+1, z);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Right depth neighbour
    //             if (x + 1 < numX && z + 1 < numZ){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x+1, y, z+1);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Bottom depth neighbour
    //             if (y + 1 < numY && z + 1 < numZ){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x, y+1, z+1);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }

    //             // Diagonal depth neighbour
    //             if (x + 1 < numX && y + 1 < numY && z + 1 < numZ){
    //                 Link newLink;
    //                 newLink.index[0] = i;
    //                 newLink.index[1] = getFlatIndex(x+1, y+1, z+1);
    //                 newLink.distance = vectorLength(objects[newLink.index[0]].midP - objects[newLink.index[1]].midP);

    //                 links.push_back(newLink);
    //             }
    //         }
    //     }
    // }

    {
    int numX = 10;
    int numY = 10;
    int radius = 1;
    float springStrength = 0.005;

    int baseIndex = objects.size();

    for (int i = 0; i < numX * numY; i++){
        createObj(objects, 
                  {radius + (float)(i % numX) * (380 / 9 - radius), -200, 
                   radius + (float)(i / numY) * (380 / 9 - radius)},
                  {0, 127, 255}, 
                  radius, 
                  2, 
                  false, 
                  (i == 0 || i == numX - 1 || i == numX * numY - numX || i == numX * numY - 1),
                  false,
                  false);
    }

    for (int i = baseIndex; i < baseIndex + numX * numY; i++) {
        int localIndex = i - baseIndex;
        int row = localIndex / numX;
        int col = localIndex % numX;

        // Skip last row/column for faces
        if (row < numY - 1 && col < numX - 1) {
            createFace(faces, i + numX, i + 1, i, {0, 127, 255}, false);

            createFace(faces, i + numX, i + numX + 1, i + 1, {0, 127, 255}, false);
        }

        // Right neighbor spring
        if (col < numX - 1) {
            createSpring(objects, springs, i, i+1, springStrength);
        }

        // Down neighbor spring
        if (row < numY - 1) {
            createSpring(objects, springs, i, i + numX, springStrength);
        }

        // Down-right diagonal spring
        if (col < numX - 1 && row < numY - 1) {
            createSpring(objects, springs, i + 1, i + numX, springStrength);
        }
    }
    }

    // SPRING

    // {
    // int radius = 1;

    // vec3d positions[4] = {
    //     {50, -1, 50, 1},
    //     {150, -1, 50, 1},
    //     {150, -1, 150, 1},
    //     {50, -1, 150, 1}
    // };

    // // Bottom layer
    // for (int i = 0; i < 4; i++){
    //     createObj(objects, 
    //               positions[i],
    //               {0, 127, 255},
    //               radius,
    //               2,
    //               true,
    //               true,
    //               false,
    //               true
    //     );
    // }

    // createFace(faces, 2, 1, 0, {255, 255, 255}, false);

    // createFace(faces, 3, 2, 0, {255, 255, 255}, false);

    // // Connect them
    // for (int i = 0; i < 4; i++){
    //     createLink(objects, links, i, (i + 1) % 4, false);
    // }

    // // Top layer
    // for (int i = 0; i < 4; i++){
    //     createObj(objects, 
    //               {positions[i].x, positions[i].y - 50, positions[i].z},
    //               {0, 127, 255},
    //               radius,
    //               2,
    //               true,
    //               false,
    //               false,
    //               true
    //     );
    // }

    // createFace(faces, 6, 5, 4, {255, 255, 255}, false);

    // createFace(faces, 7, 6, 4, {255, 255, 255}, false);

    // // Connect them
    // for (int i = 0; i < 4; i++){
    //     createLink(objects, links, 4 + i, 4 + (i + 1) % 4, false);
    // }

    // // Connect the bottom and top layer
    // for (int i = 0; i < 4; i++){
    //     createSpring(objects, springs, 4 + i, i, 1e-5);
    // }

    // // Connect the bottom and top layer with diagonals
    // for (int i = 0; i < 4; i++){
    //     createSpring(objects, springs, i, 4 + (i + 1) % 4, 1e-3);

    //     createSpring(objects, springs, i, 4 + (i + 3) % 4, 1e-3);
    // }

    // createLink(objects, links, 4, 6, false);

    // createLink(objects, links, 5, 7, false);
    // }

    // {
    // int baseIndex = objects.size();
    // int radius = 1;

    // // Top layer
    // for (int x = 0; x < 2; x++){
    //     for (int y = 0; y < 2; y++){
    //         for (int z = 0; z < 2; z++){
    //             createObj(objects, 
    //                       {50.0 + x * (radius * 2 + 10),
    //                       -(200.0 + y * (radius * 2 + 10)),
    //                       50.0 + z * (radius * 2 + 10)},
    //                       {0, 127, 255},
    //                       1,
    //                       3,
    //                       true,
    //                       false,
    //                       false,
    //                       true
    //             );
    //         }
    //     }
    // }

    // auto getFlatIndex = [baseIndex](int x, int y, int z){
    //     return baseIndex + z + y * 2 + x * 4;
    // };

    // // Connect all cube edges
    // for (int x = 0; x < 2; x++) {
    //     for (int y = 0; y < 2; y++) {
    //         for (int z = 0; z < 2; z++) {
    //             int idx = getFlatIndex(x, y, z);

    //             // Neighbor offsets along positive axes
    //             if (x + 1 < 2) {
    //                 createLink(objects, links, idx, getFlatIndex(x + 1, y, z), false);
    //             }
    //             if (y + 1 < 2) {
    //                 createLink(objects, links, idx, getFlatIndex(x, y + 1, z), false);
    //             }
    //             if (z + 1 < 2) {
    //                 createLink(objects, links, idx, getFlatIndex(x, y, z + 1), false);
    //             }
    //         }
    //     }
    // }

    // // Connect cube space diagonals (4 total)
    // createLink(objects, links, getFlatIndex(0,0,0), getFlatIndex(1,1,1), false);

    // createLink(objects, links, getFlatIndex(0,0,1), getFlatIndex(1,1,0), false);

    // createLink(objects, links, getFlatIndex(0,1,0), getFlatIndex(1,0,1), false);

    // createLink(objects, links, getFlatIndex(1,0,0), getFlatIndex(0,1,1), false);
    // }

    // Initialisation
    for (auto& face : faces){
        face.startIndex = meshObj.vertices.size();

        // Vertices of the triangle
        for (int i = 0; i < 3; i++){
            meshObj.vertices.push_back(objects[face.index[i]].midP);
        }

        face.endIndex = meshObj.vertices.size();

        vec3d side1 = objects[face.index[0]].midP - objects[face.index[1]].midP;
        vec3d side2 = objects[face.index[0]].midP - objects[face.index[2]].midP;

        vec3d normal = calculateNormal(side1, side2);

        // Create objects, because we can only draw links between two objects
        // The midpoint of the triangle face
        vec3d midP = vectorDiv(objects[face.index[0]].midP + objects[face.index[1]].midP + objects[face.index[2]].midP, 3.0f);
        createObj(objects, midP, {0, 0, 255}, 1, 2, false, false, false, false);

        // The unit normal point
        vec3d normalP = midP + normal;

        face.normal = objects.size();
        createObj(objects, normalP, {0, 0, 255}, 1, 2, false, false, false, false);

        // Create a link
        createLink(objects, links, face.normal - 1, face.normal, true);

        triangle newTri;
        newTri.p[0] = face.startIndex;
        newTri.p[1] = face.startIndex + 1;
        newTri.p[2] = face.startIndex + 2;
        newTri.cullBackface = face.cullBackface;
        newTri.mat.diffuse = {face.color.r / 255, face.color.g / 255, face.color.b / 255};
        newTri.mat.specular = {0.5, 0.5, 0.5};
        newTri.mat.shininess = 50;

        meshObj.tris.push_back(newTri);
    }

    initDrawLinks(objects, links, 0);
    initDrawSprings(objects, springs, 0);

    size_t maxVertices = meshObj.tris.size() * 18;
    DynamicMesh meshTri = CreateDynamicMesh(maxVertices);
    DynamicMesh meshLine = CreateDynamicMesh(maxVertices);

    // 3 vertices per triangle
    // 6 floats per vertex
    // 18 floats per triangle
    vector<float> drawArray(meshObj.tris.size() * 18);

    // 3 lines per triangle
    // 3 vertices per line
    // 5 floats per vertex
    vector<float> lineArray(meshObj.tris.size() * 45);

    // Projection matrix
    float near = 0.1f;
    float far = 10000.0f;
    float FOV = 90.0f;
    float aspectRatio = (float)WINDOW_HEIGHT / (float) WINDOW_WIDTH;

    double FPS = 120;
    double speed = 1;

    mat4x4 matProj;
    matProj = matrixMakeProjection(FOV, aspectRatio, near, far);

    //cout << rotation.x << ' ' << rotation.y << ' ' << rotation.z << endl;

    mat4x4 matRotX, matRotY, matRotZ;
    matRotX = matrixMakeRotationX(rotation.x * PI / 180.0);
    matRotY = matrixMakeRotationY(rotation.y * PI / 180.0);
    matRotZ = matrixMakeRotationZ(rotation.z * PI / 180.0);

    mat4x4 matTrans;
    matTrans = matrixMakeTranslation(0.0f, 0.0f, 0.0f);

    mat4x4 matWorld;
    matWorld = matrixMakeIdentity();
    matWorld = matrixMultiplyMatrix(matWorld, matRotX);
    matWorld = matrixMultiplyMatrix(matWorld, matRotY);
    matWorld = matrixMultiplyMatrix(matWorld, matRotZ);
    matWorld = matrixMultiplyMatrix(matWorld, matTrans);

    vector<vec3d> vertexNormals(meshObj.vertices.size(), {0.0f, 0.0f, 0.0f});

    // Update vertex normals
    for (const triangle& tri : meshObj.tris) {
        // Grab the three vertex positions
        vec3d v0 = meshObj.vertices[tri.p[0]];
        vec3d v1 = meshObj.vertices[tri.p[1]];
        vec3d v2 = meshObj.vertices[tri.p[2]];

        // Compute edges
        vec3d edge1 = vectorSub(v1, v0);
        vec3d edge2 = vectorSub(v2, v0);

        // Un‐normalized face normal (magnitude == parallelogram area)
        vec3d faceCross = vectorCross(edge1, edge2);

        // Convert to triangle area and unit normal
        float triArea = 0.5f * vectorLength(faceCross);
        vec3d faceNormal = vectorNormalise(faceCross);

        // Weighted face normal (area factor is optional but smooths better on irregular meshes)
        vec3d weightedNormal = vectorMul(faceNormal, triArea);

        // Accumulate into each of the three vertex slots
        vertexNormals[tri.p[0]] = vectorAdd(vertexNormals[tri.p[0]], weightedNormal);
        vertexNormals[tri.p[1]] = vectorAdd(vertexNormals[tri.p[1]], weightedNormal);
        vertexNormals[tri.p[2]] = vectorAdd(vertexNormals[tri.p[2]], weightedNormal);
    }

    for (size_t i = 0; i < vertexNormals.size(); ++i) {
        vertexNormals[i] = vectorNormalise(vertexNormals[i]);
    }

    for (triangle& tri : meshObj.tris) {
        tri.normals[0] = vertexNormals[tri.p[0]];
        tri.normals[1] = vertexNormals[tri.p[1]];
        tri.normals[2] = vertexNormals[tri.p[2]];
    }

    int detail = 10;

    // Create the enclosing rectangle to start off with
    rectangle encloseRect;
    createEnclosingRectangle(meshObj, encloseRect);
    boundingBs.push_back(encloseRect);

    vec3d initMin = encloseRect.minVert;
    vec3d initMax = encloseRect.maxVert;

    // Subdivide the bounding box
    boundingBs = makeGrid(encloseRect, detail);

    // // Add our enclosing rect to the mesh so we can visualise it
    // addBoundingBoxToMesh(meshObj, boundingBs);

    SDL_Color colors[9] = {
        {255, 0, 0, 255},   // Red
        {180, 255, 180, 255},   // Green
        {0, 0, 255, 255},   // Blue
        {255, 255, 0, 255}, // Yellow
        {255, 165, 0, 255}, // Orange
        {128, 0, 128, 255}, // Purple
        {0, 255, 255, 255}, // Cyan
        {255, 192, 203, 255}, // Pink
        {255, 255, 255, 255}, // White
    };


    SDL_Color white = {255, 255, 255, 255};
    SDL_Color black = {0, 0, 0, 255};

    bool rightMouseDown = false;
    bool doRotation = true;
    SDL_FPoint startPan;

    float yaw = 0;
    float pitch = 0;

    vec3d camera;

    clock_t startTime = clock();
    clock_t endTime = clock();
    double deltaTime = 0;
    double elapsedSeconds = 0;
    double fixedStep = 1000.0f / 120.0f;
    double forceMagnitude = 1;

    int lastFPS = FPS;
    int frames = 0;
    int numClippedTris = 0;
    int numTrianglesDrawn = 0;
    int numIterations = 0;

    int previousTime = 0;
    int nextTime = 0;

    vec3d gravity = {0, 0.001, 0};

    bool step = false;

    while (run) {

        bool stepFrame = false;

        while (step && !stepFrame){
            while (SDL_PollEvent(&event)){
                switch (event.type){
                    case SDL_EVENT_KEY_DOWN: {
                        if (event.key.key == SDLK_RIGHT){
                            stepFrame = true;
                        }
                        // allow the user to switch between stepping frames and running the program normally
                        if (event.key.key == SDLK_K){
                            step = !step;
                        }
                        break;
                    }
                    case SDL_EVENT_QUIT: {
                        stepFrame = true;
                        run = false;
                        break;
                    }
                }
            }
        }

        startTime = clock();

        // Update the bounding boxes because the scene is dynamic
        for (auto& box : boundingBs){
            box.trisIndex.clear();
        }

        // Put the triangles into the bounding boxes
        for (size_t i = 0; i < meshObj.tris.size(); ++i) {
            putTriangleInBoundingBox(meshObj.tris[i], i, boundingBs, initMin, initMax, detail);
        }

        size_t quarter = boundingBs.size() / 4;
        int remainder = boundingBs.size() % 4;

        // First quarter  
        vector<rectangle> boxesA;
        vector<triangleP> computeA;

        // Second quarter  
        vector<rectangle> boxesB;
        vector<triangleP> computeB;

        // Third quarter  
        vector<rectangle> boxesC;
        vector<triangleP> computeC;

        // Fourth quarter  
        vector<rectangle> boxesD;
        vector<triangleP> computeD;

        boxesA.insert(boxesA.begin(), boundingBs.begin(), boundingBs.begin() + quarter);
        boxesB.insert(boxesB.begin(), boundingBs.begin() + quarter, boundingBs.begin() + quarter * 2);
        boxesC.insert(boxesC.begin(), boundingBs.begin() + quarter * 2, boundingBs.begin() + quarter * 3);
        boxesD.insert(boxesD.begin(), boundingBs.begin() + quarter * 3, boundingBs.begin() + quarter * 4 + remainder);

        for (auto& object : objects) {
            // Verlet integration
            if (!object.anchored){
                vec3d oldPos = object.midP;
                object.midP = object.midP * 2.0f - object.prevMid + object.acceleration * fixedStep * fixedStep;
                object.prevMid = oldPos;

                // Gravity
                object.acceleration = gravity;

                // Update the velocity vector
                vec3d velocity = object.midP - object.prevMid;

                // Put the point the proper distance away from the current midPoint of our ball
                if (object.velocityVec != -1){
                    vec3d vector = object.midP;
                    vec3d radiusVec = vectorNormalise(velocity) * object.radius;
                    vec3d scale = velocity * 5;

                    vector = vector + radiusVec + scale;

                    objects[object.velocityVec].midP = vector;
                    objects[object.velocityVec].prevMid = vector;
                }
            }
        }

        for (int iter = 0; iter < 30; iter++){
            // Apply spring forces
            for (auto& spring : springs) {
                springForce(objects, spring);
            }

            // Apply link forces
            for (auto& link : links) {
                // Don't apply physics to this link, it only exists to visualize some link
                if (link.visual) continue;

                linkForce(objects, link);
            }

            // Apply collisions (pairwise)
            for (int i = 0; i < objects.size(); i++) {
                for (int j = i+1; j < objects.size(); j++) {
                    // Prevent collision between two balls that are part of the same object
                    if (objects[i].objIndex == objects[j].objIndex && !objects[i].selfCollision) continue;

                    // Prevent collision between any objects
                    if (!objects[i].doCollide || !objects[j].doCollide) continue;

                    cout << i << ' ' << j << endl;

                    collisionForce(objects[i], objects[j]);
                }
            }

            // Apply face collisions
            for (auto& object : objects){
                collisionForceFace(objects, object, faces);
            }

            // Apply bounds
            for (auto& object : objects) {
                boundForce(object);
            }
        }

        // 3. Apply displacement to mesh AFTER physics is stable
        for (auto& object : objects){
            if (!object.anchored){
                // Compute current geometric center of this sphere’s vertices
                vec3d center = {0,0,0};
                int count = object.endIndex - object.startIndex;
                for (int i = object.startIndex; i < object.endIndex; i++) {
                    center = vectorAdd(center, meshObj.vertices[i]);
                }
                center = center * (1.0f / count);

                // Compute offset needed to align sphere with physics midP
                vec3d correction = object.midP - center;

                // Apply correction to all vertices of this sphere
                for (int i = object.startIndex; i < object.endIndex; i++) {
                    meshObj.vertices[i] = vectorAdd(meshObj.vertices[i], correction);
                }
            }
        }
        
        // Update the triangles too
        for (auto& face : faces){
            object v0 = objects[face.index[0]];
            object v1 = objects[face.index[1]];
            object v2 = objects[face.index[2]];

            // Update normals
            vec3d side1 = v0.midP - v1.midP;
            vec3d side2 = v0.midP - v2.midP;

            vec3d normal = calculateNormal(side1, side2);

            vec3d midP = vectorDiv(v0.midP + v1.midP + v2.midP, 3.0f);
            vec3d planeVelocity = vectorDiv((v0.midP - v0.prevMid + v1.midP - v1.prevMid + v2.midP - v2.prevMid), 3.0f);

            // The unit normal point (to visualize the normal)
            vec3d normalP = midP - normal * -planeVelocity * 5;

            objects[face.normal - 1].midP = midP;
            objects[face.normal - 1].prevMid = midP;

            objects[face.normal].midP = normalP;
            objects[face.normal].prevMid = objects[face.normal].midP;

            vec3d displacements[3];
            
            displacements[0] = objects[face.index[0]].midP - meshObj.vertices[face.startIndex];
            displacements[1] = objects[face.index[1]].midP - meshObj.vertices[face.startIndex + 1];
            displacements[2] = objects[face.index[2]].midP - meshObj.vertices[face.startIndex + 2];

            for (int i = 0; i < 3; i++){
                vec3d& updateVec = meshObj.vertices[face.startIndex + i];

                updateVec = vectorAdd(updateVec, displacements[i]);
            }
        }

        // Set the colors of the balls based on if theyre anchored or not
        for (auto& object: objects){
            for (int i = object.startTri; i < object.endTri; i++){
                meshObj.tris.at(i).mat.diffuse = object.anchored ? vec3d{255, 0, 0} : object.color;
            }
        }

        drawLinks(objects, links);
        drawSprings(objects, springs);

        // How many times to print the FPS and number of drawn triangles every second
        const int drawPerSecond = 24;
        previousTime = (int)elapsedSeconds / (1000 / drawPerSecond);
        elapsedSeconds += max(1000 * deltaTime, 1000.0 / FPS);
        frames++;

        nextTime = (int)elapsedSeconds / (1000 / drawPerSecond);

        if (nextTime - previousTime > 0){
            std::cout << "FPS: " << lastFPS << ' ' << "Drawn triangles: " << numTrianglesDrawn << std::endl;
            //std::cout << camera.x << ' ' << camera.y << ' ' << camera.z << endl;
        }

        // Update the rotation matrices
        matRotX = matrixMakeRotationX(rotation.x * PI / 180.0);
        matRotY = matrixMakeRotationY(rotation.y * PI / 180.0);
        matRotZ = matrixMakeRotationZ(rotation.z * PI / 180.0);

        matWorld = matrixMakeIdentity();
        matWorld = matrixMultiplyMatrix(matWorld, matRotX);
        matWorld = matrixMultiplyMatrix(matWorld, matRotY);
        matWorld = matrixMultiplyMatrix(matWorld, matRotZ);
        matWorld = matrixMultiplyMatrix(matWorld, matTrans);

        vector<vec3d> vertexNormals(meshObj.vertices.size(), {0.0f, 0.0f, 0.0f});

        // Update vertex normals
        for (const triangle& tri : meshObj.tris) {
            // Grab the three vertex positions
            vec3d v0 = meshObj.vertices[tri.p[0]];
            vec3d v1 = meshObj.vertices[tri.p[1]];
            vec3d v2 = meshObj.vertices[tri.p[2]];

            // Compute edges
            vec3d edge1 = vectorSub(v1, v0);
            vec3d edge2 = vectorSub(v2, v0);

            // Un‐normalized face normal (magnitude == parallelogram area)
            vec3d faceCross = vectorCross(edge1, edge2);

            // Convert to triangle area and unit normal
            float triArea = 0.5f * vectorLength(faceCross);
            vec3d faceNormal = vectorNormalise(faceCross);

            // Weighted face normal (area factor is optional but smooths better on irregular meshes)
            vec3d weightedNormal = vectorMul(faceNormal, triArea);

            // Accumulate into each of the three vertex slots
            vertexNormals[tri.p[0]] = vectorAdd(vertexNormals[tri.p[0]], weightedNormal);
            vertexNormals[tri.p[1]] = vectorAdd(vertexNormals[tri.p[1]], weightedNormal);
            vertexNormals[tri.p[2]] = vectorAdd(vertexNormals[tri.p[2]], weightedNormal);
        }

        for (size_t i = 0; i < vertexNormals.size(); ++i) {
            vertexNormals[i] = vectorNormalise(vertexNormals[i]);
        }

        for (triangle& tri : meshObj.tris) {
            tri.normals[0] = vertexNormals[tri.p[0]];
            tri.normals[1] = vertexNormals[tri.p[1]];
            tri.normals[2] = vertexNormals[tri.p[2]];
        }

        // Rotate the camera
        // yaw and pitch in radians
        // Convert cordinates on a sphere to a vector
        float cordX = cosf(pitch) * sinf(yaw);
        float cordY = sinf(pitch);
        float cordZ = cosf(pitch) * cosf(yaw);
        vec3d lookDir = {cordX, cordY, cordZ};

        // Use world up for a standard FPS camera
        vec3d worldUp = {0, 1, 0};

        // Recompute camera right and up
        vec3d rightDir = vectorNormalise(vectorCross(lookDir, worldUp));
        vec3d camUp = vectorCross(rightDir, lookDir);

        // Build view matrix
        vec3d target = vectorAdd(camera, lookDir);
        mat4x4 matCamera = matrixPointAt(camera, target, camUp);
        mat4x4 matView = matrixQuickInverse(matCamera);

        // Get mouse state
        SDL_PumpEvents();

        float x, y;
        Uint32 buttons = SDL_GetMouseState(&x, &y);

        vec3d right = vectorMul(vectorCross(camUp, lookDir), speed);
        vec3d forward = vectorMul(lookDir, speed);
        vec3d up = vectorMul(camUp, speed);

        if (buttons & SDL_BUTTON_MASK(SDL_BUTTON_LEFT)){
            // Spawn a ball
            int index1 = objects.size();
            createObj(objects, vec3d{camera.x, camera.y, camera.z} + forward, {0, 0, 255}, 10, 0, true, false, true, true);

            // Velocity vector point
            int index2 = objects.size();
            createObj(objects, vec3d{camera.x, camera.y, camera.z} + forward, {255, 0, -255}, 10, 0, false, false, true, false);

            // The place where the velocity vector point is stored
            objects[index1].velocityVec = index2;

            int startIndex = links.size();
            // Create a visual link between them
            createLink(objects, links, index1, index2, true);

            initDrawLinks(objects, links, startIndex);
        }

        const bool* keyboardState = SDL_GetKeyboardState(NULL);
        while (SDL_PollEvent(&event)){
            switch (event.type){
                case SDL_EVENT_QUIT:
                    run = false;
                    break;
                case SDL_EVENT_MOUSE_WHEEL:
                    speed -= event.wheel.y * speed / 10;
                    break;
                case SDL_EVENT_MOUSE_BUTTON_DOWN:
                    if (event.button.button == SDL_BUTTON_RIGHT){
                        rightMouseDown = true;
                        startPan = {event.button.x, event.button.y};
                    }
                    break;
                case SDL_EVENT_MOUSE_BUTTON_UP:
                    if (event.button.button == SDL_BUTTON_RIGHT){
                        rightMouseDown = false;
                    }
                    break;
                case SDL_EVENT_KEY_DOWN:
                    // Anchor the object the user is looking at
                    if (event.key.key == SDLK_R){
                        int nearest = nearestSphere(objects, forward, camera);
                        cout << nearest << endl;

                        objects.at(nearest).anchored = !objects.at(nearest).anchored;
                    }
                    break;
            }
        }

        if (rightMouseDown){
            SDL_FPoint cur = {x, y};

            yaw -= (startPan.x - cur.x) / 100;
            pitch -= (startPan.y - cur.y) / 100;

            if (pitch < -PI / 2){
                pitch = -PI / 2;
            }
            if (pitch > PI / 2){
                pitch = PI / 2;
            }

            startPan = cur;
        }

        if (keyboardState[SDL_SCANCODE_Q]){
            camera = vectorAdd(camera, up);
        }
        if (keyboardState[SDL_SCANCODE_E]){
            camera = vectorSub(camera, up);
        }
        if (keyboardState[SDL_SCANCODE_A]){
            camera = vectorSub(camera, right);
        }
        if (keyboardState[SDL_SCANCODE_D]){
            camera = vectorAdd(camera, right);
        }
        if (keyboardState[SDL_SCANCODE_W]){
            camera = vectorAdd(camera, forward);
        }
        if (keyboardState[SDL_SCANCODE_S]){
            camera = vectorSub(camera, forward);
        }


        if (keyboardState[SDL_SCANCODE_F]){
            for (auto& object : objects){
                vec3d force = camera - object.midP;
                double distance = vectorLength(force);

                if (distance < 30 && !object.anchored){
                    object.midP = object.prevMid + vectorMul(forward, forceMagnitude);
                }
            }
        }

        computeA.clear();
        computeB.clear();
        computeC.clear();
        computeD.clear();

        // Create four threads (max amount of concurrency on my machine)
        thread t1(
            render3D,
            std::ref(computeA),
            std::cref(boxesA),
            std::cref(meshObj), std::cref(matProj),
            std::cref(matWorld), std::cref(matView),
            std::cref(camera), std::cref(matCamera));

        thread t2(
            render3D,
            std::ref(computeB),
            std::cref(boxesB),
            std::cref(meshObj), std::cref(matProj),
            std::cref(matWorld), std::cref(matView),
            std::cref(camera), std::cref(matCamera));

        thread t3(
            render3D,
            std::ref(computeC),
            std::cref(boxesC),
            std::cref(meshObj), std::cref(matProj),
            std::cref(matWorld), std::cref(matView),
            std::cref(camera),  std::cref(matCamera));

        thread t4(
            render3D,
            std::ref(computeD),
            std::cref(boxesD),
            std::cref(meshObj), std::cref(matProj),
            std::cref(matWorld), std::cref(matView),
            std::cref(camera),  std::cref(matCamera));

        t1.join();
        t2.join();
        t3.join();
        t4.join();

        // Just put the arrays together and let the depth buffer handle things
        std::vector<triangleP> vecTrianglesToRaster;

        vecTrianglesToRaster.reserve(computeA.size() + computeB.size() + computeC.size() + computeD.size());

        vecTrianglesToRaster.insert(vecTrianglesToRaster.end(), computeA.begin(), computeA.end());
        vecTrianglesToRaster.insert(vecTrianglesToRaster.end(), computeB.begin(), computeB.end());
        vecTrianglesToRaster.insert(vecTrianglesToRaster.end(), computeC.begin(), computeC.end());
        vecTrianglesToRaster.insert(vecTrianglesToRaster.end(), computeD.begin(), computeD.end());

        numTrianglesDrawn = 0;

        float avg = 0;
        float numAvg = 0;

        for (auto& triToRaster : vecTrianglesToRaster){
            float points[3][3] = { 0 };
            numTrianglesDrawn++;
            for (int i = 0; i < 3; i++){
                float newX = (triToRaster.p[i].x / (float)WINDOW_WIDTH) * 2.0f - 1.0f;
                float newY = 1.0f - (triToRaster.p[i].y / (float)WINDOW_HEIGHT) * 2.0f;

                // Triangles
                drawArray.push_back(newX);
                drawArray.push_back(newY);
                drawArray.push_back(triToRaster.p[i].z);
                drawArray.push_back(triToRaster.vertCols[i].r / 255);
                drawArray.push_back(triToRaster.vertCols[i].g / 255);
                drawArray.push_back(triToRaster.vertCols[i].b / 255);

                if (i == 0 || i == 1) {
                    // First two vertices: just store them and also emit line vertices
                    lineArray.push_back(newX);
                    lineArray.push_back(newY);
                    lineArray.push_back(triToRaster.p[i].z - 0.0001);                // z for lines
                    lineArray.push_back(black.r / 255.0f);
                    lineArray.push_back(black.g / 255.0f);
                    lineArray.push_back(black.b / 255.0f);

                    points[i][0] = newX;
                    points[i][1] = newY;
                    points[i][2] = triToRaster.p[i].z - 0.0001;
                } else {
                    // Close edges: (v0, v2) and (v1, v2)
                    // Edge v0 -> v2
                    lineArray.push_back(points[0][0]); // x0
                    lineArray.push_back(points[0][1]); // y0
                    lineArray.push_back(points[0][2]);         // z
                    lineArray.push_back(black.r / 255.0f);
                    lineArray.push_back(black.g / 255.0f);
                    lineArray.push_back(black.b / 255.0f);

                    lineArray.push_back(newX);         // x2
                    lineArray.push_back(newY);         // y2
                    lineArray.push_back(triToRaster.p[i].z - 0.0001);         // z
                    lineArray.push_back(black.r / 255.0f);
                    lineArray.push_back(black.g / 255.0f);
                    lineArray.push_back(black.b / 255.0f);

                    // Edge v1 -> v2
                    lineArray.push_back(points[1][0]); // x1
                    lineArray.push_back(points[1][1]); // y1
                    lineArray.push_back(points[1][2]);         // z
                    lineArray.push_back(black.r / 255.0f);
                    lineArray.push_back(black.g / 255.0f);
                    lineArray.push_back(black.b / 255.0f);

                    lineArray.push_back(newX);         // x2
                    lineArray.push_back(newY);         // y2
                    lineArray.push_back(triToRaster.p[i].z - 0.0001);         // z
                    lineArray.push_back(black.r / 255.0f);
                    lineArray.push_back(black.g / 255.0f);
                    lineArray.push_back(black.b / 255.0f);
                }
            }
        }

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT);

        // Triangles
        UpdateDynamicMesh(meshTri, drawArray);
        GLsizei triCount = static_cast<GLsizei>(drawArray.size() / 6); // 6 floats per vertex

        // Lines
        UpdateDynamicMesh(meshLine, lineArray);
        GLsizei lineCount = static_cast<GLsizei>(lineArray.size() / 6);

        // Draw both
        DrawDynamicMeshTriangle(meshTri, shaderProgram, triCount);
        //DrawDynamicMeshLine(meshLine, shaderProgram, lineCount);

        drawArray.clear();
        lineArray.clear();

        SDL_GL_SwapWindow(window);

        if (elapsedSeconds >= 1000){
            lastFPS = frames;

            frames = 0;
            elapsedSeconds = 0;
        }

        // Cap the FPS rate
        endTime = clock();
        deltaTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;

        if (deltaTime < 1000.0 / FPS){
            SDL_Delay(1000.0 / FPS - deltaTime);
        }
    }

    return 0;
}


