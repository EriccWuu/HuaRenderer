#include <iostream>
#include <sstream>
#include "Model.h"

Model::Model(const std::string filename) {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Can't open file" << filename << '\n';
        return;
    }

    std::vector<vec3> verts{};
    std::vector<vec3> norms{};
    std::vector<vec2> texcoords{};
    std::vector<Index> faces;
    std::vector<int> face_vert{};
    std::vector<int> face_tex{};
    std::vector<int> face_norm{};

    std:: string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            for (int i = 0; i < 3; i ++) iss >> v[i];
            verts.push_back(v);
        }
        else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            vec3 n;
            for (int i = 0; i < 3; i ++) iss >> n[i];
            norms.push_back(n.normalize());
        }
        else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            vec2 uv;
            for (int i = 0; i < 2; i ++) iss >> uv[i];
            uv.y = 1 - uv.y;    // flip in y = 0.5
            texcoords.push_back(uv);
        }
        else if (!line.compare(0, 2, "f ")) {
            int f, t, n;
            iss >> trash;
            int cnt = 0;
            while (iss >> f >> trash >> t >> trash >> n) {
                faces.push_back({--f, --t, --n});
                face_vert.push_back(f);
                face_tex.push_back(t);
                face_norm.push_back(n);
                cnt ++;
            }
            if (cnt != 3) {
                std::cerr << "Error: the obj file is supposed to be triangulated." << std::endl;
                return ;
            }
        }
    }
    attribute.vertices.swap(verts);
    attribute.normals.swap(norms);
    attribute.texcoords.swap(texcoords);
    mesh.indices.swap(faces);

    makeVerticesUnique();

    std::string texfile;
    texfile = texture_path(filename, "_diffuse.tga");
    // texDiffuse = Texture(texfile, "diffuse");
    // diffuseTexPtr = std::make_shared<Texture>(texDiffuse);
    diffuseTexPtr = std::make_shared<Texture>(texfile, "diffuse");
    texfile = texture_path(filename, "_spec.tga");
    // texSpecular = Texture(texfile, "specular");
    // specularTexPtr = std::make_shared<Texture>(texSpecular);
    specularTexPtr = std::make_shared<Texture>(texfile, "specular");
    texfile = texture_path(filename, "_nm_tangent.tga");
    // texNormal = Texture(texfile, "normal");
    // normalTexPtr = std::make_shared<Texture>(texNormal);
    normalTexPtr = std::make_shared<Texture>(texfile, "normal");

    std::cerr << "- Model Info:" << " #v#" << nverts() << " #f#" << nfaces() << " #vt#" << attribute.texcoords.size() << " #vn#" << attribute.normals.size() << std::endl;
}

std::string Model::texture_path(const std::string &filename, const std::string &suffix) {
    size_t dot = filename.find_last_of(".");
    if (dot == std::string::npos) return "";
    std::string texfile = filename.substr(0, dot) + suffix;
    return texfile;
}

void Model::load_texture(const std::string filename, const std::string suffix, TGAImage &img) { 
    auto texfile = texture_path(filename, suffix);
    if (texfile == "") return;
    std::cerr << "Loading texture file \"" << texfile << "\"...\n";
    std::cerr << (img.read_tga_file(texfile.c_str()) ? "- Load successfully" : "- Fail to load.") << std::endl;
}

void Model::makeVerticesUnique() {
    std::vector<Vertex> vertices;
    std::vector<size_t> indices;
    std::unordered_map<Vertex, int> uniqueVert;
    int n = mesh.indices.size() / 3;
    for (int i = 0; i < n; ++ i) {
        auto i0 = mesh.indices[3*i];
        auto i1 = mesh.indices[3*i + 1];
        auto i2 = mesh.indices[3*i + 2];
        vec3 v0 = attribute.vertices[i0.vertex_index];
        vec3 v1 = attribute.vertices[i1.vertex_index];
        vec3 v2 = attribute.vertices[i2.vertex_index];
        vec2 uv0 = attribute.texcoords[i0.texcoord_index];
        vec2 uv1 = attribute.texcoords[i1.texcoord_index];
        vec2 uv2 = attribute.texcoords[i2.texcoord_index];
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec2 duv1 = uv1 - uv0;
        vec2 duv2 = uv2 - uv0;
        double k = 1.f / (duv1.x * duv2.y - duv1.y * duv2.x);
        vec3 tangent = k * (edge1 * duv2.y - edge2 * duv1.y);

        for (int j = 0; j < 3; ++ j) {
            auto idx = mesh.indices[3 * i + j];
            Vertex vertex;
            vertex.pos = vec4(attribute.vertices[idx.vertex_index], 1.f);
            vertex.normal = vec4(attribute.normals[idx.normal_index], 0.f);
            vertex.texcoord = attribute.texcoords[idx.texcoord_index];
            vertex.color = ONE_VEC3;
            vertex.tangent = vertex.tangent + tangent;
            if (uniqueVert.count(vertex) == 0) {
                uniqueVert[vertex] = (size_t)vertices.size();
                vertices.push_back(vertex);
            }
            indices.push_back(uniqueVert[vertex]);
        }
    }
    uniqueVerticesPtr = std::make_shared<std::vector<Vertex>>(vertices);
    uniqueIndicesPtr = std::make_shared<std::vector<size_t>>(indices);
}

std::shared_ptr<std::vector<Vertex>> Model::uniqueVertices() {
    return uniqueVerticesPtr;
}

std::shared_ptr<std::vector<size_t>> Model::uniqueIndices() {
    return uniqueIndicesPtr;
}

std::shared_ptr<Texture> Model::diffuseTexture() {
    return diffuseTexPtr;
}

std::shared_ptr<Texture> Model::specularTexture() {
    return specularTexPtr;
}

std::shared_ptr<Texture> Model::normalTexture() {
    return normalTexPtr;
}

int Model::nverts() const { 
    return attribute.vertices.size();
}

int Model::nfaces() const { 
    return mesh.indices.size() / 3;
}

Index Model::face(const int i) const { 
    return mesh.indices[i];
}

vec3 Model::normal(const int iface, const int nthvert) const { 
    return attribute.normals[mesh.indices[3 * iface + nthvert].normal_index];
}

vec3 Model::vert(const int i) const { 
    return attribute.vertices[i]; 
}

vec3 Model::vert(const int iface, const int nthvert) const { 
    return attribute.vertices[mesh.indices[3 * iface + nthvert].vertex_index];
}

vec2 Model::uv(const int iface, const int nthvert) const { 
    return attribute.texcoords[mesh.indices[3 * iface + nthvert].texcoord_index];
}

vec4 Model::diffuse(vec2 uv) const {
    return diffuseTexPtr->value(uv);
}

vec4 Model::specular(vec2 uv) const {
    return specularTexPtr->value(uv);
}

vec4 Model::normal(const vec2 &uv) const { 
    return normalTexPtr->value(uv);
}