#include "simulation.h"

#include <iostream>
#include <unordered_map>

#include "graphics/MeshLoader.h"
#include "tiny_obj_loader.h"
#include <QFileInfo>
#include <QString>
#include "Eigen/Geometry"

static std::string filename = "example-meshes/cube.mesh";

using namespace Eigen;
using namespace std;

Simulation::Simulation()
{
}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    std::vector<Vector3f> vertices;
    std::vector<Vector4i> tetras;
    if(MeshLoader::loadTetMesh(filename, vertices, tetras)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...
        std::vector<Vector3i> faces;

        // Get surface tets
        std::unordered_map<std::string,int> faceCount;
        // For every tetrahedron
        for (Vector4i tet : tetras) {
            // For every face on this tetrahedron
            for (int i = 0; i < 4; i++) {
                std::vector<int> verts = {tet[i],tet[(i+1)%4],tet[(i+2)%4]};
                std::sort(verts.begin(),verts.end());
                std::string key = std::to_string(verts[0])+" "+std::to_string(verts[1])+" "+std::to_string(verts[2]);
                // Add to count
                if(!faceCount.count(key)) {
                    faceCount[key] = 0;
                }
                faceCount[key] += 1;
            }
        }

        mass = std::vector<float>(vertices.size());
        for (Vector4i tet : tetras) {
            // Store the tet
            tets.push_back(tet);

            // Calculate the beta
            std::vector<Vector3f> v = {vertices[tet[0]],vertices[tet[1]],
                    vertices[tet[2]],vertices[tet[3]]};
            betas.push_back(getBarycentric(v).inverse());

            // Calculate volume of v
            float V = abs((v[1]-v[0]).cross(v[2]-v[0])
                        .dot(v[3]-v[0]))/6.f;
            float m = _rho*V;

            // Add quarter of mass to each vertex
            for (int i = 0; i < 4; i++) {
                mass[tet[i]] += m/4.f;
            }

            // For every face on this tetrahedron
            for (int i = 0; i < 4; i++) {
                std::vector<int> verts = {tet[i],tet[(i+1)%4],tet[(i+2)%4]};
                std::sort(verts.begin(),verts.end());
                std::string key = std::to_string(verts[0])+" "+std::to_string(verts[1])+" "+std::to_string(verts[2]);
                // Emplace if there's only 1
                if (faceCount[key] == 1) {
                    // Get counterclockwise face
                    Vector3f a = vertices[verts[1]] - vertices[verts[0]];
                    Vector3f b = vertices[verts[2]] - vertices[verts[0]];

                    // Compare cross product to vector from centroid of this face to v4
                    Vector3f inNorm = vertices[tet[(i+3)%4]] - vertices[verts[0]];

                    // a cross b points in, choose opposite direction
                    if (a.cross(b).dot(inNorm) > 0.f) {
                        faces.emplace_back(verts[0],verts[2],verts[1]);
                    } else {
                        faces.emplace_back(verts[0],verts[1],verts[2]);
                    }
                }
            }
        }

        // Store position and velocity
        for (Vector3f vec : vertices) {
            pos.push_back(vec);
            vel.emplace_back(0.f,0.f,0.f);
        }

        // Store num verts
        numVerts = vertices.size();

        m_shape.init(vertices, faces, tetras);
        surface = faces;
    }

    groundHeight = -2.f;
    m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, -groundHeight, 0)));

    initGround();
    initObject();
}

void Simulation::update(float seconds)
{
    // STUDENTS: This method should contain all the time-stepping logic for your simulation.
    //   Specifically, the code you write here should compute new, updated vertex positions for your
    //   simulation mesh, and it should then call m_shape.setVertices to update the display with those
    //   newly-updated vertices.

    // STUDENTS: As currently written, the program will just continually compute simulation timesteps as long
    //    as the program is running (see View::tick in view.cpp) . You might want to e.g. add a hotkey for pausing
    //    the simulation, and perhaps start the simulation out in a paused state.

    // Calculate initial derivative position and velocity
    std::vector<Vector3f> forces = getForces(pos,vel);

    // Calculate initial vel derivative for every vertex
    std::vector<Vector3f> d_vel = std::vector<Vector3f>(numVerts);
    for (int i = 0; i < numVerts; i++) {
        d_vel[i] = forces[i] / mass[i] + _gravity;

        // Add ground force if needed
        Vector3f penalty = groundCollision(pos[i]);//+ sphereCollision(pos[i]);
        if (penalty != Vector3f::Zero()) {
            d_vel[i] += penalty;
        }
    }
    // Calculate pos derivative for every vertex
    std::vector<Vector3f> d_pos = vel;

    // Calculate midpoint position and velocity
    std::vector<Vector3f> midPos = std::vector<Vector3f>(numVerts);
    std::vector<Vector3f> midVel = std::vector<Vector3f>(numVerts);
    for (int i = 0; i < numVerts; i++) {
        midVel[i] = vel[i] + d_vel[i]*seconds/2.f;
        midPos[i] = pos[i] + d_pos[i]*seconds/2.f;
    }

    // Calculate midpoint derivative position and velocity
    std::vector<Vector3f> midForces = getForces(midPos,midVel);

    // Calculate midpoint vel derivative for every vertex
    std::vector<Vector3f> midD_vel = std::vector<Vector3f>(numVerts);
    for (int i = 0; i < numVerts; i++) {
        midD_vel[i] = midForces[i] / mass[i] + _gravity;

        // Add ground force if needed
        Vector3f penalty = groundCollision(midPos[i]);// + sphereCollision(midPos[i]);
        if (penalty != Vector3f::Zero()) {
            midD_vel[i] += penalty;
        }
    }
    // Calculate pos derivative for every vertex
    std::vector<Vector3f> midD_pos = midVel;

    // Calculate final position and velocity
    for (int i = 0; i < numVerts; i++) {
        vel[i] = midVel[i] + midD_vel[i]*seconds/2.f;
        pos[i] = midPos[i] + midD_pos[i]*seconds/2.f;
    }

    // Update shape
    m_shape.setVertices(pos);
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
//    m_obj.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

std::vector<Vector3f> Simulation::getForces(std::vector<Vector3f> positions, std::vector<Vector3f> velocities) {
    std::vector<Vector3f> forces = std::vector<Vector3f>();
    for (int i = 0; i < numVerts; i++) {
        forces.push_back(Vector3f::Zero());
    }

    for (int i = 0; i < tets.size(); i++) {
        Vector4i tet = tets[i];

        std::vector<Vector3f> p = {positions[tet[0]],positions[tet[1]],
                positions[tet[2]],positions[tet[3]]};
        Matrix3f P = getBarycentric(p);
        Matrix3f F_P = P*betas[i];
        Matrix3f strain_P = F_P.transpose()*F_P - Eigen::Matrix3f::Identity();
        Matrix3f stress = _lambda*strain_P.trace()*Eigen::Matrix3f::Identity() + 2.f*_myu*strain_P;

        std::vector<Vector3f> v = {velocities[tet[0]],velocities[tet[1]],
                velocities[tet[2]],velocities[tet[3]]};
        Matrix3f V = getBarycentric(v);
        Matrix3f F_V = V*betas[i];
        Matrix3f strain_V = F_P.transpose()*F_V + F_V.transpose()*F_P;
        stress += _phi*strain_V.trace()*Eigen::Matrix3f::Identity() + 2.f*_psi*strain_V;

        for (int j = 0; j < 4; j++) {
            Vector3f a = positions[tet[(j+2)%4]] - positions[tet[(j+1)%4]];
            Vector3f b = positions[tet[(j+3)%4]] - positions[tet[(j+1)%4]];

            // Compare cross product to vector from centroid of this face to v4
            Vector3f in = positions[tet[j]] - positions[tet[(j+1)%4]];

            Vector3f norm = 0.5f*a.cross(b);
            norm = norm.dot(in) > 0 ? norm : -norm;

            // Vertex j's force
            Vector3f force = -stress*norm;

            forces[tet[j]] += force;
        }
    }
    return forces;
}

Matrix3f Simulation::getBarycentric(std::vector<Vector3f> v) {
    Matrix3f mat;
    mat << v[0]-v[3],
        v[1]-v[3],
        v[2]-v[3];
    return mat;
}

void Simulation::initGround()
{
    std::vector<Vector3f> groundVerts;
    std::vector<Vector3i> groundFaces;

    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);

    m_ground.init(groundVerts, groundFaces);
}

void Simulation::initObject() {
    std::vector<Vector3f> _vertexArray;
    std::vector<Vector3i> _faceArray;

    tinyobj::attrib_t attrib;
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    QFileInfo info(QString("sphere-obj/spherical.obj"));
    string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                                info.absoluteFilePath().toStdString().c_str(), (info.absolutePath().toStdString() + "/").c_str(), true);
    if(!err.empty()) {
        cerr << err << endl;
    }

    if(!ret) {
        cerr << "Failed to load/parse .obj file" << endl;
        return;
    }

    for(size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            Vector3i face;
            for(size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

        face[v] = idx.vertex_index;

            }
            _faceArray.push_back(face);

            index_offset += fv;
        }
    }

    float scaleFactor = 1.f/Vector3f(attrib.vertices[0], attrib.vertices[1], attrib.vertices[2]).norm();
    Matrix4f scale = Affine3f(Scaling<float>(scaleFactor,scaleFactor,scaleFactor)).matrix();

    for(size_t i = 0; i < attrib.vertices.size(); i += 3) {
        Vector4f vec = scale*Vector4f(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2], 1.f);
        _vertexArray.emplace_back(vec[0], vec[1], vec[2]);
    }

    m_obj.init(_vertexArray,_faceArray);
}

Vector3f Simulation::groundCollision(Vector3f vec) {
    float dist = Vector4f(vec[0],vec[1],vec[2],-1.f).dot(Vector4f(0,1,0,groundHeight));

    if (dist < 0.f) {
        return Vector3f(0,1,0)*_penalty(abs(dist));
    }
    return Vector3f::Zero();
}

Vector3f Simulation::sphereCollision(Vector3f vec) {
    float currRad = (vec-Vector3f(0,groundHeight,0)).norm();

    if (currRad < 1.f) {
        Vector3f norm = vec-Vector3f(0,groundHeight,0);
        norm.normalize();
        return norm*_penalty(vec.norm());
    }
    return Vector3f::Zero();
}

bool Simulation::drag(Vector3f d, Vector3f vec, Vector3f d_v) {
    vec += Vector3f(0,groundHeight,0);

    float min_t = INFINITY;
    Vector3f min_coplanar = Vector3f::Zero();
    Vector3i min_face = Vector3i::Zero();

    for (Vector3i face : surface) {
        std::vector<Vector3f> tri = {
            pos[face[0]], pos[face[1]],pos[face[2]]
        };

        Vector3f norm = (tri[1]-tri[0]).cross(tri[2]-tri[0]);

        if (norm.dot(d) == 0.f) {
            continue;
        } else if (norm.dot(d) > 0.f) {
            norm = -norm;
        }

        // Find plane intersection
        float t = (norm.dot(tri[0]-vec))/(norm.dot(d));
        Vector3f coplanar = vec+t*d;

        Vector3f u = tri[1]-tri[0];
        Vector3f v = tri[2]-tri[0];
        Vector3f v_ind = v - v.dot(u)*u;

        Vector3f x = coplanar-tri[0];

        // Find u and v components
        float u_comp = u.dot(x)/pow(u.norm(),2.0);
        float v_comp = v_ind.dot(x)/pow(v_ind.norm(),2.0);

        // Check for intersection
        if (u_comp > 0.f && v_comp > 0.f && u_comp + v_comp < 1.f) {
            if (t < min_t) {
                min_t = t;
                min_coplanar = coplanar;
                min_face = face;
            }
        }
    }

    if (isinf(min_t)) {
        return false;
    }

    int index = 0;
    for (int i = 0; i < 3; i++) {
        float dist = (min_coplanar-pos[min_face[index]]).norm();
        if ((min_coplanar-pos[min_face[i]]).norm() < dist) {
            index = i;
        }
    }

    pos[min_face[index]] += d_v*(1.f*min_t);

    m_shape.setVertices(pos);
    return true;
}
