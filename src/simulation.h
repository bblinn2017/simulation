#ifndef SIMULATION_H
#define SIMULATION_H

#include "graphics/shape.h"
#include <memory>

class Shader;

struct Tet;

class Simulation
{
public:
    Simulation();

    void init();

    void update(float seconds);

    void draw(Shader *shader);

    void toggleWire();

    bool drag(Eigen::Vector3f d, Eigen::Vector3f v, Eigen::Vector3f d_v);

    bool paused = true;
private:
    Shape m_shape;
    Shape m_ground;
    Shape m_obj;
    std::vector<Eigen::Vector3i> surface;

    void initGround();
    void initObject();

    std::vector<Eigen::Vector3f> getForces(std::vector<Eigen::Vector3f> positions, std::vector<Eigen::Vector3f> velocities);
    Eigen::Matrix3f getBarycentric(std::vector<Eigen::Vector3f> v);

    Eigen::Vector3f sphereCollision(Eigen::Vector3f vec);
    Eigen::Vector3f groundCollision(Eigen::Vector3f vec);

    std::vector<Eigen::Vector3f> pos;
    std::vector<Eigen::Vector3f> vel;
    std::vector<float> mass;

    std::vector<Eigen::Matrix3f> betas;
    std::vector<Eigen::Vector4i> tets;

    int numVerts;
    float groundHeight;

    const Eigen::Vector3f _gravity = Eigen::Vector3f(0.f, -1.f, 0.f);
    const float _lambda = 1e0f; //incompressibility for the whole material
    const float _myu = 1e1f; //rigidity for the whole material
    const float _phi = 1.5f; //coefficients of viscosity
    const float _psi = 1.5f;
    const float _rho = 5.f; //density
    const float _penalty(float dist) {
        return pow(dist,0.1)*150;
    } // penalty function

};

#endif // SIMULATION_H
