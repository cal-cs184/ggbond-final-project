#ifndef COLLISIONOBJECT_SPHERE_H
#define COLLISIONOBJECT_SPHERE_H

#include "../misc/sphere_drawing.h"
#include "collisionObject.h"

using namespace CGL;
using namespace std;

struct Sphere : public CollisionObject {
  public:
    Sphere(const Vector3D& origin, double radius, double friction, int num_lat = 40, int num_lon = 40)
        : origin(origin), radius(radius), radius2(radius * radius), friction(friction),
          m_sphere_mesh(Misc::SphereMesh(num_lat, num_lon)) {}

    void render(GLShader& shader);
    void collide(PointMass& pm);

    void set_use_sdf(bool enable) { use_sdf = enable; }
    bool get_use_sdf() const { return use_sdf; }
    void set_use_ccd(bool enable) { use_ccd = enable; }
    bool get_use_ccd() const { return use_ccd; }
    void set_use_ray_marching(bool enable) { use_ray_marching = enable; }
    bool get_use_ray_marching() const { return use_ray_marching; }
    void set_hide_when_dynamic_sdf(bool hide) { hide_when_dynamic_sdf = hide; }
    bool get_hide_when_dynamic_sdf() const { return hide_when_dynamic_sdf; }

  private:
    Vector3D origin;
    double radius;
    double radius2;

    double friction;

    Misc::SphereMesh m_sphere_mesh;

    bool use_sdf = false;
    bool use_ccd = false;
    bool use_ray_marching = false;
    bool hide_when_dynamic_sdf = false;
};

#endif /* COLLISIONOBJECT_SPHERE_H */
