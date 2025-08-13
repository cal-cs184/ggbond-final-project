#ifndef COLLISIONOBJECT_DYNAMIC_SDF_H
#define COLLISIONOBJECT_DYNAMIC_SDF_H

#include "../misc/sphere_drawing.h"
#include "collisionObject.h"

// A dynamic collision object based on SDF, used to demonstrate the advantages of CPU-side "ray marching / SDF"
// collision Shape example: A sphere oscillating on the X-axis over time + an oscillating box, taking the union
struct DynamicSDFObject : public CollisionObject {
  public:
    explicit DynamicSDFObject(double friction = 0.2)
        : m_friction(friction), m_time(0.0), m_enabled(false), m_debug_sphere(Misc::SphereMesh(24, 24)) {}

    // For external per-frame updates
    void set_time(double t) { m_time = t; }
    void set_enabled(bool e) { m_enabled = e; }
    bool enabled() const { return m_enabled; }

    // Required interfaces
    void render(GLShader& shader) override;
    void collide(PointMass& pm) override;

  private:
    // SDF definitions (simplified version consistent with shaders)
    double sdf_sphere(const Vector3D& p, double r) const { return p.norm() - r; }

    double sdf_box(const Vector3D& p, const Vector3D& b) const {
        Vector3D q = Vector3D(fabs(p.x) - b.x, fabs(p.y) - b.y, fabs(p.z) - b.z);
        double outside = Vector3D(std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0)).norm();
        double inside = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0);
        return outside + inside;
    }

    // Dynamic scene: union(sphere, box)
    double sdf_scene(const Vector3D& p, double t) const;

    // SDF normal (for placement)
    Vector3D sdf_normal(const Vector3D& p, double t) const;

    // SDF point cloud rendering (more accurate SDF visualization)
    void renderSDFPointCloud(GLShader& shader);

  private:
    double m_friction; // [0,1]
    double m_time;
    bool m_enabled;
    Misc::SphereMesh m_debug_sphere; // Only for visualization
};

#endif // COLLISIONOBJECT_DYNAMIC_SDF_H
