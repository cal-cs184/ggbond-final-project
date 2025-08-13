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

    // Configurable parameters
    void set_num_samples(int samples) { m_num_samples = samples; }
    void set_surface_threshold(double threshold) { m_surface_threshold = threshold; }
    void set_max_distance(double distance) { m_max_distance = distance; }
    void set_min_sphere_size(double size) { m_min_sphere_size = size; }
    void set_max_sphere_size(double size) { m_max_sphere_size = size; }
    void set_motion_speed(double speed) { m_motion_speed = speed; }
    void set_sphere_radius(double radius) { m_sphere_radius = radius; }
    void set_box_size(double size) { m_box_size = size; }

    // Getters
    int get_num_samples() const { return m_num_samples; }
    double get_surface_threshold() const { return m_surface_threshold; }
    double get_max_distance() const { return m_max_distance; }
    double get_min_sphere_size() const { return m_min_sphere_size; }
    double get_max_sphere_size() const { return m_max_sphere_size; }
    double get_motion_speed() const { return m_motion_speed; }
    double get_sphere_radius() const { return m_sphere_radius; }
    double get_box_size() const { return m_box_size; }

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

    // Configurable parameters
    int m_num_samples = 8000;
    double m_surface_threshold = 0.02;
    double m_max_distance = 2.5;
    double m_min_sphere_size = 0.003;
    double m_max_sphere_size = 0.04;
    double m_motion_speed = 2.0;
    double m_sphere_radius = 1.0;
    double m_box_size = 0.8;
};

#endif // COLLISIONOBJECT_DYNAMIC_SDF_H
