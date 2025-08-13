#ifndef COLLISIONOBJECT_DYNAMIC_SDF_H
#define COLLISIONOBJECT_DYNAMIC_SDF_H

#include "../misc/sphere_drawing.h"
#include "collisionObject.h"

// 一个基于SDF的动态碰撞体，用于展示 CPU 端“ray marching / SDF”碰撞优势
// 形状示例：随时间在 X 轴上往复的球体 + 往复的盒子，取并集
struct DynamicSDFObject : public CollisionObject {
  public:
    explicit DynamicSDFObject(double friction = 0.2)
        : m_friction(friction), m_time(0.0), m_enabled(false), m_debug_sphere(Misc::SphereMesh(24, 24)) {}

    // 供外部每帧更新
    void set_time(double t) { m_time = t; }
    void set_enabled(bool e) { m_enabled = e; }
    bool enabled() const { return m_enabled; }

    // 必要接口
    void render(GLShader& shader) override;
    void collide(PointMass& pm) override;

  private:
    // SDF 定义（与着色器中保持一致的简化版本）
    double sdf_sphere(const Vector3D& p, double r) const { return p.norm() - r; }

    double sdf_box(const Vector3D& p, const Vector3D& b) const {
        Vector3D q = Vector3D(fabs(p.x) - b.x, fabs(p.y) - b.y, fabs(p.z) - b.z);
        double outside = Vector3D(std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0)).norm();
        double inside = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0);
        return outside + inside;
    }

    // 动态场景：并集(球, 盒子)
    double sdf_scene(const Vector3D& p, double t) const;

    // SDF 法线（用于放置）
    Vector3D sdf_normal(const Vector3D& p, double t) const;

  private:
    double m_friction; // [0,1]
    double m_time;
    bool m_enabled;
    Misc::SphereMesh m_debug_sphere; // 仅用于可视化
};

#endif // COLLISIONOBJECT_DYNAMIC_SDF_H
