#include <algorithm>

#include "dynamicSDF.h"

using namespace CGL;

// 简单渲染：在体素/模型空间里画两个辅助小球便于定位（可按需删去）
void DynamicSDFObject::render(GLShader& shader) {
    // 由屏幕上的 ray-marching shader 负责显示真实体；
    // 这里保持空，避免与画面认知冲突。
    (void)shader;
}

// 场景 SDF：两个随时间移动的原始体并集
double DynamicSDFObject::sdf_scene(const Vector3D& p, double t) const {
    // 动态球体
    Vector3D sphere_center = Vector3D(sin(t * 0.5) * 2.0, 0.0, 0.0);
    double d_sphere = sdf_sphere(p - sphere_center, 1.0);

    // 动态盒子
    Vector3D box_center = Vector3D(cos(t * 0.5) * 2.0, 0.0, 0.0);
    double d_box = sdf_box(p - box_center, Vector3D(0.8, 0.8, 0.8));

    return std::min(d_sphere, d_box); // 并集
}

// 数值法线
Vector3D DynamicSDFObject::sdf_normal(const Vector3D& p, double t) const {
    const double h = 1e-4;
    const Vector3D kx(h, 0, 0), ky(0, h, 0), kz(0, 0, h);
    double dx = sdf_scene(p + kx, t) - sdf_scene(p - kx, t);
    double dy = sdf_scene(p + ky, t) - sdf_scene(p - ky, t);
    double dz = sdf_scene(p + kz, t) - sdf_scene(p - kz, t);
    Vector3D n(dx, dy, dz);
    double len = n.norm();
    if (len < 1e-12)
        return Vector3D(0, 1, 0);
    return n / len;
}

void DynamicSDFObject::collide(PointMass& pm) {
    if (!m_enabled)
        return;

    const double mu = std::max(0.0, std::min(1.0, m_friction));

    // 线段 CCD + SDF marching（保守步进）
    Vector3D x0 = pm.last_position;
    Vector3D x1 = pm.position;
    Vector3D seg = x1 - x0;
    double L = seg.norm();
    if (L < 1e-12)
        return;
    Vector3D dir = seg / L;

    // marching 参数
    const int max_steps = 64;
    const double eps = 1e-5; // 命中阈值

    Vector3D p = x0;
    double traveled = 0.0;
    bool hit = false;
    Vector3D hit_p;

    for (int i = 0; i < max_steps; i++) {
        double d = sdf_scene(p, m_time);
        if (d <= eps) {
            hit = true;
            hit_p = p;
            break;
        }
        double step = std::min(std::abs(d), L - traveled);
        if (step <= 1e-8)
            break;
        p = p + dir * step;
        traveled += step;
        if (traveled >= L - 1e-9)
            break;
    }

    if (!hit) {
        // 兜底：若终点处在内部也算命中
        if (sdf_scene(x1, m_time) <= 0.0) {
            hit = true;
            hit_p = x1;
        }
    }

    if (hit) {
        // 投影到表面
        double d = sdf_scene(hit_p, m_time);
        Vector3D n = sdf_normal(hit_p, m_time);
        Vector3D place = hit_p - d * n; // 更靠近表面的点

        // 摩擦响应（与现有实现保持一致的手感）
        pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
    }
}
