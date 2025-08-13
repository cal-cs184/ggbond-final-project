#include <algorithm>
#include <nanogui/nanogui.h>

#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

// 关键参数（先走保守值，观察FPS后再放宽）
static const int kMaxMarchSteps = 48;          // marching 兜底步数
static const double kCCDTriggerRatio = 0.3;    // 仅当 L > 0.3R 时跑CCD
static const double kNearBandMultiplier = 3.0; // 距离表面 < 3*thickness 才跑CCD
static const double kDegenerateLenEps = 1e-12;
static const double kDegenerateDir2Eps = 1e-16;

// 解析线段-球相交（半径已含 thickness），返回是否命中 & 命中点x_hit（限制在段内）
static bool segmentSphereHit(const Vector3D& x0, const Vector3D& x1, const Vector3D& center, double radius,
                             Vector3D& x_hit) {
    Vector3D d = x1 - x0;
    double L = d.norm();
    if (L < kDegenerateLenEps)
        return false;
    Vector3D dir = d / L;

    Vector3D oc = x0 - center;
    double b = 2.0 * dot(oc, dir);
    double c = dot(oc, oc) - radius * radius;
    double disc = b * b - 4.0 * c; // a=1
    if (disc < 0.0)
        return false;

    double sqrtD = sqrt(std::max(0.0, disc));
    double t0 = (-b - sqrtD) * 0.5;
    double t1 = (-b + sqrtD) * 0.5;

    // 取第一次进入的解（段内）
    if (t0 >= 0.0 && t0 <= L) {
        x_hit = x0 + t0 * dir;
        return true;
    }
    if (t1 >= 0.0 && t1 <= L) {
        x_hit = x0 + t1 * dir;
        return true;
    }
    return false;
}

void Sphere::collide(PointMass& pm) {
    // ---- 公共参数 ----
    const double mu = clamp(friction, 0.0, 1.0); // (friction)
    // 建议用粒子间距: thickness = 0.5 * particle_spacing; 这里给个保底
    // const double thickness = std::max(1e-4, 1e-3 * radius);
    // const double Rexp = radius + thickness; // 扩展半径（检测与放置一致）
    // const double sdf_epsilon = std::max(1e-6, 0.5 * thickness);

    auto safeDir = [&](const Vector3D& v) -> Vector3D {
        double l = v.norm();
        if (l > kDegenerateLenEps)
            return v / l;
        Vector3D last = pm.last_position - origin;
        if (last.norm2() > kDegenerateDir2Eps)
            return last.unit();
        return Vector3D(0, 1, 0);
    };

    // ---- Legacy 点-球（最便宜）----
    if (!use_sdf && !use_ccd && !use_ray_marching) {
        Vector3D to = pm.position - origin;
        // Vector3D to = pm.last_position - origin;
        if (to.norm2() <= radius2) {
            Vector3D n = safeDir(to);
            Vector3D tgt = origin + radius * n;
            pm.position = pm.last_position + (1.0 - mu) * (tgt - pm.last_position);
        }
        return;
    }

    // ---- 优先解析CCD（仅在"可能撞"时触发）----
    if (use_ccd) {
        Vector3D x0 = pm.last_position;
        Vector3D x1 = pm.position; // 预测位置
        Vector3D seg = x1 - x0;
        double L = seg.norm();

        // 触发条件1：步长足够大
        if (L > kCCDTriggerRatio * radius) {
            // 粗判定A：起点离表面太远且位移不够覆盖 => 不可能命中（triangle inequality）
            double d0 = (x0 - origin).norm() - radius;
            if (d0 <= L * kNearBandMultiplier) { // 放宽一点，避免漏边界
                // 解析交：命中则直接放置（极快）
                Vector3D xhit;
                if (segmentSphereHit(x0, x1, origin, radius, xhit)) {
                    Vector3D n = safeDir(xhit - origin);
                    Vector3D place = origin + radius * n;

                    pm.position = pm.position + (1.0 - mu) * (place - pm.position);
                    return;
                }
            }
        }
        // 走到这里：解析CCD判否或 miss ——> 兜底：SDF DCD
    }

    // ---- 纯Ray Marching碰撞检测 ----
    if (use_ray_marching) {
        Vector3D ray_origin = pm.last_position;
        Vector3D ray_direction = pm.position - pm.last_position;
        double ray_length = ray_direction.norm();

        if (ray_length < kDegenerateLenEps) {
            return; // 质点没有移动
        }

        ray_direction = ray_direction / ray_length; // 归一化

        // Ray marching参数
        const int max_steps = kMaxMarchSteps;
        const double step_size = ray_length / max_steps;
        const double epsilon = 1e-6;

        Vector3D current_pos = ray_origin;
        bool hit_found = false;
        Vector3D hit_point;

        // Ray marching循环
        for (int i = 0; i < max_steps; i++) {
            // 计算当前点到球心的距离
            Vector3D to_center = current_pos - origin;
            double distance_to_center = to_center.norm();

            // 计算有向距离场（SDF）
            double sdf = distance_to_center - radius;

            // 检查是否命中球体表面
            if (sdf <= epsilon) {
                hit_found = true;
                hit_point = current_pos;
                break;
            }

            // 计算下一步的步长（使用SDF值进行自适应步长）
            double step = std::min(step_size, std::abs(sdf));

            // 更新位置
            current_pos = current_pos + ray_direction * step;

            // 检查是否超出射线长度
            if ((current_pos - ray_origin).norm() > ray_length) {
                break;
            }
        }

        // 如果找到碰撞点，进行碰撞响应
        if (hit_found) {
            Vector3D n = safeDir(hit_point - origin);
            Vector3D place = origin + radius * n;

            // 应用摩擦的碰撞响应
            pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
            return;
        }

        // 如果没有找到碰撞点，检查终点是否在球体内部（兜底检测）
        Vector3D end_to_center = pm.position - origin;
        if (end_to_center.norm2() <= radius2) {
            Vector3D n = safeDir(end_to_center);
            Vector3D place = origin + radius * n;
            pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
        }
    }

    // ---- SDF 离散投影（DCD，便宜且稳）----
    {
        Vector3D to = pm.position - origin;
        double r_len = to.norm();
        double sd = r_len - radius; // 一致半径
        if (sd <= 0.0) {
            Vector3D n = safeDir(to);
            Vector3D tgt = origin + radius * n;
            pm.position = pm.last_position + (1.0 - mu) * (tgt - pm.last_position);
        }
    }
}

void Sphere::render(GLShader& shader) {
    // 渲染半径略减，避免扁平三角形视觉穿插
    m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}