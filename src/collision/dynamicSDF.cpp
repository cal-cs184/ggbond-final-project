#include <algorithm>
#include <cmath>
#include <nanogui/nanogui.h>

#include "dynamicSDF.h"

using namespace CGL;

// 渲染动态SDF对象：显示两个组成部分的近似可视化
void DynamicSDFObject::render(GLShader& shader) {
    if (!m_enabled) {
        return; // 如果未启用，不渲染
    }

    // // 方法1：Mesh近似渲染（简单但不够精确）
    // // 计算当前时间下两个对象的位置
    // Vector3D sphere_center = Vector3D(sin(m_time * 0.5) * 2.0, 0.0, 0.0);
    // Vector3D box_center = Vector3D(cos(m_time * 0.5) * 2.0, 0.0, 0.0);

    // // 渲染球体部分（使用debug_sphere）
    // nanogui::Color sphere_color(0.8f, 0.3f, 0.3f, 1.0f); // 红色
    // shader.setUniform("u_color", sphere_color);
    // m_debug_sphere.draw_sphere(shader, sphere_center, 1.0);

    // // 渲染盒子部分（用球体近似，因为盒子mesh需要额外实现）
    // // 使用稍小的球体来近似盒子
    // nanogui::Color box_color(0.3f, 0.8f, 0.3f, 1.0f); // 绿色
    // shader.setUniform("u_color", box_color);
    // m_debug_sphere.draw_sphere(shader, box_center, 0.8);

    // 方法2：SDF点云采样渲染（更精确的SDF可视化）
    renderSDFPointCloud(shader);
}

// 新增：SDF点云采样渲染函数
void DynamicSDFObject::renderSDFPointCloud(GLShader& shader) {
    // 精确的SDF点云采样：严格限制在表面附近
    const int num_samples = 8000;          // 增加采样数量
    const double surface_threshold = 0.02; // 更严格的阈值
    const double max_distance = 2.5;       // 限制采样范围

    std::vector<Vector3D> surface_points;
    std::vector<Vector3D> surface_normals;

    // 计算当前时间下两个对象的位置
    Vector3D sphere_center = Vector3D(sin(m_time * 0.5) * 2.0, 0.0, 0.0);
    Vector3D box_center = Vector3D(cos(m_time * 0.5) * 2.0, 0.0, 0.0);

    // 策略1：精确的球体表面采样
    for (int i = 0; i < num_samples / 3; i++) {
        // 在球体表面精确采样
        double theta = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
        double phi = (rand() / (double)RAND_MAX) * M_PI;

        // 从球体表面开始，向内或向外偏移
        double r = 1.0 + (rand() / (double)RAND_MAX - 0.5) * 0.1; // 更小的偏移范围

        Vector3D sample_point;
        sample_point.x = sphere_center.x + r * sin(phi) * cos(theta);
        sample_point.y = sphere_center.y + r * cos(phi);
        sample_point.z = sphere_center.z + r * sin(phi) * sin(theta);

        // 严格检查SDF值
        double sdf_value = sdf_scene(sample_point, m_time);
        if (std::abs(sdf_value) < surface_threshold && sample_point.norm() < max_distance) {
            surface_points.push_back(sample_point);
            surface_normals.push_back(sdf_normal(sample_point, m_time));
        }
    }

    // 策略2：精确的盒子表面采样
    for (int i = 0; i < num_samples / 3; i++) {
        Vector3D box_size(0.8, 0.8, 0.8);

        // 随机选择盒子的一个面
        int face = rand() % 6; // 6个面
        Vector3D sample_point;

        switch (face) {
        case 0: // +X面
            sample_point.x = box_center.x + box_size.x + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 1: // -X面
            sample_point.x = box_center.x - box_size.x + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 2: // +Y面
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + box_size.y + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 3: // -Y面
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y - box_size.y + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 4: // +Z面
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + box_size.z + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            break;
        case 5: // -Z面
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z - box_size.z + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            break;
        }

        // 严格检查SDF值
        double sdf_value = sdf_scene(sample_point, m_time);
        if (std::abs(sdf_value) < surface_threshold && sample_point.norm() < max_distance) {
            surface_points.push_back(sample_point);
            surface_normals.push_back(sdf_normal(sample_point, m_time));
        }
    }

    // 策略3：SDF表面精确投影
    for (int i = 0; i < num_samples / 3; i++) {
        // 在合理范围内随机采样
        Vector3D sample_point;
        sample_point.x = (rand() / (double)RAND_MAX - 0.5) * 4.0;
        sample_point.y = (rand() / (double)RAND_MAX - 0.5) * 2.0;
        sample_point.z = (rand() / (double)RAND_MAX - 0.5) * 2.0;

        if (sample_point.norm() < max_distance) {
            double sdf_value = sdf_scene(sample_point, m_time);

            // 如果接近表面，进行精确投影
            if (std::abs(sdf_value) < surface_threshold * 2.0) {
                // 使用牛顿法精确投影到表面
                Vector3D projected_point = sample_point;
                for (int iter = 0; iter < 5; iter++) {
                    double current_sdf = sdf_scene(projected_point, m_time);
                    if (std::abs(current_sdf) < surface_threshold * 0.1)
                        break;

                    Vector3D normal = sdf_normal(projected_point, m_time);
                    projected_point = projected_point - normal * current_sdf;
                }

                // 最终检查
                double final_sdf = sdf_scene(projected_point, m_time);
                if (std::abs(final_sdf) < surface_threshold && projected_point.norm() < max_distance) {
                    surface_points.push_back(projected_point);
                    surface_normals.push_back(sdf_normal(projected_point, m_time));
                }
            }
        }
    }

    // 渲染SDF表面点云
    if (!surface_points.empty()) {
        // 根据SDF值设置不同颜色
        for (size_t i = 0; i < surface_points.size(); i++) {
            double sdf_val = sdf_scene(surface_points[i], m_time);

            // 根据SDF值设置颜色（负值=内部，正值=外部）
            nanogui::Color sdf_color;
            if (sdf_val < 0) {
                // 内部：蓝色
                sdf_color = nanogui::Color(0.2f, 0.6f, 1.0f, 0.9f);
            } else {
                // 外部：绿色
                sdf_color = nanogui::Color(0.2f, 1.0f, 0.6f, 0.7f);
            }

            shader.setUniform("u_color", sdf_color);

            // 根据SDF值调整球体大小（越接近表面越小）
            double sphere_size = 0.02 * (1.0 - std::abs(sdf_val) / surface_threshold);
            sphere_size = std::max(sphere_size, 0.003); // 更小的最小尺寸
            sphere_size = std::min(sphere_size, 0.04);  // 更小的最大尺寸

            m_debug_sphere.draw_sphere(shader, surface_points[i], sphere_size);
        }
    }
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
