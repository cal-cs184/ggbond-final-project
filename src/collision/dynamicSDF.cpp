#include <algorithm>
#include <cmath>
#include <nanogui/nanogui.h>

#include "dynamicSDF.h"

using namespace CGL;

// Render dynamic SDF object: display approximate visualization of two components
void DynamicSDFObject::render(GLShader& shader) {
    if (!m_enabled) {
        return; // If not enabled, don't render
    }

    // // Method 1: Mesh approximation rendering (simple but not precise enough)
    // // Calculate positions of two objects at current time
    // Vector3D sphere_center = Vector3D(sin(m_time * 0.5) * 2.0, 0.0, 0.0);
    // Vector3D box_center = Vector3D(cos(m_time * 0.5) * 2.0, 0.0, 0.0);

    // // Render sphere part (using debug_sphere)
    // nanogui::Color sphere_color(0.8f, 0.3f, 0.3f, 1.0f); // Red
    // shader.setUniform("u_color", sphere_color);
    // m_debug_sphere.draw_sphere(shader, sphere_center, 1.0);

    // // Render box part (approximated with sphere, since box mesh requires additional implementation)
    // // Use slightly smaller sphere to approximate box
    // nanogui::Color box_color(0.3f, 0.8f, 0.3f, 1.0f); // Green
    // shader.setUniform("u_color", box_color);
    // m_debug_sphere.draw_sphere(shader, box_center, 0.8);

    // Method 2: SDF point cloud sampling rendering (more precise SDF visualization)
    renderSDFPointCloud(shader);
}

// New: SDF point cloud sampling rendering function
void DynamicSDFObject::renderSDFPointCloud(GLShader& shader) {
    // Precise SDF point cloud sampling: strictly limited to near surface
    const int num_samples = 8000;          // Increase sampling count
    const double surface_threshold = 0.02; // Stricter threshold
    const double max_distance = 2.5;       // Limit sampling range

    std::vector<Vector3D> surface_points;
    std::vector<Vector3D> surface_normals;

    // Calculate positions of two objects at current time
    Vector3D sphere_center = Vector3D(sin(m_time * 0.5) * 2.0, 0.0, 0.0);
    Vector3D box_center = Vector3D(cos(m_time * 0.5) * 2.0, 0.0, 0.0);

    // Strategy 1: Precise sphere surface sampling
    for (int i = 0; i < num_samples / 3; i++) {
        // Precise sampling on sphere surface
        double theta = (rand() / (double)RAND_MAX) * 2.0 * M_PI;
        double phi = (rand() / (double)RAND_MAX) * M_PI;

        // Start from sphere surface, offset inward or outward
        double r = 1.0 + (rand() / (double)RAND_MAX - 0.5) * 0.1; // Smaller offset range

        Vector3D sample_point;
        sample_point.x = sphere_center.x + r * sin(phi) * cos(theta);
        sample_point.y = sphere_center.y + r * cos(phi);
        sample_point.z = sphere_center.z + r * sin(phi) * sin(theta);

        // Strictly check SDF value
        double sdf_value = sdf_scene(sample_point, m_time);
        if (std::abs(sdf_value) < surface_threshold && sample_point.norm() < max_distance) {
            surface_points.push_back(sample_point);
            surface_normals.push_back(sdf_normal(sample_point, m_time));
        }
    }

    // Strategy 2: Precise box surface sampling
    for (int i = 0; i < num_samples / 3; i++) {
        Vector3D box_size(0.8, 0.8, 0.8);

        // Randomly select one face of the box
        int face = rand() % 6; // 6 faces
        Vector3D sample_point;

        switch (face) {
        case 0: // +X face
            sample_point.x = box_center.x + box_size.x + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 1: // -X face
            sample_point.x = box_center.x - box_size.x + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 2: // +Y face
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + box_size.y + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 3: // -Y face
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y - box_size.y + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            sample_point.z = box_center.z + (rand() / (double)RAND_MAX - 0.5) * box_size.z * 2.0;
            break;
        case 4: // +Z face
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z + box_size.z + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            break;
        case 5: // -Z face
            sample_point.x = box_center.x + (rand() / (double)RAND_MAX - 0.5) * box_size.x * 2.0;
            sample_point.y = box_center.y + (rand() / (double)RAND_MAX - 0.5) * box_size.y * 2.0;
            sample_point.z = box_center.z - box_size.z + (rand() / (double)RAND_MAX - 0.5) * 0.1;
            break;
        }

        // Strictly check SDF value
        double sdf_value = sdf_scene(sample_point, m_time);
        if (std::abs(sdf_value) < surface_threshold && sample_point.norm() < max_distance) {
            surface_points.push_back(sample_point);
            surface_normals.push_back(sdf_normal(sample_point, m_time));
        }
    }

    // Strategy 3: SDF surface precise projection
    for (int i = 0; i < num_samples / 3; i++) {
        // Random sampling within reasonable range
        Vector3D sample_point;
        sample_point.x = (rand() / (double)RAND_MAX - 0.5) * 4.0;
        sample_point.y = (rand() / (double)RAND_MAX - 0.5) * 2.0;
        sample_point.z = (rand() / (double)RAND_MAX - 0.5) * 2.0;

        if (sample_point.norm() < max_distance) {
            double sdf_value = sdf_scene(sample_point, m_time);

            // If close to surface, perform precise projection
            if (std::abs(sdf_value) < surface_threshold * 2.0) {
                // Use Newton's method to precisely project to surface
                Vector3D projected_point = sample_point;
                for (int iter = 0; iter < 5; iter++) {
                    double current_sdf = sdf_scene(projected_point, m_time);
                    if (std::abs(current_sdf) < surface_threshold * 0.1)
                        break;

                    Vector3D normal = sdf_normal(projected_point, m_time);
                    projected_point = projected_point - normal * current_sdf;
                }

                // Final check
                double final_sdf = sdf_scene(projected_point, m_time);
                if (std::abs(final_sdf) < surface_threshold && projected_point.norm() < max_distance) {
                    surface_points.push_back(projected_point);
                    surface_normals.push_back(sdf_normal(projected_point, m_time));
                }
            }
        }
    }

    // Render SDF surface point cloud
    if (!surface_points.empty()) {
        // Set different colors based on SDF value
        for (size_t i = 0; i < surface_points.size(); i++) {
            double sdf_val = sdf_scene(surface_points[i], m_time);

            // Set color based on SDF value (negative = inside, positive = outside)
            nanogui::Color sdf_color;
            if (sdf_val < 0) {
                // Inside: blue
                sdf_color = nanogui::Color(0.2f, 0.6f, 1.0f, 0.9f);
            } else {
                // Outside: green
                sdf_color = nanogui::Color(0.2f, 1.0f, 0.6f, 0.7f);
            }

            shader.setUniform("u_color", sdf_color);

            // Adjust sphere size based on SDF value (smaller when closer to surface)
            double sphere_size = 0.02 * (1.0 - std::abs(sdf_val) / surface_threshold);
            sphere_size = std::max(sphere_size, 0.003); // Smaller minimum size
            sphere_size = std::min(sphere_size, 0.04);  // Smaller maximum size

            m_debug_sphere.draw_sphere(shader, surface_points[i], sphere_size);
        }
    }
}

// Scene SDF: union of two primitives moving over time
double DynamicSDFObject::sdf_scene(const Vector3D& p, double t) const {
    // Dynamic sphere
    Vector3D sphere_center = Vector3D(sin(t * 0.5) * 2.0, 0.0, 0.0);
    double d_sphere = sdf_sphere(p - sphere_center, 1.0);

    // Dynamic box
    Vector3D box_center = Vector3D(cos(t * 0.5) * 2.0, 0.0, 0.0);
    double d_box = sdf_box(p - box_center, Vector3D(0.8, 0.8, 0.8));

    return std::min(d_sphere, d_box); // Union
}

// Numerical normal
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

    // Line segment CCD + SDF marching (conservative stepping)
    Vector3D x0 = pm.last_position;
    Vector3D x1 = pm.position;
    Vector3D seg = x1 - x0;
    double L = seg.norm();
    if (L < 1e-12)
        return;
    Vector3D dir = seg / L;

    // Marching parameters
    const int max_steps = 64;
    const double eps = 1e-5; // Hit threshold

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
        // Fallback: if endpoint is inside, also count as hit
        if (sdf_scene(x1, m_time) <= 0.0) {
            hit = true;
            hit_p = x1;
        }
    }

    if (hit) {
        // Project to surface
        double d = sdf_scene(hit_p, m_time);
        Vector3D n = sdf_normal(hit_p, m_time);
        Vector3D place = hit_p - d * n; // Point closer to surface

        // Friction response (maintain consistent feel with existing implementation)
        pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
    }
}
