#include <algorithm>
#include <nanogui/nanogui.h>

#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

// Key parameters (start with conservative values, observe FPS before relaxing)
static const int kMaxMarchSteps = 48;          // marching fallback steps
static const double kCCDTriggerRatio = 0.3;    // Only run CCD when L > 0.3R
static const double kNearBandMultiplier = 3.0; // Only run CCD when distance to surface < 3*thickness
static const double kDegenerateLenEps = 1e-12;
static const double kDegenerateDir2Eps = 1e-16;

// Solve line segment-sphere intersection (radius already includes thickness), returns whether hit & hit point x_hit
// (constrained within segment)
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

    // Take the first entry solution (within segment)
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
    // ---- Common parameters ----
    const double mu = clamp(friction, 0.0, 1.0); // (friction)
    // Recommended to use particle spacing: thickness = 0.5 * particle_spacing; here gives a fallback
    // const double thickness = std::max(1e-4, 1e-3 * radius);
    // const double Rexp = radius + thickness; // Extended radius (detection and placement consistent)
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

    // ---- Legacy point-sphere (cheapest) ----
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

    // ---- Prioritize analytical CCD (only triggered when "likely to collide") ----
    if (use_ccd) {
        Vector3D x0 = pm.last_position;
        Vector3D x1 = pm.position; // Predicted position
        Vector3D seg = x1 - x0;
        double L = seg.norm();

        // Trigger condition 1: step size is large enough
        if (L > kCCDTriggerRatio * radius) {
            // Rough check A: starting point too far from surface and displacement not enough to cover => impossible hit
            // (triangle inequality)
            double d0 = (x0 - origin).norm() - radius;
            if (d0 <= L * kNearBandMultiplier) { // Relax a bit to avoid missing boundaries
                // Analytical intersection: if hit, place directly (very fast)
                Vector3D xhit;
                if (segmentSphereHit(x0, x1, origin, radius, xhit)) {
                    Vector3D n = safeDir(xhit - origin);
                    Vector3D place = origin + radius * n;

                    pm.position = pm.position + (1.0 - mu) * (place - pm.position);
                    return;
                }
            }
        }
        // If we reach here: analytical CCD failed or missed ——> fallback: SDF DCD
    }

    // ---- Pure Ray Marching collision detection ----
    if (use_ray_marching) {
        Vector3D ray_origin = pm.last_position;
        Vector3D ray_direction = pm.position - pm.last_position;
        double ray_length = ray_direction.norm();

        if (ray_length < kDegenerateLenEps) {
            return; // Point mass hasn't moved
        }

        ray_direction = ray_direction / ray_length; // Normalize

        // Ray marching parameters
        const int max_steps = kMaxMarchSteps;
        const double step_size = ray_length / max_steps;
        const double epsilon = 1e-6;

        Vector3D current_pos = ray_origin;
        bool hit_found = false;
        Vector3D hit_point;

        // Ray marching loop
        for (int i = 0; i < max_steps; i++) {
            // Calculate distance from current point to sphere center
            Vector3D to_center = current_pos - origin;
            double distance_to_center = to_center.norm();

            // Calculate signed distance field (SDF)
            double sdf = distance_to_center - radius;

            // Check if hit sphere surface
            if (sdf <= epsilon) {
                hit_found = true;
                hit_point = current_pos;
                break;
            }

            // Calculate next step size (use SDF value for adaptive step size)
            double step = std::min(step_size, std::abs(sdf));

            // Update position
            current_pos = current_pos + ray_direction * step;

            // Check if exceeded ray length
            if ((current_pos - ray_origin).norm() > ray_length) {
                break;
            }
        }

        // If collision point found, perform collision response
        if (hit_found) {
            Vector3D n = safeDir(hit_point - origin);
            Vector3D place = origin + radius * n;

            // Apply frictional collision response
            pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
            return;
        }

        // If no collision point found, check if endpoint is inside sphere (fallback detection)
        Vector3D end_to_center = pm.position - origin;
        if (end_to_center.norm2() <= radius2) {
            Vector3D n = safeDir(end_to_center);
            Vector3D place = origin + radius * n;
            pm.position = pm.last_position + (1.0 - mu) * (place - pm.last_position);
        }
    }

    // ---- SDF discrete projection (DCD, cheap and stable) ----
    {
        Vector3D to = pm.position - origin;
        double r_len = to.norm();
        double sd = r_len - radius; // Consistent radius
        if (sd <= 0.0) {
            Vector3D n = safeDir(to);
            Vector3D tgt = origin + radius * n;
            pm.position = pm.last_position + (1.0 - mu) * (tgt - pm.last_position);
        }
    }
}

void Sphere::render(GLShader& shader) {
    // Don't render if hidden when dynamic SDF is enabled
    if (hide_when_dynamic_sdf) {
        return;
    }

    // Slightly reduce rendering radius to avoid visual intersection of flat triangles
    m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}