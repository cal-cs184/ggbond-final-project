#include "iostream"
#include <cmath>
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass& pm) {
    // TODO (Part 3): Handle collisions with planes.

    // Calculate signed distance from point to plane
    Vector3D to_point = pm.position - point;
    double current_distance = dot(to_point, normal);

    // Calculate signed distance from last position to plane
    Vector3D to_last_point = pm.last_position - point;
    double last_distance = dot(to_last_point, normal);

    // If the point is penetrating the plane or crossed the plane between the last and current
    if (current_distance <= 0 || (current_distance > 0 && last_distance < 0) ||
        (current_distance < 0 && last_distance > 0)) {
        // Compute the intersection point (tangent point) of the segment with the plane
        Vector3D direction = pm.position - pm.last_position;
        double denom = dot(direction, normal);
        if (fabs(denom) < 1e-10) {
            return; // Direction parallel to plane, nothing to do
        }
        double t = -last_distance / denom; // Could be >1 or <0 if already inside but still works
        if (t < 0.0)
            t = 0.0; // Clamp to last_position for penetrating cases
        if (t > 1.0)
            t = 1.0;
        Vector3D tangent_point = pm.last_position + t * direction;

        // Move the point slightly above the plane along the normal
        Vector3D corrected_point = tangent_point + SURFACE_OFFSET * normal;

        // Compute correction vector
        Vector3D correction = corrected_point - pm.last_position;

        // Apply correction with friction factor
        pm.position = pm.last_position + (1.0 - friction) * correction;
    }
}

void Plane::render(GLShader& shader) {
    nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

    Vector3f sPoint(point.x, point.y, point.z);
    Vector3f sNormal(normal.x, normal.y, normal.z);
    Vector3f sParallel(normal.y - normal.z, normal.z - normal.x, normal.x - normal.y);
    sParallel.normalize();
    Vector3f sCross = sNormal.cross(sParallel);

    MatrixXf positions(3, 4);
    MatrixXf normals(3, 4);

    positions.col(0) << sPoint + 2 * (sCross + sParallel);
    positions.col(1) << sPoint + 2 * (sCross - sParallel);
    positions.col(2) << sPoint + 2 * (-sCross + sParallel);
    positions.col(3) << sPoint + 2 * (-sCross - sParallel);

    normals.col(0) << sNormal;
    normals.col(1) << sNormal;
    normals.col(2) << sNormal;
    normals.col(3) << sNormal;

    if (shader.uniform("u_color", false) != -1) {
        shader.setUniform("u_color", color);
    }
    shader.uploadAttrib("in_position", positions);
    if (shader.attrib("in_normal", false) != -1) {
        shader.uploadAttrib("in_normal", normals);
    }

    shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
}
