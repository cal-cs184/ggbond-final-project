#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass& pm) {
    // TODO (Part 3): Handle collisions with spheres.

    // Calculate the distance from point mass to sphere center
    Vector3D to_point = pm.position - origin;
    double distance_squared = to_point.norm2();

    // Check if point mass is inside or on the sphere surface
    if (distance_squared <= radius2) {
        // Point mass is inside the sphere, need to move it to surface

        // Calculate the tangent point on sphere surface
        Vector3D direction = to_point.unit();
        Vector3D tangent_point = origin + radius * direction;

        // Calculate correction vector from last_position to tangent point
        Vector3D correction = tangent_point - pm.last_position;

        // Apply correction scaled by friction
        pm.position = pm.last_position + (1.0 - friction) * correction;
    }
}

void Sphere::render(GLShader& shader) {
    // We decrease the radius here so flat triangles don't behave strangely
    // and intersect with the sphere when rendered
    m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
