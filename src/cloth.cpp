#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points, int num_height_points, float thickness) {
    this->width = width;
    this->height = height;
    this->num_width_points = num_width_points;
    this->num_height_points = num_height_points;
    this->thickness = thickness;

    buildGrid();
    buildClothMesh();
}

Cloth::~Cloth() {
    point_masses.clear();
    springs.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void Cloth::buildGrid() {
    // Build a grid of masses and springs.
    double width_step = width / (num_width_points - 1);
    double height_step = height / (num_height_points - 1);

    // Create point masses
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            // Calculate position
            double x = i * width_step;
            double y = 0;
            double z = j * height_step;
            Vector3D position = Vector3D(x, y, z);

            // Check if point mass is pinned
            bool isPinned = false;
            for (auto pin : pinned) {
                if (pin.size() == 2 && pin[0] == i && pin[1] == j) {
                    isPinned = true;
                    break;
                }
            }

            // Create point mass and add to vector
            PointMass pm = PointMass(position, isPinned);
            point_masses.push_back(pm);
        }
    }

    // Create springs
    // Structural springs (adjacent points)
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            // Current point mass index
            int idx = j * num_width_points + i;

            // Horizontal structural springs
            if (i < num_width_points - 1) {
                int right_idx = j * num_width_points + (i + 1);
                Spring spring = Spring(&point_masses[idx], &point_masses[right_idx], STRUCTURAL);
                springs.push_back(spring);
            }

            // Vertical structural springs
            if (j < num_height_points - 1) {
                int down_idx = (j + 1) * num_width_points + i;
                Spring spring = Spring(&point_masses[idx], &point_masses[down_idx], STRUCTURAL);
                springs.push_back(spring);
            }

            // Shearing springs (diagonal connections)
            if (i < num_width_points - 1 && j < num_height_points - 1) {
                // Diagonal down-right
                int diag_idx = (j + 1) * num_width_points + (i + 1);
                Spring spring = Spring(&point_masses[idx], &point_masses[diag_idx], SHEARING);
                springs.push_back(spring);

                // Diagonal down-left (from right neighbor to down neighbor)
                int right_idx = j * num_width_points + (i + 1);
                int down_idx = (j + 1) * num_width_points + i;
                Spring spring2 = Spring(&point_masses[right_idx], &point_masses[down_idx], SHEARING);
                springs.push_back(spring2);
            }

            // Bending springs (connections with distance 2)
            // Horizontal bending
            if (i < num_width_points - 2) {
                int bend_idx = j * num_width_points + (i + 2);
                Spring spring = Spring(&point_masses[idx], &point_masses[bend_idx], BENDING);
                springs.push_back(spring);
            }

            // Vertical bending
            if (j < num_height_points - 2) {
                int bend_idx = (j + 2) * num_width_points + i;
                Spring spring = Spring(&point_masses[idx], &point_masses[bend_idx], BENDING);
                springs.push_back(spring);
            }
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters* cp,
                     vector<Vector3D> external_accelerations, vector<CollisionObject*>* collision_objects) {
    double mass = width * height * cp->density / num_width_points / num_height_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    // Compute total force acting on each point mass
    for (auto& pm : point_masses) {
        // Reset forces
        pm.forces = Vector3D(0, 0, 0);

        // Add all external forces
        for (auto& acc : external_accelerations) {
            pm.forces += mass * acc;
        }
    }

    // Add spring forces
    for (auto& spring : springs) {
        // Skip if spring constraint is disabled
        if ((spring.spring_type == STRUCTURAL && !cp->enable_structural_constraints) ||
            (spring.spring_type == SHEARING && !cp->enable_shearing_constraints) ||
            (spring.spring_type == BENDING && !cp->enable_bending_constraints)) {
            continue;
        }

        // Calculate spring force using Hooke's law: F = -k(|x_i - x_j| - l) * (x_i - x_j)/|x_i - x_j|
        Vector3D p1 = spring.pm_a->position;
        Vector3D p2 = spring.pm_b->position;

        Vector3D dir = p1 - p2;
        double length = dir.norm();

        // Avoid division by zero
        if (length < 0.0001)
            continue;

        dir = dir.unit();
        double displacement = length - spring.rest_length;
        Vector3D force = -cp->ks * displacement * dir;

        // Apply force to both points (equal and opposite)
        spring.pm_a->forces += force;
        spring.pm_b->forces += -force;
    }

    // Use Verlet integration to compute new point mass positions
    for (auto& pm : point_masses) {
        if (pm.pinned)
            continue; // Skip pinned points

        Vector3D acceleration = pm.forces / mass;

        // Calculate new position using Verlet integration
        Vector3D new_position =
            pm.position + (1 - cp->damping) * (pm.position - pm.last_position) + acceleration * delta_t * delta_t;

        // Update positions
        pm.last_position = pm.position;
        pm.position = new_position;
    }

    // Handle self-collisions
    build_spatial_map();
    for (auto& pm : point_masses) {
        self_collide(pm, simulation_steps);
    }

    // Handle collisions with other primitives
    for (auto& pm : point_masses) {
        // Skip pinned points
        if (pm.pinned)
            continue;

        // Check collision with each collision object
        for (auto& co : *collision_objects) {
            co->collide(pm);
        }
    }

    // Constrain spring lengths (Provot's constraint)
    // Prevent springs from stretching more than 10% per timestep
    for (auto& spring : springs) {
        Vector3D p1 = spring.pm_a->position;
        Vector3D p2 = spring.pm_b->position;

        Vector3D dir = p1 - p2;
        double length = dir.norm();

        // Check if spring is stretched more than 10%
        double max_length = spring.rest_length * 1.1;

        if (length > max_length) {
            // Calculate correction
            dir = dir.unit();
            double correction = length - max_length;

            // Apply correction based on whether points are pinned
            if (spring.pm_a->pinned && spring.pm_b->pinned) {
                // Both pinned, do nothing
            } else if (spring.pm_a->pinned) {
                // Only pm_a is pinned
                spring.pm_b->position += correction * dir;
            } else if (spring.pm_b->pinned) {
                // Only pm_b is pinned
                spring.pm_a->position -= correction * dir;
            } else {
                // Neither is pinned, move both
                spring.pm_a->position -= 0.5 * correction * dir;
                spring.pm_b->position += 0.5 * correction * dir;
            }
        }
    }
}

void Cloth::build_spatial_map() {
    for (const auto& entry : map) {
        delete (entry.second);
    }
    map.clear();

    // Build a spatial map out of all of the point masses
    for (auto& pm : point_masses) {
        float hash = hash_position(pm.position);

        if (map.find(hash) == map.end()) {
            // Create a new vector for this hash
            map[hash] = new vector<PointMass*>();
        }

        map[hash]->push_back(&pm);
    }
}

void Cloth::self_collide(PointMass& pm, double simulation_steps) {
    // Handle self-collision for a given point mass
    float hash = hash_position(pm.position);

    // Skip if pinned
    if (pm.pinned)
        return;

    // Check for nearby point masses in the same hash bucket
    Vector3D correction = Vector3D(0, 0, 0);
    int correction_count = 0;

    // Check point masses in the current bucket
    if (map.find(hash) != map.end()) {
        vector<PointMass*>* nearby_points = map[hash];

        for (PointMass* nearby : *nearby_points) {
            // Skip self-comparison
            if (&pm == nearby)
                continue;

            // Calculate distance between points
            Vector3D dir = pm.position - nearby->position;
            double dist = dir.norm();

            // Check if points are too close (less than 2x thickness)
            if (dist < 2 * thickness) {
                // Normalize direction
                if (dist < 0.0001)
                    continue; // Avoid division by zero

                Vector3D correction_dir = dir.unit();

                // Calculate correction amount
                double correction_amount = 2 * thickness - dist;

                // Add to total correction
                correction += correction_dir * correction_amount;
                correction_count++;
            }
        }
    }

    // Apply average correction
    if (correction_count > 0) {
        pm.position += correction / correction_count / simulation_steps;
    }
}

float Cloth::hash_position(Vector3D pos) {
    // Hash a 3D position into a unique float identifier
    // Use a grid-based spatial hashing approach

    // Define grid cell size (based on cloth thickness)
    float w = 3 * thickness;

    // Calculate grid cell indices
    int i = floor(pos.x / w);
    int j = floor(pos.y / w);
    int k = floor(pos.z / w);

    // Use a simple hash function to convert 3D indices to a float
    // This is a simple spatial hash function that should work for our purposes
    float hash = i * 73856093 + j * 19349663 + k * 83492791;

    // Convert to a float in a reasonable range
    return hash;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    PointMass* pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->last_position = pm->start_position;
        pm++;
    }
}

void Cloth::buildClothMesh() {
    if (point_masses.size() == 0)
        return;

    ClothMesh* clothMesh = new ClothMesh();
    vector<Triangle*> triangles;

    // Create vector of triangles
    for (int y = 0; y < num_height_points - 1; y++) {
        for (int x = 0; x < num_width_points - 1; x++) {
            PointMass* pm = &point_masses[y * num_width_points + x];
            // Get neighboring point masses:
            /*                      *
             * pm_A -------- pm_B   *
             *             /        *
             *  |         /   |     *
             *  |        /    |     *
             *  |       /     |     *
             *  |      /      |     *
             *  |     /       |     *
             *  |    /        |     *
             *      /               *
             * pm_C -------- pm_D   *
             *                      *
             */

            float u_min = x;
            u_min /= num_width_points - 1;
            float u_max = x + 1;
            u_max /= num_width_points - 1;
            float v_min = y;
            v_min /= num_height_points - 1;
            float v_max = y + 1;
            v_max /= num_height_points - 1;

            PointMass* pm_A = pm;
            PointMass* pm_B = pm + 1;
            PointMass* pm_C = pm + num_width_points;
            PointMass* pm_D = pm + num_width_points + 1;

            Vector3D uv_A = Vector3D(u_min, v_min, 0);
            Vector3D uv_B = Vector3D(u_max, v_min, 0);
            Vector3D uv_C = Vector3D(u_min, v_max, 0);
            Vector3D uv_D = Vector3D(u_max, v_max, 0);

            // Both triangles defined by vertices in counter-clockwise orientation
            triangles.push_back(new Triangle(pm_A, pm_C, pm_B, uv_A, uv_C, uv_B));
            triangles.push_back(new Triangle(pm_B, pm_C, pm_D, uv_B, uv_C, uv_D));
        }
    }

    // For each triangle in row-order, create 3 edges and 3 internal halfedges
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge* h1 = new Halfedge();
        Halfedge* h2 = new Halfedge();
        Halfedge* h3 = new Halfedge();

        // Allocate new edges on heap
        Edge* e1 = new Edge();
        Edge* e2 = new Edge();
        Edge* e3 = new Edge();

        // Assign a halfedge pointer to the triangle
        t->halfedge = h1;

        // Assign halfedge pointers to point masses
        t->pm1->halfedge = h1;
        t->pm2->halfedge = h2;
        t->pm3->halfedge = h3;

        // Update all halfedge pointers
        h1->edge = e1;
        h1->next = h2;
        h1->pm = t->pm1;
        h1->triangle = t;

        h2->edge = e2;
        h2->next = h3;
        h2->pm = t->pm2;
        h2->triangle = t;

        h3->edge = e3;
        h3->next = h1;
        h3->pm = t->pm3;
        h3->triangle = t;
    }

    // Go back through the cloth mesh and link triangles together using halfedge
    // twin pointers

    // Convenient variables for math
    int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points - 1) * 2;

    bool topLeft = true;
    for (int i = 0; i < triangles.size(); i++) {
        Triangle* t = triangles[i];

        if (topLeft) {
            // Get left triangle, if it exists
            if (i % num_width_tris != 0) { // Not a left-most triangle
                Triangle* temp = triangles[i - 1];
                t->pm1->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm1->halfedge->twin = nullptr;
            }

            // Get triangle above, if it exists
            if (i >= num_width_tris) { // Not a top-most triangle
                Triangle* temp = triangles[i - num_width_tris + 1];
                t->pm3->halfedge->twin = temp->pm2->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle to bottom right; guaranteed to exist
            Triangle* temp = triangles[i + 1];
            t->pm2->halfedge->twin = temp->pm1->halfedge;
        } else {
            // Get right triangle, if it exists
            if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
                Triangle* temp = triangles[i + 1];
                t->pm3->halfedge->twin = temp->pm1->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle below, if it exists
            if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
                Triangle* temp = triangles[i + num_width_tris - 1];
                t->pm2->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm2->halfedge->twin = nullptr;
            }

            // Get triangle to top left; guaranteed to exist
            Triangle* temp = triangles[i - 1];
            t->pm1->halfedge->twin = temp->pm2->halfedge;
        }

        topLeft = !topLeft;
    }

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}
