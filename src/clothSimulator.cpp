#include <cmath>
#include <glad/glad.h>

#include <CGL/vector3D.h>
#include <nanogui/nanogui.h>

#include "clothSimulator.h"

#include "camera.h"
#include "cloth.h"
#include "collision/dynamicSDF.h"
#include "collision/sphere.h"
#include "misc/camera_info.h"
#include "misc/file_utils.h"
// Needed to generate stb_image binaries. Should only define in exactly one source file importing stb_image.h.
#define STB_IMAGE_IMPLEMENTATION
#include "misc/stb_image.h"

using namespace nanogui;
using namespace std;

Vector3D load_texture(int frame_idx, GLuint handle, const char* where) {
    Vector3D size_retval;

    if (strlen(where) == 0)
        return size_retval;

    glActiveTexture(GL_TEXTURE0 + frame_idx);
    glBindTexture(GL_TEXTURE_2D, handle);

    int img_x, img_y, img_n;
    unsigned char* img_data = stbi_load(where, &img_x, &img_y, &img_n, 3);
    size_retval.x = img_x;
    size_retval.y = img_y;
    size_retval.z = img_n;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img_x, img_y, 0, GL_RGB, GL_UNSIGNED_BYTE, img_data);
    stbi_image_free(img_data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    return size_retval;
}

void load_cubemap(int frame_idx, GLuint handle, const std::vector<std::string>& file_locs) {
    glActiveTexture(GL_TEXTURE0 + frame_idx);
    glBindTexture(GL_TEXTURE_CUBE_MAP, handle);
    for (int side_idx = 0; side_idx < 6; ++side_idx) {

        int img_x, img_y, img_n;
        unsigned char* img_data = stbi_load(file_locs[side_idx].c_str(), &img_x, &img_y, &img_n, 3);
        glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + side_idx, 0, GL_RGB, img_x, img_y, 0, GL_RGB, GL_UNSIGNED_BYTE,
                     img_data);
        stbi_image_free(img_data);
        std::cout << "Side " << side_idx << " has dimensions " << img_x << ", " << img_y << std::endl;

        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    }
}

void ClothSimulator::load_textures() {
    glGenTextures(1, &m_gl_texture_1);
    glGenTextures(1, &m_gl_texture_2);
    glGenTextures(1, &m_gl_texture_3);
    glGenTextures(1, &m_gl_texture_4);
    glGenTextures(1, &m_gl_cubemap_tex);
    glGenTextures(1, &m_cloth_particles_tex);

    m_gl_texture_1_size = load_texture(1, m_gl_texture_1, (m_project_root + "/textures/texture_1.png").c_str());
    m_gl_texture_2_size = load_texture(2, m_gl_texture_2, (m_project_root + "/textures/texture_2.png").c_str());
    m_gl_texture_3_size = load_texture(3, m_gl_texture_3, (m_project_root + "/textures/texture_3.png").c_str());
    m_gl_texture_4_size = load_texture(4, m_gl_texture_4, (m_project_root + "/textures/texture_4.png").c_str());

    std::cout << "Texture 1 loaded with size: " << m_gl_texture_1_size << std::endl;
    std::cout << "Texture 2 loaded with size: " << m_gl_texture_2_size << std::endl;
    std::cout << "Texture 3 loaded with size: " << m_gl_texture_3_size << std::endl;
    std::cout << "Texture 4 loaded with size: " << m_gl_texture_4_size << std::endl;

    // Initialize cloth particle texture
    // We'll set the size based on the cloth dimensions when the cloth is loaded
    m_cloth_tex_width = 32;  // Default size, will be updated later
    m_cloth_tex_height = 32; // Default size, will be updated later

    // Initialize the cloth particles texture with empty data
    glActiveTexture(GL_TEXTURE6); // Use texture unit 6 for cloth particles
    glBindTexture(GL_TEXTURE_2D, m_cloth_particles_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Allocate empty texture (will be filled later)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, m_cloth_tex_width, m_cloth_tex_height, 0, GL_RGB, GL_FLOAT, nullptr);

    std::cout << "Initialized cloth particles texture with size: " << m_cloth_tex_width << "x" << m_cloth_tex_height
              << std::endl;

    std::vector<std::string> cubemap_fnames = {
        m_project_root + "/textures/cube/posx.jpg", m_project_root + "/textures/cube/negx.jpg",
        m_project_root + "/textures/cube/posy.jpg", m_project_root + "/textures/cube/negy.jpg",
        m_project_root + "/textures/cube/posz.jpg", m_project_root + "/textures/cube/negz.jpg"};

    load_cubemap(5, m_gl_cubemap_tex, cubemap_fnames);
    std::cout << "Loaded cubemap texture" << std::endl;
}

void ClothSimulator::load_shaders() {
    std::set<std::string> shader_folder_contents;
    bool success = FileUtils::list_files_in_directory(m_project_root + "/shaders", shader_folder_contents);
    if (!success) {
        std::cout << "Error: Could not find the shaders folder!" << std::endl;
    }

    std::string std_vert_shader = m_project_root + "/shaders/Default.vert";

    for (const std::string& shader_fname : shader_folder_contents) {
        std::string file_extension;
        std::string shader_name;

        FileUtils::split_filename(shader_fname, shader_name, file_extension);

        if (file_extension != "frag") {
            std::cout << "Skipping non-shader file: " << shader_fname << std::endl;
            continue;
        }

        std::cout << "Found shader file: " << shader_fname << std::endl;

        // Check if there is a proper .vert shader or not for it
        std::string vert_shader = std_vert_shader;
        std::string associated_vert_shader_path = m_project_root + "/shaders/" + shader_name + ".vert";
        if (FileUtils::file_exists(associated_vert_shader_path)) {
            vert_shader = associated_vert_shader_path;
        }

        std::shared_ptr<GLShader> nanogui_shader = make_shared<GLShader>();
        nanogui_shader->initFromFiles(shader_name, vert_shader, m_project_root + "/shaders/" + shader_fname);

        // Special filenames are treated a bit differently
        ShaderTypeHint hint;
        if (shader_name == "Wireframe") {
            hint = ShaderTypeHint::WIREFRAME;
            std::cout << "Type: Wireframe" << std::endl;
        } else if (shader_name == "Normal") {
            hint = ShaderTypeHint::NORMALS;
            std::cout << "Type: Normal" << std::endl;
        } else if (shader_name == "NewRayMarching") {
            hint = ShaderTypeHint::VOLUME_RENDERING;
            std::cout << "Type: Volume Rendering" << std::endl;
        } else if (shader_name == "DynamicRayMarching") {
            hint = ShaderTypeHint::VOLUME_RENDERING;
            std::cout << "Type: Dynamic Ray Marching" << std::endl;
        } else {
            hint = ShaderTypeHint::PHONG;
            std::cout << "Type: Custom" << std::endl;
        }

        UserShader user_shader(shader_name, nanogui_shader, hint);

        shaders.push_back(user_shader);
        shaders_combobox_names.push_back(shader_name);
    }

    // Assuming that it's there, use "Wireframe" by default
    for (size_t i = 0; i < shaders_combobox_names.size(); ++i) {
        if (shaders_combobox_names[i] == "Wireframe") {
            active_shader_idx = i;
            break;
        }
    }
}

ClothSimulator::ClothSimulator(std::string project_root, Screen* screen) : m_project_root(project_root) {
    this->screen = screen;

    this->load_shaders();
    this->load_textures();

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);
}

ClothSimulator::~ClothSimulator() {
    for (auto shader : shaders) {
        shader.nanogui_shader->free();
    }
    glDeleteTextures(1, &m_gl_texture_1);
    glDeleteTextures(1, &m_gl_texture_2);
    glDeleteTextures(1, &m_gl_texture_3);
    glDeleteTextures(1, &m_gl_texture_4);
    glDeleteTextures(1, &m_gl_cubemap_tex);
    glDeleteTextures(1, &m_cloth_particles_tex);

    if (cloth)
        delete cloth;
    if (cp)
        delete cp;
    if (collision_objects)
        delete collision_objects;
}

void ClothSimulator::loadCloth(Cloth* cloth) {
    this->cloth = cloth;

    // Update cloth particles texture size based on cloth dimensions
    m_cloth_tex_width = cloth->num_width_points;
    m_cloth_tex_height = cloth->num_height_points;

    // Resize the cloth particles texture
    glActiveTexture(GL_TEXTURE6);
    glBindTexture(GL_TEXTURE_2D, m_cloth_particles_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, m_cloth_tex_width, m_cloth_tex_height, 0, GL_RGB, GL_FLOAT, nullptr);

    std::cout << "Resized cloth particles texture to: " << m_cloth_tex_width << "x" << m_cloth_tex_height << std::endl;
}

void ClothSimulator::loadClothParameters(ClothParameters* cp) { this->cp = cp; }

void ClothSimulator::loadCollisionObjects(vector<CollisionObject*>* objects) { this->collision_objects = objects; }

/**
 * Initializes the cloth simulation and spawns a new thread to separate
 * rendering from simulation.
 */
void ClothSimulator::init() {

    // Initialize GUI
    screen->setSize(default_window_size);
    initGUI(screen);

    // Initialize camera

    CGL::Collada::CameraInfo camera_info;
    camera_info.hFov = 50;
    camera_info.vFov = 35;
    camera_info.nClip = 0.01;
    camera_info.fClip = 10000;

    // Try to intelligently figure out the camera target

    Vector3D avg_pm_position(0, 0, 0);

    for (auto& pm : cloth->point_masses) {
        avg_pm_position += pm.position / cloth->point_masses.size();
    }

    CGL::Vector3D target(avg_pm_position.x, avg_pm_position.y / 2, avg_pm_position.z);
    CGL::Vector3D c_dir(0., 0., 0.);
    canonical_view_distance = max(cloth->width, cloth->height) * 0.9;
    scroll_rate = canonical_view_distance / 10;

    view_distance = canonical_view_distance * 2;
    min_view_distance = canonical_view_distance / 10.0;
    max_view_distance = canonical_view_distance * 20.0;

    // canonicalCamera is a copy used for view resets

    camera.place(target, acos(c_dir.y), atan2(c_dir.x, c_dir.z), view_distance, min_view_distance, max_view_distance);
    canonicalCamera.place(target, acos(c_dir.y), atan2(c_dir.x, c_dir.z), view_distance, min_view_distance,
                          max_view_distance);

    screen_w = default_window_size(0);
    screen_h = default_window_size(1);

    camera.configure(camera_info, screen_w, screen_h);
    canonicalCamera.configure(camera_info, screen_w, screen_h);
}

bool ClothSimulator::isAlive() { return is_alive; }

void ClothSimulator::drawContents() {
    glEnable(GL_DEPTH_TEST);

    if (!is_paused) {
        vector<Vector3D> external_accelerations = {gravity};

        for (int i = 0; i < simulation_steps; i++) {
            cloth->simulate(frames_per_sec, simulation_steps, cp, external_accelerations, collision_objects);
        }

        // Update simulation time for animations
        simulation_time += 1.0 / frames_per_sec;
    }

    // Update cloth particles texture with current positions
    if (cloth && cloth->point_masses.size() > 0) {
        // Create a buffer to hold particle positions
        std::vector<float> particle_positions(m_cloth_tex_width * m_cloth_tex_height * 3);

        // Fill the buffer with particle positions
        for (int j = 0; j < cloth->num_height_points; j++) {
            for (int i = 0; i < cloth->num_width_points; i++) {
                int idx = j * cloth->num_width_points + i;
                int tex_idx = (j * m_cloth_tex_width + i) * 3;

                if (idx < cloth->point_masses.size()) {
                    Vector3D pos = cloth->point_masses[idx].position;
                    particle_positions[tex_idx] = pos.x;
                    particle_positions[tex_idx + 1] = pos.y;
                    particle_positions[tex_idx + 2] = pos.z;
                }
            }
        }

        // Update the texture with new particle positions
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, m_cloth_particles_tex);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_cloth_tex_width, m_cloth_tex_height, GL_RGB, GL_FLOAT,
                        particle_positions.data());
    }

    // Bind the active shader

    const UserShader& active_shader = shaders[active_shader_idx];

    GLShader& shader = *active_shader.nanogui_shader;
    shader.bind();

    // Prepare the camera projection matrix

    Matrix4f model;
    model.setIdentity();

    Matrix4f view = getViewMatrix();
    Matrix4f projection = getProjectionMatrix();

    Matrix4f viewProjection = projection * view;

    shader.setUniform("u_model", model);
    shader.setUniform("u_view_projection", viewProjection);

    // Get camera position for all shaders
    Vector3D cam_pos = camera.position();

    switch (active_shader.type_hint) {
    case WIREFRAME:
        shader.setUniform("u_color", color, false);
        drawWireframe(shader);
        break;
    case NORMALS:
        drawNormals(shader);
        break;
    case PHONG:
        // Others
        shader.setUniform("u_color", color, false);
        shader.setUniform("u_cam_pos", Vector3f(cam_pos.x, cam_pos.y, cam_pos.z), false);
        shader.setUniform("u_light_pos", Vector3f(0.5, 2, 2), false);
        shader.setUniform("u_light_intensity", Vector3f(3, 3, 3), false);
        shader.setUniform("u_texture_1_size", Vector2f(m_gl_texture_1_size.x, m_gl_texture_1_size.y), false);
        shader.setUniform("u_texture_2_size", Vector2f(m_gl_texture_2_size.x, m_gl_texture_2_size.y), false);
        shader.setUniform("u_texture_3_size", Vector2f(m_gl_texture_3_size.x, m_gl_texture_3_size.y), false);
        shader.setUniform("u_texture_4_size", Vector2f(m_gl_texture_4_size.x, m_gl_texture_4_size.y), false);
        // Textures
        shader.setUniform("u_texture_1", 1, false);
        shader.setUniform("u_texture_2", 2, false);
        shader.setUniform("u_texture_3", 3, false);
        shader.setUniform("u_texture_4", 4, false);

        shader.setUniform("u_normal_scaling", m_normal_scaling, false);
        shader.setUniform("u_height_scaling", m_height_scaling, false);

        shader.setUniform("u_texture_cubemap", 5, false);

        // Pass cloth particles texture and dimensions
        shader.setUniform("u_cloth_particles_tex", 6, false);
        shader.setUniform("u_cloth_tex_width", m_cloth_tex_width, false);
        shader.setUniform("u_cloth_tex_height", m_cloth_tex_height, false);
        drawPhong(shader);
        break;
    case VOLUME_RENDERING:
        // Set uniforms for volume rendering
        shader.setUniform("u_model", model, false);
        shader.setUniform("u_view_projection", viewProjection, false);
        // Pass inverse view-projection for correct ray directions in volume rendering
        Matrix4f invViewProjection = viewProjection.inverse();
        shader.setUniform("u_inv_view_projection", invViewProjection, false);
        shader.setUniform("u_cam_pos", Vector3f(cam_pos.x, cam_pos.y, cam_pos.z), false);
        shader.setUniform("u_light_pos", Vector3f(0.5, 2, 2), false);
        shader.setUniform("u_light_intensity", Vector3f(3, 3, 3), false);
        shader.setUniform("u_time", (float)simulation_time, false);
        shader.setUniform("u_texture_1", 1, false);
        shader.setUniform("u_cloth_particles_tex", 6, false);
        shader.setUniform("u_cloth_tex_width", m_cloth_tex_width, false);
        shader.setUniform("u_cloth_tex_height", m_cloth_tex_height, false);
        drawVolumeRendering(shader);
        // In volume rendering mode, don't render collision objects, only show ray marching scene
        break;
    }

    for (CollisionObject* co : *collision_objects) {
        // Try to pass the SDF toggle to spheres (safe dynamic_cast)
        if (auto* sphere = dynamic_cast<Sphere*>(co)) {
            sphere->set_use_sdf(use_sdf_collision);
            sphere->set_use_ccd(use_ccd_collision);
            sphere->set_use_ray_marching(use_ray_marching_collision);
            sphere->set_hide_when_dynamic_sdf(use_dynamic_sdf_object);
        }
        // If it's a dynamic SDF object, synchronize time and toggle
        if (auto* sdfObj = dynamic_cast<DynamicSDFObject*>(co)) {
            sdfObj->set_time(simulation_time);
            sdfObj->set_enabled(use_dynamic_sdf_object);

            // Update configurable parameters
            sdfObj->set_num_samples(dynamic_sdf_num_samples);
            sdfObj->set_surface_threshold(dynamic_sdf_surface_threshold);
            sdfObj->set_max_distance(dynamic_sdf_max_distance);
            sdfObj->set_min_sphere_size(dynamic_sdf_min_sphere_size);
            sdfObj->set_max_sphere_size(dynamic_sdf_max_sphere_size);
            sdfObj->set_motion_speed(dynamic_sdf_motion_speed);
            sdfObj->set_sphere_radius(dynamic_sdf_sphere_radius);
            sdfObj->set_box_size(dynamic_sdf_box_size);
        }
        co->render(shader);
    }
}

void ClothSimulator::drawWireframe(GLShader& shader) {
    int num_structural_springs =
        2 * cloth->num_width_points * cloth->num_height_points - cloth->num_width_points - cloth->num_height_points;
    int num_shear_springs = 2 * (cloth->num_width_points - 1) * (cloth->num_height_points - 1);
    int num_bending_springs = num_structural_springs - cloth->num_width_points - cloth->num_height_points;

    int num_springs = cp->enable_structural_constraints * num_structural_springs +
                      cp->enable_shearing_constraints * num_shear_springs +
                      cp->enable_bending_constraints * num_bending_springs;

    MatrixXf positions(4, num_springs * 2);
    MatrixXf normals(4, num_springs * 2);

    // Draw springs as lines

    int si = 0;

    for (int i = 0; i < cloth->springs.size(); i++) {
        Spring s = cloth->springs[i];

        if ((s.spring_type == STRUCTURAL && !cp->enable_structural_constraints) ||
            (s.spring_type == SHEARING && !cp->enable_shearing_constraints) ||
            (s.spring_type == BENDING && !cp->enable_bending_constraints)) {
            continue;
        }

        Vector3D pa = s.pm_a->position;
        Vector3D pb = s.pm_b->position;

        Vector3D na = s.pm_a->normal();
        Vector3D nb = s.pm_b->normal();

        positions.col(si) << pa.x, pa.y, pa.z, 1.0;
        positions.col(si + 1) << pb.x, pb.y, pb.z, 1.0;

        normals.col(si) << na.x, na.y, na.z, 0.0;
        normals.col(si + 1) << nb.x, nb.y, nb.z, 0.0;

        si += 2;
    }

    // shader.setUniform("u_color", nanogui::Color(1.0f, 1.0f, 1.0f, 1.0f), false);
    shader.uploadAttrib("in_position", positions, false);
    // Commented out: the wireframe shader does not have this attribute
    // shader.uploadAttrib("in_normal", normals);

    shader.drawArray(GL_LINES, 0, num_springs * 2);
}

void ClothSimulator::drawNormals(GLShader& shader) {
    int num_tris = cloth->clothMesh->triangles.size();

    MatrixXf positions(4, num_tris * 3);
    MatrixXf normals(4, num_tris * 3);

    for (int i = 0; i < num_tris; i++) {
        Triangle* tri = cloth->clothMesh->triangles[i];

        Vector3D p1 = tri->pm1->position;
        Vector3D p2 = tri->pm2->position;
        Vector3D p3 = tri->pm3->position;

        Vector3D n1 = tri->pm1->normal();
        Vector3D n2 = tri->pm2->normal();
        Vector3D n3 = tri->pm3->normal();

        positions.col(i * 3) << p1.x, p1.y, p1.z, 1.0;
        positions.col(i * 3 + 1) << p2.x, p2.y, p2.z, 1.0;
        positions.col(i * 3 + 2) << p3.x, p3.y, p3.z, 1.0;

        normals.col(i * 3) << n1.x, n1.y, n1.z, 0.0;
        normals.col(i * 3 + 1) << n2.x, n2.y, n2.z, 0.0;
        normals.col(i * 3 + 2) << n3.x, n3.y, n3.z, 0.0;
    }

    shader.uploadAttrib("in_position", positions, false);
    shader.uploadAttrib("in_normal", normals, false);

    shader.drawArray(GL_TRIANGLES, 0, num_tris * 3);
}

void ClothSimulator::drawPhong(GLShader& shader) {
    int num_tris = cloth->clothMesh->triangles.size();

    MatrixXf positions(4, num_tris * 3);
    MatrixXf normals(4, num_tris * 3);
    MatrixXf uvs(2, num_tris * 3);
    MatrixXf tangents(4, num_tris * 3);

    for (int i = 0; i < num_tris; i++) {
        Triangle* tri = cloth->clothMesh->triangles[i];

        Vector3D p1 = tri->pm1->position;
        Vector3D p2 = tri->pm2->position;
        Vector3D p3 = tri->pm3->position;

        Vector3D n1 = tri->pm1->normal();
        Vector3D n2 = tri->pm2->normal();
        Vector3D n3 = tri->pm3->normal();

        positions.col(i * 3) << p1.x, p1.y, p1.z, 1.0;
        positions.col(i * 3 + 1) << p2.x, p2.y, p2.z, 1.0;
        positions.col(i * 3 + 2) << p3.x, p3.y, p3.z, 1.0;

        normals.col(i * 3) << n1.x, n1.y, n1.z, 0.0;
        normals.col(i * 3 + 1) << n2.x, n2.y, n2.z, 0.0;
        normals.col(i * 3 + 2) << n3.x, n3.y, n3.z, 0.0;

        uvs.col(i * 3) << tri->uv1.x, tri->uv1.y;
        uvs.col(i * 3 + 1) << tri->uv2.x, tri->uv2.y;
        uvs.col(i * 3 + 2) << tri->uv3.x, tri->uv3.y;

        tangents.col(i * 3) << 1.0, 0.0, 0.0, 1.0;
        tangents.col(i * 3 + 1) << 1.0, 0.0, 0.0, 1.0;
        tangents.col(i * 3 + 2) << 1.0, 0.0, 0.0, 1.0;
    }

    shader.uploadAttrib("in_position", positions, false);
    shader.uploadAttrib("in_normal", normals, false);
    shader.uploadAttrib("in_uv", uvs, false);
    shader.uploadAttrib("in_tangent", tangents, false);

    shader.drawArray(GL_TRIANGLES, 0, num_tris * 3);
}

// ----------------------------------------------------------------------------
// CAMERA CALCULATIONS
//
// OpenGL 3.1 deprecated the fixed pipeline, so we lose a lot of useful OpenGL
// functions that have to be recreated here.
// ----------------------------------------------------------------------------

void ClothSimulator::resetCamera() { camera.copy_placement(canonicalCamera); }

Matrix4f ClothSimulator::getProjectionMatrix() {
    Matrix4f perspective;
    perspective.setZero();

    double cam_near = camera.near_clip();
    double cam_far = camera.far_clip();

    double theta = camera.v_fov() * PI / 360;
    double range = cam_far - cam_near;
    double invtan = 1. / tanf(theta);

    perspective(0, 0) = invtan / camera.aspect_ratio();
    perspective(1, 1) = invtan;
    perspective(2, 2) = -(cam_near + cam_far) / range;
    perspective(3, 2) = -1;
    perspective(2, 3) = -2 * cam_near * cam_far / range;
    perspective(3, 3) = 0;

    return perspective;
}

Matrix4f ClothSimulator::getViewMatrix() {
    Matrix4f lookAt;
    Matrix3f R;

    lookAt.setZero();

    // Convert CGL vectors to Eigen vectors
    // TODO: Find a better way to do this!

    CGL::Vector3D c_pos = camera.position();
    CGL::Vector3D c_udir = camera.up_dir();
    CGL::Vector3D c_target = camera.view_point();

    Vector3f eye(c_pos.x, c_pos.y, c_pos.z);
    Vector3f up(c_udir.x, c_udir.y, c_udir.z);
    Vector3f target(c_target.x, c_target.y, c_target.z);

    R.col(2) = (eye - target).normalized();
    R.col(0) = up.cross(R.col(2)).normalized();
    R.col(1) = R.col(2).cross(R.col(0));

    lookAt.topLeftCorner<3, 3>() = R.transpose();
    lookAt.topRightCorner<3, 1>() = -R.transpose() * eye;
    lookAt(3, 3) = 1.0f;

    return lookAt;
}

// ----------------------------------------------------------------------------
// EVENT HANDLING
// ----------------------------------------------------------------------------

bool ClothSimulator::cursorPosCallbackEvent(double x, double y) {
    if (left_down && !middle_down && !right_down) {
        if (ctrl_down) {
            mouseRightDragged(x, y);
        } else {
            mouseLeftDragged(x, y);
        }
    } else if (!left_down && !middle_down && right_down) {
        mouseRightDragged(x, y);
    } else if (!left_down && !middle_down && !right_down) {
        mouseMoved(x, y);
    }

    mouse_x = x;
    mouse_y = y;

    return true;
}

bool ClothSimulator::mouseButtonCallbackEvent(int button, int action, int modifiers) {
    switch (action) {
    case GLFW_PRESS:
        switch (button) {
        case GLFW_MOUSE_BUTTON_LEFT:
            left_down = true;
            break;
        case GLFW_MOUSE_BUTTON_MIDDLE:
            middle_down = true;
            break;
        case GLFW_MOUSE_BUTTON_RIGHT:
            right_down = true;
            break;
        }
        return true;

    case GLFW_RELEASE:
        switch (button) {
        case GLFW_MOUSE_BUTTON_LEFT:
            left_down = false;
            break;
        case GLFW_MOUSE_BUTTON_MIDDLE:
            middle_down = false;
            break;
        case GLFW_MOUSE_BUTTON_RIGHT:
            right_down = false;
            break;
        }
        return true;
    }

    return false;
}

void ClothSimulator::mouseMoved(double x, double y) { y = screen_h - y; }

void ClothSimulator::mouseLeftDragged(double x, double y) {
    float dx = x - mouse_x;
    float dy = y - mouse_y;

    camera.rotate_by(-dy * (PI / screen_h), -dx * (PI / screen_w));
}

void ClothSimulator::mouseRightDragged(double x, double y) {
    camera.move_by(mouse_x - x, y - mouse_y, canonical_view_distance);
}

bool ClothSimulator::keyCallbackEvent(int key, int scancode, int action, int mods) {
    ctrl_down = (bool)(mods & GLFW_MOD_CONTROL);

    if (action == GLFW_PRESS) {
        switch (key) {
        case GLFW_KEY_ESCAPE:
            is_alive = false;
            break;
        case 'r':
        case 'R':
            cloth->reset();
            break;
        case ' ':
            resetCamera();
            break;
        case 'p':
        case 'P':
            is_paused = !is_paused;
            break;
        case 'n':
        case 'N':
            if (is_paused) {
                is_paused = false;
                drawContents();
                is_paused = true;
            }
            break;
        }
    }

    return true;
}

bool ClothSimulator::dropCallbackEvent(int count, const char** filenames) { return true; }

bool ClothSimulator::scrollCallbackEvent(double x, double y) {
    camera.move_forward(y * scroll_rate);
    return true;
}

bool ClothSimulator::resizeCallbackEvent(int width, int height) {
    screen_w = width;
    screen_h = height;

    camera.set_screen_size(screen_w, screen_h);
    return true;
}

void ClothSimulator::initGUI(Screen* screen) {
    Window* window;

    window = new Window(screen, "Simulation");
    window->setPosition(Vector2i(default_window_size(0) - 245, 15));
    window->setLayout(new GroupLayout(15, 6, 14, 5));

    // Spring types

    new Label(window, "Spring types", "sans-bold");

    {
        Button* b = new Button(window, "structural");
        b->setFlags(Button::ToggleButton);
        b->setPushed(cp->enable_structural_constraints);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { cp->enable_structural_constraints = state; });

        b = new Button(window, "shearing");
        b->setFlags(Button::ToggleButton);
        b->setPushed(cp->enable_shearing_constraints);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { cp->enable_shearing_constraints = state; });

        b = new Button(window, "bending");
        b->setFlags(Button::ToggleButton);
        b->setPushed(cp->enable_bending_constraints);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { cp->enable_bending_constraints = state; });
    }

    // Mass-spring parameters

    new Label(window, "Parameters", "sans-bold");

    {
        Widget* panel = new Widget(window);
        GridLayout* layout = new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
        layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
        layout->setSpacing(0, 10);
        panel->setLayout(layout);

        new Label(panel, "density :", "sans-bold");

        FloatBox<double>* fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(cp->density / 10);
        fb->setUnits("g/cm^2");
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { cp->density = (double)(value * 10); });

        new Label(panel, "ks :", "sans-bold");

        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(cp->ks);
        fb->setUnits("N/m");
        fb->setSpinnable(true);
        fb->setMinValue(0);
        fb->setCallback([this](float value) { cp->ks = value; });
    }

    // Simulation constants

    new Label(window, "Simulation", "sans-bold");

    {
        Widget* panel = new Widget(window);
        GridLayout* layout = new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
        layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
        layout->setSpacing(0, 10);
        panel->setLayout(layout);

        new Label(panel, "frames/s :", "sans-bold");

        IntBox<int>* fsec = new IntBox<int>(panel);
        fsec->setEditable(true);
        fsec->setFixedSize(Vector2i(100, 20));
        fsec->setFontSize(14);
        fsec->setValue(frames_per_sec);
        fsec->setSpinnable(true);
        fsec->setCallback([this](int value) { frames_per_sec = value; });

        new Label(panel, "steps/frame :", "sans-bold");

        IntBox<int>* num_steps = new IntBox<int>(panel);
        num_steps->setEditable(true);
        num_steps->setFixedSize(Vector2i(100, 20));
        num_steps->setFontSize(14);
        num_steps->setValue(simulation_steps);
        num_steps->setSpinnable(true);
        num_steps->setMinValue(0);
        num_steps->setCallback([this](int value) { simulation_steps = value; });
    }

    // Damping slider and textbox

    new Label(window, "Damping", "sans-bold");

    {
        Widget* panel = new Widget(window);
        panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 5));

        Slider* slider = new Slider(panel);
        slider->setValue(cp->damping);
        slider->setFixedWidth(105);

        TextBox* percentage = new TextBox(panel);
        percentage->setFixedWidth(75);
        percentage->setValue(to_string(cp->damping));
        percentage->setUnits("%");
        percentage->setFontSize(14);

        slider->setCallback([percentage](float value) { percentage->setValue(std::to_string(value)); });
        slider->setFinalCallback([&](float value) {
            cp->damping = (double)value;
            // cout << "Final slider value: " << (int)(value * 100) << endl;
        });
    }

    // Collision options
    new Label(window, "Sphere Collision", "sans-bold");
    // {
    //     Button* b = new Button(window, "Use SDF for sphere");
    //     b->setFlags(Button::ToggleButton);
    //     b->setPushed(use_sdf_collision);
    //     b->setFontSize(14);
    //     b->setChangeCallback([this](bool state) { use_sdf_collision = state; });
    // }
    {
        Button* b = new Button(window, "Use CCD");
        b->setFlags(Button::ToggleButton);
        b->setPushed(use_ccd_collision);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { use_ccd_collision = state; });
    }
    {
        Button* b = new Button(window, "Use Pure Ray Marching");
        b->setFlags(Button::ToggleButton);
        b->setPushed(use_ray_marching_collision);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { use_ray_marching_collision = state; });
    }

    new Label(window, "Dynamic SDF Object", "sans-bold");
    {
        Button* b = new Button(window, "Enable (collision with ray marching)");
        b->setFlags(Button::ToggleButton);
        b->setPushed(use_dynamic_sdf_object);
        b->setFontSize(14);
        b->setChangeCallback([this](bool state) { use_dynamic_sdf_object = state; });
    }

    // Dynamic SDF Object Parameters
    {
        Widget* panel = new Widget(window);
        GridLayout* layout = new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
        layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
        layout->setSpacing(0, 10);
        panel->setLayout(layout);

        new Label(panel, "Particle Count :", "sans-bold");
        IntBox<int>* ib = new IntBox<int>(panel);
        ib->setEditable(true);
        ib->setFixedSize(Vector2i(100, 20));
        ib->setFontSize(14);
        ib->setValue(dynamic_sdf_num_samples);
        ib->setSpinnable(true);
        ib->setCallback([this](int value) { dynamic_sdf_num_samples = value; });

        new Label(panel, "Surface Threshold :", "sans-bold");
        FloatBox<double>* fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_surface_threshold);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_surface_threshold = value; });

        new Label(panel, "Max Distance :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_max_distance);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_max_distance = value; });

        new Label(panel, "Min Sphere Size :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_min_sphere_size);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_min_sphere_size = value; });

        new Label(panel, "Max Sphere Size :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_max_sphere_size);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_max_sphere_size = value; });

        new Label(panel, "Motion Speed :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_motion_speed);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_motion_speed = value; });

        new Label(panel, "Sphere Radius :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_sphere_radius);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_sphere_radius = value; });

        new Label(panel, "Box Size :", "sans-bold");
        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(dynamic_sdf_box_size);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { dynamic_sdf_box_size = value; });
    }

    // Gravity

    new Label(window, "Gravity", "sans-bold");

    {
        Widget* panel = new Widget(window);
        GridLayout* layout = new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
        layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
        layout->setSpacing(0, 10);
        panel->setLayout(layout);

        new Label(panel, "x :", "sans-bold");

        FloatBox<double>* fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(gravity.x);
        fb->setUnits("m/s^2");
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { gravity.x = value; });

        new Label(panel, "y :", "sans-bold");

        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(gravity.y);
        fb->setUnits("m/s^2");
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { gravity.y = value; });

        new Label(panel, "z :", "sans-bold");

        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(gravity.z);
        fb->setUnits("m/s^2");
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { gravity.z = value; });
    }

    window = new Window(screen, "Appearance");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout(15, 6, 14, 5));

    // Appearance

    {

        ComboBox* cb = new ComboBox(window, shaders_combobox_names);
        cb->setFontSize(14);
        cb->setCallback([this, screen](int idx) { active_shader_idx = idx; });
        cb->setSelectedIndex(active_shader_idx);
    }

    // Shader Parameters

    new Label(window, "Color", "sans-bold");

    {
        ColorWheel* cw = new ColorWheel(window, color);
        cw->setColor(this->color);
        cw->setCallback([this](const nanogui::Color& color) { this->color = color; });
    }

    new Label(window, "Parameters", "sans-bold");

    {
        Widget* panel = new Widget(window);
        GridLayout* layout = new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
        layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
        layout->setSpacing(0, 10);
        panel->setLayout(layout);

        new Label(panel, "Normal :", "sans-bold");

        FloatBox<double>* fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(this->m_normal_scaling);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { this->m_normal_scaling = value; });

        new Label(panel, "Height :", "sans-bold");

        fb = new FloatBox<double>(panel);
        fb->setEditable(true);
        fb->setFixedSize(Vector2i(100, 20));
        fb->setFontSize(14);
        fb->setValue(this->m_height_scaling);
        fb->setSpinnable(true);
        fb->setCallback([this](float value) { this->m_height_scaling = value; });
    }
}

void ClothSimulator::drawVolumeRendering(GLShader& shader) {
    // Create a full-screen quad for volume rendering
    // This quad covers the entire screen in clip space
    MatrixXf positions(4, 4);
    MatrixXf normals(4, 4);
    MatrixXf tangents(4, 4);
    MatrixXf uvs(2, 4);

    // Full-screen quad vertices in clip space
    positions.col(0) << -1.0f, -1.0f, 0.0f, 1.0f; // Bottom-left
    positions.col(1) << 1.0f, -1.0f, 0.0f, 1.0f;  // Bottom-right
    positions.col(2) << -1.0f, 1.0f, 0.0f, 1.0f;  // Top-left
    positions.col(3) << 1.0f, 1.0f, 0.0f, 1.0f;   // Top-right

    // Normals (not used in ray marching but required by vertex shader)
    normals.col(0) << 0.0f, 0.0f, 1.0f, 0.0f;
    normals.col(1) << 0.0f, 0.0f, 1.0f, 0.0f;
    normals.col(2) << 0.0f, 0.0f, 1.0f, 0.0f;
    normals.col(3) << 0.0f, 0.0f, 1.0f, 0.0f;

    // Tangents (not used in ray marching but required by vertex shader)
    tangents.col(0) << 1.0f, 0.0f, 0.0f, 0.0f;
    tangents.col(1) << 1.0f, 0.0f, 0.0f, 0.0f;
    tangents.col(2) << 1.0f, 0.0f, 0.0f, 0.0f;
    tangents.col(3) << 1.0f, 0.0f, 0.0f, 0.0f;

    // UV coordinates for the quad
    uvs.col(0) << 0.0f, 0.0f; // Bottom-left
    uvs.col(1) << 1.0f, 0.0f; // Bottom-right
    uvs.col(2) << 0.0f, 1.0f; // Top-left
    uvs.col(3) << 1.0f, 1.0f; // Top-right

    // Upload attributes
    shader.uploadAttrib("in_position", positions, false);
    shader.uploadAttrib("in_normal", normals, false);
    shader.uploadAttrib("in_tangent", tangents, false);
    shader.uploadAttrib("in_uv", uvs, false);

    // Draw the quad using triangle strip
    shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
}
