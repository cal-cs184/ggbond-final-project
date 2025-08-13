#ifndef CGL_CLOTH_SIMULATOR_H
#define CGL_CLOTH_SIMULATOR_H

#include <memory>
#include <nanogui/nanogui.h>

#include "camera.h"
#include "cloth.h"
#include "collision/collisionObject.h"

using namespace nanogui;

struct UserShader;
enum ShaderTypeHint { WIREFRAME = 0, NORMALS = 1, PHONG = 2, VOLUME_RENDERING = 3 };

class ClothSimulator {
  public:
    ClothSimulator(std::string project_root, Screen* screen);
    ~ClothSimulator();

    void init();

    void loadCloth(Cloth* cloth);
    void loadClothParameters(ClothParameters* cp);
    void loadCollisionObjects(vector<CollisionObject*>* objects);
    virtual bool isAlive();
    virtual void drawContents();

    // Screen events

    virtual bool cursorPosCallbackEvent(double x, double y);
    virtual bool mouseButtonCallbackEvent(int button, int action, int modifiers);
    virtual bool keyCallbackEvent(int key, int scancode, int action, int mods);
    virtual bool dropCallbackEvent(int count, const char** filenames);
    virtual bool scrollCallbackEvent(double x, double y);
    virtual bool resizeCallbackEvent(int width, int height);

  private:
    virtual void initGUI(Screen* screen);
    void drawWireframe(GLShader& shader);
    void drawNormals(GLShader& shader);
    void drawPhong(GLShader& shader);
    void drawVolumeRendering(GLShader& shader);

    void load_shaders();
    void load_textures();

    // File management

    std::string m_project_root;

    // Camera methods

    virtual void resetCamera();
    virtual Matrix4f getProjectionMatrix();
    virtual Matrix4f getViewMatrix();

    // Default simulation values

    int frames_per_sec = 90;
    int simulation_steps = 30;

    CGL::Vector3D gravity = CGL::Vector3D(0, -9.8, 0);
    nanogui::Color color = nanogui::Color(1.0f, 1.0f, 1.0f, 1.0f);

    Cloth* cloth;
    ClothParameters* cp;
    vector<CollisionObject*>* collision_objects;

    // OpenGL attributes

    int active_shader_idx = 0;

    vector<UserShader> shaders;
    vector<std::string> shaders_combobox_names;

    // OpenGL textures

    Vector3D m_gl_texture_1_size;
    Vector3D m_gl_texture_2_size;
    Vector3D m_gl_texture_3_size;
    Vector3D m_gl_texture_4_size;
    GLuint m_gl_texture_1;
    GLuint m_gl_texture_2;
    GLuint m_gl_texture_3;
    GLuint m_gl_texture_4;
    GLuint m_gl_cubemap_tex;

    // Ray marching related textures
    GLuint m_cloth_particles_tex; // Texture to store cloth particle positions
    int m_cloth_tex_width;        // Width of the particle texture
    int m_cloth_tex_height;       // Height of the particle texture

    // OpenGL customizable inputs

    double m_normal_scaling = 2.0;
    double m_height_scaling = 0.1;

    // Camera attributes

    CGL::Camera camera;
    CGL::Camera canonicalCamera;

    double view_distance;
    double canonical_view_distance;
    double min_view_distance;
    double max_view_distance;

    double scroll_rate;

    // Screen methods

    Screen* screen;
    void mouseLeftDragged(double x, double y);
    void mouseRightDragged(double x, double y);
    void mouseMoved(double x, double y);

    // Mouse flags

    bool left_down = false;
    bool right_down = false;
    bool middle_down = false;

    // Keyboard flags

    bool ctrl_down = false;

    // Simulation flags

    bool is_paused = true;

    // Time for animations
    double simulation_time = 0.0;

    // Collision options
    bool use_sdf_collision = false;          // toggle in GUI
    bool use_ccd_collision = false;          // toggle in GUI
    bool use_ray_marching_collision = false; // toggle in GUI
    bool use_dynamic_sdf_object = false;     // toggle in GUI for dynamic SDF object

    // Dynamic SDF Object parameters
    int dynamic_sdf_num_samples = 8000;
    double dynamic_sdf_surface_threshold = 0.02;
    double dynamic_sdf_max_distance = 2.5;
    double dynamic_sdf_min_sphere_size = 0.003;
    double dynamic_sdf_max_sphere_size = 0.04;
    double dynamic_sdf_motion_speed = 2.0;
    double dynamic_sdf_sphere_radius = 1.0;
    double dynamic_sdf_box_size = 0.8;

    // Screen attributes

    int mouse_x;
    int mouse_y;

    int screen_w;
    int screen_h;

    bool is_alive = true;

    Vector2i default_window_size = Vector2i(1024, 800);
};

struct UserShader {
    UserShader(std::string display_name, std::shared_ptr<GLShader> nanogui_shader, ShaderTypeHint type_hint)
        : display_name(display_name), nanogui_shader(nanogui_shader), type_hint(type_hint) {}

    std::shared_ptr<GLShader> nanogui_shader;
    std::string display_name;
    ShaderTypeHint type_hint;
};

#endif // CGL_CLOTH_SIM_H
