#version 330

// Uniform variables
uniform mat4 u_model;
uniform mat4 u_view_projection;

// Input attributes
in vec4 in_position;
in vec4 in_normal;
in vec4 in_tangent;
in vec2 in_uv;

// Output variables
out vec4 v_position;
out vec4 v_normal;
out vec2 v_uv;
out vec4 v_tangent;

void main() {
    // Transform vertex data
    v_position = u_model * in_position;
    v_normal = normalize(u_model * in_normal);
    v_uv = in_uv;
    v_tangent = normalize(u_model * in_tangent);
    
    // Set final position
    gl_Position = u_view_projection * u_model * in_position;
} 