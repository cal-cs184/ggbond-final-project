#version 330

// Uniform variables (declared but not used for full-screen quad)
uniform mat4 u_model;
uniform mat4 u_view_projection;

// Input vertex attributes for full-screen quad
in vec4 in_position;
in vec2 in_uv;

// Output variables to fragment shader
out vec2 v_uv;

void main() {
    // Pass UV coordinates directly
    v_uv = in_uv;
    
    // For full-screen quad, position is already in clip space
    // But we need to use the uniforms to prevent compiler optimization
    vec4 transformed_pos = u_model * in_position;
    gl_Position = u_view_projection * transformed_pos;
} 