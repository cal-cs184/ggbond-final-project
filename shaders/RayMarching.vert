#version 330

// Uniform variables
uniform mat4 u_model;
uniform mat4 u_view_projection;

// Input vertex attributes
in vec4 in_position;
in vec4 in_normal;
in vec4 in_tangent;
in vec2 in_uv;

// Output variables to fragment shader
out vec4 v_position;
out vec4 v_normal;
out vec2 v_uv;
out vec4 v_tangent;
out vec3 v_world_position;
out vec3 v_view_dir;

// Camera position for view direction calculation
uniform vec3 u_cam_pos;

void main() {
  // Transform position to world space
  v_position = u_model * in_position;
  v_world_position = v_position.xyz;
  
  // Transform normal and tangent to world space
  v_normal = normalize(u_model * in_normal);
  v_tangent = normalize(u_model * in_tangent);
  
  // Pass UV coordinates
  v_uv = in_uv;
  
  // Calculate view direction (from camera to fragment)
  v_view_dir = normalize(v_world_position - u_cam_pos);
  
  // Set final position in clip space
  gl_Position = u_view_projection * u_model * in_position;
}