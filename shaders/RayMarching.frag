#version 330

// Uniform variables
uniform mat4 u_view_projection;
uniform mat4 u_model;

uniform float u_normal_scaling;
uniform float u_height_scaling;

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// Textures
uniform sampler2D u_texture_1; // Used for noise texture
uniform sampler2D u_texture_2; // Used for color gradient
uniform sampler2D u_texture_3; // Reserved for future use
uniform sampler2D u_texture_4; // Reserved for future use
uniform samplerCube u_texture_cubemap;

// Cloth particles texture
uniform sampler2D u_cloth_particles_tex;
uniform int u_cloth_tex_width;
uniform int u_cloth_tex_height;

// Input from vertex shader
in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;
in vec3 v_world_position;
in vec3 v_view_dir;

// Output color
out vec4 out_color;

// Ray marching parameters
const int MAX_STEPS = 100;
const float MAX_DIST = 100.0;
const float EPSILON = 0.001;

// Volume rendering parameters
const float DENSITY_MULTIPLIER = 0.5;
const float ABSORPTION = 0.5;
const vec3 SCATTER_COLOR = vec3(0.7, 0.3, 1.0); // Purple-ish color

// 2D density field parameters
const int DENSITY_FIELD_RESOLUTION = 64;
const float DENSITY_FIELD_SIZE = 10.0;
const float PARTICLE_INFLUENCE_RADIUS = 1.5;

// Noise functions for procedural effects
float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x) {
    // A simple 3D noise function
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    
    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

// Function to get cloth particle position from texture
vec3 getClothParticlePosition(int x, int y) {
    // Calculate UV coordinates for texture lookup
    float u = (float(x) + 0.5) / float(u_cloth_tex_width);
    float v = (float(y) + 0.5) / float(u_cloth_tex_height);
    
    // Sample the texture to get particle position
    return texture(u_cloth_particles_tex, vec2(u, v)).xyz;
}

// Function to calculate 2D density field from cloth particles
float calculate2DDensity(vec2 pos_xz) {
    float density = 0.0;
    
    // Calculate approximate grid position for optimization
    float grid_cell_size = PARTICLE_INFLUENCE_RADIUS * 0.5;
    vec2 grid_pos = floor(pos_xz / grid_cell_size);
    
    // Define a search window (only check nearby particles)
    int search_radius = 3;
    
    // Iterate through nearby cloth particles only
    for (int dy = -search_radius; dy <= search_radius; dy++) {
        // Calculate y index in cloth grid
        int base_y = int(grid_pos.y) + dy;
        if (base_y < 0 || base_y >= u_cloth_tex_height) continue;
        
        for (int dx = -search_radius; dx <= search_radius; dx++) {
            // Calculate x index in cloth grid
            int base_x = int(grid_pos.x) + dx;
            if (base_x < 0 || base_x >= u_cloth_tex_width) continue;
            
            // Get particle position
            vec3 particle_pos = getClothParticlePosition(base_x, base_y);
            
            // Project particle position to xz-plane
            vec2 particle_pos_xz = vec2(particle_pos.x, particle_pos.z);
            
            // Calculate distance in xz-plane
            float dist = length(pos_xz - particle_pos_xz);
            
            // If within influence radius, add to density
            if (dist < PARTICLE_INFLUENCE_RADIUS) {
                // Use smooth falloff function
                float falloff = 1.0 - smoothstep(0.0, PARTICLE_INFLUENCE_RADIUS, dist);
                density += falloff * falloff; // Squared for smoother falloff
            }
        }
    }
    
    // Normalize density
    return min(density * 0.2, 1.0);
}

// Function to sample density at a point in space
float sampleDensity(vec3 pos) {
    // Get 2D density from cloth projection onto xz-plane
    float density_2d = calculate2DDensity(vec2(pos.x, pos.z));
    
    // Calculate vertical falloff based on height
    float y_min = 0.0;  // Minimum height where density starts
    float y_max = 2.0;  // Maximum height where density ends
    float y_falloff = 1.0 - smoothstep(y_min, y_max, pos.y);
    
    // Combine 2D density with height falloff
    float base_density = density_2d * y_falloff;
    
    // Add some noise for interesting effects
    float time_factor = u_height_scaling * 10.0; // Use height scaling as time
    float noise_val = noise(pos * 3.0 + vec3(0.0, time_factor * 0.1, time_factor * 0.2));
    float noise_val2 = noise(pos * 7.0 - vec3(time_factor * 0.05, 0.0, time_factor * 0.15));
    
    // Create swirling effect
    vec2 swirl = vec2(
        sin(pos.y * 2.0 + time_factor * 0.2) * 0.5,
        cos(pos.x * 2.0 + time_factor * 0.3) * 0.5
    );
    pos.xz += swirl * noise_val * base_density;
    
    // Combine base density with noise for final effect
    float final_density = base_density * (0.7 + 0.3 * noise_val) + 0.2 * noise_val2 * base_density;
    
    return final_density * DENSITY_MULTIPLIER;
}

// Get color from gradient texture based on density and height
vec3 getDensityColor(float density, float height, float noise) {
    // Use texture2 as a color gradient
    float gradient_pos = density * 0.7 + height * 0.2 + noise * 0.1;
    vec3 color = texture(u_texture_2, vec2(gradient_pos, 0.5)).rgb;
    
    // If no texture is available, use a default color gradient
    if (color == vec3(0.0)) {
        // Create a colorful gradient
        color = mix(
            mix(vec3(0.2, 0.1, 0.7), vec3(0.7, 0.0, 1.0), density), // Purple to violet
            mix(vec3(1.0, 0.3, 0.5), vec3(0.9, 0.6, 0.1), density), // Pink to orange
            height
        );
    }
    
    return color;
}

// Ray marching function with optimizations
vec4 rayMarch(vec3 ro, vec3 rd) {
    float t = 0.0;
    vec3 color = vec3(0.0);
    float transmittance = 1.0;
    
    // Early exit if ray is pointing away from the volume of interest
    if (rd.y > 0.0 && ro.y > 3.0) {
        return vec4(0.0, 0.0, 0.0, 0.0);
    }
    
    // Calculate entry and exit points for the volume
    float y_min = 0.0;
    float y_max = 3.0;
    
    // Skip empty space before entering the volume
    if (ro.y > y_max && rd.y < 0.0) {
        // Ray starts above the volume, calculate intersection with top plane
        float t_entry = (y_max - ro.y) / rd.y;
        t = max(0.0, t_entry);
    }
    
    // Adaptive step size based on viewing angle
    float base_step = 0.1;
    float angle_factor = abs(dot(rd, vec3(0.0, 1.0, 0.0)));
    float view_step_factor = mix(0.5, 1.5, angle_factor);
    
    // Step through the volume
    for(int i = 0; i < MAX_STEPS; i++) {
        if(t > MAX_DIST || transmittance < 0.01) break;
        
        // Current position along the ray
        vec3 pos = ro + rd * t;
        
        // Exit if we're below the volume
        if (pos.y < y_min) break;
        
        // Exit if we're above the volume and moving upward
        if (pos.y > y_max && rd.y > 0.0) break;
        
        // Sample density at current position
        float density = sampleDensity(pos);
        
        if(density > EPSILON) {
            // Calculate lighting
            vec3 light_dir = normalize(u_light_pos - pos);
            float diffuse = max(0.0, dot(v_normal.xyz, light_dir));
            
            // Add some ambient light
            float ambient = 0.3;
            
            // Get noise value for color variation (only calculate if needed)
            float time_factor = u_height_scaling * 10.0;
            float noise_val = noise(pos * 2.5 + vec3(time_factor * 0.1));
            
            // Get color based on density, height and noise
            float normalized_height = smoothstep(0.0, 2.0, pos.y);
            vec3 volume_color = getDensityColor(density, normalized_height, noise_val);
            
            // Calculate color based on density, lighting and volume color
            vec3 sample_color = volume_color * (ambient + diffuse) * density;
            
            // Add some emissive glow
            float glow = density * density * 0.5;
            sample_color += volume_color * glow;
            
            // Accumulate color with transmittance
            color += sample_color * transmittance;
            
            // Update transmittance (light absorption)
            transmittance *= exp(-density * ABSORPTION);
            
            // Reduce step size in dense regions
            base_step = mix(0.05, 0.2, exp(-density * 2.0));
        }
        
        // Adaptive step size based on density and viewing angle
        float step_size = base_step * view_step_factor;
        t += step_size;
    }
    
    return vec4(color, 1.0 - transmittance);
}

void main() {
    // Ray origin (camera position)
    vec3 ro = u_cam_pos;
    
    // Ray direction (from vertex shader)
    vec3 rd = normalize(v_world_position - u_cam_pos);
    
    // Perform ray marching
    vec4 volume_color = rayMarch(ro, rd);
    
    // Combine with surface color (make surface more transparent)
    vec4 surface_color = (vec4(1.0, 1.0, 1.0, 0.0) + v_normal) / 2.0;
    surface_color.a = 0.3; // Make cloth surface semi-transparent
    
    // Mix volume and surface color based on volume alpha
    out_color = mix(surface_color, volume_color, volume_color.a * 0.9);
    
    // Add some environment reflection
    vec3 reflection = reflect(rd, normalize(v_normal.xyz));
    vec4 env_color = texture(u_texture_cubemap, reflection);
    
    // Blend with environment reflection
    out_color = mix(out_color, env_color, 0.1 * (1.0 - volume_color.a));
    
    // Ensure alpha is set (more transparent where volume is thin)
    out_color.a = max(0.3, volume_color.a);
}