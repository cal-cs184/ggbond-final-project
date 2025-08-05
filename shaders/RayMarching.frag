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

// Simplified parameters
const int   MAX_STEPS                = 40;   // Reduced step count for better performance
const float MAX_DIST                 = 30.0;
const float EPSILON                  = 0.002;
const float DENSITY_MULTIPLIER       = 0.2;
const float ABSORPTION               = 0.3;
const float PARTICLE_INFLUENCE_RADIUS = 0.8;
const float STEP_SIZE                = 0.06; // Larger step for performance, acceptable quality
const float GRID_SIZE                = 0.05;

// Simplified noise function
float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

// Get cloth particle position
vec3 getClothParticlePosition(int x, int y) {
    float u = (float(x) + 0.5) / float(u_cloth_tex_width);
    float v = (float(y) + 0.5) / float(u_cloth_tex_height);
    return texture(u_cloth_particles_tex, vec2(u, v)).xyz;
}

// Simplified density calculation – only checks nearby particles
float calculateDensity(vec3 pos) {
    float density = 0.0;

    // Compute grid position
    int grid_x = int(pos.x / GRID_SIZE);
    int grid_z = int(pos.z / GRID_SIZE);

    // Only check a 3×3 grid neighborhood
    for (int dz = -1; dz <= 1; dz++) {
        for (int dx = -1; dx <= 1; dx++) {
            int x = grid_x + dx;
            int z = grid_z + dz;

            if (x < 0 || x >= u_cloth_tex_width || z < 0 || z >= u_cloth_tex_height) continue;

            vec3 particle_pos = getClothParticlePosition(x, z);
            float dist = length(pos - particle_pos);

            if (dist < PARTICLE_INFLUENCE_RADIUS) {
                float falloff = 1.0 - smoothstep(0.0, PARTICLE_INFLUENCE_RADIUS, dist);
                density += falloff;
            }
        }
    }

    return min(density * 0.3, 1.0);
}

// Enhanced volume color with rich gradients
vec3 getVolumeColor(float density, float height) {
    // Rich color palette based on density and height
    vec3 lowDensity = mix(vec3(0.1, 0.3, 0.8), vec3(0.8, 0.2, 0.9), height);
    vec3 highDensity = mix(vec3(0.2, 0.8, 0.6), vec3(0.9, 0.7, 0.3), height);
    vec3 color = mix(lowDensity, highDensity, density);
    
    // Add subtle noise for texture variation
    float noise = fract(sin(dot(vec2(height, density), vec2(12.9898, 78.233))) * 43758.5453);
    color += noise * 0.1;
    
    return color * (0.6 + 0.4 * density);
}

// Simplified ray marching
vec4 rayMarch(vec3 ro, vec3 rd) {
    float t             = 0.0;
    vec3  color         = vec3(0.0);
    float transmittance = 1.0;

    // // Early exit if the ray points upward above the volume
    // if (rd.y > 0.0 && ro.y > 2.0) {
    //     return vec4(0.0, 0.0, 0.0, 0.0);
    // }

    for (int i = 0; i < MAX_STEPS; i++) {
        // if (t > MAX_DIST || transmittance < 0.1) break;

        vec3 pos = ro + rd * t;

        // if (pos.y < 0.0 || pos.y > 2.0) break;

        float density = calculateDensity(pos);

        // if (density > EPSILON) {
            // Simplified lighting
            vec3  light_dir = normalize(u_light_pos - pos);
            float diffuse   = max(0.0, dot(v_normal.xyz, light_dir));
            float ambient   = 0.25;

            float height        = smoothstep(0.0, 2.0, pos.y);
            vec3  volume_color  = getVolumeColor(density, height);

            vec3 sample_color = volume_color * (ambient + diffuse * 0.5) * density;

            // Add subtle glow effect
            float glow = density * 0.1;
            sample_color += volume_color * glow;

            color += sample_color * transmittance;
            transmittance *= exp(-density * ABSORPTION);
        // }

        t += STEP_SIZE;
    }

    return vec4(color, 1.0 - transmittance);
}

void main() {
    vec3 ro = u_cam_pos;
    vec3 rd = normalize(v_world_position - u_cam_pos);

    vec4 volume_color = rayMarch(ro, rd);

    // Improve surface color
    vec4 surface_color = vec4(0.8, 0.8, 0.8, 0.2);
    surface_color.rgb += v_normal.xyz * 0.2;

    // Better blending between surface and volume
    // out_color = mix(surface_color, volume_color, 0.8); // Fixed blend weight

    // // Add environment reflection
    // vec3 reflection = reflect(rd, normalize(v_normal.xyz));
    // vec4 env_color  = texture(u_texture_cubemap, reflection);
    // out_color = mix(out_color, env_color, 0.1 * (1.0 - volume_color.a));

    // Ensure reasonable opacity
    // out_color.a = max(1.0, volume_color.a);
    // out_color.a = clamp(volume_color.a, 0.3, 1.0);

    // Enhanced final output with better volume effect
    if (volume_color.a < 0.01) {
        // Fallback to surface color for empty areas
        out_color = surface_color;
    } else {
        // Rich volume rendering with proper alpha
        out_color = volume_color;
        out_color.a = clamp(volume_color.a, 0.2, 0.9); // Semi-transparent volume
    }
}