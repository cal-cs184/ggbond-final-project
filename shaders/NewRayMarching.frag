#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform sampler2D u_texture_1;

// Cloth particles texture for dynamic density calculation
uniform sampler2D u_cloth_particles_tex;
uniform int u_cloth_tex_width;
uniform int u_cloth_tex_height;

// For full-screen quad rendering
in vec2 v_uv;

out vec4 out_color;

// Smart ray marching parameters
const int MAX_STEPS = 24;
const float MAX_DIST = 6.0;
const float STEP_SIZE = 0.15;

// Optimized density function with early exit
float calculateDensity(vec3 pos) {
    float density = 0.0;
    int samples = 0;
    
    // Sample particles with early exit
    for (int i = 0; i < u_cloth_tex_width && samples < 16; i += 3) {
        for (int j = 0; j < u_cloth_tex_height && samples < 16; j += 3) {
            samples++;
            
            float u = (float(i) + 0.5) / float(u_cloth_tex_width);
            float v = (float(j) + 0.5) / float(u_cloth_tex_height);
            vec3 particle_pos = texture(u_cloth_particles_tex, vec2(u, v)).xyz;
            
            float dist = length(pos - particle_pos);
            if (dist < 1.2) {
                float falloff = 1.0 - dist / 1.2;
                density += falloff * 0.4;
                
                // Early exit if we have enough density
                if (density > 0.8) break;
            }
        }
        if (density > 0.8) break;
    }
    
    return min(density, 1.0);
}

void main() {
    // TRUE RAY MARCHING: Convert UV to ray direction
    vec2 ndc = v_uv * 2.0 - 1.0;
    vec3 ray_direction = normalize(vec3(ndc.x, ndc.y, -1.0));
    
    // Ray marching from camera position
    vec3 ray_origin = u_cam_pos;
    vec3 color = vec3(0.0);
    float transmittance = 1.0;
    
    // TRUE RAY MARCHING LOOP
    for (int i = 0; i < MAX_STEPS; i++) {
        vec3 sample_pos = ray_origin + ray_direction * (i * STEP_SIZE);
        
        // Early exit conditions
        if (i * STEP_SIZE > MAX_DIST || transmittance < 0.05) break;
        
        // Calculate density at this sample point
        float density = calculateDensity(sample_pos);
        
        if (density > 0.05) {
            // Simple lighting
            vec3 light_dir = normalize(u_light_pos - sample_pos);
            float diffuse = max(0.0, dot(normalize(vec3(0, 1, 0)), light_dir));
            vec3 sample_color = vec3(0.2, 0.6, 0.9) * (0.3 + 0.7 * diffuse) * density;
            
            // Accumulate color with proper volume rendering equation
            color += sample_color * transmittance * STEP_SIZE;
            transmittance *= (1.0 - density * 0.15);
        }
    }
    
    // Add background
    vec3 background = vec3(0.05, 0.05, 0.1);
    vec3 final_color = mix(background, color, 1.0 - transmittance);
    
    out_color = vec4(final_color, 1.0);
}
