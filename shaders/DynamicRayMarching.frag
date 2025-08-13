#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;
uniform mat4 u_inv_view_projection;
uniform float u_time; // Time parameter for animation

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;
in vec4 v_tangent;

out vec4 out_color;

// Ray marching parameters
const int MAX_STEPS = 64;
const float MAX_DIST = 10.0;
const float EPSILON = 0.001;

// Noise function - for creating dynamic shapes
float noise(vec3 p) {
    return fract(sin(dot(p, vec3(12.9898, 78.233, 45.164))) * 43758.5453);
}

float smoothNoise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    f = f * f * (3.0 - 2.0 * f); // smoothstep
    
    float a = noise(i);
    float b = noise(i + vec3(1.0, 0.0, 0.0));
    float c = noise(i + vec3(0.0, 1.0, 0.0));
    float d = noise(i + vec3(1.0, 1.0, 0.0));
    float e = noise(i + vec3(0.0, 0.0, 1.0));
    float f1 = noise(i + vec3(1.0, 0.0, 1.0));
    float g = noise(i + vec3(0.0, 1.0, 1.0));
    float h = noise(i + vec3(1.0, 1.0, 1.0));
    
    return mix(mix(mix(a, b, f.x), mix(c, d, f.x), f.y),
               mix(mix(e, f1, f.x), mix(g, h, f.x), f.y), f.z);
}

// Fractal noise - creates more complex shapes
float fractalNoise(vec3 p, int octaves) {
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;
    
    for (int i = 0; i < octaves; i++) {
        value += amplitude * smoothNoise(p * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }
    
    return value;
}

// Basic shape SDF functions
float sphereSDF(vec3 p, float r) {
    return length(p) - r;
}

float boxSDF(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float torusSDF(vec3 p, vec2 t) {
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}

// Boolean operations
float unionSDF(float a, float b) {
    return min(a, b);
}

float intersectionSDF(float a, float b) {
    return max(a, b);
}

float differenceSDF(float a, float b) {
    return max(a, -b);
}

// Deformation functions
vec3 twist(vec3 p, float k) {
    float c = cos(k * p.y);
    float s = sin(k * p.y);
    mat2 m = mat2(c, -s, s, c);
    return vec3(m * p.xz, p.y);
}

vec3 bend(vec3 p, float k) {
    float c = cos(k * p.x);
    float s = sin(k * p.x);
    mat2 m = mat2(c, -s, s, c);
    return vec3(m * p.xy, p.z);
}

// Main SDF function - completely consistent with CPU-side DynamicSDFObject: union(sphere, box)
float sceneSDF(vec3 p) {
    float t = u_time * 0.5;
    vec3 sphere_center = vec3(sin(t) * 2.0, 0.0, 0.0);
    vec3 box_center = vec3(cos(t) * 2.0, 0.0, 0.0);

    float d_sphere = sphereSDF(p - sphere_center, 1.0);
    float d_box = boxSDF(p - box_center, vec3(0.8));
    return min(d_sphere, d_box);
}

// Calculate normal
vec3 calculateNormal(vec3 p) {
    const float h = 0.001;
    const vec2 k = vec2(1.0, -1.0);
    return normalize(k.xyy * sceneSDF(p + k.xyy * h) +
                    k.yyx * sceneSDF(p + k.yyx * h) +
                    k.yxy * sceneSDF(p + k.yxy * h) +
                    k.xxx * sceneSDF(p + k.xxx * h));
}

// Ray marching main function
float rayMarch(vec3 ro, vec3 rd) {
    float depth = 0.0;
    
    for (int i = 0; i < MAX_STEPS; i++) {
        vec3 p = ro + rd * depth;
        float dist = sceneSDF(p);
        
        if (dist < EPSILON) {
            return depth;
        }
        
        depth += dist;
        
        if (depth > MAX_DIST) {
            break;
        }
    }
    
    return -1.0;
}

// Lighting calculation
vec3 calculateLighting(vec3 p, vec3 normal, vec3 ro, vec3 rd) {
    vec3 light_dir = normalize(u_light_pos - p);
    vec3 view_dir = normalize(ro - p);
    vec3 half_dir = normalize(light_dir + view_dir);
    
    // Diffuse reflection
    float diffuse = max(0.0, dot(normal, light_dir));
    
    // Specular reflection
    float specular = pow(max(0.0, dot(normal, half_dir)), 32.0);
    
    // Ambient light
    float ambient = 0.2;
    
    // Position-based dynamic color
    vec3 base_color = vec3(0.6, 0.8, 1.0);
    float noise_color = fractalNoise(p * 3.0 + u_time * 0.2, 2);
    vec3 color = mix(base_color, vec3(1.0, 0.5, 0.2), noise_color * 0.5);
    
    return color * (ambient + diffuse) + vec3(1.0) * specular * 0.5;
}

void main() {
    // Calculate ray direction
    vec2 ndc = v_uv * 2.0 - 1.0;
    vec4 clip = vec4(ndc, -1.0, 1.0);
    vec4 world = u_inv_view_projection * clip;
    world /= world.w;
    vec3 ray_direction = normalize(world.xyz - u_cam_pos);
    
    vec3 ray_origin = u_cam_pos;
    
    // Execute ray marching
    float dist = rayMarch(ray_origin, ray_direction);
    
    if (dist > 0.0) {
        // Calculate collision point
        vec3 hit_point = ray_origin + ray_direction * dist;
        vec3 normal = calculateNormal(hit_point);
        
        // Calculate lighting
        vec3 color = calculateLighting(hit_point, normal, ray_origin, ray_direction);
        
        // Add fog effect
        float fog = 1.0 - exp(-dist * 0.1);
        color = mix(color, vec3(0.5, 0.7, 1.0), fog);
        
        out_color = vec4(color, 1.0);
    } else {
        // Background
        vec3 background = vec3(0.05, 0.05, 0.1);
        out_color = vec4(background, 1.0);
    }
} 