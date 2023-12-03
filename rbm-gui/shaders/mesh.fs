#version 330 core

in vec4 vertexColor;
in vec3 Normal_worldspace;
out vec4 color;

void main(){
    float cosTheta = clamp(dot(Normal_worldspace, vec3(0, 0, 1)), 0, 1);
    float factor = 0.5 + 0.5 * cosTheta;
    color = vertexColor * factor;
}
