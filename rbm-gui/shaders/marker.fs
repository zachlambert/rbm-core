#version 330 core

in vec4 vertexColor;
in vec3 Position_cameraspace;
in vec3 Normal_cameraspace;
out vec4 color;

void main(){
    // Light as if the light comes from the camera view
    float cosTheta = clamp(dot(Normal_cameraspace, vec3(0, 0, 1)), 0, 1);
    color = vertexColor * cosTheta;
}
