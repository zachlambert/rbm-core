#version 330 core

// Input vertex data: position and normal
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexNormal_modelspace;
layout(location = 2) in vec3 vertexColor_in;

out vec4 vertexColor;
out vec3 Normal_worldspace;

uniform mat4 MVP;
uniform mat4 M;

void main(){
    gl_Position = MVP * vec4(vertexPosition_modelspace, 1);
    Normal_worldspace = normalize((M * vec4(vertexNormal_modelspace, 0)).xyz);
    vertexColor = vec4(vertexColor_in, 1);
}
