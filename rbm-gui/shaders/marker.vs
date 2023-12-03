#version 330 core

// Input vertex data: position and normal
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexNormal_modelspace;

out vec4 vertexColor;
out vec3 Position_cameraspace;
out vec3 Normal_cameraspace;

uniform mat4 MVP;
uniform mat4 MV;
uniform vec4 color;

void main(){
    gl_Position = MVP * vec4(vertexPosition_modelspace, 1);
    Position_cameraspace = (MV * vec4(vertexPosition_modelspace, 1)).xyz;
    Normal_cameraspace = normalize((MV * vec4(vertexNormal_modelspace, 0)).xyz);
    vertexColor = color;
}
