#version 330 core

// Input vertex data: position and normal
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec4 vertexColor_in;
layout(location = 2) in float vertexNormal_modelspace;

out vec4 vertexColor;
uniform mat4 MVP;

void main(){
    gl_Position = MVP * vec4(vertexPosition_modelspace, 1);
    vertexColor = vertexColor_in;
}
