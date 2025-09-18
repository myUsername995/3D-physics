#version 330 core

layout(location = 0) in vec2 aPos;
layout(location = 1) in vec3 aInstanceXYR;
layout(location = 2) in vec3 aInstanceColor;

uniform vec2 uScreenSize;

out vec3 vColor;

void main() {
    vec2 pos = aPos * aInstanceXYR.z + aInstanceXYR.xy;
    pos = (pos / uScreenSize) * 2.0 - 1.0;
    pos.y *= -1.0;
    gl_Position = vec4(pos, 0.0, 1.0);

    vColor = aInstanceColor;
}