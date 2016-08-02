#version 400
 
layout (location = 0) 
in vec3 Pos;
in uint gl_VertexID;
out vec3 vPos;

void main()
{
    vPos = Pos;
	
	if(gl_VertexID == 8)
	vPos.y += 1.0;
}