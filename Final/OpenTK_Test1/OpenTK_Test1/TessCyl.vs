#version 400
 
in uint gl_VertexID;
in uint gl_InstanceID;
out vec3 vPos;

out VS_OUT
{
	vec3 Pos;
} Out;

void main()
{
	uint x = gl_VertexID & 7;
	uint y = gl_VertexID >> 3;
	
	const float pi = 3.141592653589793238462643383279502884197169399375105820974944592;
	
	float lx = float(x) / 8.0 * 2.0 * pi;
	float ly = float(y);
	
	Out.Pos.y = ly;
	Out.Pos.x = cos(lx);
	Out.Pos.z = sin(lx);
	// sincos(lx, vPos.z, vPos.x);
	
	if(gl_InstanceID != 0)
	Out.Pos.y -= 17.0;
}