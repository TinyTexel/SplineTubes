#version 430

#define float2 vec2
#define float3 vec3
#define float4 vec4

float saturate(float x) { return clamp(x, 0.0, 1.0); }
vec2 saturate(vec2 x) { return clamp(x, vec2(0.0), vec2(1.0)); }
vec3 saturate(vec3 x) { return clamp(x, vec3(0.0), vec3(1.0)); }
vec4 saturate(vec4 x) { return clamp(x, vec4(0.0), vec4(1.0)); }

#define lerp(a, b, l) mix(a, b, l)
#define mul(a, b) (a * b)
// #define Pos gl_Position

float SinDot(float2 a, float2 b)
{
	return a.x * b.y - a.y * b.x;
}

int ToInt(bool cond) { return cond ? 1 : 0; }
float ToFloat(bool cond) { return cond ? 1.0 : 0.0; }

layout (std140) uniform StaticBuffer
{ 
	vec2 PixelCount;
};

layout (std140) uniform PerFrameBuffer
{ 
    mat4 ViewProjMat;
	vec3 CamPos;
};

layout(triangles) in;
layout (triangle_strip, max_vertices = 3) out;
 
in TE_OUT
{
	float3 WPos; 		
	float3 Normal;			
	float3 UVW0;
} In[3];

out GS_OUT
{
	float3 WPos;
	float3 Normal;
	float3 UVW;
	float3 UVW0;
} Out;
 
 
void main()
{
	// float3 normal = normalize(cross(In[1].WPos - In[0].WPos, In[2].WPos - In[0].WPos));
	
	for(uint i = 0; i < gl_in.length(); i++)
	{		
		gl_Position = gl_in[i].gl_Position;

		Out.WPos   = In[i].WPos;
		Out.Normal = In[i].Normal;
		Out.UVW0   = In[i].UVW0;
		Out.UVW  = float3(ToFloat(i == 0), ToFloat(i == 1), ToFloat(i == 2));

		EmitVertex();
	}
}