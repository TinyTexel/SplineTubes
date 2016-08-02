#version 430
 
#define float2 vec2
#define float3 vec3
#define float4 vec4

#define saturate(x) clamp(x, 0.0, 1.0)
#define lerp(a, b, l) mix(a, b, l)
 #define mul(a, b) (a * b)
 

const float Pi = 3.141592653589793238462643383279502884197169399375105820974944592;
const float Rcp3 = 1.0 / 3.0;

const float RcpMaxLId = 1.0 / 31.0;

const uint SegmentSize = 1 + 16;
const uint NodeSize = 4;

layout (std140) uniform PerFrameBuffer
{ 
    mat4 ViewProjMat;
	vec3 CamPos;
};

layout (std140) uniform StaticBuffer
{ 
	vec2 PixelCount;
};


#define TUBES_DATA_TYPE_IS_SSBO

#ifdef TUBES_DATA_TYPE_IS_SSBO
	layout (std430) buffer TubesData
	{
		float4 TubesBuffer[];
	};
#else
	layout (std140) uniform TubesData
	{
		float4 TubesBuffer[8];
	};
#endif


float3 CompleteVec(float2 xy)
{
	return float3(xy, sqrt(1.0 - saturate(dot(xy, xy))));
}

in uint gl_VertexID;
in uint gl_InstanceID;

#define VId gl_VertexID
#define IId gl_InstanceID

out VS_OUT
{
	float3 Pos; 		
	float3 LPos; 		
	float3 Normal;		
	float3 OVec1;		
} Out;


void main()
{		
	float lx = float(VId - 1) / 8.0 * 2.0 * Pi;
	
	float3 Pos;
	Pos.y = 0.0;
	Pos.x = cos(lx);
	Pos.z = sin(lx);
		
	if(VId == 0) Pos = float3(0.0);	
		
	Out.Pos = Pos.xyz;
	Out.LPos = Pos.xyz;
	Out.Normal = float3(0.0, 1.0, 0.0);
	// Out.OVec1 = ovec1;	
	Out.OVec1 = float3(1.0, 0.0, 0.0);	
	
	gl_Position = mul(float4(Pos, 1.0), ViewProjMat);
}