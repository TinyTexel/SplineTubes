
// struct VS_IN
// {
	// float4 pos : Input0;
// };

// struct PS_IN
// {
	// float4 pos : SV_POSITION;
	
	// float3 wpos : TEXCOORD0;
// };

struct VS_IN
{
	float4 pos : Input0;
};

struct HS_IN
{
	float3 pos 		: POSITION;
	float3 lPos 	: LPOS;
	float3 normal 	: NORMAL;
	float tessF : TESSFAC;
};

struct HS_OUT
{
	float3 pos 		: POSITION;
	float3 lPos 	: LPOS;
	float3 normal 	: NORMAL;
	float tessF : TESSFAC;
};

struct PCF_OUT
{
    float edges[3] : SV_TessFactor;
    float inside : SV_InsideTessFactor;
};

struct PS_IN
{
	float4 pos : SV_POSITION;
	
	float3 wpos : TEXCOORD0;
	float3 normal 	: NORMAL;
	float tessF : TESSFAC;
};

cbuffer buffer0 : register(b0)
{
	float4x4 matViewProj 	: packoffset(c0);
	float4x4 matView 		: packoffset(c4);
	float3 eyePos 			: packoffset(c8);
};

cbuffer pixelInfoCB : register(b2)
{
	float2 PixelCount : packoffset(c0);
	float2 PixelSize : packoffset(c0.z);
};

static const float pi = 3.141592653589793238462643383279502884197169399375105820974944592f;

// mip_map_level(in vec2 texture_coordinate)
// {
    // // The OpenGL Graphics System: A Specification 4.2
    // //  - chapter 3.9.11, equation 3.21
 
 
    // vec2  dx_vtc        = dFdx(texture_coordinate);
    // vec2  dy_vtc        = dFdy(texture_coordinate);
    // float delta_max_sqr = max(dot(dx_vtc, dx_vtc), dot(dy_vtc, dy_vtc));
 
 
    // //return max(0.0, 0.5 * log2(delta_max_sqr) - 1.0); // == log2(sqrt(delta_max_sqr));
    // return 0.5 * log2(delta_max_sqr); // == log2(sqrt(delta_max_sqr));
// }

float3 EvalCSpline(float3 p1, float3 t1, float3 p2, float3 t2, float l)
{
	float3 h1 = p1 + t1 * 0.33333f;
	float3 h2 = p2 - t2 * 0.33333f;
	
	float3 a1 = lerp(p1, h1, l);
	float3 a2 = lerp(h1, h2, l);
	float3 a3 = lerp(h2, p2, l);
	
	float3 b1 = lerp(a1, a2, l);
	float3 b2 = lerp(a2, a3, l);
	
	// return lerp(p1, p2, l);
	return lerp(b1, b2, l);
}

float3 EvalCSplineTangent(float3 p1, float3 t1, float3 p2, float3 t2, float l)
{
	// 
	return 	+ 6.f * p1 * (l - 1.f) * l 
			- 6.f * p2 * (l - 1.f) * l 
			+ t1 
			- 4.f * l * t1 
			+ 3.f * l * l * t1 
			- 2.f * l * t2 
			+ 3.f * l * l * t2;
}


float3 EvalCSplineBitangent(float3 p1, float3 tt1, float3 p2, float3 tt2, float l)
{
//  2 t^2 	(-3 c1z t1y +3 c2z t1y +3 c1y t1z -3 c2y t1z +3 t1z t2y -3 t1y t2z  +3 c1z t2y -3 c2z t2y -3 c1y t2z +3 c2y t2z) + 
//  2 t 	(+6 c1z t1y -6 c2z t1y -6 c1y t1z +6 c2y t1z -3 t1z t2y +3 t1y t2z) +
//  2 		(-3 c1z t1y +3 c2z t1y +3 c1y t1z -3 c2y t1z +1 t1z t2y -1 t1y t2z)  
//  
//  2 t^2 	(+3 c1z t1x -3 c2z t1x -3 c1x t1z +3 c2x t1z -3 t1z t2x +3 t1x t2z  -3 c1z t2x +3 c2z t2x +3 c1x t2z -3 c2x t2z) +
//  2 t 	(-6 c1z t1x +6 c2z t1x +6 c1x t1z -6 c2x t1z +3 t1z t2x -3 t1x t2z) +
//  2 		(+3 c1z t1x -3 c2z t1x -3 c1x t1z +3 c2x t1z -1 t1z t2x +1 t1x t2z) 
//  
//  2 t^2 	(-3 c1y t1x +3 c2y t1x +3 c1x t1y -3 c2x t1y +3 t1y t2x -3 t1x t2y  +3 c1y t2x -3 c2y t2x -3 c1x t2y +3 c2x t2y) + 
//  2 t 	(+6 c1y t1x -6 c2y t1x -6 c1x t1y +6 c2x t1y -3 t1y t2x +3 t1x t2y) +
//  2 		(-3 c1y t1x +3 c2y t1x +3 c1x t1y -3 c2x t1y +1 t1y t2x -1 t1x t2y) 
	
	float3 t1 = p1.zzy * tt1.yxx;
	float3 t2 = p2.zzy * tt1.yxx;
	float3 t3 = p1.yxx * tt1.zzy;
	float3 t4 = p2.yxx * tt1.zzy;
	float3 t5 = tt1.zzy * tt2.yxx;
	float3 t6 = tt1.yxx * tt2.zzy;

	float3 t7 = p1.zzy * tt2.yxx;
	float3 t8 = p2.zzy * tt2.yxx;
	float3 t9 = p1.yxx * tt2.zzy;
	float3 t0 = p2.yxx * tt2.zzy;
	
	float cxA = 	   -t1.x +t2.x +t3.x -t4.x    +t5.x -t6.x +t7.x -t8.x -t9.x +t0.x; 
	float cxB = 2.f * (+t1.x -t2.x -t3.x +t4.x)   -t5.x +t6.x; 
	float cxC = 	   -t1.x +t2.x +t3.x -t4.x + (+t5.x -t6.x) * 0.333333f; 
	
	float cyA = 	   +t1.y -t2.y -t3.y +t4.y    -t5.y +t6.y -t7.y +t8.y +t9.y -t0.y; 
	float cyB = 2.f * (-t1.y +t2.y +t3.y -t4.y)   +t5.y -t6.y; 
	float cyC = 	   +t1.y -t2.y -t3.y +t4.y + (-t5.y +t6.y) * 0.333333f; 
	
	float czA = 	   -t1.z +t2.z +t3.z -t4.z    +t5.z -t6.z +t7.z -t8.z -t9.z +t0.z; 
	float czB = 2.f * (+t1.z -t2.z -t3.z +t4.z)   -t5.z +t6.z; 
	float czC = 	   -t1.z +t2.z +t3.z -t4.z + (+t5.z -t6.z) * 0.333333f; 
	
	return l * (l * float3(cxA, cyA, czA) + float3(cxB, cyB, czB)) + float3(cxC, cyC, czC);
}

/* Always works if the input is non-zero.
 * Doesn’t require the input to be normalised.
 * Doesn’t normalise the output. 
 http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts */
float3 GetOrthoVec(float3 v)
{
    return abs(v.x) > abs(v.z) ? float3(-v.y, v.x, 0.0f) : float3(0.0f, -v.z, v.y);
}

float3x3 CSysFromAxis(float3 axis)
{
	float3x3 mat;
	mat[2] = axis;
	mat[1] = normalize(GetOrthoVec(axis));
	mat[0] = cross(axis, -mat[1]);
	
	return mat;
}

struct Quat
{
	float4 Val;
	
	static Quat New(float4 xyzw)
	{
		Quat quat;
		
		quat.Val = xyzw;
		
		return quat;
	}
	
	static Quat New(float w, float3 xyz = 0.f)
	{
		Quat quat;
		
		quat.Val.w = w;
		quat.Val.xyz = xyz;
		
		return quat;
	}
	
	void Conj()
	{
		Val.xyz = -Val.xyz;
	}
	
	void Normalize()
	{
		Val = normalize(Val);
	}
	
	Quat Normalized()
	{
		return Quat::New(normalize(Val));
	}
		
	static Quat New(float3 from, float3 to)
	{
		Quat quat;
		
		quat.Val.xyz = cross(from, to);
		quat.Val.w = dot(from, to) + sqrt(dot(from, from) * dot(to, to));
		// quat.Val.w = dot(from, to) + 1.f;
		
		return quat.Normalized();
	}
};

// Source: https://github.com/sharpdx/SharpDX/blob/master/Source/SharpDX.Mathematics/Vector3.cs
float3 Rotate(float3 vec, Quat rotQ)
{
	float4 rot = rotQ.Val;
	
	float x = rot.x + rot.x;
	float y = rot.y + rot.y;
	float z = rot.z + rot.z;
	float wx = rot.w * x;
	float wy = rot.w * y;
	float wz = rot.w * z;
	float xx = rot.x * x;
	float xy = rot.x * y;
	float xz = rot.x * z;
	float yy = rot.y * y;
	float yz = rot.y * z;
	float zz = rot.z * z;

	return	float3(((vec.x * ((1.0f - yy) - zz)) + (vec.y * (xy - wz))) + (vec.z * (xz + wy)),
				  ((vec.x * (xy + wz)) + (vec.y * ((1.0f - xx) - zz))) + (vec.z * (yz - wx)),
				  ((vec.x * (xz - wy)) + (vec.y * (yz + wx))) + (vec.z * ((1.0f - xx) - yy)));
}

// static const float hc0 = 12.9898;
// static const float hc1 = 27.3972;
// static const float hc2 = 33.7311;


// float Hash(float v)
// {
    // return frac(sin(v) * 43758.5453);
// }

// float Hash(float2 v)
// {
	// return Hash(v.y + v.x * hc0);
// }

// float Hash(float3 v)
// {
	// return Hash(v.y + v.x * hc0 + v.z * hc2);
	// return Hash(v.y + Hash(v.x * hc0 + Hash(v.z)) * hc2);
// }

// float2 Hash2(float2 v)
// {
	// return float2(Hash(v.y + v.x * hc0), Hash(v.x + v.y * hc2));
// }

// float2 Hash2(float3 v)
// {
	// return float2(Hash(v.y + v.x * hc0 + v.z * hc2), 
                  // Hash(v.z + v.x * hc1 + v.y * hc0));
// }

// float3 Hash3(float3 v)
// {
	// return float3(Hash(v.y + v.x * hc0 + v.z * hc2), 
                  // Hash(v.x + v.y * hc2 + v.z * hc1),
                  // Hash(v.z + v.x * hc1 + v.y * hc0));
// }

// A single iteration of Bob Jenkins' One-At-A-Time hashing algorithm.
uint Hash0(uint x) 
{
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
	
    return x;
}

// Construct a float with in interval (-1:1) using low 23 bits.
// All zeroes yields 0.0, all ones yields the next smallest representable value below 1.0.
float ConstructFloat(uint m) 
{
    // const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    // const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    // m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    // m |= ieeeOne;                          // Add fractional part to 1.0

	float flt = asfloat(m & 0x007FFFFFu | 0x3F800000u);
    // return flt - 1.f;       
        
    return m >> 31 ? flt - 1.f : flt - 2.f;             
	// return m >> 31 ? 1.f : 0.f;
}

float2 ConstructFloat(uint2 m) { return float2(ConstructFloat(m.x), ConstructFloat(m.y)); }
float3 ConstructFloat(uint3 m) { return float3(ConstructFloat(m.xy), ConstructFloat(m.z)); }
float4 ConstructFloat(uint4 m) { return float4(ConstructFloat(m.xyz), ConstructFloat(m.w)); }


uint Hash(uint v, uint r1) 								{ return Hash0(v.x ^ r1); }
uint Hash(uint2 v, uint r1, uint r2) 					{ return Hash0(Hash(v.x, r1) ^ (v.y ^ r2)); }
uint Hash(uint3 v, uint r1, uint r2, uint r3) 			{ return Hash0(Hash(v.xy, r1, r2) ^ (v.z ^ r3)); }

static const uint RndNum1 = 0xc36a1b77;
static const uint RndNum2 = 0x3fe36a19;
static const uint RndNum3 = 0x11d33555;


uint Hash(uint3 v)   { return Hash(v, RndNum1, RndNum2, RndNum3); }
float Hash(float3 v)   { return ConstructFloat(Hash(asuint(v))); }


void EvalTube(float l, float2 xy, out float3 pos, out float3 normal)
{
	float3 p1 = float3(0.f, -1.f, 0.f);
	float3 p2 = float3(4.f, 0.f, 0.f);
	
	// float3 p2 = float3(0.f, 1.f, 0.f);
	float3 t1 = float3(2.f, 4.f, 4.f)*2.f;
	float3 t2 = -float3(2.f, -4.f, 4.f)*2.f;
	
	p1 = float3(0.f, 1.f, 0.f);
	p2 = float3(0.f, -1.f, 0.f);

	t1 = float3(2.f, 4.f, 4.f)*2.f;
	t2 = -float3(2.f, -4.f, 4.f)*2.f;
	
	// t1 = float3(0.f, -1.f, 0.f);
	t2 = float3(0.f, -1.f, 0.f);
	
	float seed = 5.f;
	t1 = float3(Hash(seed * 11.333f + 0.2345f), 
				Hash(seed * 7.3242f + 0.654f),
				Hash(seed * 21.675f + 3.4567f)) * 2.f - 1.f;
	
	t2 = float3(Hash(seed * 3.456f + 5.2435f), 
				Hash(seed * 13.546 + 11.5655f),
				Hash(seed * 17.3457f + 1.8976f)) * 2.f - 1.f;	
	t1 *= 4.f;			
	t2 *= 4.f;			
	
	p1 = float3(0.f, -1.f, 0.f);
	p2 = float3(0.f, 1.f, 0.f);
	
	t1 = float3(1.f, 1.f, 0.f) * 4.f;
	t2 = float3(-1.f, 1.f, 0.f) * 4.f;
	
	p1 *= 3.f; p2 *= 3.f;
	t1 *= 3.f; t2 *= 3.f;
	
	float o = 0.001f;
	float3 s1 = EvalCSpline(p1, t1, p2, t2, l);
	float3 s2 = EvalCSpline(p1, t1, p2, t2, l + o);
	
	float3 t = (s2 - s1) / o;
	t = EvalCSplineTangent(p1, t1, p2, t2, l);
	float3 nt = normalize(t);
	
	float3 bi = normalize(EvalCSplineBitangent(p1, t1, p2, t2, l));
	
	Quat q = Quat::New(float3(0.f, 1.f, 0.f), nt);
	// Quat::New(float3(0.f, 1.f, 0.f), float3(nt.x, abs(nt.y), nt.z));
	// q = Quat::New(float4(0.f, 0.f, 0.f, 1.f));
	float3 ovec1 = Rotate(float3(1.f, 0.f, 0.f), q);
	ovec1 = bi;
	// ovec1 = float3(1.f, 0.f, 0.f);
	// float3 ovec2 = Rotate(float3(0.f, 0.f, 1.f), q);
	float3 ovec2 = cross(nt, ovec1);
	// ovec1 = cross(nt, ovec2);
	
	ovec1 = normalize(ovec1);
	ovec2 = normalize(ovec2);
	
	float3x3 mat = CSysFromAxis(nt);
	
	float r = 0.3f;
	// Out.wpos.xyz = s1;
	// Out.wpos.xyz = s1 + (mat[1] * In.pos.x + mat[0] * In.pos.z) * r;
	float e = 1.f;
	normal = mat[1] * pow(xy.x, e) + mat[0] * pow(xy.y, e);
	// normal = -ovec1 * pow(xy.x, e) + ovec2 * pow(xy.y, e);
	pos = s1 + normal * r;
	// Out.pos.xyz = s1 + (-ovec1 * pow(In.pos.x, e) + ovec2 * pow(In.pos.z, e)) * r;
}

float EvalTessFac(float3 id0, float3 id1)
{
	return lerp(1.f, 21.f, max(Hash(id0) * 0.5f + 0.5f, Hash(id1) * 0.5f + 0.5f));
	return max(saturate(Hash(id0) * 0.5f + 0.5f) * 20.f + 1.f, saturate(Hash(id1) * 0.5f + 0.5f) * 20.f + 1.f);
}

HS_IN VS( VS_IN In )
{
	HS_IN Out = (HS_IN)0;
	
	Out.pos = In.pos.xyz;
	
	
	Out.lPos = In.pos;
	EvalTube(In.pos.y, In.pos.xz, Out.pos, Out.normal);
	
	Out.tessF = saturate(Hash(Out.pos.xyz) * 0.5f + 0.5f) * 6.f + 1.f;
	// Out.tessF = saturate(Out.normal.z) * 5.f +1.f;
	
	Out.tessF = pow(1.f - saturate(dot(normalize(eyePos - Out.pos), Out.normal)), 2.f) * 8.f + 1.f;
	
	// Out.tessF = 1.f;
	// if(In.pos.y > 0.5f)
	// Out.wpos.x += 1.f;
	
	// Out.pos = mul(float4(Out.wpos.xyz, 1.f), matViewProj);
	
	return Out;
}

float2 EvalScreenSpacePosition(float3 pos)
{
	float4 pPos = mul(float4(pos, 1.f), matViewProj);
	
	return pPos.xy / pPos.w;
}

float EvalSSTessFac(float3 pos0, float3 pos1)
{
	float2 sPos0 = EvalScreenSpacePosition(pos0);
	float2 sPos1 = EvalScreenSpacePosition(pos1);
	
	float len = length((sPos1 - sPos0) * PixelCount.xy * 0.5f);
	// float len = length((sPos0 - sPos1) * PixelCount.xy * 0.5f);
	
	// return 1.f;
	// return max(3.f, len / 32.f);
	return 1.f + min(len / 64.f, 7.f);
}

PCF_OUT PCF(InputPatch<HS_IN, 3> patch, uint patchId : SV_PrimitiveID)
{    
	PCF_OUT Out;
	
	float tessellationAmount = 1.0f;
	
	// InputPatch p = patch[patchId];
	
	float3 pos[3];
	float3 normal[3];
	float eyeFac[3];
	float maxTF = 1.f;
	
	[unroll]
	for(uint i = 0; i < 3; ++i)
	{
		pos[i] = patch[i].pos;
		normal[i] = patch[i].normal;
		
		// eyeFac[i] = 1.f - saturate(dot(normalize(eyePos - pos[i]), normal[i]));
		eyeFac[i] = saturate(normal[i].z);
		eyeFac[i] = saturate(patch[i].normal.z);
		// float tf = lerp(1.f, 4.f, eyeFac[i]);
		// maxTF = max(maxTF, tf);
		
		// Out.edges[i] = tf;
	}
	

	
	for(uint i = 0; i < 3; ++i)
	{
		// float tf = lerp(3.f, 9.f, max(eyeFac[i], eyeFac[(i + 1) % 3]));
		float tf = lerp(1.f, 4.f, max(eyeFac[(i + 1) % 3], eyeFac[(i + 2) % 3]));
		// maxTF = max(maxTF, tf);
		
		Out.edges[i] = tf;
	}
	
	// Out.edges[0] = (max(saturate(patch[1].normal.z), saturate(patch[2].normal.z)))*5+1;
	// Out.edges[1] = (max(saturate(patch[2].normal.z), saturate(patch[0].normal.z)))*5+1;
	// Out.edges[2] = (max(saturate(patch[0].normal.z), saturate(patch[1].normal.z)))*5+1;
	
	// Out.edges[0] = lerp(1.f, 5.f, max(saturate(patch[1].normal.z), saturate(patch[2].normal.z)));
	// Out.edges[1] = lerp(1.f, 5.f, max(saturate(patch[2].normal.z), saturate(patch[0].normal.z)));
	// Out.edges[2] = lerp(1.f, 5.f, max(saturate(patch[0].normal.z), saturate(patch[1].normal.z)));
	
	// eyeFac[0] = saturate(patch[0].normal.z);
	// eyeFac[1] = saturate(patch[1].normal.z);
	// eyeFac[2] = saturate(patch[2].normal.z);
	
	Out.edges[0] = lerp(1.f, 5.f, max(eyeFac[1], eyeFac[2]));
	Out.edges[1] = lerp(1.f, 5.f, max(eyeFac[2], eyeFac[0]));
	Out.edges[2] = lerp(1.f, 5.f, max(eyeFac[0], eyeFac[1]));
	
	// float eyeFac0 = saturate(patch[0].normal.z);
	// float eyeFac1 = saturate(patch[1].normal.z);
	// float eyeFac2 = saturate(patch[2].normal.z);
	
	// Out.edges[0] = lerp(1.f, 5.f, max(eyeFac1, eyeFac2));
	// Out.edges[1] = lerp(1.f, 5.f, max(eyeFac2, eyeFac0));
	// Out.edges[2] = lerp(1.f, 5.f, max(eyeFac0, eyeFac1));
	
	// Out.edges[0] = ((saturate(patch[1].normal.z) * 5.f +1.f) + (saturate(patch[2].normal.z) * 5.f +1.f)) * 0.5f;
	// Out.edges[1] = ((saturate(patch[2].normal.z) * 5.f +1.f) + (saturate(patch[0].normal.z) * 5.f +1.f)) * 0.5f;
	// Out.edges[2] = ((saturate(patch[0].normal.z) * 5.f +1.f) + (saturate(patch[1].normal.z) * 5.f +1.f)) * 0.5f;
	
	// Out.edges[0] = EvalTessFac(patch[1].pos, patch[2].pos);	
	// Out.edges[1] = EvalTessFac(patch[2].pos, patch[0].pos);	
	// Out.edges[2] = EvalTessFac(patch[0].pos, patch[1].pos);	
	
	// Out.edges[0] = max(patch[1].tessF, patch[2].tessF);	
	// Out.edges[1] = max(patch[2].tessF, patch[0].tessF);	
	// Out.edges[2] = max(patch[0].tessF, patch[1].tessF);
	
	Out.edges[0] = (patch[1].tessF + patch[2].tessF) * 0.5f;	
	Out.edges[1] = (patch[2].tessF + patch[0].tessF) * 0.5f;	
	Out.edges[2] = (patch[0].tessF + patch[1].tessF) * 0.5f;
	
	// Out.edges[0] = EvalSSTessFac(patch[1].pos, patch[2].pos);
	// Out.edges[1] = EvalSSTessFac(patch[2].pos, patch[0].pos);
	// Out.edges[2] = EvalSSTessFac(patch[0].pos, patch[1].pos);
	
	Out.edges[0] = max(Out.edges[0], EvalSSTessFac(patch[1].pos, patch[2].pos));
	Out.edges[1] = max(Out.edges[1], EvalSSTessFac(patch[2].pos, patch[0].pos));
	Out.edges[2] = max(Out.edges[2], EvalSSTessFac(patch[0].pos, patch[1].pos));
	
	// float3 normalF = normalize(cross(pos[1] - pos[0], pos[2] - pos[0]));
	
	// tessellationAmount = lerp(1.0f, 4.f, saturate(normalF.z));	
	

    // Out.edges[0] = tessellationAmount;
    // Out.edges[1] = tessellationAmount;
    // Out.edges[2] = tessellationAmount;

    // Out.inside = tessellationAmount;
    Out.inside = max(Out.edges[0], max(Out.edges[1], Out.edges[2]));
    Out.inside = (Out.edges[0] + Out.edges[1] + Out.edges[2]) * 0.33333f;
	// Out.edges[0] = Out.edges[1] = Out.edges[2] = 1.f;
    // Out.inside = max(3.0f, Out.inside);
    // Out.inside = 3.0f;

    return Out;
}

[domain("tri")]
[partitioning("fractional_odd")]
[outputtopology("triangle_cw")]
[outputcontrolpoints(3)]
[patchconstantfunc("PCF")]
// [maxtessfactor(3)]
HS_OUT HS(InputPatch<HS_IN, 3> patch, uint pointId : SV_OutputControlPointID, uint patchId : SV_PrimitiveID)
{
    HS_OUT Out;

    // Set the position for this control point as the output position.
    Out.pos = patch[pointId].pos;
    Out.lPos = patch[pointId].lPos;
    Out.normal = patch[pointId].normal;
    Out.tessF = patch[pointId].tessF;

    return Out;
}

[domain("tri")]
PS_IN DS(PCF_OUT In, float3 uvw : SV_DomainLocation, const OutputPatch<HS_OUT, 3> patch)
{
    PS_IN Out;

    // Determine the position of the new vertex.
    Out.wpos = uvw.x * patch[0].pos + uvw.y * patch[1].pos + uvw.z * patch[2].pos;
    float3 lPos = uvw.x * patch[0].lPos + uvw.y * patch[1].lPos + uvw.z * patch[2].lPos;
	
	EvalTube(lPos.y, normalize(lPos.xz), Out.wpos, Out.normal);
	
	Out.tessF = uvw.x * patch[0].tessF + uvw.y * patch[1].tessF + uvw.z * patch[2].tessF;
	// vPos.xyz = normalize(vPos.xyz);
	
    // Calculate the position of the new vertex against the world, view, and projection matrices.
    // Out.pos = mul(float4(vPos, 1.0f), matViewProj);
	Out.pos = mul(float4(Out.wpos.xyz, 1.f), matViewProj);

    return Out;
}

float3 PS( PS_IN In ) : SV_Target
{	
	float3 hnormals0 = normalize(-cross(ddx(In.wpos), ddy(In.wpos)));
	
	float3 eyeVec = normalize(eyePos - In.wpos);
	
	float br = (In.tessF - 1.f) / 20.f;

	// return 1.f - dot(eyeVec, hnormals0);
	// return 1.f - mul(matView, hnormals0).z;
	
	// float3 absPos = abs(In.wpos.xyz);
	// float3 signPos = sign(In.wpos.xyz);
	// float3 hnormals = float3(0.f, signPos.y, 0.f);
	// float3 uvw = In.wpos.xyz * 0.5f + 0.5f;
	// float2 uv = uvw.xz;
	
	// if(absPos.x > absPos.y)
	// {
		// if(absPos.x > absPos.z)
		// {
			// hnormals = float3(signPos.x, 0.f, 0.f);
			// uv.xy = uvw.yz;
		// }
		// else
		// {
			// hnormals = float3(0.f, 0.f, signPos.z);
			// uv.xy = uvw.xy;
		// }
	// }
	// else
	// {
		// if(absPos.z > absPos.y)
		// {
			// hnormals = float3(0.f, 0.f, signPos.z);
			// uv.xy = uvw.xy;
		// }
	// }
	
	// return 1.f;
	// return float4(uv.xy, 0, 0);
	// return normalize(In.normal.xyz);// * 0.5f + 0.5f;
	return hnormals0.xyz;// * 0.5f + 0.5f;
	// return float4(1.0f, 0.0f, 0.0f, 1.0f);
}
