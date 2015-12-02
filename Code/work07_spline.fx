
struct VS_IN
{
	float2 pos : Input0;
};

struct PS_IN
{
	float4 pos : SV_POSITION;

	float2 uv : TEXCOORD0;
};

cbuffer pixelInfoCB : register(b0)
{
	float2 PixelCount : packoffset(c0);
	float2 PixelSize : packoffset(c0.z);
};

cbuffer generalInfoCB : register(b1)
{
	float2 CamPos : packoffset(c0);
	float Zoom : packoffset(c0.z);
	float Time : packoffset(c0.w);
	float3 Mouse : packoffset(c1);
};

PS_IN VS( VS_IN In )
{
	PS_IN Out = (PS_IN)0;
	
	Out.pos = float4(In.pos, 0, 1.f);

	Out.uv = In.pos * float2(0.5f, -0.5f) + 0.5f;
	
	Out.uv *= PixelCount;
	
	return Out;
}


#include<Shaders/const.fx>
#include<Shaders/utils.fx>

#include<Shaders/complex.fx>

#include<Shaders/hash.fx>
#include<Shaders/noise.fx>
#include<Shaders/q2.fx>
#include<Shaders/q2_noise.fx>
#include<Shaders/meta.fx>
#include<Shaders/wei_du_ke.fx>

#include<Shaders/color.fx>
#include<Shaders/draw.fx>

float EvalChecker(float2 pos)
{
	float2 iPos = floor(pos);
	
	float2 stripes = frac(iPos * 0.5f) * 2.f;
	float checker = (1.f - stripes.x) * (1.f - stripes.y) + stripes.x * stripes.y;
	
	return saturate(checker);
}

float EvalChecker(float3 pos)
{
	float3 iPos = floor(pos);
	
	float3 stripes = frac(iPos * 0.5f) * 2.f;
	float checker = saturate((1.f - stripes.x) * (1.f - stripes.y) + stripes.x * stripes.y);
	
	if(stripes.z) checker = 1.f - checker;
	
	return checker;
}

float3 AngToDir(float2 ang)
{
	float sinPhi, cosPhi; 
	float sinTheta, cosTheta; 
	
	sincos(ang.x, sinPhi, cosPhi);
	sincos(ang.y, sinTheta, cosTheta);
	
	return float3(cosPhi * sinTheta, cosTheta, sinPhi * sinTheta);
	// return float3(sinPhi * cosTheta, sinTheta, -cosPhi * cosTheta);
}

float3 VecToAng(float3 vec)
{
	float len = length(vec);
	
	if(len < 0.0001f) return 0.f;
	
	float3 dir = vec / len;
	float2 dir2 = normalize(vec.xz);
	
	return float3(atan2(dir2.y, dir2.x), acos(dir.y), len);
}

void AngToSys(float2 ang, out float3 front, out float3 right, out float3 top)
{
	float2 scPhi; sincos(ang.x, scPhi.x, scPhi.y);
	float2 scTheta; sincos(ang.y, scTheta.x, scTheta.y);
	
	 front = float3(scPhi.y * scTheta.x, scTheta.y, scPhi.x * scTheta.x);
	 right = float3(scPhi.x, 0.f, -scPhi.y);
	 top = cross(front, right);
	// return 
}



float2 RotLeft(float2 v)
{
	return float2(-v.y, v.x);
}

float2 RotRight(float2 v)
{
	return float2(v.y, -v.x);
}


class IntrsectRes
{
	float mask;
	float t1, t2;

	static IntrsectRes New()
	{
		IntrsectRes res;
		res.mask = 0.f;
		res.t1 = res.t2 = 9999999999.f;
		
		return res;
	}
};

struct Intrsect
{
	float2 LineLine(float2 line1Start, float2 line1Offset, float2 line2Start, float2 line2Offset)
	{
		float div = (line2Offset.x * line1Offset.y) - (line1Offset.x * line2Offset.y);

		// if (abs(div) < 0.00001f)
			// return 0.f;

		div = 1.f / div;

		float n1 = (line1Offset.y * line1Start.x) - (line1Offset.x * line1Start.y);
		float n2 = (line2Offset.x * line2Start.y) - (line2Offset.y * line2Start.x);

		return float2
		(
			((line2Offset.x * n1) + (line1Offset.x * n2)) * div,
			((line2Offset.y * n1) + (line1Offset.y * n2)) * div
		); 
	}
	
	static IntrsectRes RaySphere(float3 rayStart, float3 rayDir, float3 sphPos, float sphRadius)
	{
		rayStart -= sphPos;
		
		float a = dot(rayDir, rayDir);
		
		float b = 2.f * dot(rayStart, rayDir);
		
		float c = dot(rayStart, rayStart) - (sphRadius * sphRadius);
		
		float discr = b * b - 4.f * a * c;
		
		IntrsectRes res = IntrsectRes::New();
		// [branch]
		if(discr < 0.f) return res;
		
		float d = sqrt(discr);
		
		float t1 = (-b + d) / (2.f * a);
		float t2 = (-b - d) / (2.f * a);
		
		res.mask = 1.f;
		res.t1 = min(t1, t2);
		res.t2 = max(t1, t2);
		
		return res;
	}
};

float2 SampleDisk(float2 h)
{
	float sinTheta, cosTheta; 
		
	sincos(Pi * h.x, sinTheta, cosTheta);
 
    return sqrt(h.y * 0.5f + 0.5f) * float2(cosTheta, sinTheta);
}


struct Ray
{
	float3 Start;
	float3 Dir;
};

struct Cam
{
	float3 Front, Right, Up, Pos;
	float2 SensorSize;
	float AxeLen;
	
	static Cam New(float2 ang, float fov, float2 pxCount, float sensorWidth = 1.f)
	{
		Cam cam;
		
		float sinPhi, cosPhi; 
		float sinTheta, cosTheta; 
		
		sincos(ang.x, sinPhi, cosPhi);
		sincos(ang.y, sinTheta, cosTheta);
		
		float3 back  = float3(sinPhi * cosTheta, sinTheta, -cosPhi * cosTheta);
		
		cam.Front = -back;
		cam.Right = float3(cosPhi, 0.f, sinPhi);
		cam.Up   = cross(back, cam.Right);
		
		cam.Pos = 0.f;
		cam.SensorSize = float2(sensorWidth, sensorWidth * pxCount.y / pxCount.x);
		cam.DetAxeLen(fov);
		
		return cam;
	}
	
	void DetAxeLen(float fov)
	{
		this.AxeLen = tan(Pi05 - fov * 0.5f) * SensorSize * 0.5f;
	}
	
	Ray CreateRay(float2 uv, float2 pxCount)
	{
		Ray ray;
		
		uv -= pxCount * 0.5f;
		uv += 0.5f;
		uv /= pxCount.x;
		uv *= this.SensorSize.x;
		
		ray.Start = this.Pos;
		ray.Dir = normalize(this.Front * this.AxeLen + (this.Right * uv.x + this.Up * uv.y));
	
		return ray;
	}
	
	Ray CreateRay(float2 uv, float rayId, float time, float2 pxCount)
	{
		Ray ray;
		
		uv -= pxCount * 0.5f;
		uv += 0.5f;
		uv /= pxCount.x;
		uv *= this.SensorSize.x;
		
		float4 h = float4(uv, rayId, time);
		float2 lensPos = SampleDisk(float2(HashV1(h), HashV2(h))) * 0.9f;
		float3 lensOff = this.Right * lensPos.x + this.Up * lensPos.y;
		
		ray.Start = this.Pos;
		ray.Start += lensOff;
		float m = 16.f;
		ray.Dir = normalize((this.Front * this.AxeLen + (this.Right * uv.x + this.Up * uv.y)) * m - lensOff);		

		return ray;
	}
};


/* Always works if the input is non-zero.
 * Doesn’t require the input to be normalised.
 * Doesn’t normalise the output. 
 http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts */
float3 GetOrthoVec(float3 v)
{
    return abs(v.x) > abs(v.z) ? float3(-v.y, v.x, 0.0f) : float3(0.0f, -v.z, v.y);
}


float3 PS( PS_IN In ) : SV_Target
{	
	IntrsectRes resi;
	
	float2 uv = floor(In.uv.xy);
	float2 fuv = uv * PixelSize.xy;

	//region CGIII Ü-2.1.b
	if(0)
	{
		float x = fuv.x;
		float y = 1.f - fuv.y;
		
		float f = (3.f / 4.f) * (x*x * y + 4.f * x * y*y + 2.f * x * y);
		float fm = (3.f / 4.f) * (1.f + 4.f + 2.f);

		OFUNC(sampleGen1, INull,
		float2 Eval(float2 r)
		{
			return r;
		})
		
		OFUNC(sampleGen2, INull,
		float2 Eval(float2 r)
		{
			float x = sqrt(24.f * r.x + 49.f - 7.f) / 3.f;
			
			float y = sqrt(6.f * r.y + 1.f - 1.f) / 3.f;
			
			return float2(x, y);
		})
		
		OFUNC(sampleGen3, INull,
		float2 Eval(float2 r)
		{
			//   49 (1 + I Sqrt[3]))/(343 - 108 r1 + 6 Sqrt[6] Sqrt[r1 (-343 + 54 r1)])^(1/3) + 
			//	  (1 - I Sqrt[3])  (343 - 108 r1 + 6 Sqrt[6] Sqrt[r1 (-343 + 54 r1)])^(1/3))
				  
			Cmplx tx = Pow((Pow(cmplx(r.x * (-343.f + 54.f * r.x)), 0.5f).Mul(6.f * sqrt(6.f))).Sub(108.f * r.x).Add(343.f), 1.f / 3.f);
			// 1/3 (-7 - 49/(343 - 108 r1 + 6 Sqrt[6] Sqrt[r1 (-343 + 54 r1)])^(1/3)
			//            - (343 - 108 r1 + 6 Sqrt[6] Sqrt[r1 (-343 + 54 r1)])^(1/3))
			float x = (-14.f + ((cmplx(1.f, sqrt(3.f)).Mul(49.f).Div(tx)).Add(cmplx(1.f, -sqrt(3.f)).Mul(tx)).R())) / 6.f;
			
			
			// 1/32 (-4 (2 + x) + 
			// (2 (1 + I Sqrt[3]) (2 + x)^2)/
			// (8 + 12 x + 6 x^2 + x^3 - 32 r2 (14 + 3 x) + 8 Sqrt[r2 (14 + 3 x) (-(2 + x)^3 + 16 r2 (14 + 3 x))])^(1/3) + 
			// 2 (1 - I Sqrt[3]) *
			// (8 + 12 x + 6 x^2 + x^3 - 32 r2 (14 + 3 x) + 8 Sqrt[r2 (14 + 3 x) (-(2 + x)^3 + 16 r2 (14 + 3 x))])^(1/3))
	  
	  		Cmplx ty = Pow( Pow(cmplx(r.y * (14 + 3 * x) * (-Pow3(2 + x) + 16 * r.y * (14 + 3 * x))), 0.5f).Mul(8.f).Add(
							8 + 12 * x + 6 * Pow2(x) + Pow3(x) - 32 * r.y * (14 + 3 * x)), 1.f / 3.f);
			// float ty = pow(8 + 12 * x + 6 * Pow2(x) + Pow3(x) - 32 * r.y * (14 + 3 * x) + 
					   // 8 * sqrt(r.y * (14 + 3 * x) * (-Pow3(2 + x) + 16 * r.y * (14 + 3 * x))), 1.f / 3.f);
			
			float y = (-4.f * (2.f + x) + (cmplx(1.f, sqrt(3.f)).Mul(2.f * Pow2(2.f + x)).Div(ty).Add(cmplx(1.f, -sqrt(3.f)).Mul(2.f).Mul(ty))).R()) / 32.f;
			
			return float2(x, y);
		})
		
		float dots = 0.f;
		
		for(float i = 0.f; i < 50.f; ++i)
		{
			dots = max(dots, Draw::Dot(float2(x,y), sampleGen3.Eval(float2(HashV1(i), HashV2(i)) * 0.5f + 0.5f), 0.002f));
			// dots = max(dots, Draw::Dot(fuv, float2(sampleGen3.Eval(float2(HashV1(0.f), HashV2(i)) * 0.5f + 0.5f).x, 0.5f), 0.002f));
			// dots = max(dots, Draw::Dot(fuv, float2(sampleGen3.Eval(float2(0.9f, 0.5f)).y, 0.5f), 0.002f));
		}
		
		//return Draw::Line(float2(x,y) * 2.f - 1.f, 0.f, Pow(cmplx(float2(-1.f, -1.f)), 0.99f).RI(), 0.01f);
		return lerp(f / fm, float3(1.f, 0.f, 0.f), dots);
		//return Draw::Isolines(f, 3.f);
		return f / fm;
		return float3(fuv, 0.f);
	}
	//endregion
	
	// Simplex stuffs
	{
		float2 pos0 = floor(In.uv.xy);
		float2 pos = ((pos0 - PixelCount * 0.5f) * Zoom + PixelCount * 0.5f - CamPos);
	// float2 pos = (pos0 - camPos) * zoom;
	
		// return Simplex(pos);// * 0.5 + 0.5;
	}
	
	float2 uv2 = fuv - 0.5f;
	
	float2 cAng = float2(0.f, Pi * 0.1f);
	
	cAng.x = -Mouse.x * 0.01;
	cAng.y = Mouse.y * 0.01;

	// float sinPhi, cosPhi; 
	// float sinTheta, cosTheta; 
	
	// sincos(cAng.x, sinPhi, cosPhi);
	// sincos(cAng.y, sinTheta, cosTheta);
	
	// float3 back  = float3(sinPhi * cosTheta, sinTheta, -cosPhi * cosTheta);
	// float3 right = float3(cosPhi, 0.f, sinPhi);
	// float3 top   = cross(back, right);
	
	Cam cam = Cam::New(cAng, Pi * 0.25f, PixelCount);	
	
	cam.Pos = float3(0.f, 0.f, 0.f) - cam.Front * (1.f + Mouse.z * 0.8);

	
	// Ray ray = Ray::New(cam, uv, PixelCount);
	
	
	// float aspect = PixelCount.x * PixelSize.y;
	// ray.Dir = normalize(cam.Front + 1.5f * (cam.Right * uv2.x * aspect + cam.Up * uv2.y));
	// float dist = rcp(abs(ray.y)) * abs(cCenter.y);
	// float3 hit = cCenter + ray * dist;
	// float mask = abs(hit.y) < 0.1 &&  SqrLen(hit.xz) < Pow2(10);
	// dist += 999999.f * (1.f - mask);
	
	IntrsectRes res = IntrsectRes::New();
	float radius = 2.f;
	
	float rayCount = 1.f;
	float3 fCol = 0.f;
	float time = 0.f;
	for(float i = 0.f; i < rayCount; ++i)
	{
		float tOff = (Hash(float4(i, time, uv)) * 0.5f + 0.5f) * 0;
		// float3 sphPos = float3(0.f, sin((Time + tOff) * 0.02f) * 1.f, 0.f);
		float3 sphPos = float3(0.f, ((time + tOff) * 0.0f) * 1.f, 0.f);
		sphPos.y -= 2.f;
		
		Ray ray = cam.CreateRay(uv, PixelCount);
		// Ray ray = cam.CreateRay(uv, i, time, PixelCount);
		
		float t0 = 1000000.f;
		// if(0)
		{
			float3x3 rayMat;
			rayMat[0] = ray.Dir;
			rayMat[1] = normalize(GetOrthoVec(ray.Dir));
			rayMat[2] = cross(ray.Dir, rayMat[1]);
			

			
			float3 s = float3(2.f, -3.f, 0.f); 
			float3 d = float3(1.f, 3.f, 1.f); 
			float r = 0.75f;
			
			// rayMat = transpose(rayMat);
			
			s -= ray.Start; 
			s = mul(rayMat, s);
			d = mul(rayMat, d);
			
			float3 s2 = s * s;
			float3 d2 = d * d;
			float r2 = r * r;
			
			float discr = d2.y * r2 + d2.z * r2 - d2.z * s2.y + 2.f * d.y * d.z * s.y * s.z - d2.y * s2.z;
			// discr = 0.f;
			float2 l2 = (-d.y * s.y - d.z * s.z + sqrt(discr) * float2(-1.f, 1.f)) / (d2.y + d2.z);
			
			float t = 1000000.f;
			float l = 0.f;
			
			if(0)
			{
				float count = 10.f;
				float rCount = rcp(count);
				
				float discr2 = 0;
				
				for(float j = 0; j <= count; ++j)
				{
					float t_l = lerp(l2.y, l2.x, j * rCount);
					
					discr2 = r2 - Pow2(d.y * t_l + s.y) - Pow2(d.z * t_l + s.z);
					discr2 = max(0.f, discr2);
					float t_t = d.x * t_l + s.x - sqrt(discr2);
					
					if(t_t < t)
					{
						l = t_l;
						t = t_t;
					}
				}
			}
			
			{
				float count = 9.f;
				float rCount = rcp(count);
				
				float discr2 = 0;
				
				// l2 = float2(min(l2.x, l2.y), max(l2.x, l2.y));
				
				for(float j = 0; j < count; ++j)
				{
					float t_l = (l2.x + l2.y) * 0.5f;
					
					discr2 = r2 - Pow2(d.y * t_l + s.y) - Pow2(d.z * t_l + s.z);
					discr2 = max(0.f, discr2);
					float t_t = d.x - (-2.f * d.y * (d.y * t_l + s.y) - 2.f * d.z * (d.z * t_l + s.z)) * 0.5f * rsqrt(discr2);
					
					if(t_t > 0)
					{
						l2.y = t_l;
					}
					else
					{
						l2.x = t_l;
					}
				}
				
				l = (l2.x + l2.y) * 0.5f;
				
				discr2 = r2 - Pow2(d.y * l + s.y) - Pow2(d.z * l + s.z);
				discr2 = max(0.f, discr2);
				t = d.x * l + s.x - sqrt(discr2);
			}
			
			if(0)
			{
				OFUNC(CylThing, INull,
				float2 Eval(float3 s, float3 d, float r)
				{
					float3 s2 = s * s;
					float3 d2 = d * d;
					float  r2 = r * r;
					
					float discr = d2.x * dot(d, d) * (d2.z * (r2 - s2.y) + 2.f * d.y * d.z * s.y * s.z + d2.y * (r2 - s2.z));
					
					float fTerm = (d2.y * d.y * s.y) + (d.y * d2.z * s.y) + (d2.y * d.z * s.z) + (d2.z * d.z * s.z) + d2.x * (d.y * s.y + d.z * s.z);
					
					float divTerm =  (d2.y + d2.z) * (d2.x + d2.y + d2.z);
					
					return -(fTerm + float2(-sqrt(discr), sqrt(discr))) / divTerm;
				})
			
				float2 t_l = CylThing.Eval(s, d, r);
				float2 discr2 = r2 - Pow2(d.y * t_l + s.y) - Pow2(d.z * t_l + s.z);
				discr2 = max(0.f, discr2);
				float2 t_t = d.x * t_l + s.x - sqrt(discr2);
				
				if(t_t.x < t_t.y)
				{
					l = t_l.x;
					t = t_t.x;
				}
				else
				{
					l = t_l.y;
					t = t_t.y;
				}
			}
				
			float3 cp = s + d * l;
			float3 ip = ray.Start + ray.Dir * t;
				
			cp = mul(cp, rayMat); cp += ray.Start;
			// ip = mul(ip, rayMat); ip += ray.Start;
			
			float3 n = normalize(ip - cp);
			// n = mul(rayMat, n);
			// n = mul(n, rayMat);
			// return t > 0;
			// return l > 0 && l < 1;
			// return ray.Start;
			// if(discr >= 0.f) fCol = 1.f;
			if(discr >= 0.f && t > 0 && l > 0 && l < 1) 
			{
				fCol = n;// * 0.5f + 0.5f;
			
				t0 = t;
			}
			// return fCol;
		}
		
				// if(0)
		{
			float3x3 rayMat;
			rayMat[0] = ray.Dir;
			rayMat[1] = normalize(GetOrthoVec(ray.Dir));
			rayMat[2] = cross(ray.Dir, rayMat[1]);
			

			float3 c1 = float3(-1.f, -2.f, 0.f); 
			float3 c2 = float3(3.f, -2.f, 0.f); 
			float3 h = float3(2.f, sin(Time * 0.00f) * 8.f + 8.f, 0.f); 
			// float3 c1 = float3(-1.f, 0.f, 0.f); 
			// float3 c2 = float3(3.f, 0.f, 0.f); 
			// float3 h = float3(2.5f, 3.f, 2.f); 
			float r = 0.5f;
			
			// rayMat = transpose(rayMat);
			
			c1 -= ray.Start; c1 = mul(rayMat, c1);
			c2 -= ray.Start; c2 = mul(rayMat, c2);
			h -= ray.Start; h = mul(rayMat, h);
			
			float3 c1_2 = c1 * c1;
			float3 c2_2 = c2 * c2;
			float3 h2 = h * h;
			float r2 = r * r;
			
			OFUNC(qSpline, INull,
			float3 Eval(float l)
			{
				return lerp(lerp(c1, h, l), lerp(h, c2, l), l);
			})
			
			OFUNC(discrPoly2, INull,
			float Eval(float a, float b, float c)
			{
				return b * b - 4.f * a * c;
			})
			
			OFUNC(discrPoly4, INull,
			float Eval(float a, float b, float c, float d, float e)
			{
				float a2 = a * a; float a3 = a2 * a;
				float b2 = b * b; float b3 = b2 * b;
				float c2 = c * c; float c3 = c2 * c;
				float d2 = d * d; float d3 = d2 * d;
				float e2 = e * e; float e3 = e2 * e;
				
				return 	256.f * a3 * e3 - 192.f * a2 * b * d * e2 - 128.f * a2 * c2 * e2 + 144.f * a2 * c * d2 * e - 27.f * a2 * d2 * d2
						+ 144.f * a * b2 * c * e2 - 6.f * a * b2 * d2 * e - 80.f * a * b * c2 * d * e + 18.f * a * b * c * d3 + 16.f * a * c2 * c2 * e
						- 4.f * a * c3 * d2 - 27.f * b2 * b2 * e2 + 18.f * b3 * c * d * e - 4.f * b3 * d3 - 4.f * b2 * c3 * e + b2 * c2 * d2;
			})
	
			OFUNC(hasRealRootsPoly4, INull,
			float Eval(float a, float b, float c, float d, float e)
			{
				float P = 8.f * a * c - 3.f * b * b;
				float D = 64.f * a*a*a * e - 16.f * a*a * c*c + 16.f * a * b*b * c - 16.f * a*a * b * d - 3.f * b*b*b*b;
				
				float discr = discrPoly4.Eval(a, b, c, d, e);
				
				return discr < 0.f || (P < 0.f && D < 0.f);// && p < 0.f;
			})
			
			// -c1y^2 - c1z^2 + r^2
			float discr_e = -c1_2.y - c1_2.z + r2;
			// (4 c1y^2 + 4 c1z^2 - 4 c1y hy - 4 c1z hz) t
			float discr_d = 4.f * c1_2.y + 4.f * c1_2.z - 4.f * c1.y * h.y - 4.f * c1.z * h.z;
			// (-6 c1y^2 - 6 c1z^2 - 2 c1y c2y - 2 c1z c2z + 12 c1y hy - 4 hy^2 + 12 c1z hz - 4 hz^2) t^2
			float discr_c = -6.f * c1_2.y - 6.f * c1_2.z - 2.f * c1.y * c2.y - 2.f * c1.z * c2.z + 12.f * c1.y * h.y - 4.f * h2.y + 12.f * c1.z * h.z - 4.f * h2.z;
			// (4 c1y^2 + 4 c1z^2 + 4 c1y c2y + 4 c1z c2z - 12 c1y hy - 4 c2y hy + 8 hy^2 - 12 c1z hz - 4 c2z hz + 8 hz^2) t^3
			float discr_b = 4.f * c1_2.y + 4.f * c1_2.z + 4.f * c1.y * c2.y + 4.f * c1.z * c2.z - 12.f * c1.y * h.y - 4.f * c2.y * h.y + 8.f * h2.y - 12.f * c1.z * h.z - 4.f * c2.z * h.z + 8.f * h2.z;
			// (-c1y^2 - c1z^2 - 2 c1y c2y - c2y^2 - 2 c1z c2z - c2z^2 + 4 c1y hy + 4 c2y hy - 4 hy^2 + 4 c1z hz + 4 c2z hz - 4 hz^2) t^4
			float discr_a = -c1_2.y - c1_2.z - 2.f * c1.y * c2.y - c2_2.y - 2.f * c1.z * c2.z - c2_2.z + 4.f * c1.y * h.y + 4.f * c2.y * h.y - 4.f * h2.y + 4.f * c1.z * h.z + 4.f * c2.z * h.z - 4.f * h2.z;
			
			
			OFUNC(discrSpline, IEVAL(float, (float), 1),
			float Eval(float l)
			{				
				float l1 = l - 1.f;
				float2 term = Pow2(l1 * l1 * c1.yz + l * (l * c2.yz - 2.f * l1 * h.yz));
				
				return r2 - term.x - term.y;
			})
			
			OFUNC(qSplineIDist, IEVAL(float, (float), 1),
			float Eval(float l)
			{
				float term = l*l * (c1.x + c2.x - 2.f * h.x) + l * (-2.f * c1.x + 2.f * h.x) + c1.x;
				
				float discr = discrSpline.Eval(l);
				// discr = max(0.f, discr);
				
				if(discr < 0.f)
				return 999999.f;
				else
				return term - sqrt(discr);
			})
			
			// OFUNC(qSplineD1, IEVAL(float, (float), 1),
			// float Eval(float t)
			// {
			// -2 c1[[1]] + 2 t c1[[1]] + 2 t c2[[1]] + 2 h[[1]] - 4 t h[[1]] - 
			// ( 
				// (
				// -4 
				// ((-1 + t)   c1[[2]] + t    c2[[2]] +  (1 - 2 t) h[[2]]) 
				// ((-1 + t)^2 c1[[2]] + t (t c2[[2]] - 2 (-1 + t) h[[2]])) - 
				// 4 
				// ((-1 + t)   c1[[3]] + t    c2[[3]] +  (1 - 2 t) h[[3]]) 
				// ((-1 + t)^2 c1[[3]] + t (t c2[[3]] - 2 (-1 + t) h[[3]]))
				// )
			// )
			// /(2 sqrt(discr))
		 
			// })
			
			OFUNC(discrSplineNorm, IEVAL(Cmplx, (Cmplx), 1),
			Cmplx Eval(Cmplx l)
			{				
				Cmplx l1 = l.Sub(1.f);
				Cmplx term1 = l1.Mul(l1).Mul(c1.y).Add(l.Mul(l.Mul(c2.y).Sub(l1.Mul(h.y).Mul(2.f)))); term1 = term1.Mul(term1);
				Cmplx term2 = l1.Mul(l1).Mul(c1.z).Add(l.Mul(l.Mul(c2.z).Sub(l1.Mul(h.z).Mul(2.f)))); term2 = term2.Mul(term2);
				
				return (cmplx(r2).Sub(term1).Sub(term2)).Div(discr_a);
			})
			
			OFUNC(discrSplineD1, INull,
			float Eval(float l)
			{				
				float l1 = l - 1.f;
				// ((-1 + t) c1[[2]] + t c2[[2]] + (1 - 2 t) h[[2]]) ((-1 + t)^2 c1[[2]] + t (t c2[[2]] - 2 (-1 + t) h[[2]]))
				float2 term = -4.f * (l1 * c1.yz * l * c2.yz + (1.f - 2.f * l) * h.yz) 
								   * (l1 * l1 * c1.yz + l * (l * c2.yz - 2.f * l1 * h.yz));
		
				return term.x + term.y;
			})
			
			float isDpos = 0.f;
			
			//region linear roots
			{
				float count = 20.f;
				float rCount = rcp(count);
				
				for(float j = 0; j <= count; ++j)
				// isDpos += discrSplineNorm.Eval(cmplx(j * rCount)).r > 0.f;
				isDpos += discrSpline.Eval(j * rCount) > 0.f;
			}
			//endregion
			
			//region Durand-Kerner roots
			isDpos = 0.f;
			if(0)
			{
				float mi = 0.1f;
				
				// Cmplx root1 = cmplx(1.f, 0.f);
				// Cmplx root2 = cmplx(0.4f, 0.9f);
				// Cmplx root3 = root2.Mul(root2);
				// Cmplx root4 = root3.Mul(root2);
				
				Cmplx root1 = cmplx(0.2f, mi);
				Cmplx root2 = cmplx(0.4f, mi);
				Cmplx root3 = cmplx(0.6f, mi);
				Cmplx root4 = cmplx(0.8f, mi);

				EvalWeiDuKe(discrSplineNorm, 10, root1, root2, root3, root4);
				
				float err = 0.001f;
				
				if(root1.r > 0.f && root1.r < 1.f && abs(discrSpline.Eval(root1.r)) < err && abs(root1.i) < err) isDpos = 1.f;
				if(root2.r > 0.f && root2.r < 1.f && abs(discrSpline.Eval(root2.r)) < err && abs(root2.i) < err) isDpos = 1.f;
				if(root3.r > 0.f && root3.r < 1.f && abs(discrSpline.Eval(root3.r)) < err && abs(root3.i) < err) isDpos = 1.f;
				if(root4.r > 0.f && root4.r < 1.f && abs(discrSpline.Eval(root4.r)) < err && abs(root4.i) < err) isDpos = 1.f;
			}
			//endregion
			
			//region binary roots
			// if(0)
			{
				isDpos = 0.f;

				float2 ival = float2(0.f, 1.f);
				// ival = float2(-r, 1.f);
				float c0 = discr_e;
				float c1 = discr_d;
				float c2 = discr_c;
				float c3 = discr_b;
				float c4 = discr_a;
				
				OFUNC(poly4, IEVAL(float, (float), 1),
				float Eval(float x)
				{
					return (c4 * Pow4(x)) + (c3 * Pow3(x)) + (c2 * Pow2(x)) + (c1 * x) + c0;
				})
				
				OFUNC(poly4D1, IEVAL(float, (float), 1),
				float Eval(float x)
				{
					return (4.f * c4 * Pow3(x)) + (3.f * c3 * Pow2(x)) + (2.f * c2 * x) + c1;
				})
				
				OFUNC(poly4D2, IEVAL(float, (float), 1),
				float Eval(float x)
				{
					return (12.f * c4 * Pow2(x)) + (6.f * c3 * x) + (2.f * c2);
				})
				
				OFUNC(poly4D2Roots, INull,
				float2 Eval()
				{
					float a = 12.f * c4;
					float b =  6.f * c3;
					float c =  2.f * c2;
				
					float discr = b * b - 4.f * a * c;
					discr = max(0.f, discr);
					
					return (-b + float2(1.f, -1.f) * sqrt(discr)) / (2.f * a);
					return -2.f * c / (b + float2(1.f, -1.f) * sqrt(discr));
				})
				
				OFUNC(binRootFinder, INull,
				float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
				{
					// float n = func.Eval(x1) < 0.f ? x1 : x2;
					// float p = x1 + x2 - n;
					
					if(func.Eval(p) < 0.f) return p;
					if(func.Eval(n) > 0.f) return n;
					
					for(uint i = 0; i < iCount; ++i)
					{
						float m = (n + p) * 0.5f;
						
						float f = func.Eval(m);
						
						if(f < 0.f)
						n = m;
						else
						p = m;
					}
					
					return (n + p) * 0.5f;
				})
				
				float2 d2Roots = poly4D2Roots.Eval();
				
				d2Roots = saturate(d2Roots);
				
				uint iCount1 = 10;
				uint iCount2 = 20;
				
				float3 d1Roots = 0.f;
				d1Roots.x = binRootFinder.Eval(d2Roots.x, ival.x, poly4D1, iCount1);
				d1Roots.y = binRootFinder.Eval(d2Roots.x, d2Roots.y, poly4D1, iCount1);
				d1Roots.z = binRootFinder.Eval(ival.y, d2Roots.y, poly4D1, iCount1);
				
				float3 fd;
				fd.x = discrSpline.Eval(d1Roots.x);
				fd.y = discrSpline.Eval(d1Roots.y);
				fd.z = discrSpline.Eval(d1Roots.z);
				
				float2 fd0;
				fd0.x = discrSpline.Eval(ival.x);
				fd0.y = discrSpline.Eval(ival.y);
				
				float2 fdT;
				fdT.x = qSplineIDist.Eval(d1Roots.x);
				fdT.y = qSplineIDist.Eval(d1Roots.z);
				// d1Roots = saturate(d1Roots);
				
				float4 roots = -1.f;
				{
					float n1 = ival.x;
					float p1 = d1Roots.x;
					float p2 = d1Roots.x;
					float n2 = d1Roots.y;
					
					if(fd.x > 0 || fd.z > 0)
					{
						if(fd.x > 0 && fd.z > 0)
						{
							if(fd.y < 0)
							{
								if(fdT.x > fdT.y)
								{
									n1 = d1Roots.y;
									p1 = p2 = d1Roots.z;
									n2 = ival.y;
								}
							}
							else
							{
								n1 = ival.x;
								p1 = d1Roots.x;
								
								p2 = d1Roots.z;
								n2 = ival.y;
							}
						}
						else if(fd.x < 0)
						{
							n1 = d1Roots.y;
							p1 = p2 = d1Roots.z;
							n2 = ival.y;
						}
						
						// if(n1 == ival.x && fd0.x > 0.f) roots.x = n1; else
						roots.x = binRootFinder.Eval(n1, p1, poly4, iCount2);
						// roots.x = binRootFinder.Eval(lerp(n1, p1, 0.05f), lerp(p1, n1, 0.05f), poly4, iCount2);
						
						// if(n2 == ival.y && fd0.y > 0.f) roots.z = n2; else
						roots.z = binRootFinder.Eval(n2, p2, poly4, iCount2);
						// roots.z = binRootFinder.Eval(lerp(n2, p2, 0.05f), lerp(p2, n2, 0.05f), poly4, iCount2);
					}
					
					// if(discrSpline.Eval(d1Roots.x) < 0.f)// || d1Roots.x < 0.f
					// {
						// n = d1Roots.y;
						// p = d1Roots.z;
					// }
					
					// roots.x = binRootFinder.Eval(n, p, poly4, iCount2);
					// roots.x = binRootFinder.Eval(ival.x, d1Roots.x, poly4, iCount2);
					// roots.z = binRootFinder.Eval(d1Roots.y, d1Roots.z, poly4, iCount2);
					
					// roots.z = binRootFinder.Eval(ival.y, d1Roots.z, poly4, iCount2);
					// roots.z = binRootFinder.Eval(d1Roots.y, d1Roots.x, poly4, iCount2);
				}
				// roots.z = binRootFinder.Eval(ival.x, d1Roots.z, poly4);
				
				roots = saturate(roots);
				
				float err = 0.001f;
				
				// if(roots.x >= 0.f && roots.x <= 1.f && abs(discrSpline.Eval(roots.x)) < err) isDpos = 1.f;
				// if(roots.z >= 0.f && roots.z <= 1.f && abs(discrSpline.Eval(roots.z)) < err) isDpos = 1.f;
				
				// if(abs(discrSpline.Eval(roots.x)) < err) isDpos = 1.f;
				// if(abs(discrSpline.Eval(roots.z)) < err) isDpos = 1.f;
				if(discrSpline.Eval(roots.x) >= -err) isDpos = 1.f;
				if(discrSpline.Eval(roots.z) >= -err) isDpos = 1.f;
				if(qSplineIDist.Eval(roots.x) < 0.f) isDpos = 0.f;
				// if(roots.x > 0.9f || roots.x < 0.f) isDpos = 0.f;
				fCol = roots.z;
				
				//region finalize
				{
				
					OFUNC(minFinder, INull,
					float Eval(float a, float b, IEVAL(float, (float), 1) func, uint count, float count2 = 3.f)
					{
						float steps = 1.f / count2;
						float t = 10000.f;
						float l0 = 0.f;
						float l = 0.f;
						float l1 = 0.f;
						
						
						for(uint j = 0; j < count; ++j)
						{
							bool blob = false;
							t = 10000.f;
							l = l0;
							l1 = l0;
							
							for(float i = 0.f; i <= count2; ++i)
							// for(float i = count; i >= 0.f; --i)
							{
								float tempL = lerp(a, b, steps * i);
								// tempL = steps * i;
								float tempT = func.Eval(tempL);
								
								if(tempT < t)
								{
									t = tempT;
									
									l0 = l;
									l = tempL;
									l1 = tempL;
									
									blob = true;
								}
								else
								{
									if(blob) l1 = tempL;
									blob = false;
								}
							}
							
							// 
							a = l0;
							b = l1;
						}
						
						return l;
					})
					
					float l = minFinder.Eval(roots.x, roots.z, qSplineIDist, 10, 3);
					float t = qSplineIDist.Eval(l);
					// if(0)
					// for(float i = 0.f; i <= count; ++i)
					// {
						// float tempL = lerp(l0, l1, steps * i);

						// float tempT = qSplineIDist.Eval(tempL);
						
						// if(tempT < t)
						// {
							// t = tempT;
							// l = tempL;
						// }
					// }
					
					// l = (roots.x + roots.z) * 0.5f;
					// l = 0.f;
					// t = qSplineIDist.Eval(l);
					
					float3 sp = qSpline.Eval(l);
					sp = mul(sp, rayMat); sp += ray.Start;
					
					float3 ip = ray.Start + ray.Dir * t;
						
					float3 n = normalize(ip - sp);
					fCol = n * 0.5f + 0.5f;
				}
				//endregion
				
				// if(fd.x < 0.f && fd.x > -1.f || fd.z < 0.f && fd.z > -1.f) 
				// if(fd.x < 0.f || fd.z < 0.f) 
				// {
					// fCol = -max(fd.x, fd.z);
					// isDpos = 1.f;
				// }
				// if(fd.x < 0.f || fd.z < 0.f) 
				// if(fd.x < 0.f) 
				{
					// fCol = -max(fd.x, fd.z);
					
					// if(fd.x < 0.f && fd.z > 0.f)
					// if(fd.x > -1.f)
					// {fCol = -fd.x; isDpos = 1.f;}
					
					// if(fd.x > fd.z && fd.x < 0.f || fd.z > 0.f) {fCol = -fd.x; isDpos = 1.f;}
					// if(fd.z > fd.x && fd.z < 0.f || fd.x > 0.f) {fCol = -fd.z; isDpos = 1.f;}
					
					// fCol = -fd.x;
					// fCol = -fd.z;
					// isDpos = 1.f;
					// return sqrt(fCol + r2);
					// return -Draw::ApplyDistScale(sqrt(fCol + r2) - r);
					// return Draw::Isolines(sqrt(fCol + r2) - r);
				}
				// fCol = 1.f;
				
			}
			//endregion
			
			// float discr = d2.y * r2 + d2.z * r2 - d2.z * s2.y + 2.f * d.y * d.z * s.y * s.z - d2.y * s2.z;
			// discr = 0.f;
			// float2 l2 = (-d.y * s.y - d.z * s.z + sqrt(discr) * float2(-1.f, 1.f)) / (d2.y + d2.z);
			
			// isDpos = discrSpline.Eval(0) >= 0.f;
			
			if(!isDpos)
			// if(hasRealRootsPoly4.Eval(discr_a, discr_b, discr_c, discr_d, discr_e))
			// if(discrPoly2.Eval(discr_c, discr_d, discr_e) > 0.f)
			// if(abs(discr_b) > 0.001f)
			{
				fCol = 0.f;
				t0 = 0.f;
			}
			
			// fCol = 0.f;
			t0 = 0.f;
			// fCol = lerp(fCol, 0.5f, discrSpline.Eval(1+r) >= 0.f);
			
			// fCol = lerp(fCol, saturate(abs(discr_d)), 0.5f);
			// fCol = lerp(fCol, saturate(abs(discr_c)), 0.5f);
			// fCol = lerp(fCol, saturate(abs(discr_b)), 0.5f);
			// if(-discr_a < 1.f) {fCol = lerp(fCol, -discr_a , 0.5f); t0 = 0.f;}
			// if(discr_b < 0) {fCol = 0.5f; t0 = 0.f;}
			// if(discr_d < 0) {fCol = 0.5f; t0 = 0.f;}
			
			float t = 1000000.f;
			float l = 0.f;
			
			// if(0)
			// {
				// float count = 10.f;
				// float rCount = rcp(count);
				
				// float discr2 = 0;
				
				// for(float j = 0; j <= count; ++j)
				// {
					// float t_l = lerp(l2.y, l2.x, j * rCount);
					
					// discr2 = r2 - Pow2(d.y * t_l + s.y) - Pow2(d.z * t_l + s.z);
					// discr2 = max(0.f, discr2);
					// float t_t = d.x * t_l + s.x - sqrt(discr2);
					
					// if(t_t < t)
					// {
						// l = t_l;
						// t = t_t;
					// }
				// }
			// }
			
			// {
				// float count = 9.f;
				// float rCount = rcp(count);
				
				// float discr2 = 0;
				
				// // l2 = float2(min(l2.x, l2.y), max(l2.x, l2.y));
				
				// for(float j = 0; j < count; ++j)
				// {
					// float t_l = (l2.x + l2.y) * 0.5f;
					
					// discr2 = r2 - Pow2(d.y * t_l + s.y) - Pow2(d.z * t_l + s.z);
					// discr2 = max(0.f, discr2);
					// float t_t = d.x - (-2.f * d.y * (d.y * t_l + s.y) - 2.f * d.z * (d.z * t_l + s.z)) * 0.5f * rsqrt(discr2);
					
					// if(t_t > 0)
					// {
						// l2.y = t_l;
					// }
					// else
					// {
						// l2.x = t_l;
					// }
				// }
				
				// l = (l2.x + l2.y) * 0.5f;
				
				// discr2 = r2 - Pow2(d.y * l + s.y) - Pow2(d.z * l + s.z);
				// discr2 = max(0.f, discr2);
				// t = d.x * l + s.x - sqrt(discr2);
			// }

				// OFUNC(CylThing, INull,
				// float2 Eval(float3 s, float3 d, float r)
				// {

				// })

				
			float3 cp = 0.f; //s + d * l;
			float3 ip = ray.Start + ray.Dir * t;
				
			cp = mul(cp, rayMat); cp += ray.Start;
			// ip = mul(ip, rayMat); ip += ray.Start;
			
			float3 n = normalize(ip - cp);

			// if(discr >= 0.f && t > 0 && l > 0 && l < 1) 
			// {
				// fCol = n;// * 0.5f + 0.5f;
			
				// t0 = t;
			// }
			// return fCol;
		}
		
		res = Intrsect::RaySphere(ray.Start, ray.Dir, sphPos, radius);

		float3 hit = ray.Start + ray.Dir * res.t1;

		// if(!res.mask) return 0.f;

		float3 n = normalize(hit);
		float2 np = n.xz / abs(n.y);

		float checker = EvalChecker(hit);

		// if(dot(fCol, 1.f) == 0.f)
		if(res.t1 < t0)
		fCol = saturate(hit - sphPos);
	}
	
	//region WeiDuKe root finding
	if(0)
	{
		fCol = 0.f;
		
		float2 uv = fuv; uv.y = 1.f - uv.y;
		// uv.x += 
		uv.y = MapTo(uv.y, -1.f, 1.f);
		uv.x = MapTo(uv.x, -1.f, 1.f);
		uv *= 2.f;
		
		// fCol = Pow2(uv.x) - uv.y;
		// fCol = 1.f - saturate(abs(fCol.x) * RcpLen(float2(ddx(fCol.x), ddy(fCol.y))));
		fCol = lerp(saturate(fCol), 1.f, Draw::CoordSys(uv));
		
		float c0 = 0.25f;
		float c1 = 0.9f;
		float c2 = -0.2f;
		float c3 = -2.f;
		float c4 = -3.f;
		
		fCol = lerp(fCol, Col3::G, Draw::Poly(uv, c0, c1, c2, c3, c4));
		
		c0 /= c4;
		c1 /= c4;
		c2 /= c4;
		c3 /= c4;
		c4 = 1.f;
		
		// Cmplx root1 = cmplx(1.f, 0.f);
		// Cmplx root2 = cmplx(0.4f, 0.9f);
		// Cmplx root3 = root2.Mul(root2);
		// Cmplx root4 = root3.Mul(root2);
		
		float mi = 0.02f;
		Cmplx root1 = cmplx(-2.f, mi);
		Cmplx root2 = cmplx(-1.f, mi);
		Cmplx root3 = cmplx(0.f, mi);
		Cmplx root4 = cmplx(1.f, mi);
		
		Cmplx cc0 = cmplx(c0);
		Cmplx cc1 = cmplx(c1);
		Cmplx cc2 = cmplx(c2);
		Cmplx cc3 = cmplx(c3);
		Cmplx cc4 = cmplx(c4);
		
		OFUNC(poly4, IEVAL(Cmplx, (Cmplx), 1),
		Cmplx Eval(Cmplx x)
		{
			return EvalPoly(x, cc0, cc1, cc2, cc3, cc4);
		})

		EvalWeiDuKe(poly4, 4, root1, root2, root3, root4);
		
		float dotR = 0.04f;
		fCol = lerp(fCol, Col3::R, Draw::Dot(uv, root1.RI(), dotR));
		fCol = lerp(fCol, Col3::R, Draw::Dot(uv, root2.RI(), dotR));
		fCol = lerp(fCol, Col3::R, Draw::Dot(uv, root3.RI(), dotR));
		fCol = lerp(fCol, Col3::R, Draw::Dot(uv, root4.RI(), dotR));
		
		float ply = EvalPoly(cmplx(uv.x), cc0, cc1, cc2, cc3, cc4).r - uv.y;
		ply = 1.f - saturate(abs(ply) / length(float2(ddx(ply), ddy(ply))));
		
		fCol = lerp(fCol, Col3::R, ply);
		
		// return ply;
		// Cmplx cn1 = cmplx(-1.f, -1.f);
		// Cmplx cn2 = cmplx(-1.2f, 0.25f);
		// Cmplx cn3 = cn1.Mul(cn2);
		// Cmplx cn4 = cn3.Mul(cn2.Rcp());
		
		// fCol = lerp(fCol, Col3::R, Draw::Line(uv, 0.f, cn1.RI()));
		// fCol = lerp(fCol, Col3::Y, Draw::Line(uv, 0.f, cn2.RI()));
		// fCol = lerp(fCol, Col3::C, Draw::Line(uv, 0.f, cn3.RI()));
		// fCol = lerp(fCol, Col3::M, Draw::Line(uv, 0.f, cn4.RI()));
		// fCol = lerp(fCol, 1.f, Draw::Line(uv, float2(qs, float2(2.f, 0.f)));
	}
	//endregion
	
	//region binary root finding
	 if(0)
	{
		// fCol = 0.f;
		
		float2 uv = fuv; uv.y = 1.f - uv.y;
		// uv.x += 
		uv.y = MapTo(uv.y, -1.f, 1.f);
		uv.x = MapTo(uv.x, -1.f, 1.f);
		uv *= 2.f;
		
		// fCol = Pow2(uv.x) - uv.y;
		// fCol = 1.f - saturate(abs(fCol.x) * RcpLen(float2(ddx(fCol.x), ddy(fCol.y))));
		fCol = lerp(saturate(fCol), 1.f, Draw::CoordSys(uv));
		
		float2 ival = float2(-1.2f, 0.75f);
		float c0 = 0.1f;
		float c1 = 0.5f;
		float c2 = 0.f;
		float c3 = 0.f;
		float c4 = -4.f;
		
		OFUNC(poly4, IEVAL(float, (float), 1),
		float Eval(float x)
		{
			return (c4 * Pow4(x)) + (c3 * Pow3(x)) + (c2 * Pow2(x)) + (c1 * x) + c0;
		})
		
		OFUNC(poly4D1, IEVAL(float, (float), 1),
		float Eval(float x)
		{
			return (4.f * c4 * Pow3(x)) + (3.f * c3 * Pow2(x)) + (2.f * c2 * x) + c1;
		})
		
		OFUNC(poly4D2, IEVAL(float, (float), 1),
		float Eval(float x)
		{
			return (12.f * c4 * Pow2(x)) + (6.f * c3 * x) + (2.f * c2);
		})
		
		OFUNC(poly4D2Roots, INull,
		float2 Eval()
		{
			float a = 12.f * c4;
			float b =  6.f * c3;
			float c =  2.f * c2;
		
			float discr = b * b - 4.f * a * c;
			discr = max(0.f, discr);
			
			return (-b + float2(1.f, -1.f) * sqrt(discr)) / (2.f * a);
			return -2.f * c / (b + float2(1.f, -1.f) * sqrt(discr));
		})
		
		OFUNC(binRootFinder, INull,
		float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
		{
			// float n = func.Eval(x1) < 0.f ? x1 : x2;
			// float p = x1 + x2 - n;
			
			for(uint i = 0; i < iCount; ++i)
			{
				float m = (n + p) * 0.5f;
				
				float f = func.Eval(m);
				
				if(f < 0.f)
				n = m;
				else
				p = m;
			}
			
			return (n + p) * 0.5f;
		})
		
		float2 d2Roots = poly4D2Roots.Eval();
		
		uint iCount1 = 2;
		uint iCount2 = 6;
		
		float3 d1Roots = 0.f;
		d1Roots.x = binRootFinder.Eval(d2Roots.x, ival.x, poly4D1, iCount1);
		d1Roots.y = binRootFinder.Eval(d2Roots.x, d2Roots.y, poly4D1, iCount1);
		d1Roots.z = binRootFinder.Eval(ival.y, d2Roots.y, poly4D1, iCount1);
		
		float4 roots = 0.f;
		roots.x = binRootFinder.Eval(ival.x, d1Roots.x, poly4, iCount2);
		roots.z = binRootFinder.Eval(d1Roots.y, d1Roots.z, poly4, iCount2);
		// roots.z = binRootFinder.Eval(ival.x, d1Roots.z, poly4);
		
		{
			float dotR = 0.02f;
			
			fCol = lerp(fCol, Col3::G, Draw::Poly(uv, c0, c1, c2, c3, c4));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(ival.x, poly4.Eval(ival.x)), dotR));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(ival.y, poly4.Eval(ival.y)), dotR));
			
			fCol = lerp(fCol, Col3::R, Draw::Curve(uv, poly4D1));
			fCol = lerp(fCol, Col3::B, Draw::Curve(uv, poly4D2));
			
			
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(d2Roots.x, 0.f), dotR));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(d2Roots.y, 0.f), dotR));
			
			// fCol = lerp(fCol, Col3::B, Draw::Dot(uv, float2(d2Roots.x, poly4D1.Eval(d2Roots.x)), dotR));
			// fCol = lerp(fCol, Col3::B, Draw::Dot(uv, float2(d2Roots.y, poly4D1.Eval(d2Roots.y)), dotR));
			
			fCol = lerp(fCol, Col3::C, Draw::Dot(uv, float2(d1Roots.x, poly4.Eval(d1Roots.x)), dotR));
			fCol = lerp(fCol, Col3::C, Draw::Dot(uv, float2(d1Roots.y, poly4.Eval(d1Roots.y)), dotR));
			fCol = lerp(fCol, Col3::C, Draw::Dot(uv, float2(d1Roots.z, poly4.Eval(d1Roots.z)), dotR));
			
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(roots.x, 0.f), dotR));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(roots.z, 0.f), dotR));
		}
	}
	//endregion
	
	//region min finder
	if(0)
	{
		fCol = 0.f;
		
		float2 uv = fuv; uv.y = 1.f - uv.y;
		// uv.x += 
		uv.y = MapTo(uv.y, -1.f, 1.f);
		uv.x = MapTo(uv.x, -1.f, 1.f);
		uv *= 2.f;
		
		OFUNC(minFinder, INull,
		float Eval(float a, float b, IEVAL(float, (float), 1) func, uint count, float count2 = 3.f)
		{
			uint noo = 1;
			
			float steps = 1.f / count2;
			float t = 10000.f;
			float l0 = 0.f;
			float l = 0.f;
			float l1 = 0.f;
			
			
			for(uint j = 0; j < count; ++j)
			{
				bool blob = false;
				t = 10000.f;
				l = l0;
				l1 = l0;
				
				for(float i = 0.f; i <= count2; ++i)
				// for(float i = count; i >= 0.f; --i)
				{
					float tempL = lerp(a,b,  steps * i);
					// tempL = steps * i;
					float tempT = func.Eval(tempL);
					
					if(j == noo)
					fCol = lerp(fCol, Col3::C, Draw::Dot(uv, float2(tempL, tempT), 0.025f));
					
					if(tempT < t)
					{
						t = tempT;
						
						l0 = l;
						l = tempL;
						l1 = tempL;
						
						blob = true;
					}
					else
					{
						if(blob) l1 = tempL;
						blob = false;
					}
				}
				
				// 
				a = l0;
				b = l1;
				
				if(j == noo)fCol = lerp(fCol, Col3::Y, Draw::Circle(uv, float2(a, func.Eval(a)), 0.05f));
				if(j == noo)fCol = lerp(fCol, Col3::Y, Draw::Circle(uv, float2(b, func.Eval(b)), 0.05f));
			}
			
			return l;
		})
		
		
		// fCol = Pow2(uv.x) - uv.y;
		// fCol = 1.f - saturate(abs(fCol.x) * RcpLen(float2(ddx(fCol.x), ddy(fCol.y))));
		fCol = lerp(saturate(fCol), 1.f, Draw::CoordSys(uv));
		
		float2 ival = float2(-0.7f, 0.75f);
		float c0 = 0.1f;
		float c1 = 0.5f;
		float c2 = 0.f;
		float c3 = 0.f;
		float c4 = -4.f;
		
		OFUNC(poly4, IEVAL(float, (float), 1),
		float Eval(float x)
		{
			return -((c4 * Pow4(x)) + (c3 * Pow3(x)) + (c2 * Pow2(x)) + (c1 * x) + c0);
		})
		
		float minl = minFinder.Eval(ival.x, ival.y, poly4, 5, 3);
		float mint = poly4.Eval(minl);
		
		{
			float dotR = 0.02f;
			
			fCol = lerp(fCol, Col3::G, Draw::Curve(uv, poly4));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(ival.x, poly4.Eval(ival.x)), dotR));
			fCol = lerp(fCol, Col3::R, Draw::Dot(uv, float2(ival.y, poly4.Eval(ival.y)), dotR));
			
			fCol = lerp(fCol, Col3::B, Draw::Dot(uv, float2(minl, poly4.Eval(minl)), dotR));
			
			// fCol = lerp(fCol, Col3::R, Draw::Curve(uv, poly4D1));
			// fCol = lerp(fCol, Col3::B, Draw::Curve(uv, poly4D2));
		}
	}
	//endregion
	
	// return fCol / rayCount;
	return GammaEncode(fCol / rayCount);
	// return pos.y;
	// return hit;
	// return checker;
	// return hit * 0.125f;
	
	// return 0.f;
	return res.mask;	
}
