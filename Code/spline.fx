#ifndef SPLINE_FX
#define SPLINE_FX

// #ifndef COMPLEX_FX
	// #error [SPLINE_FX depends on: COMPLEX_FX]
// #endif

// #ifndef META_FX
	// #error [SPLINE_FX depends on: META_FX]
// #endif


// #ifndef DEF_IEvalCmplxCmplx
// #define DEF_IEvalCmplxCmplx
		// DEF_IEVAL(Cmplx, (Cmplx), 1)
// #endif

#define DEF_EVALPOLY(N)                                    \
float EvalPoly(float x, float c[N+1], uint d = 0)            \
{                                                          \
	float xn = 1.0;                                        \
	float r = 0.0;                                         \
	                                                       \
	[unroll]                                               \
	for(uint i = d; i < N+1; ++i)                            \
	{                                                      \
		float n = 1.0;                                     \
		                                                   \
		[unroll]                                           \
		for(uint j = 0; j < d; ++j)                        \
		n *= i - j;                                        \
		                                                   \
		r += n * c[i] * xn;                                \
		xn *= x;                                           \
	}                                                      \
                                                           \
	return r;                                              \
}

DEF_EVALPOLY(0)
DEF_EVALPOLY(1)
DEF_EVALPOLY(2)
DEF_EVALPOLY(3)
DEF_EVALPOLY(4)
DEF_EVALPOLY(5)
DEF_EVALPOLY(7)
DEF_EVALPOLY(8)
DEF_EVALPOLY(9)
DEF_EVALPOLY(10)

// #ifndef DEF_IEvalCmplxCmplx
// #define DEF_IEvalCmplxCmplx
		// DEF_IEVAL(Cmplx, (Cmplx), 1)
// #endif


#define DEF_POLY(N)                                         \
struct Poly##N : IEVAL(float, (float), 1)                   \
{                                                           \
	float C[N+1];                                             \
	                                                        \
	static Poly##N New(float c[N+1])                          \
	{                                                       \
		Poly##N poly;                                       \
		                                                    \
		poly.C = c;                                         \
		                                                    \
		return poly;                                        \
	}                                                       \
                                                            \
	float Eval(float x, uint d)                    	        \
	{                                                       \
		return EvalPoly(x, C, d);                           \
	}                                                       \
	                                                        \
    float Eval(float x)                     				\
	{                                                       \
		return EvalPoly(x, C, 0);                           \
	}                                                       \
	                                                        \
	IEVAL(float, (float), 1) New(uint d = 0)                \
	{                                                       \
		float l[N+1] = C;                                     \
		                                                    \
		OFUNC(poly, IEVAL(float, (float), 1),               \
		float Eval(float x)                                 \
		{                                                   \
			return EvalPoly(x, l, d);                       \
		})                                                  \
		                                                    \
		return poly;                                        \
	}                                                       \
};

DEF_POLY(0)
DEF_POLY(1)
DEF_POLY(2)
DEF_POLY(3)
DEF_POLY(4)
DEF_POLY(5)

struct Hit
{
	bool Mask;
	float Depth;
	float3 Col;
	float3 Normal;
	
	static Hit New()
	{
		Hit hit;

		hit.Mask = false;
		hit.Depth = 99999999999.0;
		hit.Col = 0.0;
		hit.Normal = 0.0;
		
		return hit;
	}
};

float2 ItervalAdd(float2 a, float2 b) {return a + b;}
float2 ItervalSub(float2 a, float2 b) {return a - b.yx;}
float2 ItervalMul(float2 a, float2 b)
{
	float v0 = a.x * b.x;
	float v1 = a.y * b.x;
	float v2 = a.x * b.y;
	float v3 = a.y * b.y;
	
	return float2(min(min(v0, v1), min(v2, v3)), max(max(v0, v1), max(v2, v3)));
}

float2 IntervalEvalPoly(float2 x, float c[5])
{
	float2 y = c[4].xx;
	
	for(int i = 3; i >= 0; --i)
	y = ItervalMul(y, x) + c[i].xx; 
	
	return y;
}

void Pow2(in float c[3], out float o_c[5])
{
	o_c[0] = c[0] * c[0]; 
	o_c[1] =  2.0 * c[0] * c[1];
	
	o_c[2] = c[1] * c[1] +  2.0 * c[0] * c[2];
	
	o_c[3] =  2.0 * c[2] * c[1];
	o_c[4] = c[2] * c[2];
}

void Pow2(in float c[4], out float o_c[7])
{
	o_c[0] = c[0] * c[0]; 
	o_c[1] =  2.0 * c[0] * c[1]; 
	o_c[2] =  2.0 * c[0] * c[2] + c[1] * c[1]; 
	
	o_c[3] =  2.0 * (c[1] * c[2] + c[0] * c[3]); 
	
	o_c[4] =  2.0 * c[3] * c[1] + c[2] * c[2]; 
	o_c[5] =  2.0 * c[3] * c[2]; 
	o_c[6] = c[3] * c[3];
}

#define DEF_POLYADD(N) void Add(in float a[N], in float b[N], out float o_c[N]) { [unroll]for(uint i = 0; i < N; ++i){ o_c[i] = a[i] + b[i]; } }
#define DEF_POLYADD(N) void Sub(in float a[N], in float b[N], out float o_c[N]) { [unroll]for(uint i = 0; i < N; ++i){ o_c[i] = a[i] - b[i]; } }

DEF_POLYADD(5)

OFUNC(FindRoots_BinRootFinder, INull,
float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
{
	float ny = func.Eval(n);
	float py = func.Eval(p);
	
	if(ny > 0.0 && py < 0.0)
	{
		float t = n;
		n = p; p = t;
	}
	
	[loop]
	for(uint i = 0; i < iCount; ++i)
	{
		float m = (n + p) * 0.5;
		
		float f = func.Eval(m);
		
		if(f < 0.0)
		n = m;
		else
		p = m;
	}
	
	return (n + p) * 0.5;
})

#define DEF_FINDROOTS(n)                                                                                                           \
void FindRoots(IEVAL(float, (float), 1) func, float x_i[n], float m_i[n], out float x_o[n+1], out float m_o[n+1], uint iCount)       \
{	                                                                                                                               \
    m_o[0] = m_o[n] = 1.0;                                                                                                       \
	                                                                                                                               \
	x_o[0] = x_i[0];                                                                                                               \
	                                                                                                                               \
	uint j = 0;                                                                                                                    \
	                                                                                                                               \
	float x_l = x_i[0];                                                                                                            \
	float y_l = func.Eval(x_l);                                                                                                    \
	float sy_l = sign(y_l);                                                                                                        \
                                                                                                                                   \
	for(uint i = 1; i < n; ++i)                                                                                                    \
	{                                                                                                                              \
		float x_r = x_i[i];                                                                                                          \
		float y_r = func.Eval(x_r);                                                                                                \
		float sy_r = sign(y_r);                                                                                                    \
		                                                                                                                           \
		x_o[i] = 0.0;																												\
		                                                                                                                           \
		if(m_i[i] == 1.0)                                                                                                                 \
		{                                                                                                                          \
			if(sy_l != sy_r)                                                                                                       \
			{			                                                                                                           \
				x_o[i] = FindRoots_BinRootFinder.Eval(x_l, x_r, func, iCount);                                                     \
				m_o[i] = 1.0;                                                                                                     \
			}			                                                                                                           \
			else                                                                                                                   \
			m_o[i] = 0.0;                                                                                                          \
			                                                                                                                       \
			x_l = x_r;                                                                                                             \
			y_l = y_r;                                                                                                             \
			sy_l = sy_r;                                                                                                           \
		}                                                                                                                          \
		else                                                                                                                       \
		m_o[i] = -1.0;                                                                                                             \
	}                                                                                                                              \
	                                                                                                                               \
	x_o[n] = x_i[n - 1];                                                                                                             \
}

DEF_FINDROOTS(2)
DEF_FINDROOTS(3)
DEF_FINDROOTS(4)
DEF_FINDROOTS(5)
DEF_FINDROOTS(6)

float EvalPolyD0(float x, float c[3]) { return EvalPoly(x, c[0], c[1], c[2]); }
float EvalPolyD1(float x, float c[3]) { return EvalPoly(x, c[1], c[2] * 2.0); }
float EvalPolyD2(float x, float c[3]) { return EvalPoly(x, c[2] * 2.0);       }

float EvalPolyD0(float x, float c[5]) { return EvalPoly(x, c[0], c[1], c[2], c[3], c[4]);             }
float EvalPolyD1(float x, float c[5]) { return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0, c[4] * 4.0); }
float EvalPolyD2(float x, float c[5]) { return EvalPoly(x, c[2] * 2.0, c[3] * 6.0, c[4] * 12.0);      }
float EvalPolyD3(float x, float c[5]) { return EvalPoly(x, c[3] * 6.0, c[4] * 24.0);                  }

namespace Poly
{
	float EvalD0(float x, float c[3]) { return EvalPoly(x, c[0], c[1], c[2]); }
	float EvalD1(float x, float c[3]) { return EvalPoly(x, c[1], c[2] * 2.0); }
	float EvalD2(float x, float c[3]) { return EvalPoly(x, c[2] * 2.0);       }

	float EvalD0(float x, float c[5]) { return EvalPoly(x, c[0], c[1], c[2], c[3], c[4]);             }
	float EvalD1(float x, float c[5]) { return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0, c[4] * 4.0); }
	float EvalD2(float x, float c[5]) { return EvalPoly(x, c[2] * 2.0, c[3] * 6.0, c[4] * 12.0);      }
	float EvalD3(float x, float c[5]) { return EvalPoly(x, c[3] * 6.0, c[4] * 24.0);                  }

	float2 EvalD01(float x, float c[5])
	{
		float2 r;
		
		float c00 =           c[4];
		float c01 = c00 * x + c[3];
		float c02 = c01 * x + c[2];
		float c03 = c02 * x + c[1];
			  r.x = c03 * x + c[0];
			  
		float c10 =           c00;
		float c11 = c10 * x + c01;
		float c12 = c11 * x + c02;
		      r.y = c12 * x + c03;
		
		return r;
	}
	
	float3 EvalD012(float x, float c[5])
	{
		float3 r;
		
		float c00 =           c[4];
		float c01 = c00 * x + c[3];
		float c02 = c01 * x + c[2];
		float c03 = c02 * x + c[1];
			  r.x = c03 * x + c[0];
			  
		float c10 =           c00;
		float c11 = c10 * x + c01;
		float c12 = c11 * x + c02;
		      r.y = c12 * x + c03;
		
		float c20 =           c10;
		float c21 = c20 * x + c11;
		      r.z = c21 * x + c12;
			  
		return r;
	}
	
	float4 EvalD0123(float x, float c[5])
	{
		float4 r;
		
		float c00 =           c[4];
		float c01 = c00 * x + c[3];
		float c02 = c01 * x + c[2];
		float c03 = c02 * x + c[1];
			  r.x = c03 * x + c[0];
			  
		float c10 =           c00;
		float c11 = c10 * x + c01;
		float c12 = c11 * x + c02;
		      r.y = c12 * x + c03;
		
		float c20 =           c10;
		float c21 = c20 * x + c11;
		      r.z = c21 * x + c12;
		
		float c30 =           c20;
		      r.w = c30 * x + c21;
			  
		return r;
	}
	
	IEVAL(float, (float), 1) NewD0(float c[3]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD0(x, c);}) return poly; }
	IEVAL(float, (float), 1) NewD1(float c[3]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD1(x, c);}) return poly; }
	IEVAL(float, (float), 1) NewD2(float c[3]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD2(x, c);}) return poly; }
	
	IEVAL(float, (float), 1) NewD0(float c[5]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD0(x, c);}) return poly; }
	IEVAL(float, (float), 1) NewD1(float c[5]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD1(x, c);}) return poly; }
	IEVAL(float, (float), 1) NewD2(float c[5]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD2(x, c);}) return poly; }
	IEVAL(float, (float), 1) NewD3(float c[5]) { OFUNC(poly, IEVAL(float, (float), 1), float Eval(float x) { return EvalD3(x, c);}) return poly; }
};

void SplinePointsToPolyCoeffs(float p0, float h, float p1, out float o_c[3])
{
	o_c[0] = p0;
	o_c[1] = -2.0 * p0 + 2.0 * h;
	o_c[2] =   p0 + p1 - 2.0 * h;
}

void SplinePointsToPolyCoeffs(float p0, float h0, float h1, float p1, out float o_c[4])
{
	o_c[0] = p0; 
	o_c[1] =  3.0 * (h0 - p0);
	o_c[2] = -6.0 * h0 + 3.0 * h1 + 3.0 * p0;
	o_c[3] =  3.0 * h0 - 3.0 * h1 - p0 + p1;
}

float MyPow2(float x) { return sign(x) * (x * x); }
float MySqrt(float x) { return sign(x) * sqrt(abs(x)); }

Hit EvalSplineISect(float3 rayStart, float3 rayDir, float3 s, float3 t, float3 h, float r = 0.25, uint mode = 0, float2 uv = 0.0)
{
	Hit hit = Hit::New();
	
	float rs = 0.25;
	float rh = 2.0;
	float rt = 0.125;
	
	// float rs = 0.25;
	// float rh = 0.8;
	// float rt = 0.33;
	
	s = 0.0;
	t = float3(2.0, 0.0, 0.0);
	h = float3(1.0, 3.0, 0.0);
	
	// { float3 o = s; s = t; t = o; }
	// { float ro = rs; rs = rt; rt = ro; }
	
	float3x3 rayMat;
	rayMat[0] = rayDir;
	rayMat[1] = normalize(GetOrthoVec(rayDir));
	rayMat[2] = cross(rayDir, rayMat[1]);
	

	// float3 c1 = float3(-1.0, -2.0, 0.0); 
	// float3 c2 = float3(3.0, -3.0, 0.0); 
	// float3 h = float3(-8.0, sin(Time * 0.0f) * 8.0 + 8.0, -2.0); 
	
	 // c1 = float3(-2.0, 0.0, 0.0); 
	 // c2 = float3( 2.0, 0.0, 0.0); 
	 // h  = float3( 0.0, 8.0, 0.0); 
	
	float3 fCol = 0.0;
	// float3 tc = c2; c2 = c1; c1 = tc;
	
	// h.xyz *= 4.0;
	// float3 c1 = float3(-1.0, 0.0, 0.0); 
	// float3 c2 = float3(3.0, 0.0, 0.0); 
	// float3 h = float3(2.5f, 3.0, 2.0); 
	// float r = 0.25f;
	
	// rayMat = transpose(rayMat);
	
	s -= rayStart; s = mul(rayMat, s);
	t -= rayStart; t = mul(rayMat, t);
	h -= rayStart; h = mul(rayMat, h);
	
	// float3 c1_2 = c1 * c1;
	// float3 c2_2 = c2 * c2;
	// float3 h2 = h * h;
	// float r2 = r * r
	
	
	float curveX[3];
	float curveY[3];
	float curveZ[3];

	SplinePointsToPolyCoeffs(s.x, h.x, t.x, curveX);
	SplinePointsToPolyCoeffs(s.y, h.y, t.y, curveY);
	SplinePointsToPolyCoeffs(s.z, h.z, t.z, curveZ);
	
	
	float polyB_C[5];
	
	float rcurve[3];
	SplinePointsToPolyCoeffs(rs, rh, rt, rcurve);
	
	Pow2(rcurve, polyB_C);

	{
		float c[5];
		Pow2(curveY, c);
		
		Sub(polyB_C, c, polyB_C);
	}
	
	{
		float c[5];
		Pow2(curveZ, c);
		
		Sub(polyB_C, c, polyB_C);
	}
	
	
	OFUNC(qSpline, INull,
	float3 Eval(float l)
	{
		return float3(EvalPolyD0(l, curveX), 
					  EvalPolyD0(l, curveY), 
					  EvalPolyD0(l, curveZ));
	})
	
	
	OFUNC(qSplineIDist, IEVAL(float, (float), 1),
	float Eval(float t)
	{		
		float term  = Poly::EvalD0(t, curveX);
		float discr = Poly::EvalD0(t, polyB_C);
		
		if(discr < 0.0)
		return 999999.0;
		else
		return term - sqrt(discr);
	})
	
	OFUNC(qSplineD1, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 		
		float f1D1 = Poly::EvalD1(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		
		// float2 r = Poly::EvalD01(t, polyB_C);
		
		// f2D0 = r.x;
		// f2D1 = r.y;
		
		// return f1D1 * 2.0 * sqrt(max(0.0, f2D0)) - f2D1;
		return f1D1 - f2D1 * 0.5 * rsqrt(max(0.0, f2D0));
	})
	
	OFUNC(qSplineD1_2, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 		
		float f1D1 = Poly::EvalD1(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		
		return 4.0 * Pow2(f1D1) * f2D0 - Pow2(f2D1);
	})
	
	OFUNC(qSplineD1_31, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 		
		float f1D1 = Poly::EvalD1(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		
		// return SqrLen((cmplx(f1D1).Mul(SqrtR(f2D0))).Mul(2.0).Sub(f2D1).RI());
		// return 2.0 * (cmplx(f1D1).Mul(SqrtR(f2D0))).r - f2D1;
		return 2.0 * f1D1 * MySqrt(f2D0) - f2D1;
		return 2.0 * f1D1 * sqrt(max(0.0, f2D0)) - f2D1;
	})
	
	OFUNC(qSplineD1_32, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 		
		float f1D1 = Poly::EvalD1(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		
		// return 2.0 * MyPow2(f1D1) * f2D0 - MyPow2(f2D1);
		return Pow3(2.0 * f1D1) * (f2D0 * sqrt(abs(f2D0))) - Pow3(f2D1);
		return MyPow2(2.0 * f1D1 * MySqrt(f2D0)) - MyPow2(f2D1);
		return 2.0 * f1D1 * sqrt(max(0.0, f2D0)) - f2D1;
	})
	
	OFUNC(qSplineD2, IEVAL(float, (float), 1),
	float Eval(float t)
	{		
		float f1D1 = Poly::EvalD1(t, curveX);
		float f1D2 = Poly::EvalD2(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		float f2D2 = Poly::EvalD2(t, polyB_C);
		
		// float3 r = Poly::EvalD012(t, polyB_C);
		
		// f2D0 = r.x;
		// f2D1 = r.y;
		// f2D2 = r.z;

		f2D0 = max(0.0, f2D0);
		
		// return Pow2(f2D1) / (4.0 * pow(f2D0, 1.5f)) + f1D2 - f2D2 / (2.0 * sqrt(f2D0));
		// return Pow2(f2D1) / (4.0 * sqrt(f2D0) * f2D0) + f1D2 - f2D2 / (2.0 * sqrt(f2D0));

		return (Pow2(f2D1) / f2D0 * 0.25 - f2D2 * 0.5) * rsqrt(f2D0) + f1D2;	
	})
	
	OFUNC(qSplineD2_2, IEVAL(float, (float), 1),
	float Eval(float t)
	{		
		float f1D1 = Poly::EvalD1(t, curveX);
		float f1D2 = Poly::EvalD2(t, curveX);
		
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		float f2D2 = Poly::EvalD2(t, polyB_C);
		
		return f1D1 * f2D1 + 2.0 * f2D0 * f1D2 - f2D2 * MySqrt(f2D0);
		return Pow2(f2D1) * 0.25 * MySqrt(f2D0) - f2D2 * 0.5 * MySqrt(f2D0) * f2D0 + f1D2 * f2D0;	
		return (Pow2(f2D1) / f2D0 * 0.25 - f2D2 * 0.5) * MySqrt(f2D0) + f1D2;	
	})
	
	OFUNC(qSplineD3, IEVAL(float, (float), 1),
	float Eval(float t)
	{	
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		float f2D2 = Poly::EvalD2(t, polyB_C);
		float f2D3 = Poly::EvalD3(t, polyB_C);
		
		// float4 r = Poly::EvalD0123(t, polyB_C);
		
		// f2D0 = r.x;
		// f2D1 = r.y;
		// f2D2 = r.z;
		// f2D3 = r.w;

		f2D0 = max(0.0, f2D0);
		
		// return rcp(8.0 * Pow2(f2D0) * sqrt(f2D0)) * (-3.0 * Pow3(f2D1) + 6.0 * f2D0 * f2D1 * f2D2 - 4.0 * Pow2(f2D0) * f2D3);
		return (-3.0 * Pow3(f2D1) + 6.0 * f2D0 * f2D1 * f2D2 - 4.0 * Pow2(f2D0) * f2D3) / Pow2(f2D0) * rsqrt(f2D0);// * 0.125;
		return (-3.0 * Pow3(f2D1) / Pow2(f2D0) + 6.0 * (f2D0 / Pow2(f2D0)) * f2D1 * f2D2 - 4.0 * f2D3) * rsqrt(f2D0) * 0.125;
		return (-3.0 * Pow3(f2D1) * Pow2(rcp(f2D0)) + 6.0 * rcp(f2D0) * f2D1 * f2D2 - 4.0 * f2D3) * rsqrt(f2D0);// * 0.125;
	})
	
	OFUNC(qSplineD3_2, IEVAL(float, (float), 1),
	float Eval(float t)
	{	
		float f2D0 = Poly::EvalD0(t, polyB_C);
		float f2D1 = Poly::EvalD1(t, polyB_C);
		float f2D2 = Poly::EvalD2(t, polyB_C);
		float f2D3 = Poly::EvalD3(t, polyB_C);
		
		return (-3.0 * Pow3(f2D1) + 6.0 * f2D0 * f2D1 * f2D2 - 4.0 * Pow2(f2D0) * f2D3) * f2D0;
		return (-3.0 * Pow3(f2D1) + 6.0 * f2D0 * f2D1 * f2D2 - 4.0 * Pow2(f2D0) * f2D3) / Pow2(f2D0) * f2D0;
	})
	
	
	float isDpos = 0.0;
		
	OFUNC(quadPolyRoots, INull,
	float2 Eval(float c, float b, float a)
	{		
		float discr = b * b - 4.0 * a * c;
		
		return (-b + float2(1.0, -1.0) * sqrt(max(0.0, discr))) / (2.0 * a);
		// return -2.0 * c / (b + float2(1.0, -1.0) * sqrt(discr));
	})
	
	OFUNC(binRootFinder, INull,
	float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
	{
		// float ny = func.Eval(n);
		// float py = func.Eval(p);
		
		// if(sign(ny) == sign(py))
		// {
			// if(abs(ny) < abs(py)) 
			// return n; 
			// else 
			// return p;
		// }
		// else
		// if(ny > 0.0 && py < 0.0)
		// {
			// float t = n;
			// n = p; p = t;
		// }
		
		// float n = func.Eval(x1) < 0.0 ? x1 : x2;
		// float p = x1 + x2 - n;
		
		if(func.Eval(n) > 0.0) return n;
		if(func.Eval(p) < 0.0) return p;		
		
		[loop]
		for(uint i = 0; i < iCount; ++i)
		{
			float m = (n + p) * 0.5;
			
			float f = func.Eval(m);
			
			if(f < 0.0)
			n = m;
			else
			p = m;
		}
		
		return (n + p) * 0.5;
	})
	
	
	float l1 = 0.0;
	float l2 = 0.0;
	
	float rootType = 0.0;
	float4 roots = -1.0;
	
	uint iCount = 10;

	float x0[2] = {0.0, 1.0};
	float m0[2] = {1.0, 1.0};
	
	float x1[3]; float m1[3];
	FindRoots(Poly::NewD3(polyB_C), x0, m0, x1, m1, iCount);
	
	float x2[4]; float m2[4];	
	FindRoots(Poly::NewD2(polyB_C), x1, m1, x2, m2, iCount);
	
	float x3[5]; float m3[5];	
	FindRoots(Poly::NewD1(polyB_C), x2, m2, x3, m3, iCount);
	
	float x4[6]; float m4[6];	
	for(uint i = 0; i < 6; ++i)m4[i] = 4.0;
	// if((m3[1] && Poly::EvalD0(x3[1], polyB_C) > 0.0) || 
	   // (m3[2] && Poly::EvalD0(x3[2], polyB_C) > 0.0) ||
	   // (m3[3] && Poly::EvalD0(x3[3], polyB_C) > 0.0))
	{		
		FindRoots(Poly::NewD0(polyB_C), x3, m3, x4, m4, iCount);
		
		float rn = 0.0;
		
		if(Poly::EvalD0(0.0, polyB_C) >= 0.0) 
		{ 
								  roots.x =  0.0; rn = 1.0; rootType = 1.5;
		}
		
		if(m4[1] == 1.0) 
		{
			if(rn == 0.0) 		{ roots.x = x4[1]; rn = 1.0; rootType = 1.5;						}
			else				{ roots.y = x4[1]; rn = 2.0; rootType = 1.0; 						}	
		}
		// else if(m4[1] == 0.0 && rootType == 1.5) rootType = 0.0;
		else if(rootType == 1.5) rootType = 0.0;
		
		if(m4[2] == 1.0) 
		{
			if     (rn == 0.0) 	{ roots.x = x4[2]; rn = 1.0; rootType = 1.5;						 }
			else if(rn == 1.0)	{ roots.y = x4[2]; rn = 2.0; rootType = rootType == 0.0 ? 3.0 : 1.0; }	
			else 				{ roots.z = x4[2]; rn = 3.0; rootType = 2.0; 						 }	
		}
		// else if(m4[2] == 0.0 && rootType == 1.5) rootType = 0.0;
		else if(rootType == 1.5) rootType = 0.0;
		
		if(m4[3] == 1.0) 
		{
			if     (rn == 0.0) 	{ roots.x = x4[3]; rn = 1.0; rootType = 1.5;						 }
			else if(rn == 1.0)	{ roots.y = x4[3]; rn = 2.0; rootType = rootType == 0.0 ? 3.0 : 1.0; }	
			else if(rn == 2.0)	{ roots.z = x4[3]; rn = 3.0; rootType = 2.0; 						 }	
			else 				{ roots.w = x4[3]; rn = 4.0; rootType = 2.0; 						 }	
		}
		// else if(m4[3] == 0.0 && rootType == 1.5) rootType = 0.0;
		else if(rootType == 1.5) rootType = 0.0;
		
		if(m4[4] == 1.0) 
		{
			if     (rn == 0.0) 	{ roots.x = x4[4]; rn = 1.0; rootType = 1.5;						 }
			else if(rn == 1.0)	{ roots.y = x4[4]; rn = 2.0; rootType = rootType == 0.0 ? 3.0 : 1.0; }	
			else if(rn == 2.0)	{ roots.z = x4[4]; rn = 3.0; rootType = 2.0; 						 }	
			else 				{ roots.w = x4[4]; rn = 4.0; rootType = 2.0; 						 }	
		} 
		// else if(m4[4] == 0.0 && rootType == 1.5) rootType = 0.0;
		else if(rootType == 1.5) rootType = 0.0;
		
		if(Poly::EvalD0(1.0, polyB_C) > 0.0) 
		{ 
			if(rn == 1.0) 		{ roots.y =   1.0; rn = 1.0; rootType = rootType == 0.0 ? 3.0 : 1.0; }
			else				{ roots.w =   1.0; rn = 2.0; rootType = 2.0; 						 }	
		}
		
		if(rootType == 1.5) rootType = 1.0;
		
		if(rootType == 1.0) rootType = 3.0;
		
			
		//region finalize
		if(rootType > 0.0)
		{		
			const uint count = 10;
			
			if(rootType == 3.0)
			{
				// uint iCount = 10;

				// float x0[2] = {roots.x, roots.y};
				// float m0[2] = {1.0, 1.0};
				
				// float x1[3]; float m1[3];
				// FindRoots(qSplineD3, x0, m0, x1, m1, iCount);
				
				// float x2[4]; float m2[4];	
				// FindRoots(qSplineD2, x1, m1, x2, m2, iCount);
				
				// float x3[5]; float m3[5];	
				// FindRoots(qSplineD1, x2, m2, x3, m3, iCount);
	
	
				float rootD3 = binRootFinder.Eval(roots.x, roots.y, qSplineD3, count);
				
				float2 rootsD2;
				rootsD2.x = binRootFinder.Eval(rootD3, roots.x, qSplineD2, count);
				rootsD2.y = binRootFinder.Eval(rootD3, roots.y, qSplineD2, count);
				
				l1 = binRootFinder.Eval(roots.x, rootsD2.x, qSplineD1, count);
				l2 = binRootFinder.Eval(rootsD2.y, roots.y, qSplineD1, count);
				// l1 = rootD3;
			}
			else
			{
				l1 = binRootFinder.Eval(roots.x, roots.y, qSplineD1, count);
				
				if(rootType == 2.0)
				l2 = binRootFinder.Eval(roots.z, roots.w, qSplineD1, count);
				else
				l2 = l1;
			}
			
			float t1 = qSplineIDist.Eval(l1);
			float t2 = qSplineIDist.Eval(l2);				
			
			float t = t1 < t2 ? t1 : t2;
			float l = t1 < t2 ? l1 : l2;
			
			// if(t2 < t1)
			// {
				// t1 = t2; 
				// l1 = l2;
			// }
			
			
			float3 sp = qSpline.Eval(l);
			sp = mul(sp, rayMat); sp += rayStart;
			
			float3 ip = rayStart + rayDir * t;
				
			float3 n = normalize(ip - sp);
			
			float3 viewVec = normalize(rayStart - ip);
			
			fCol = (n * 0.5 + 0.5);
			// fCol = float3(rootType == 1.0, rootType == 3.0, rootType == 2.0);
			
			fCol = pow(fCol, 2.2);
			// fCol *= saturate(dot(viewVec, n));

			// fCol *= 1.0 - saturate(Draw::ApplyDistScale(sqrt(min(-fd.x + r2, -fd.z + r2)) - r) + 1.0);
			
			hit.Depth = t;
			hit.Col = fCol;	
		}
		//endregion
	}
	
	OFUNC(ivalTest, INull,
	float Eval(float a, float b)
	{
		float2 v = IntervalEvalPoly(float2(a, b), polyB_C);
		// float2 v = IntervalEvalPoly(float2(a, a), polyB_C);
		
		
		// return Ternary(EvalPoly(a, polyB_C) > 0.0, 1.0, 0.0);
		return Ternary(v.y > 0.0, 1.0, 0.0);
		// return v.x < 0.0 && v.y > 0.0 ? 1.0 : 0.0;
	})
	
	// float2 poo = IntervalEvalPoly(float2(0.0, 1.0), polyB_C);
	// hit.Col = 0.5;
	// hit.Col = poo.x < 0.0 && poo.y > 0.0 ? 1.0 : 0.0;
	
	{
		float mask = 0.0;
		
		float tCount = 16.0;
		float stepS = 1.0 / tCount;
		
		for(float i = 0.0; i < tCount-1.0; ++i)
		mask += ivalTest.Eval(stepS * i, saturate(stepS * i + stepS));
		
		hit.Col = lerp(hit.Col, 1.0, saturate(mask) * 0.5);
	}
	
	
	//region Draw
	if(mode)
	{
		Hit hit = Hit::New();
		
		
		float3 col = 0.0;//lerp(saturate(0.0), 1.0, Draw::CoordSys(uv));
		// col = lerp(col, Col3::R, Draw::Curve(uv, poly));
		// col = lerp(col, Col3::C, Draw::Poly(uv, 0.0, 0.0, 1.0));
		col = lerp(col, Col3::G, Draw::Curve(uv, qSplineIDist) * 0.5);
		
		col = lerp(col, Col3::G, Draw::CircleS(uv, float2(l1, qSplineIDist.Eval(l1)), 4.0));
		// if(rootType > 1.0) 
		col = lerp(col, Col3::G, Draw::CircleS(uv, float2(l2, qSplineIDist.Eval(l2)), 6.0));
		
		col = lerp(col, Col3::O, Draw::Curve(uv, qSplineD1) * 0.5);
		// col = lerp(col, Col3::R, Draw::Curve(uv, qSplineD1_2) * 0.5);
		col = lerp(col, Col3::B, Draw::Curve(uv, qSplineD1_31) * 0.9);
		col = lerp(col, Col3::Y, Draw::Curve(uv, qSplineD2_2) * 0.9);
		col = lerp(col, Col3::R, Draw::Curve(uv, qSplineD3_2, 0.1) * 0.9);
		// col = lerp(col, Col3::C, Draw::Curve(uv, qSplineD1_32) * 0.9);
		// col = lerp(col, Col3::C, Draw::Curve(uv, qSplineD2, 0.1) * 0.5);
		// col = lerp(col, Col3::Y, Draw::Curve(uv, qSplineD3, 0.1) * 0.5);
		
		// col = lerp(col, Col3::R, Draw::Curve(uv, Poly::NewD0(polyB_C)) * 0.5);
		// col = lerp(col, Col3::B, Draw::Curve(uv, Poly::NewD1(polyB_C)) * 0.25);
		// col = lerp(col, Col3::O, Draw::Curve(uv, Poly::NewD1B(polyB_C)) * 0.5);
		// col = lerp(col, Col3::B, Draw::Curve(uv, polyBD1) * 0.5);
		// col = lerp(col, Col3::C, Draw::Curve(uv, polyBD2) * 0.5);
		// col = lerp(col, Col3::B, Draw::Curve(uv, polyD2));
		// col = lerp(col, Col3::C, Draw::Curve(uv, polyD3));
		// col = lerp(col, Col3::Y, Draw::Curve(uv, polyD4));
		// col = lerp(col, Col3::M, Draw::Curve(uv, polyD5));
		
		col = lerp(col, Col3::Y, Draw::CircleS(uv, float2(roots.x, 0.0), 2.0));
		col = lerp(col, Col3::Y, Draw::CircleS(uv, float2(roots.y, 0.0), 4.0));
		
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(d1Roots.x, 0.0), 6.0));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(d1Roots.y, 0.0), 8.0));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(d1Roots.z, 0.0), 10.0));
		
		col = lerp(col, Col3::O, Draw::CircleS(uv, float2(roots.z, 0.0), 6.0));
		col = lerp(col, Col3::O, Draw::CircleS(uv, float2(roots.w, 0.0), 8.0));
		
		if(rootType > 0.0) col = lerp(col, Col3::R, Draw::CircleS(uv, float2(0.0, 0.5), 4.0));
		if(rootType > 1.0) col = lerp(col, Col3::R, Draw::CircleS(uv, float2(0.0, 0.5), 6.0));
		if(rootType > 2.0) col = lerp(col, Col3::R, Draw::CircleS(uv, float2(0.0, 0.5), 8.0));
		
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(0.0 , 4.0),  4.0 + m4[1] * 2.0));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(0.33, 4.0),  4.0 + m4[2] * 2.0));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(0.66 , 4.0), 4.0 + m4[3] * 2.0));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(1.0 , 4.0),  4.0 + m4[4] * 2.0));
		
		// col = lerp(col, Col3::G, Draw::CircleS(uv, float2(x1*4.0-2.0, poly.Eval(x1*4.0-2.0)), 4));
		
		// col = lerp(col, Col3::O, Draw::CircleS(uv, float2(x2*4.0-2.0, poly.Eval(x2*4.0-2.0)), 6));
		// col = lerp(col, Col3::Y, Draw::CircleS(uv, float2(x3*4.0-2.0, poly.Eval(x3*4.0-2.0)), 8));
		
		// col = lerp(col, Col3::B, Draw::CircleS(uv, float2(m0_0, 0.0), 2));
		// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(m0_1, 0.0), 4));
		
		// col = lerp(col, Col3::B, Draw::CircleS(uv, float2(x2*4.0-2.0, 0.0), 2));
		
		// hit.Mask = Draw::Curve(uv, poly);
		// hit.Mask = Draw::CoordSys(uv);
		// hit.Col = Col3::R;
		
		hit.Col = col;
		
		return hit;
	}
	//endregion
	
	return hit;
}

#endif