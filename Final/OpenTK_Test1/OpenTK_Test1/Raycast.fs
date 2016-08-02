
#ifndef RAYCAST_DEFS
	#error [Raycast.fs depends on: Raycast.defs]
#endif

// source: https://gist.github.com/paniq/c86df99358c76906a6c4
// solve_quadratic0 & solve_quadratic by Leonard Ritter (paniq)
//-----------------------------------------------------------//
float solve_quadratic0(float3 fa) 
{
    float a = fa.z;
    float b = fa.y;
    float c = fa.x;
	
    // the quadratic solve doesn't work for a=0
    // so we need a branch here.
    if (a == 0.0) {
        return -c / b;
    } else { 
        // (-b +- sqrt(b*b - 4.0*a*c)) / 2.0*a
        float k = -0.5*b/a;
        float q = sqrt(k*k - c/a);
        // pick the closest positive root
	return k + ((k <= q)?q:-q);
    }
}

float solve_quadratic(float3 fa, float x) 
{
    float a = fa.z;
    float b = fa.y;
    float c = fa.x;

    // the quadratic solve doesn't work for a=0
    // so we need a branch here.
    if (a == 0.0) {
        return -c / b;
    } else { 
        // (-b +- sqrt(b*b - 4.0*a*c)) / 2.0*a
        float k = -0.5*b/a;
        float q = sqrt(k*k - c/a);
        float q0 = k - q;
        float q1 = k + q;
        // pick the closest root right of x
        return (q0 <= x)?q1:q0;
    }
}
//-----------------------------------------------------------//

float2 SolveQuadratic(float a, float b, float c) 
{
    if (a == 0.0) 
	{
		if(b == 0.0)
		return float2(-2.0, 2.0);
		else
        return float2(-c / b, 2.0);
    } 
	else 
	{ 
		float discr = b * b - 4.0 * a * c;
		
		if(discr == 0.0)
		return float2(-b / (2.0 * a), 2.0);
		
		if(discr < 0.0)
		return float2(-2.0, 2.0);
		
		float2 r = (-float2(b) + float2(-1.0, 1.0) * float2(sqrt(discr))) / (2.0 * a);
		
		return r.x < r.y ? r.xy : r.yx;
    }
}


struct Hit
{
	bool Mask;
	float Depth;
	float3 Col;
	float3 Normal;
	float L;
};

Hit Hit_New()
{
	Hit hit;

	hit.Mask = false;
	hit.Depth = 99999999999.0;
	hit.Col = float3(0.0);
	hit.Normal = float3(0.0);
	hit.L = 0.0;
	
	return hit;
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

#define DEF_POLYADD(N)                                                                                                     \
void Add(in float a[N], in float b[N], out float o_c[N]) { for(uint i = 0; i < N; ++i){ o_c[i] = a[i] + b[i]; } }  \
void Sub(in float a[N], in float b[N], out float o_c[N]) { for(uint i = 0; i < N; ++i){ o_c[i] = a[i] - b[i]; } }

DEF_POLYADD(5)
DEF_POLYADD(7)

float EvalPoly(float x, float c0, float c1 = 0, float c2 = 0, float c3 = 0, float c4 = 0, float c5 = 0, float c6 = 0)
{
	return x * (x * (x * (x * (x * (x * c6 + c5) + c4) + c3) + c2) + c1) + c0;
}

float EvalPolyD0(float x, float c[3]) { return EvalPoly(x, c[0], c[1], c[2]); }
float EvalPolyD1(float x, float c[3]) { return EvalPoly(x, c[1], c[2] * 2.0); }
float EvalPolyD2(float x, float c[3]) { return EvalPoly(x, c[2] * 2.0);       }

	float EvalPolyD0(float x, float c[4]) { return EvalPoly(x, c[0], c[1], c[2], c[3]);       }
	float EvalPolyD1(float x, float c[4]) { return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0); }
	float EvalPolyD2(float x, float c[4]) { return EvalPoly(x, c[2] * 2.0, c[3] * 6.0);       }
	float EvalPolyD3(float x, float c[4]) { return EvalPoly(x, c[3] * 6.0);                   }
	
float EvalPolyD0(float x, float c[5]) { return EvalPoly(x, c[0], c[1], c[2], c[3], c[4]);             }
float EvalPolyD1(float x, float c[5]) { return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0, c[4] * 4.0); }
float EvalPolyD2(float x, float c[5]) { return EvalPoly(x, c[2] * 2.0, c[3] * 6.0, c[4] * 12.0);      }
float EvalPolyD3(float x, float c[5]) { return EvalPoly(x, c[3] * 6.0, c[4] * 24.0);                  }

	float EvalPolyD0(float x, float c[7]) { return EvalPoly(x, c[0], c[1], c[2], c[3], c[4], c[5], c[6]);                         }
	float EvalPolyD1(float x, float c[7]) { return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0, c[4] * 4.0, c[5] * 5.0, c[6] * 6.0); }
	float EvalPolyD2(float x, float c[7]) { return EvalPoly(x, c[2] * 2.0, c[3] * 6.0, c[4] * 12.0, c[5] * 20.0, c[6] * 30.0);    }
	float EvalPolyD3(float x, float c[7]) { return EvalPoly(x, c[3] * 6.0, c[4] * 24.0, c[5] * 60.0, c[6] * 120.0);               }

	// for experimental purposes only:
	// #define CUBIC

	#ifndef CUBIC
		#define N0 3
		#define N1 5
	#else
		#define N0 4
		#define N1 7
	#endif
	
#define DEF_FINDROOTS(N, D)                                                                                                       \
void FindRoots##D(float poly_C[N1], float x_i[N], float m_i[N], out float x_o[N+1], out float m_o[N+1], uint iCount) \
{	                                                                                                                           \
    m_o[0] = m_o[N] = 1.0;                                                                                                     \
	                                                                                                                           \
	x_o[0] = x_i[0];                                                                                                           \
	                                                                                                                           \
	uint j = 0;                                                                                                                \
	                                                                                                                           \
	float x_l = x_i[0];                                                                                                        \
	float y_l = EvalPoly##D(x_l, poly_C);                                                                                                \
	float sy_l = sign(y_l);                                                                                                    \
                                                                                                                               \
	for(uint i = 1; i < N; ++i)                                                                                                \
	{                                                                                                                          \
		float x_r = x_i[i];                                                                                                    \
		float y_r = EvalPoly##D(x_r, poly_C);                                                                                            \
		float sy_r = sign(y_r);                                                                                                \
		                                                                                                                       \
		x_o[i] = 0.0;																										   \
		                                                                                                                       \
		if(m_i[i] == 1.0)                                                                                                      \
		{                                                                                                                      \
			if(sy_l != sy_r)                                                                                                   \
			{			                                                                                                       \
				float n = x_l;                                                                                                 \
				float p = x_r;                                                                                                 \
				                                                                                                               \
				float ny = EvalPoly##D(n, poly_C);                                                                                       \
				float py = EvalPoly##D(p, poly_C);                                                                                       \
				                                                                                                               \
				if(ny > 0.0 && py < 0.0)                                                                                       \
				{                                                                                                              \
					float t = n;                                                                                               \
					n = p; p = t;                                                                                              \
				}                                                                                                              \
				                                                                                                               \
				for(uint j = 0; j < iCount; ++j)                                                                               \
				{                                                                                                              \
					float m = (n + p) * 0.5;                                                                                   \
					                                                                                                           \
					float f = EvalPoly##D(m, poly_C);                                                                                    \
					                                                                                                           \
					if(f < 0.0)                                                                                                \
					n = m;                                                                                                     \
					else                                                                                                       \
					p = m;                                                                                                     \
				}                                                                                                              \
				                                                                                                               \
				x_o[i] = (n + p) * 0.5;                                                                                        \
                                                                                                                               \
				m_o[i] = 1.0;                                                                                                  \
			}			                                                                                                       \
			else                                                                                                               \
			m_o[i] = 0.0;                                                                                                      \
			                                                                                                                   \
			x_l = x_r;                                                                                                         \
			y_l = y_r;                                                                                                         \
			sy_l = sy_r;                                                                                                       \
		}                                                                                                                      \
		else                                                                                                                   \
		m_o[i] = -1.0;                                                                                                         \
	}                                                                                                                          \
	                                                                                                                           \
	x_o[N] = x_i[N - 1];                                                                                                       \
}

DEF_FINDROOTS(2, D3)
DEF_FINDROOTS(3, D2)
DEF_FINDROOTS(4, D1)
DEF_FINDROOTS(5, D0)

/* Always works if the input is non-zero.
 * Doesn’t require the input to be normalised.
 * Doesn’t normalise the output. 
 http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts */
float3 GetOrthoVec(float3 v)
{
    return abs(v.x) > abs(v.z) ? float3(-v.y, v.x, 0.0f) : float3(0.0f, -v.z, v.y);
}


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

	float3 qSplineEval(float l, float curveX[N0], float curveY[N0], float curveZ[N0])
	{
		return float3(EvalPolyD0(l, curveX), 
					  EvalPolyD0(l, curveY), 
					  EvalPolyD0(l, curveZ));
	}
	
	#define qSpline_Eval(l) qSplineEval(l, curveX, curveY, curveZ)
	#define qSpline_ParasDec float curveX[N0], float curveY[N0], float curveZ[N0]
	#define qSpline_Paras curveX, curveY, curveZ
	
	float qSplineIDistEval(float t, float curveX[N0], float polyB_C[N1])
	{		
		float term  = EvalPolyD0(t, curveX);
		float discr = EvalPolyD0(t, polyB_C);
		
		if(discr < 0.0)
		return 999999.0;
		else
		return term - sqrt(discr);
	}
	
	#define qSplineIDist_Eval(t) qSplineIDistEval(t, curveX, polyB_C)
	#define qSplineIDist_ParasDec float curveX[N0], float polyB_C[N1]
	#define qSplineIDist_Paras curveX, polyB_C
	
	float qSplineD1Eval(float t, float curveX[N0], float polyB_C[N1])
	{	 		
		float f1D1 = EvalPolyD1(t, curveX);
		
		float f2D0 = EvalPolyD0(t, polyB_C);
		float f2D1 = EvalPolyD1(t, polyB_C);
		
		// float2 r = Poly::EvalD01(t, polyB_C);
		
		// f2D0 = r.x;
		// f2D1 = r.y;
		
		// return f1D1 * 2.0 * sqrt(max(0.0, f2D0)) - f2D1;
		return f1D1 - f2D1 * 0.5 * rsqrt(max(0.0, f2D0));
	}
	
	#define qSplineD1_Eval(t) qSplineD1Eval(t, curveX, polyB_C)
	#define qSplineD1_ParasDec float curveX[N0], float polyB_C[N1]
	#define qSplineD1_Paras curveX, polyB_C
	
	float qSplineD2Eval(float t, float curveX[N0], float polyB_C[N1])
	{		
		float f1D1 = EvalPolyD1(t, curveX);
		float f1D2 = EvalPolyD2(t, curveX);
		
		float f2D0 = EvalPolyD0(t, polyB_C);
		float f2D1 = EvalPolyD1(t, polyB_C);
		float f2D2 = EvalPolyD2(t, polyB_C);
		
		// float3 r = Poly::EvalD012(t, polyB_C);
		
		// f2D0 = r.x;
		// f2D1 = r.y;
		// f2D2 = r.z;

		f2D0 = max(0.0, f2D0);
		
		// return Pow2(f2D1) / (4.0 * pow(f2D0, 1.5f)) + f1D2 - f2D2 / (2.0 * sqrt(f2D0));
		// return Pow2(f2D1) / (4.0 * sqrt(f2D0) * f2D0) + f1D2 - f2D2 / (2.0 * sqrt(f2D0));

		return (Pow2(f2D1) / f2D0 * 0.25 - f2D2 * 0.5) * rsqrt(f2D0) + f1D2;	
	}
	
	#define qSplineD2_Eval(t) qSplineD2Eval(t, curveX, polyB_C)
	#define qSplineD2_ParasDec float curveX[N0], float polyB_C[N1]
	#define qSplineD2_Paras curveX, polyB_C
	
	float qSplineD3Eval(float t, float polyB_C[N1])
	{	
		float f2D0 = EvalPolyD0(t, polyB_C);
		float f2D1 = EvalPolyD1(t, polyB_C);
		float f2D2 = EvalPolyD2(t, polyB_C);
		float f2D3 = EvalPolyD3(t, polyB_C);
		
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
	}
	
	#define qSplineD3_Eval(t) qSplineD3Eval(t, polyB_C)
	#define qSplineD3_ParasDec float polyB_C[N1]
	#define qSplineD3_Paras polyB_C
	
	
	#define DEF_binRootFinder(func)                                                          \
	float func##_BinRootFinder_Eval(float n, float p, const uint iCount, func##_ParasDec)    \
	{		                                                                                 \
		if(func##Eval(n, func##_Paras) > 0.0) return n;                                                     \
		if(func##Eval(p, func##_Paras) < 0.0) return p;		                                             \
		                                                                                     \
		for(uint i = 0; i < iCount; ++i)                                                     \
		{                                                                                    \
			float m = (n + p) * 0.5;                                                         \
			                                                                                 \
			float f = func##Eval(m, func##_Paras);                                           \
			                                                                                 \
			if(f < 0.0)                                                                      \
			n = m;                                                                           \
			else                                                                             \
			p = m;                                                                           \
		}                                                                                    \
		                                                                                     \
		return (n + p) * 0.5;                                                                \
	}
	
	DEF_binRootFinder(qSplineIDist)
	DEF_binRootFinder(qSplineD1)
	DEF_binRootFinder(qSplineD2)
	DEF_binRootFinder(qSplineD3)
	
	#define binRootFinder_Eval(n, p, func, iCount) func##_BinRootFinder_Eval(n, p, iCount, func##_Paras)
	

Hit EvalSplineISect(float3 rayStart, float3 rayDir, float3 s, float3 h, float3 t, float rs, float rh, float rt)
{
	Hit hit = Hit_New();
	
	
	float3x3 rayMat;
	rayMat[0] = rayDir;
	rayMat[1] = normalize(GetOrthoVec(rayDir));
	rayMat[2] = cross(rayDir, rayMat[1]);
	
	
	float3 fCol = float3(0.0);
	
	s -= rayStart; s = mul(s, rayMat);
	t -= rayStart; t = mul(t, rayMat);
	h -= rayStart; h = mul(h, rayMat);

	float curveX[N0];
	float curveY[N0];
	float curveZ[N0];

	float rcurve[N0];
	
	
	#ifndef CUBIC
	SplinePointsToPolyCoeffs(s.x, h.x, t.x, curveX);
	SplinePointsToPolyCoeffs(s.y, h.y, t.y, curveY);
	SplinePointsToPolyCoeffs(s.z, h.z, t.z, curveZ);	
	
	SplinePointsToPolyCoeffs(rs, rh, rt, rcurve);
	#else
	
	float3 h1 = s  + (h - s)   * 0.666;
	float3 h2 = t  + (h - t)   * 0.666;
	float rh1 = rs + (rh - rs) * 0.666;
	float rh2 = rt + (rh - rt) * 0.666;
	
	// h1.z += 0.5;
	// h2.z -= 0.5;
	
	SplinePointsToPolyCoeffs(s.x, h1.x, h2.x, t.x, curveX);
	SplinePointsToPolyCoeffs(s.y, h1.y, h2.y, t.y, curveY);
	SplinePointsToPolyCoeffs(s.z, h1.z, h2.z, t.z, curveZ);	
	
	SplinePointsToPolyCoeffs(rs, rh1, rh2, rt, rcurve);
	
	// SplinePointsToPolyCoeffs(s.x, h.x, h.x, t.x, curveX);
	// SplinePointsToPolyCoeffs(s.y, h.y, h.y, t.y, curveY);
	// SplinePointsToPolyCoeffs(s.z, h.z, h.z, t.z, curveZ);	
	
	// SplinePointsToPolyCoeffs(rs, rh, rh, rt, rcurve);
	#endif
	
	
	float polyB_C[N1];
		
	Pow2(rcurve, polyB_C);

	{
		float c[N1];
		Pow2(curveY, c);
		
		Sub(polyB_C, c, polyB_C);
	}
	
	{
		float c[N1];
		Pow2(curveZ, c);
		
		Sub(polyB_C, c, polyB_C);
	}
	
	float isDpos = 0.0;
		
	// OFUNC(quadPolyRoots, INull,
	// float2 Eval(float c, float b, float a)
	// {		
		// float discr = b * b - 4.0 * a * c;
		
		// return (-b + float2(1.0, -1.0) * sqrt(max(0.0, discr))) / (2.0 * a);
		// // return -2.0 * c / (b + float2(1.0, -1.0) * sqrt(discr));
	// })
	
	
	float l1 = 0.0;
	float l2 = 0.0;
	
	float rootType = 0.0;
	float4 roots = float4(-1.0);
	
	uint iCount = 10;

	// float x0[2]; x0[0] = 0.0; x0[1] = 1.0;
	// float m0[2]; m0[0] = 1.0; m0[1] = 1.0;
	
	// float x1[3]; float m1[3];
	// FindRootsD3(polyB_C, x0, m0, x1, m1, iCount);
	
	// x1[0] = 0.0; m1[0] = 1.0; 
	// x1[1] = -(polyB_C[3] * 6.0) / (polyB_C[4] * 24.0);
	// m1[1] = (x1[1] < 0.0 || x1[1] > 1.0) ? 0.0 : 1.0;
	// x1[2] = 1.0; m1[2] = 1.0;
	
	float x2[4]; float m2[4];	
	// FindRootsD2(polyB_C, x1, m1, x2, m2, iCount);
	
	float2 r = SolveQuadratic(polyB_C[4] * 12.0, polyB_C[3] * 6.0, polyB_C[2] * 2.0);
	
	x2[0] = 0.0; m2[0] = 1.0; 
	x2[1] = r.x; m2[1] = (x2[1] <= 0.0 || x2[1] >= 1.0) ? 0.0 : 1.0;
	x2[2] = r.y; m2[2] = (x2[2] <= 0.0 || x2[2] >= 1.0) ? 0.0 : 1.0;
	x2[3] = 1.0; m2[3] = 1.0;	
	
	float x3[5]; float m3[5];	
	FindRootsD1(polyB_C, x2, m2, x3, m3, iCount);
	
	float x4[6]; float m4[6];	
	for(uint i = 0; i < 6; ++i)m4[i] = 4.0;
	// if((m3[1] && EvalPolyD0(x3[1], polyB_C) > 0.0) || 
	   // (m3[2] && EvalPolyD0(x3[2], polyB_C) > 0.0) ||
	   // (m3[3] && EvalPolyD0(x3[3], polyB_C) > 0.0))
	{		
		FindRootsD0(polyB_C, x3, m3, x4, m4, iCount);
		
		float rn = 0.0;
		
		if(EvalPolyD0(0.0, polyB_C) >= 0.0) 
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
		
		if(EvalPolyD0(1.0, polyB_C) > 0.0) 
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
				float rootD3 = binRootFinder_Eval(roots.x, roots.y, qSplineD3, count);
				
				float2 rootsD2;
				rootsD2.x = binRootFinder_Eval(rootD3, roots.x, qSplineD2, count);
				rootsD2.y = binRootFinder_Eval(rootD3, roots.y, qSplineD2, count);
				
				l1 = binRootFinder_Eval(roots.x, rootsD2.x, qSplineD1, count);
				l2 = binRootFinder_Eval(rootsD2.y, roots.y, qSplineD1, count);
				// l1 = rootD3;
			}
			else
			{
				l1 = binRootFinder_Eval(roots.x, roots.y, qSplineD1, count);
				
				if(rootType == 2.0)
				l2 = binRootFinder_Eval(roots.z, roots.w, qSplineD1, count);
				else
				l2 = l1;
			}
			
			float t1 = qSplineIDist_Eval(l1);
			float t2 = qSplineIDist_Eval(l2);				
			
			float r1 = EvalPolyD0(l1, rcurve);
			float r2 = EvalPolyD0(l2, rcurve);
			
			bool hit1 = t1 > 0.0 && r1 > 0.0;
			bool hit2 = t2 > 0.0 && r2 > 0.0;
			
			
			float t, l;
			// should do the trick
			// if(hit1 && (t1 < t2 || !hit2))
			// {
				// t = t1;
				// l = l1;
			// }
			// else
			// {
				// t = t2;
				// l = l2;
			// }
				
			if(hit1)
			{
				if(t1 < t2 || !hit2)
				{
					t = t1;
					l = l1;
				}
				else
				{
					t = t2;
					l = l2;
				}
			}
			else
			{
				t = t2;
				l = l2;
			}
			
			float3 sp = qSpline_Eval(l);
			sp = mul(rayMat, sp); sp += rayStart;//eval untransformed curve to save mat mul!
			
			float3 ip = rayStart + rayDir * t;
				
			float3 n = normalize(ip - sp);
			
			// float3 viewVec = normalize(rayStart - ip);
			
			fCol = (n * 0.5 + 0.5);
			// fCol = float3(rootType == 1.0, rootType == 3.0, rootType == 2.0);
			
			fCol = pow(fCol, float3(2.2));
			// fCol *= saturate(dot(viewVec, n));

			// if(firstHit) fCol *= 0.0;
			// if(r1 > 0.0 || r2 > 0.0) fCol *= 0.0;
			// fCol *= 1.0 - saturate(Draw::ApplyDistScale(sqrt(min(-fd.x + r2, -fd.z + r2)) - r) + 1.0);
			
			// hit.Mask = true;
			hit.Mask = hit1 || hit2;
			// hit.Mask = firstHit || secondHit;
			hit.Depth = t;
			hit.Col = fCol;
			hit.Normal = n;
			hit.L = l;
		}
		//endregion
	}

	return hit;
}


in VS_OUT
{
	float3 Pos;

#if SHADING_STYLE == COLOR_HQ || SHADING_STYLE == LIT_COLOR_HQ
	float LOff;
	float3 ColP0;
	float3 ColT0;
	float3 ColP1;
	float3 ColT1;
#endif
	
	QTube QTube;		
} In;

out vec4 FragColor;

#ifdef USE_EARLYZ
	// layout(early_fragment_tests) in;

	// #extension GL_ARB_conservative_depth : require
	layout (depth_greater) out float gl_FragDepth;
	// layout (depth_less) out float gl_FragDepth;
	// layout (depth_any) out float gl_FragDepth;
#endif
 
void main()
{
#if SHADING_STYLE == BBOX 
    float3 normal = normalize(cross(dFdx(In.Pos), dFdy(In.Pos)));
	
	FragColor = vec4(GammaDecode(normal * 0.5 + 0.5), 1.0);
	return;
#endif
	
	float3 rayDir = normalize(In.Pos - CamPos);
	
	Hit hit = EvalSplineISect(CamPos, rayDir, 
	In.QTube.Start.Pos, In.QTube.H.Pos, In.QTube.End.Pos, 
	In.QTube.Start.Rad, In.QTube.H.Rad, In.QTube.End.Rad);

	// if(!hit.Mask) 
	FragColor = float4(0.0);
	gl_FragDepth = 1.0;
	
	// if(!hit.Mask) return;
	if(!hit.Mask) discard;

	// return rayDir;
		// float3 wPos2 = CamPos + rayDir * hit.Depth;
	// float3 normal2 = normalize(cross(dFdx(wPos2), dFdy(wPos2)));
	float rim = saturate(dot(hit.Normal, -rayDir));
	// float rim = saturate(dot(normal2, -rayDir));

	// #define SHADING_STYLE LIT_COLOR_HQ

#if SHADING_STYLE == COLOR
	FragColor.rgb = EvalColCurve(In.QTube, hit.L);
	FragColor.rgb *= rim;
	
#elif SHADING_STYLE == COLOR_HQ

	#ifdef SPLIT_IN_4_SEGMENTS
		FragColor.rgb = EvalCSpline(In.ColP0, In.ColT0, In.ColP1, In.ColT1, hit.L * 0.25 + In.LOff);
		// FragColor.rgb = float3(hit.L * 0.25 + In.LOff);
	#else
		FragColor.rgb = EvalCSpline(In.ColP0, In.ColT0, In.ColP1, In.ColT1, hit.L * 0.5 + In.LOff);
		// FragColor.rgb = float3(hit.L * 0.5 + In.LOff);
	#endif

	FragColor.rgb *= rim;	
FragColor.rgb += pow(rim, 24.0);

#elif SHADING_STYLE == LIT_COLOR_HQ

	#ifdef SPLIT_IN_4_SEGMENTS
		float3 col = EvalCSpline(In.ColP0, In.ColT0, In.ColP1, In.ColT1, hit.L * 0.25 + In.LOff);
	#else
		float3 col = EvalCSpline(In.ColP0, In.ColT0, In.ColP1, In.ColT1, hit.L * 0.5 + In.LOff);
	#endif
	
	float3 lightDir = normalize(float3(1.0, 1.0, 1.0));
	float3 normal = hit.Normal;

	float3 eyeVec = -rayDir;
	float NdL = dot(lightDir, normal);
	float diff = saturate(NdL);
	float back = -NdL * 0.5 + 0.5;
	float3 H = normalize(eyeVec + normal);
	float HdN = dot(normal, H);
	float spec = pow(saturate(HdN), 100.0) * diff * 10.0;
	// float fresnel = pow(1.0 - rim, 5.0/(1.0 - abs(dot(normalize(t), eyeVec)) * 0.8));
	float fresnel = pow(1.0 - rim, 15.0);
	// float LdT = dot(normalize(t), lightDir);
	// float cfac = 1.0 / (1.0 - abs(LdT));
	
	// FragColor.rgb = float3(exp2(-(float3(1.0) - clamp(col, float3(0.01), float3(0.99))) * (-NdL * 0.5 + 0.5) * (50.00 * 1.0)));
	// FragColor.rgb = lerp(FragColor.rgb, col, diff);
	// FragColor.rgb = float3(1.0 - abs(LdT));
	// FragColor.rgb = float3(fresnel);
	fresnel *= 0.5;
	fresnel *= 1.0 - diff;
	// fresnel *= (1.0 - abs(dot(normalize(t), eyeVec)) * 0.99);
	// fresnel *= ((dot(lightDir, -eyeVec) * 0.5 + 0.5));
	// diff *= 1.0 - sqrt(rim);// * 2.0;
	// diff *= (1.0 - (1.0 - rim) * (1.0 - rim));
	FragColor.rgb = col * (diff + col * back * 0.05+ fresnel);// ;
	FragColor.rgb += pow(rim, 24.0);
	// FragColor.rgb = float3(1.0 - abs(dot(normalize(t), eyeVec)));
	// FragColor.rgb = col * (diff + back * 0.05);
	
#elif SHADING_STYLE == SMOOTH_NORMAL
	FragColor.rgb = GammaDecode(hit.Normal * 0.5 + 0.5);
	FragColor.rgb *= rim;
	FragColor.rgb += pow(rim, 24.0);

#elif SHADING_STYLE == WIP
	FragColor.rgb = lerp(float3(1.0), float3(1.0, 0.3, 0.1), rim);
	
#endif
		
	// FragColor = float4(GammaEncode(FragColor.rgb), 1.0);
	
	float3 wPos = CamPos + rayDir * hit.Depth;
	float4 pPos = mul(float4(wPos, 1.0), ViewProjMat);
	
	gl_FragDepth = pPos.z / pPos.w;
	// if(!hit.Mask) gl_FragDepth = 1.0;
	// gl_FragDepth = 1.0;
	if(gl_FragDepth < 0.0) FragColor = float4(0.0, 0.0, 0.0, 1.0);
	 // gl_FragDepth = 1.0;
}