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
	float xn = 1.f;                                        \
	float r = 0.f;                                         \
	                                                       \
	[unroll]                                               \
	for(uint i = d; i < N+1; ++i)                            \
	{                                                      \
		float n = 1.f;                                     \
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
		hit.Depth = 99999999999.f;
		hit.Col = 0.f;
		hit.Normal = 0.f;
		
		return hit;
	}
};



Hit EvalSplineISect(float3 rayStart, float3 rayDir, float3 c1, float3 c2, float3 h, float r = 0.25f, uint mode = 0)
{
	Hit hit = Hit::New();
	
	float3x3 rayMat;
	rayMat[0] = rayDir;
	rayMat[1] = normalize(GetOrthoVec(rayDir));
	rayMat[2] = cross(rayDir, rayMat[1]);
	

	// float3 c1 = float3(-1.f, -2.f, 0.f); 
	// float3 c2 = float3(3.f, -3.f, 0.f); 
	// float3 h = float3(-8.f, sin(Time * 0.0f) * 8.f + 8.f, -2.f); 
	
	 // c1 = float3(-2.f, 0.f, 0.f); 
	 // c2 = float3( 2.f, 0.f, 0.f); 
	 // h  = float3( 0.f, 8.f, 0.f); 
	
	float3 fCol = 0.f;
	// float3 tc = c2; c2 = c1; c1 = tc;
	
	// h.xyz *= 4.f;
	// float3 c1 = float3(-1.f, 0.f, 0.f); 
	// float3 c2 = float3(3.f, 0.f, 0.f); 
	// float3 h = float3(2.5f, 3.f, 2.f); 
	// float r = 0.25f;
	
	// rayMat = transpose(rayMat);
	
	c1 -= rayStart; c1 = mul(rayMat, c1);
	c2 -= rayStart; c2 = mul(rayMat, c2);
	h -= rayStart; h = mul(rayMat, h);
	
	float3 c1_2 = c1 * c1;
	float3 c2_2 = c2 * c2;
	float3 h2 = h * h;
	float r2 = r * r;
	
	OFUNC(qSpline, INull,
	float3 Eval(float l)
	{
		return lerp(lerp(c1, h, l), lerp(h, c2, l), l);
	})
	

	
	float polyA_C[3] = 
	{
		c1.x, 
		-2.f * c1.x + 2.f * h.x, 
		c1.x + c2.x - 2.f * h.x
	};
	
	Poly2 polyA = Poly2::New(polyA_C);
	
	float polyB_C[5] = 
	{
		-c1_2.y - c1_2.z + r2,
		
		4.f * c1_2.y + 4.f * c1_2.z - 4.f * c1.y * h.y - 4.f * c1.z * h.z,
		
		-6.f * c1_2.y - 6.f * c1_2.z - 2.f * c1.y * c2.y - 2.f * c1.z * c2.z + 
		12.f * c1.y * h.y - 4.f * h2.y + 12.f * c1.z * h.z - 4.f * h2.z,
		
		4.f * c1_2.y + 4.f * c1_2.z + 4.f * c1.y * c2.y + 4.f * c1.z * c2.z - 
		12.f * c1.y * h.y - 4.f * c2.y * h.y + 8.f * h2.y - 12.f * c1.z * h.z - 4.f * c2.z * h.z + 8.f * h2.z,
		
		-c1_2.y - c1_2.z - 2.f * c1.y * c2.y - c2_2.y - 2.f * c1.z * c2.z - c2_2.z + 
		4.f * c1.y * h.y + 4.f * c2.y * h.y - 4.f * h2.y + 4.f * c1.z * h.z + 4.f * c2.z * h.z - 4.f * h2.z,
	};
	
	Poly4 polyB = Poly4::New(polyB_C);
	
	
	OFUNC(qSplineIDist, IEVAL(float, (float), 1),
	float Eval(float l)
	{
		float term = l*l * (c1.x + c2.x - 2.f * h.x) + l * (-2.f * c1.x + 2.f * h.x) + c1.x;
		
		float discr = polyB.Eval(l);
		// discr = max(0.f, discr);
		
		term = polyA.Eval(l);
		discr = polyB.Eval(l);
		
		if(discr < 0.f)
		return 999999.f;
		else
		return term - sqrt(discr);
	})
	
	OFUNC(qSplineD1, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 		
		//D[1][f1][x]*2 Sqrt[f2[x]] + D[1][f2][x]
		float f1D1 = polyA.Eval(t, 1);
		
		float f2 = polyB.Eval(t);
		float f2D1 = polyB.Eval(t, 1);
		
		// return f1D1 * 2.f * sqrt(max(0.f, f2)) - f2D1;
		return f1D1 - f2D1 * 0.5f * rsqrt(max(0.f, f2));
	})
	
	OFUNC(qSplineD2, IEVAL(float, (float), 1),
	float Eval(float t)
	{		
		float f1D1 = polyA.Eval(t, 1);
		float f1D2 = polyA.Eval(t, 2);
		
		float f2 = polyB.Eval(t);
		float f2D1 = polyB.Eval(t, 1);
		float f2D2 = polyB.Eval(t, 2);
		
		f2 = max(0.f, f2);

		// return Pow2(f2D2) / (4.f * pow(f2, 1.5f)) + f1D2 - f2D2 / (2.f * sqrt(f2));
		// return Pow2(f2D2) / Pow2(f2) * 0.5f + f1D2 * sqrt(f2) * 2.f - f2D2;
		// return Pow2(f2D2 / f2) * 0.5f + f1D2 * sqrt(f2) * 2.f - f2D2;
		return (Pow2(f2D2 / f2) * 0.5f - f2D2) * rsqrt(f2) + f1D2 * 2.f ;		
	})
	
	OFUNC(qSplineD3, IEVAL(float, (float), 1),
	float Eval(float t)
	{	
		float f2 = polyB.Eval(t);
		float f2D1 = polyB.Eval(t, 1);
		float f2D2 = polyB.Eval(t, 2);
		float f2D3 = polyB.Eval(t, 3);
		
		f2 = max(0.f, f2);	   
	   
		// return 
		// f2D1 * f2D2 / (2.f * pow(f2, 1.5f)) + 
		// 0.5f * f2D1 * (-3.f * Pow2(f2D1) / (4.f * pow(f2, 2.5f)) + f2D2 / (2.f * pow(f2, 1.5f)) ) + 
		// f1D3 - 
		// f2D3 / (2.f * sqrt(f2));
		
		// return 
		// f2D1 * (f2D2 / Pow2(f2) - Pow2(f2D1) / Pow3(f2) * 0.75f + f2D2 / Pow2(f2) * 0.5f ) + 
		// f1D3 * sqrt(f2) * 2.f - 
		// f2D3;
		
		// return (f2D1 * f2D2 / Pow2(f2) - Pow3(f2D1 / f2) * 0.75f - f2D3) * rsqrt(f2) + f1D3 * 2.f;
		return f2D1 * f2D2 / Pow2(f2) - Pow3(f2D1 / f2) * 0.75f - f2D3;	
	})
	
	
	float isDpos = 0.f;
	
	
	//region binary roots
	// if(0)
	{
		isDpos = 0.f;

		float2 ival = float2(0.f, 1.f);
		// ival = float2(-r, 1.f + r);		
		
		OFUNC(quadPolyRoots, INull,
		float2 Eval(float c, float b, float a)
		{		
			float discr = b * b - 4.f * a * c;
			discr = max(0.f, discr);
			
			return (-b + float2(1.f, -1.f) * sqrt(discr)) / (2.f * a);
			// return -2.f * c / (b + float2(1.f, -1.f) * sqrt(discr));
		})
		
		OFUNC(binRootFinder, INull,
		float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
		{
			// float n = func.Eval(x1) < 0.f ? x1 : x2;
			// float p = x1 + x2 - n;
			
			if(func.Eval(n) > 0.f) return n;
			if(func.Eval(p) < 0.f) return p;		
			
			[loop]
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
		
		float2 d2Roots = quadPolyRoots.Eval(polyB.C[2], 3.f * polyB.C[3], 6.f * polyB.C[4]);
		
		// ival.x = min(ival.x, d2Roots.x);
		// ival.y = max(ival.y, d2Roots.y);
		d2Roots = saturate(d2Roots);
		// d2Roots = clamp(d2Roots, ival.x, ival.y);
		
		uint iCount1 = 10;
		uint iCount2 = 20;
		
		float3 d1Roots = 0.f;
		
		d1Roots.x = binRootFinder.Eval(d2Roots.x, ival.x, polyB.New(1), iCount1);
		d1Roots.y = binRootFinder.Eval(d2Roots.x, d2Roots.y, polyB.New(1), iCount1);
		d1Roots.z = binRootFinder.Eval(ival.y, d2Roots.y, polyB.New(1), iCount1);
		
		float3 fd;
		fd.x = polyB.Eval(d1Roots.x);
		fd.y = polyB.Eval(d1Roots.y);
		fd.z = polyB.Eval(d1Roots.z);
		
		float2 fd0;
		fd0.x = polyB.Eval(ival.x);
		fd0.y = polyB.Eval(ival.y);
		
		float2 fdT;
		fdT.x = qSplineIDist.Eval(d1Roots.x);
		fdT.y = qSplineIDist.Eval(d1Roots.z);
		// d1Roots = saturate(d1Roots);
		
		float rootType = 0.f;

		float4 rootsAlt = 0.f;
		float4 roots = 0.f;
		{
			float n1 = ival.x;
			float p1 = d1Roots.x;
			float p2 = d1Roots.x;
			float n2 = d1Roots.y;
			
			
			if(fd.x > 0 || fd.z > 0)
			{
				rootType = 1.f;
				
				if(fd.x > 0 && fd.z > 0)
				{
					if(fd.y < 0)
					{
						rootType = 2.f;
						
						if(fdT.x > fdT.y)
						{
							n1 = d1Roots.y;
							p1 = p2 = d1Roots.z;
							n2 = ival.y;
						}
					}
					else
					{
						rootType = 3.f;
						
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
				roots.x = binRootFinder.Eval(n1, p1, polyB, iCount2);
				// roots.x = binRootFinder.Eval(lerp(n1, p1, 0.05f), lerp(p1, n1, 0.05f), poly4, iCount2);
				
				// if(n2 == ival.y && fd0.y > 0.f) roots.z = n2; else
				roots.z = binRootFinder.Eval(n2, p2, polyB, iCount2);
				// roots.z = binRootFinder.Eval(lerp(n2, p2, 0.05f), lerp(p2, n2, 0.05f), poly4, iCount2);
				
				if(rootType == 2.f)
				{
					if(fdT.x > fdT.y)
					{
						n1 = d1Roots.x;
						p1 = p2 = d1Roots.x;
						n2 = ival.y;
					}
					else
					{
						n1 = d1Roots.y;
						p1 = p2 = d1Roots.z;
						n2 = ival.y;
					}
						
					rootsAlt.x = binRootFinder.Eval(n1, p1, polyB, iCount2);
					rootsAlt.z = binRootFinder.Eval(n2, p2, polyB, iCount2);
				}
			}
			
			// if(polyB.Eval(d1Roots.x) < 0.f)// || d1Roots.x < 0.f
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
		
		// roots = saturate(roots);
		
		float err = 0.001f;
		
		// if(roots.x >= 0.f && roots.x <= 1.f && abs(polyB.Eval(roots.x)) < err) isDpos = 1.f;
		// if(roots.z >= 0.f && roots.z <= 1.f && abs(polyB.Eval(roots.z)) < err) isDpos = 1.f;
		
		// if(abs(polyB.Eval(roots.x)) < err) isDpos = 1.f;
		// if(abs(polyB.Eval(roots.z)) < err) isDpos = 1.f;
		
		// if(polyB.Eval(roots.x) >= -err) isDpos = 1.f;
		// if(polyB.Eval(roots.z) >= -err) isDpos = 1.f;
		// if(qSplineIDist.Eval(roots.x) < 0.f) isDpos = 0.f;
		
		// if(roots.x > 0.9f || roots.x < 0.f) isDpos = 0.f;
		fCol = roots.z;
		
		//region finalize
		{
			float l1, l2;
			
			// roots.x = rootsAlt.x = 0.f;
			// roots.z = rootsAlt.z = 1.f;
			// rootType = 3.f;
			

			// l2 = saturate(l2);
			const uint count = 10;
			
			// if(0)
			if(rootType == 3.f)
			{
				float rootD3 = binRootFinder.Eval(roots.x, roots.z, qSplineD3, count);
				
				float2 rootsD2;
				rootsD2.x = binRootFinder.Eval(rootD3, roots.x, qSplineD2, count);
				rootsD2.y = binRootFinder.Eval(rootD3, roots.z, qSplineD2, count);
				
				l1 = binRootFinder.Eval(roots.x, rootsD2.x, qSplineD1, count);
				l2 = binRootFinder.Eval(rootsD2.y, roots.z, qSplineD1, count);
			}
			else
			{
				l1 = binRootFinder.Eval(roots.x, roots.z, qSplineD1, count);
				l2 = binRootFinder.Eval(rootsAlt.x, rootsAlt.z, qSplineD1, count);
			}
			
			float t = qSplineIDist.Eval(l1);
			float t2 = qSplineIDist.Eval(l2);				
			
			if(t2 < t)
			{
				t = t2; l1 = l2;
			}
			
			
			float3 sp = qSpline.Eval(l1);
			sp = mul(sp, rayMat); sp += rayStart;
			
			float3 ip = rayStart + rayDir * t;
				
			float3 n = normalize(ip - sp);
			fCol = n;// * 0.5f + 0.5f;
			
			float3 viewVec = normalize(rayStart - ip);
			
			fCol = (n * 0.5f + 0.5f);
			// fCol *= saturate(dot(viewVec, n));
			// return t * 0.02f;
			fCol = pow(fCol, 2.2f);
			fCol *= 1.f - saturate(Draw::ApplyDistScale(sqrt(min(-fd.x + r2, -fd.z + r2)) - r) + 1.f);
			
			hit.Depth = t;
			hit.Col = fCol;
			
		}
		//endregion
		
		// if(rootType.x != 3.f) fCol *= 0.75f;
		// fCol *= rootType * 0.333f;
		
		// if(!isDpos) return 0.f;
		
		// return fCol;				
		// return GammaEncode(fCol);
		
	}
	//endregion
	

	
	return hit;
}

#endif