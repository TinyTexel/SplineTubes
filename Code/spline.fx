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

// Hit hit()
// {

// }

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
	
	OFUNC(qSplineD1, IEVAL(float, (float), 1),
	float Eval(float t)
	{	 
		// float term1 = -2.f * c1.x + 2.f * t * c1.x + 2.f * t * c2.x + 2.f * h.x - 4.f * t * h.x;
		
		// float2 term21 = (t - 1.f) * c1.yz + t * c2.yz + (1.f - 2.f * t) * h.yz;
		// float2 term22 = Pow2(t - 1.f) * c1.yz + t * (t * c2.yz - 2.f * (t - 1.f) * h.yz);			 
	 
		// float discr = discrSpline.Eval(t);
		// discr = max(0.f, discr);
		
		// return term1 * (2.f * sqrt(discr)) - (-4.f * term21.x * term22.x - 4.f * term21.y * term22.y);
		
		float term1 = -c1.x + t * c1.x + t * c2.x + h.x - 2.f * t * h.x;
		
		float2 term21 = (t - 1.f) * c1.yz + t * c2.yz + (1.f - 2.f * t) * h.yz;
		float2 term22 = Pow2(t - 1.f) * c1.yz + t * (t * c2.yz - 2.f * (t - 1.f) * h.yz);			 
	 
		float discr = discrSpline.Eval(t);
		discr = max(0.f, discr);
		
		return term1 * sqrt(discr) + (term21.x * term22.x + term21.y * term22.y);
		// return term1 + (term21.x * term22.x + term21.y * term22.y) * rsqrt(discr);
	})
	
	OFUNC(qSplineD2, IEVAL(float, (float), 1),
	float Eval(float t)
	{
		float3 term1 = (t - 1.f) * c1 + t * c2 + (1.f - 2.f * t) * h;
		float3 term2 = 4.f * (c1 + c2 - 2.f * h);
		float3 term3; 
			term3.yz = Pow2(t - 1.f) * c1.yz + t * (t * c2.yz - 2.f * (t - 1.f) * h.yz);
  
		float2 termA = -8.f * Pow2(term1.yz) - term2.yz * term3.yz;
		float2 termB = -4.f * term1.yz * term3.yz;
		
		float discr = discrSpline.Eval(t);
			
		return -(termA.x + termA.y) + term1.x * (term1.y * term3.y + term1.z * term3.z) + term2.x * discr;
	})
	
	OFUNC(qSplineD3, IEVAL(float, (float), 1),
	float Eval(float t)
	{
		float3 term1;
			term1.yz =  t * Pow2(c2.yz) + (-1.f + t) * Pow2(c1.yz) + c2.yz * h.yz - 
						4.f * t * c2.yz * h.yz - 2.f * Pow2(h.yz) +                         		
						4.f * t * Pow2(h.yz) + c1.yz * ((-1 + 2.f * t) * c2.yz + (3.f - 4.f * t) * h.yz);
						
		float3 term2 = 	c1 + c2 - 2.f * h;
		
		float3 term3 = 	(-1.f + t) * c1 + t * c2 + (1.f - 2.f * t) * h;
		
		float3 term4;
			term4.yz = 	Pow2(-1.f + t) * c1.yz + t * (t * c2.yz - 2.f * (-1.f + t) * h.yz);
		
		
		float  termA = 24.f * (term1.y + term1.z);
		
		float2 termB = -8.f * Pow2(term3.yz) - 4.f * term2.yz * term4.yz;
		
		float2 termC = -4.f * term3.yz * term4.yz;
		
		
		return termA + 2.f * term3.x * (termB.x + termB.y) + 6.f * term2.x * (termC.x + termC.y);
		
		// return 0.f; 
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
	
	
	//region binary roots
	// if(0)
	{
		isDpos = 0.f;

		float2 ival = float2(0.f, 1.f);
		// ival = float2(-r, 1.f + r);
		
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
		
		// ival.x = min(ival.x, d2Roots.x);
		// ival.y = max(ival.y, d2Roots.y);
		// d2Roots = saturate(d2Roots);
		d2Roots = clamp(d2Roots, ival.x, ival.y);
		
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
				roots.x = binRootFinder.Eval(n1, p1, poly4, iCount2);
				// roots.x = binRootFinder.Eval(lerp(n1, p1, 0.05f), lerp(p1, n1, 0.05f), poly4, iCount2);
				
				// if(n2 == ival.y && fd0.y > 0.f) roots.z = n2; else
				roots.z = binRootFinder.Eval(n2, p2, poly4, iCount2);
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
						
					rootsAlt.x = binRootFinder.Eval(n1, p1, poly4, iCount2);
					rootsAlt.z = binRootFinder.Eval(n2, p2, poly4, iCount2);
				}
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
		
		//region linear roots
		if(0)
		{
			roots = 0.f;
			
			float count = 200.f;
			float rCount = rcp(count);
			
			{
				float tempF =  poly4.Eval(0.f);
				
				if(tempF > 0) roots.x = 0.f;
				else
				for(float j = 1; j <= count; ++j)
				{
					float f = poly4.Eval(j * rCount);
					
					if(f > 0.f)
					{
						// roots.x = j * rCount;
						roots.x = binRootFinder.Eval((j - 1) * rCount, j * rCount, poly4, 10);
						break;
					}
				}
			}
			
			{
				float tempF =  poly4.Eval(1.f);
				
				if(tempF > 0) roots.z = 1.f;
				else
				for(float j = count - 1; j >= 0; --j)
				{
					float f = poly4.Eval(j * rCount);
					
					if(f > 0.f)
					{
						// roots.z = j * rCount;
						roots.z = binRootFinder.Eval((j + 1) * rCount, j * rCount, poly4, 10);
						break;
					}
				}
			}
			
			rootType = 3.f;
		}
		//endregion
		
		// roots.z = binRootFinder.Eval(ival.x, d1Roots.z, poly4);
		
		// roots = saturate(roots);
		
		float err = 0.001f;
		
		// if(roots.x >= 0.f && roots.x <= 1.f && abs(discrSpline.Eval(roots.x)) < err) isDpos = 1.f;
		// if(roots.z >= 0.f && roots.z <= 1.f && abs(discrSpline.Eval(roots.z)) < err) isDpos = 1.f;
		
		// if(abs(discrSpline.Eval(roots.x)) < err) isDpos = 1.f;
		// if(abs(discrSpline.Eval(roots.z)) < err) isDpos = 1.f;
		
		// if(discrSpline.Eval(roots.x) >= -err) isDpos = 1.f;
		// if(discrSpline.Eval(roots.z) >= -err) isDpos = 1.f;
		// if(qSplineIDist.Eval(roots.x) < 0.f) isDpos = 0.f;
		
		// if(roots.x > 0.9f || roots.x < 0.f) isDpos = 0.f;
		fCol = roots.z;
		
		//region finalize
		{
		
			//region minFinder
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
			//endregion
			float l, l2;
			
			// roots.x = rootsAlt.x = 0.f;
			// roots.z = rootsAlt.z = 1.f;
			// rootType = 3.f;
			
			// l = minFinder.Eval(roots.x, roots.z, qSplineIDist, 2, 20);
			l = binRootFinder.Eval(roots.x, roots.z, qSplineD1, 10);
				  // l = saturate(l);
			
			l2 = binRootFinder.Eval(rootsAlt.x, rootsAlt.z, qSplineD1, 10);
			// l2 = saturate(l2);
			
			// if(0)
			if(rootType == 3.f)
			{
				const uint count = 10;
				
				float rootD3 = binRootFinder.Eval(roots.x, roots.z, qSplineD3, count);
				
				float2 rootsD2;
				rootsD2.x = binRootFinder.Eval(rootD3, roots.x, qSplineD2, count);
				rootsD2.y = binRootFinder.Eval(rootD3, roots.z, qSplineD2, count);
				
				l = binRootFinder.Eval(roots.x, rootsD2.x, qSplineD1, count);
				l2 = binRootFinder.Eval(rootsD2.y, roots.z, qSplineD1, count);
			}
			
			float t = qSplineIDist.Eval(l);
			float t2 = qSplineIDist.Eval(l2);				
			
			if(t2 < t)
			{
				t = t2; l = l2;
			}
			
			
			// t = t < t2 ? t : t2;
			
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
			// if(t < 0.f) fCol = 0.f;
			// l = (roots.x + roots.z) * 0.5f;
			// l = 0.f;
			// t = qSplineIDist.Eval(l);
			
			
			
			float3 sp = qSpline.Eval(l);
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
	
	// float discr = d2.y * r2 + d2.z * r2 - d2.z * s2.y + 2.f * d.y * d.z * s.y * s.z - d2.y * s2.z;
	// discr = 0.f;
	// float2 l2 = (-d.y * s.y - d.z * s.z + sqrt(discr) * float2(-1.f, 1.f)) / (d2.y + d2.z);
	
	// isDpos = discrSpline.Eval(0) >= 0.f;
	

	
	return hit;
}

#endif