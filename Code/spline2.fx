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



#define real float
		
//region Poyls
struct Ply1
{
	real A, B;
	
	static Ply1 New(float a, float b)
	{
		Ply1 p; p.A = a; p.B = b;
		
		return p;
	}
};

struct Ply2
{
	real A, B, C;
	
	static Ply2 New(float a, float b, float c)
	{
		Ply2 p; p.A = a; p.B = b; p.C = c;
		
		return p;
	}
};

struct Ply3
{
	real A, B, C, D;
	
	static Ply3 New(float a, float b, float c, float d)
	{
		Ply3 p; p.A = a; p.B = b; p.C = c; p.D = d;
		
		return p;
	}
};

struct Ply4
{
	real A, B, C, D, E;
	
	static Ply4 New(float a, float b, float c, float d, float e)
	{
		Ply4 p; p.A = a; p.B = b; p.C = c; p.D = d; p.E = e;
		
		return p;
	}
};

struct Ply6
{
	real A, B, C, D, E, F, G;
	
	static Ply6 New(float a, float b, float c, float d, float e, float f, float g)
	{
		Ply6 p; p.A = a; p.B = b; p.C = c; p.D = d; p.E = e; p.F = f; p.G = g;
		
		return p;
	}
	
	float Eval(real x)
	{
		// return A + B * x + C * Pow2(x) + D * Pow3(x) + E * Pow4(x) + F * Pow5(x) + G * Pow6(x);
		return (((((G * x + F) * x + E) * x + D) * x + C) * x + B) * x + A;
	}
};
//endregion
		
Ply2 Pow2i(Ply1 p)
{
	return Ply2::New(p.A * p.A, 2.0 * p.A * p.B, p.B * p.B);
}

Ply4 Pow2i(Ply2 p)
{
	return Ply4::New(p.A * p.A, 2.0 * p.A * p.B, p.B * p.B + 2.0 * p.A * p.C, 2.0 * p.B * p.C, p.C * p.C);
}

Ply6 Pow2i(Ply3 p)
{
	return Ply6::New(p.A * p.A, 2.0 * p.A * p.B, p.B * p.B + 2.0 * p.A * p.C, 2.0 * p.B * p.C + 2.0 * p.A * p.D, 
					 p.C * p.C + 2.0 * p.B * p.D, 2.0 * p.C * p.D, p.D * p.D);
}

Ply2 Muli(float v, Ply2 p)
{
	return Ply2::New(v * p.A, v * p.B, v * p.C);
}

Ply3 Muli(Ply2 p1, Ply1 p2)
{
	return Ply3::New(p1.A * p2.A, p2.A * p1.B + p1.A * p2.B, p1.B * p2.B + p2.A * p1.C, p2.B * p1.C);
}

Ply3 D0MulD1(Ply2 p)
{
	return Ply3::New(p.A * p.B, p.B * p.B + 2.0 * p.A * p.C, 3.0 * p.B * p.C, 2.0 * p.C * p.C);
}

Ply6 Muli(Ply3 p1, Ply3 p2)
{
	return Ply6::New(p1.A * p2.A,  p2.A * p1.B + p1.A * p2.B, 
					 p1.B * p2.B + p2.A * p1.C + p1.A * p2.C, 
					 p2.B * p1.C + p1.B * p2.C + p2.A * p1.D + p1.A * p2.D, 
					 p1.C * p2.C + p2.B * p1.D + p1.B * p2.D, 
					 p2.C * p1.D + p1.C * p2.D,  p1.D * p2.D);
}

Ply6 Muli(Ply4 p1, Ply2 p2)
{
	return Ply6::New(p1.A * p2.A,  p2.A * p1.B + p1.A * p2.B, 
					 p1.B * p2.B + p2.A * p1.C + p1.A * p2.C, 
					 p2.B * p1.C + p1.B * p2.C + p2.A * p1.D, 
					 p1.C * p2.C + p2.B * p1.D + p2.A * p1.E, 
					 p2.C * p1.D + p2.B * p1.E,  p2.C * p1.E);
}

Ply3 Addi(Ply3 p1, Ply3 p2)
{
	return Ply3::New(p1.A + p2.A, p1.B + p2.B, p1.C + p2.C, p1.D + p2.D);
}

Ply6 Addi(Ply6 p1, Ply6 p2)
{
	return Ply6::New(p1.A + p2.A, p1.B + p2.B, p1.C + p2.C, p1.D + p2.D, p1.E + p2.E, p1.F + p2.F, p1.G + p2.G);
}

Ply4 Addi(Ply4 p1, Ply4 p2)
{
	return Ply4::New(p1.A + p2.A, p1.B + p2.B, p1.C + p2.C, p1.D + p2.D, p1.E + p2.E);
}

Ply4 Addi(Ply4 p1, Ply2 p2)
{
	return Ply4::New(p1.A + p2.A, p1.B + p2.B, p1.C + p2.C, p1.D, p1.E);
}

Ply6 Subi(Ply6 p1, Ply6 p2)
{
	return Ply6::New(p1.A - p2.A, p1.B - p2.B, p1.C - p2.C, p1.D - p2.D, p1.E - p2.E, p1.F - p2.F, p1.G - p2.G);
}

Ply6 Subi(Ply6 p1, Ply2 p2)
{
	return Ply6::New(p1.A - p2.A, p1.B - p2.B, p1.C - p2.C, p1.D, p1.E, p1.F, p1.G);
}

Ply3 Negi(Ply3 p)
{
	return Ply3::New(-p.A, -p.B, -p.C, -p.D);
}


Hit EvalSplineISect(float3 rayStart, float3 rayDir, float3 c1, float3 c2, float3 h, float r = 0.25f, uint mode = 0, float2 uv = 0.0)
{
	Hit hit = Hit::New();
	
	float rs = 0.125;
	float rh = 2.0;
	float rt = 0.125;
	
	// float rs = 0.33;
	// float rh = 0.8;
	// float rt = 0.25;
	
	// c1 = 0.0;
	// c2 = float3(2.0, 0.0, 0.0);
	// h = float3(1.0, 2.0, 0.0);
	
	float3 c1o = c1;
	float3 c2o = c2;
	float3 ho = h;
	
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
	

	// float r2 = r * r;
	
	{		
		Ply6 hitPly;
		{
			float3 a = c1;
			float3 b = 2.0 * (h - c1);
			float3 c = c1 + c2 - 2.0 * h;
			
			float ra = rs;
			float rb = 2.0 * (rh - rs);
			float rc = rs + rt - 2.0 * rh;
			
			Ply2 px = Ply2::New(a.x, b.x, c.x);
			Ply2 py = Ply2::New(a.y, b.y, c.y);
			Ply2 pz = Ply2::New(a.z, b.z, c.z);
						
			Ply1 tx = Ply1::New(b.x, 2.0 * c.x);
			Ply1 ty = Ply1::New(b.y, 2.0 * c.y);
			Ply1 tz = Ply1::New(b.z, 2.0 * c.z);
			
			Ply2 rad = Ply2::New(ra, rb, rc);

			hitPly = 
			Addi
			(
				Addi
				(
					Pow2i(Negi(Addi(D0MulD1(py), D0MulD1(pz)))), 
					Pow2i(Muli(py, tx))
				), 
				Subi
				(
					Pow2i(Muli(pz, tx)),
					Muli(
					Pow2i(rad),// 0.33 * 0.33,
					
					Pow2i(tx))
				)
			);
			
			// hitPly = 
				// Sub
				// (
					// Add
					// (
						// Add
						// (
							// Pow2(Neg(Add(D0MulD1(py), D0MulD1(pz)))), 
							// Pow2(Mul(py, tx))
						// ), 
						// Pow2(Mul(pz, tx))
					// )
					// ,
					// Mul(
					// Pow2A(rad),
					// Pow2(tx))
				// );
				
				Ply3 pt = Addi(Addi(Muli(px, tx), Muli(py, ty)), Muli(pz, tz));
				Ply4 p2 = Addi(
				Addi(Pow2i(px), 
				Pow2i(py)), 
				Pow2i(pz));
				Ply2 tx2 = Pow2i(tx);
				
				hitPly = 
				Addi
				(
					Addi
					(
						Muli(Muli(Muli(-2.0, px), tx), pt), 
						Pow2i(pt) 
					),
					Subi
					(
						Muli(p2, tx2), 
						Muli(Pow2i(rad), tx2)
					)
				);
		}
		
		OFUNC(fooPoly, IEVAL(float, (float), 1),
		float Eval(float l)
		{
			float3 p = EvalQSpline(c1, c2, h, l);
			float3 t = EvalQSplineD1(c1, c2, h, l);
			
			float pt = dot(p, t);
			float tx2 = t.x * t.x;

			float r = EvalQSpline(rs, rt, rh, l);
			float r2 = r * r;
			
			// return Pow2(-p.y * t.y - p.z * t.z) + Pow2(p.y * t.x) + Pow2(p.z * t.x) - (r2 * tx2);

			
			// return ((-2.0 * p.x * (pt / t.x) * tx2 + Pow2(pt) + dot(p, p) * tx2) - r2 * tx2);
			return (-2.0 * p.x * pt * t.x) + Pow2(pt) + (dot(p, p) * tx2) - (r2 * tx2);

			// return (SqrLen(p - float3(pt / t.x, 0.0, 0.0)) - r2) * (t.x * t.x);
			return SqrLen(p - float3(pt / t.x, 0.0, 0.0)) - r2;
			return (-2.0 * p.x * pt + pt * pt + dot(p, p) - r2) * (t.x * t.x);
		})
		
		
		//region FFT
		OFUNC(cmplxPoly, IEVAL(Cmplx, (Cmplx), 1),
		Cmplx Eval(Cmplx l)
		{
			// l = l.Mul(0.5).Add(0.5);
			l = l.Mul(0.25).Add(0.5);
			
			float3 a = c1;
			float3 b = 2.0 * (h - c1);
			float3 c = c1 + c2 - 2.0 * h;
			
			float ra = rs;
			float rb = 2.0 * (rh - rs);
			float rc = rs + rt - 2.0 * rh;
			
			Cmplx px = EvalPoly(l, cmplx(a.x), cmplx(b.x), cmplx(c.x));
			Cmplx py = EvalPoly(l, cmplx(a.y), cmplx(b.y), cmplx(c.y));
			Cmplx pz = EvalPoly(l, cmplx(a.z), cmplx(b.z), cmplx(c.z));
			
			Cmplx tx = EvalPoly(l, cmplx(b.x), cmplx(2.0 * c.x));
			Cmplx ty = EvalPoly(l, cmplx(b.y), cmplx(2.0 * c.y));
			Cmplx tz = EvalPoly(l, cmplx(b.z), cmplx(2.0 * c.z));
			
			Cmplx r = EvalPoly(l, cmplx(ra), cmplx(rb), cmplx(rc));

			return Pow2(py.Mul(ty).Add(pz.Mul(tz)).Neg()).Add(
			Pow2(py.Mul(tx))).Add(Pow2(pz.Mul(tx))).Sub((Pow2(r).Mul(Pow2(tx))));
			// Pow2(py.Mul(tx))).Add(Pow2(pz.Mul(tx))).Sub((cmplx(r2).Mul(Pow2(tx))));
		})
		
		const uint n = 8.0;
		Cmplx a[n];
		Cmplx wn = cmplxAng(Pi2 * 0.125);
		Cmplx w = cmplx(1.0);
		
		for(uint i = 0; i < n; ++i)
		{
			a[i] = cmplxPoly.Eval(w);
			a[i] = a[i].Mul(0.125);
			w = w.Mul(wn);
		}

		Cmplx y[n];			
		EvalFFT8(a, y);
		
		float c[n];

		for(uint i = 0; i < n; ++i)
		c[i] = y[i].r;// * 0.125;
		//endregion
		
		{
			Hit hit = Hit::New();
			
			OFUNC(poly, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				// return fooPoly.Eval(x * 0.25 + 0.5);
				return EvalPoly(x, c[0], c[1], c[2], c[3], c[4], c[5], c[6]);
			})
			
			OFUNC(polyD1, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				return EvalPoly(x, c[1], c[2] * 2.0, c[3] * 3.0, c[4] * 4.0, c[5] * 5.0, c[6] * 6.0);
			})
			
			OFUNC(polyD2, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				return EvalPoly(x, c[2] * 2.0, c[3] * 6.0, c[4] * 12.0, c[5] * 20.0, c[6] * 30.0);
			})
			
			OFUNC(polyD3, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				return EvalPoly(x, c[3] * 6.0, c[4] * 24.0, c[5] * 60.0, c[6] * 120.0);
			})
			
			OFUNC(polyD4, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				return EvalPoly(x, c[4] * 24.0, c[5] * 120.0, c[6] * 360.0);
			})
			
			OFUNC(polyD5, IEVAL(float, (float), 1),
			float Eval(float x)
			{
				return EvalPoly(x, c[5] * 120.0, c[6] * 720.0);
			})
			
			OFUNC(polyD5Root, INull,
			float Eval()
			{
				return -(c[5] * 120.0) / (c[6] * 720.0);
			})
			
			OFUNC(newtonRootFinder, INull,
			float Eval(float x, IEVAL(float, (float), 1) funcD0, IEVAL(float, (float), 1) funcD1, const uint iCount)
			{
				[loop]
				for(uint i = 0; i < iCount; ++i) { x -= funcD0.Eval(x) / funcD1.Eval(x); }
				
				return x;
			})
			
			OFUNC(linRoot, INull,
			float Eval(float c0, float c1)
			{
				return -c0 / c1;
			})
			
			OFUNC(quadRoots, INull,
			float2 Eval(float c0, float c1, float c2)
			{
				float discr = c1 * c1 - 4.0 * c2 * c0;
				discr = max(0.f, discr);
				
				return (-c1 + float2(-1.0, 1.0) * sqrt(discr)) / (2.0 * c2);
				// return -2.f * c / (b + float2(1.f, -1.f) * sqrt(discr));
			})
			
			OFUNC(quadRoot, INull,
			float2 Eval(float c0, float c1, float c2, float dir)
			{
				float discr = c1 * c1 - 4.0 * c2 * c0;
				discr = max(0.f, discr);
				
				return (-c1 + dir * sqrt(discr)) / abs(c2) * 0.5;
				// return -2.f * c / (b + float2(1.f, -1.f) * sqrt(discr));
			})
			
			OFUNC(taylorO2, INull,
			float3 Eval(float a, 
			IEVAL(float, (float), 1) funcD0, 
			IEVAL(float, (float), 1) funcD1,
			IEVAL(float, (float), 1) funcD2)
			{
				float f0 = funcD0.Eval(a); // [opt] use horner scheme to eval simultaneously
				float f1 = funcD1.Eval(a);
				float f2 = funcD2.Eval(a);
				
				return float3(f0 - a * (f1 - 0.5 * (a * f2)), f1 - (a * f2), 0.5 * f2);
			})

			OFUNC(quadRootFinder, INull,
			float3 Eval(float x, 
			IEVAL(float, (float), 1) funcD0, 
			IEVAL(float, (float), 1) funcD1,
			IEVAL(float, (float), 1) funcD2, 			
			const uint iCount, float dir = -1.0)
			{
				[loop]
				for(uint i = 0; i < iCount; ++i) 
				{ 
					float3 c = taylorO2.Eval(x, funcD0, funcD1, funcD2);
					
					if(abs(c.z) < 0.0001)
					x -= funcD0.Eval(x) / funcD1.Eval(x); 
					else
					x = quadRoot.Eval(c.x, c.y, c.z, dir);
				}
				
				return x;
			})
			
			OFUNC(quadRootFinder2, INull,
			float3 Eval(float x, 
			IEVAL(float, (float), 1) funcD0, 
			IEVAL(float, (float), 1) funcD1,
			IEVAL(float, (float), 1) funcD2, 			
			const uint iCount, float dir = -1.0, float curve = 1.0)
			{
				[loop]
				for(uint i = 0; i < iCount; ++i) 
				{ 
					float3 c = taylorO2.Eval(x, funcD0, funcD1, funcD2);
					
					if(c.z * curve < 0.0001)
					x -= funcD0.Eval(x) / funcD1.Eval(x); 
					else
					x = quadRoot.Eval(c.x, c.y, c.z, dir);
				}
				
				return x;
			})
			
			OFUNC(quadMinFinder, INull,
			float3 Eval(float x, 
			IEVAL(float, (float), 1) funcD0, 
			IEVAL(float, (float), 1) funcD1,
			IEVAL(float, (float), 1) funcD2, 			
			const uint iCount)
			{
				[loop]
				for(uint i = 0; i < iCount; ++i) 
				{ 
					float3 c = taylorO2.Eval(x, funcD0, funcD1, funcD2);
					
					x = linRoot.Eval(c.y, 2.0 * c.z);
				}
				
				return x;
			})
			
			OFUNC(binRootFinder, INull,
			float Eval(float n, float p, IEVAL(float, (float), 1) func, const uint iCount)
			{
				// float n = func.Eval(x1) < 0.f ? x1 : x2;
				// float p = x1 + x2 - n;
				float ny = func.Eval(n);
				float py = func.Eval(p);
				
				// [flatten]
				if(ny > 0.0 && py < 0.0)
				{
					float t = n;
					n = p; p = t;
				}
				else
				{
					if(ny > 0.f) return n;
					if(py < 0.f) return p;		
				}
				
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
			
			float left = -2.0;
			float right = 2.0;
			
			// float m5_0 = polyD5Root.Eval();
			
			// float m4_0 = binRootFinder.Eval(m5_0, left,  polyD4, 10);
			// float m4_1 = binRootFinder.Eval(m5_0, right, polyD4, 10);
			
			// float m3_0 = binRootFinder.Eval(left, m4_0, polyD3, 10);
			// float m3_1 = binRootFinder.Eval(m4_1, m4_0, polyD3, 10);
			
			// float m2_0 = binRootFinder.Eval(m3_0, left, polyD2, 10);
			// float m2_1 = binRootFinder.Eval(m3_0, m3_1, polyD2, 10);
			
			// float m1_0 = binRootFinder.Eval(left, m2_0, polyD1, 10);
			// float m1_1 = binRootFinder.Eval(m2_1, m2_0, polyD1, 10);
			
			// float m0_0 = binRootFinder.Eval(m1_0, left, poly, 10);
			// float m0_1 = binRootFinder.Eval(m1_0, m1_1, poly, 10);
			
			// float x = m0_0;
			// float x2 = m0_1;
			
			float x0, x1;
			{
				float dir, curve;
				
				float startD = polyD1.Eval(left);
				if(startD < 0.0) 
				{
					dir = -1.0;
					curve = 1.0;
				}
				else
				{
					dir = 1.0;
					curve = -1.0;
				}
				
				x0 = quadRootFinder2.Eval(left, poly, polyD1, polyD2, 8, dir, curve);
				
				
				float sig = poly.Eval(x0);
				if(abs(sig) < 0.02)
				{
					sig = 0.0;
					
					dir = -1.0;
					curve = startD < 0.0 ? 1.0 : -1.0;
				}
				else
				{
					dir = 1.0;
					curve = startD < 0.0 ? -1.0 : 1.0;
				}
				
				
				x1 = quadRootFinder2.Eval(x0, polyD1, polyD2, polyD3, 8, dir, curve);
				
				
				if(sig == 0.0)
				{
					dir = 1.0;
					
					curve = startD < 0.0 ? 1.0 : -1.0;
				}
				else
				{
					dir = sign(sig) == sign(poly.Eval(x1)) ? 1.0 : -1.0;
					
					curve = startD < 0.0 ? -1.0 : 1.0;
				}
				
				x1 = quadRootFinder2.Eval(x1, poly, polyD1, polyD2, 8, dir, curve);
			}
			
			float x2, x3;
			{
				float dir, curve;
				
				float startD = polyD1.Eval(right);
				if(startD > 0.0) 
				{
					dir = 1.0;
					curve = 1.0;
				}
				else
				{
					dir = -1.0;
					curve = -1.0;
				}
				
				x2 = quadRootFinder2.Eval(right, poly, polyD1, polyD2, 8, dir, curve);
				
				
				float sig = poly.Eval(x2);
				if(abs(sig) < 0.02)
				{
					sig = 0.0;
					
					dir = 1.0;
					curve = startD > 0.0 ? 1.0 : -1.0;
				}
				else
				{
					dir = -1.0;
					curve = startD > 0.0 ? 1.0 : -1.0;
				}			
				
				x3 = quadRootFinder2.Eval(x2, polyD1, polyD2, polyD3, 8, dir, curve);
				
				
				if(sig == 0.0)
				{
					dir = -1.0;
					
					curve = startD > 0.0 ? 1.0 : -1.0;
				}
				else
				{
					dir = sign(sig) == sign(poly.Eval(x1)) ? -1.0 : 1.0;
					
					curve = startD > 0.0 ? -1.0 : 1.0;
				}
				
				x3 = quadRootFinder2.Eval(x3, poly, polyD1, polyD2, 8, dir, curve);
			}
			// float x2, x3;
			// {
				// float dir = 1.0;
				// if(poly.Eval(right) < 0.0) dir = -dir;
				
				// x2 = quadRootFinder.Eval(right, poly, polyD1, polyD2, 8, dir);
				
				// float curve = 1.0;
				// if(poly.Eval(x2) > 0.02)
					// curve = -curve;
				// else
					// dir = -dir;
				
				// x3 = quadRootFinder.Eval(x2, polyD1, polyD2, polyD3, 8, -dir);
				// // x3 = quadRootFinder2.Eval(x3, poly, polyD1, polyD2, 8, dir, curve);
			// }
			
			// float dir2 = 1.0;
			// if(poly.Eval(right) < 0.0) dir2 = -dir2;
			
			// float x2 = quadRootFinder.Eval(right, poly, polyD1, polyD2, 8, dir2);
			// // float x3 = quadRootFinder.Eval(x2, poly, polyD1, polyD2, 8, -dir2);
			// float x3 = quadRootFinder.Eval(x2, polyD1, polyD2, polyD3, 10, -dir2);
				  // // x3 = quadRootFinderP.Eval(x3, poly, polyD1, polyD2, 8, -dir2);
			
			OFUNC(isecDist, INull,
			float Eval(float x)
			{
				x = x * 0.25 + 0.5;
				x = saturate(x);
			
				float3 p0 = EvalQSpline(c1, c2, h, x);
				float3 t0 = EvalQSplineD1(c1, c2, h, x);
				
				return abs(dot(p0, t0) / t0.x);
			})
			
			
			float mt = 0.02;
			
			float3 rm0 = float3(isecDist.Eval(x0), x0, abs(poly.Eval(x0)) < mt && x0 >= left && x0 <= right ? 1.0 : 0.0);
			float3 rm1 = float3(isecDist.Eval(x1), x1, abs(poly.Eval(x1)) < mt && x1 >= left && x1 <= right ? 1.0 : 0.0);
			float3 rm2 = float3(isecDist.Eval(x2), x2, abs(poly.Eval(x2)) < mt && x2 >= left && x2 <= right ? 1.0 : 0.0);
			float3 rm3 = float3(isecDist.Eval(x3), x3, abs(poly.Eval(x3)) < mt && x3 >= left && x3 <= right ? 1.0 : 0.0);
			
			// float rm2 = abs(poly.Eval(x2)) < mt ? 1.0 : 0.0;
			// float rm3 = abs(poly.Eval(x3)) < mt ? 1.0 : 0.0;
			// float rm4 = abs(poly.Eval(x4)) < mt ? 1.0 : 0.0;
			
			// float r0 = x;
			// float r1 = x1;
			// float r2 = x2;
			// float r3 = x3;
			
			// if(abs(poly.Eval(r0)) > mt || r0 < 0.0 || r1 > 0.0)
			// {
				// r0 = r1; r1 = r2; r2 = r3;
			// }
			
			OFUNC(MinR, INull,
			float3 Eval(float3 rm0, float3 rm1)
			{
				if((rm0.x < rm1.x && rm0.z) || !rm1.z)
				return rm0;
				else
				return rm1;
			})
			
			float3 rm = MinR.Eval(rm0, MinR.Eval(rm1, MinR.Eval(rm2, rm3)));
			
			// float x = newtonRootFinder.Eval(-2.0, poly, polyD1, 4);
			
			// float x = newtonRootFinder.Eval(-2.0, polyD1, polyD2, 15);
			// float x2 = newtonRootFinder.Eval(x, polyD2, polyD3, 20);
				  // x2 = newtonRootFinder.Eval(x2, poly, polyD1, 20);
			
			// float x2 = newtonRootFinder.Eval(x, polyD3, polyD4, 40);
				  // x2 = newtonRootFinder.Eval(x2, polyD1, polyD2, 40);
			
			// x2 = binRootFinder.Eval(x, x2, poly, 20);
			
			// if(poly.Eval(x) <= 0.0)
			// x = binRootFinder.Eval(x, -2.0, poly, 10);
			
			
			// x = newtonRootFinder.Eval(-2.0, poly, polyD1, 10);
			
			// hit.Mask = abs(EvalPoly(x0, c[0], c[1], c[2], c[3], c[4], c[5], c[6])) < 0.001;
			
			x0 = x0 * 0.25 + 0.5;
			x0 = saturate(x0);
			
			x1 = x1 * 0.25 + 0.5;
			x1 = saturate(x1);
			
			x2 = x2 * 0.25 + 0.5;
			x2 = saturate(x2);
			
			x3 = x3 * 0.25 + 0.5;
			x3 = saturate(x3);
			
			// x =(x + 2.0) / 4.0;
			// x = 0.5;
			

			
			float xf = x0;
			{
				float3 p0 = EvalQSpline(c1, c2, h, x0);
				float3 t0 = (EvalQSplineD1(c1, c2, h, x0));
				
				
				float3 p1 = EvalQSpline(c1, c2, h, x1);
				float3 t1 = (EvalQSplineD1(c1, c2, h, x1));
				
				// if(dot(p0, t0) * t1.x > dot(p1, t1) * t0.x)
				if(abs(dot(p0, t0) / t0.x) > abs(dot(p1, t1) / t1.x))
				xf = x1;
			}
			
			hit.Mask = rm.z;
			hit.Mask = rm0.z + rm1.z + rm2.z + rm3.z;
			// hit.Mask = rm0.z;
			xf = saturate(rm.y * 0.25 + 0.5);
			
			// return rm.z;
			
			float3 spos = EvalQSpline(c1o, c2o, ho, xf);
			
			IntrsectRes res = Intrsect::RaySphere(rayStart, rayDir, 
			EvalQSpline(c1o, c2o, ho, xf), EvalQSpline(rs, rt, rh, xf));
			// c1o, 0.25);
			
			// float3 spos2 = EvalQSpline(c1o, c2o, ho, x2);
			
			// IntrsectRes res2 = Intrsect::RaySphere(rayStart, rayDir, 
			// EvalQSpline(c1o, c2o, ho, x2), EvalQSpline(0.25, 0.33, 0.8, x2));
			
			// if(res2.mask && res2.t1 < res.t1)
			// {
				// res = res2;
				// spos = spos2;
			// }
			
			// res = res2;
			// spos = spos2;
			
			float3 hpos = rayStart + rayDir * res.t1;
			hit.Mask = hit.Mask && res.mask;
			
			// hit.Mask = x < -0.01;
			// IntrsectRes = Intrsect::RaySphere(float3(-2.0, 0.0, 0.0), float3(1.0, 0.0, 0.0), 
			// EvalQSpline(c1, c2, h, x), 
			// EvalQSpline(0.25, 0.33, 0.8, x));
			

			if(hit.Mask)
			// if(rm2.z)
			{
				hit.Col = 0.5;
				hit.Col = normalize(hpos - spos);
				// hit.Col = res.t1 * 0.8;
				// hit.Col = spos;
				// hit.Col = hpos;
				// hit.Col = hpos;
				hit.Normal = 1.0;
				hit.Depth = 0.0;
			}
			
			if(mode)
			{
				Hit hit = Hit::New();
				
				OFUNC(taylor, IEVAL(float, (float), 1),
				float Eval(float x)
				{
					float a = -2.0;
					a = x1*4.0-2.0;
					return poly.Eval(a) + polyD1.Eval(a) * (x - a) + 0.5 * polyD2.Eval(a) * Pow2(x - a);
				})
				
				// float3 c = taylorO2.Eval(-2.0, poly, polyD1, polyD2);
				
				float3 col = 0.0;//lerp(saturate(0.0), 1.0, Draw::CoordSys(uv));
				col = lerp(col, Col3::R, Draw::Curve(uv, poly));
				// col = lerp(col, Col3::C, Draw::Poly(uv, c.x, c.y, c.z));
				// col = lerp(col, Col3::Y, Draw::Curve(uv, taylor) * 0.25);
				// col = lerp(col, Col3::G, Draw::Curve(uv, polyD1));
				// col = lerp(col, Col3::B, Draw::Curve(uv, polyD2));
				// col = lerp(col, Col3::C, Draw::Curve(uv, polyD3));
				// col = lerp(col, Col3::Y, Draw::Curve(uv, polyD4));
				// col = lerp(col, Col3::M, Draw::Curve(uv, polyD5));
				
				col = lerp(col, Col3::C, Draw::CircleS(uv, float2(x0*4.0-2.0, poly.Eval(x0*4.0-2.0)), 2));
				col = lerp(col, Col3::G, Draw::CircleS(uv, float2(x1*4.0-2.0, poly.Eval(x1*4.0-2.0)), 4));
				
				col = lerp(col, Col3::O, Draw::CircleS(uv, float2(x2*4.0-2.0, poly.Eval(x2*4.0-2.0)), 6));
				col = lerp(col, Col3::Y, Draw::CircleS(uv, float2(x3*4.0-2.0, poly.Eval(x3*4.0-2.0)), 8));
				
				// col = lerp(col, Col3::B, Draw::CircleS(uv, float2(m0_0, 0.0), 2));
				// col = lerp(col, Col3::C, Draw::CircleS(uv, float2(m0_1, 0.0), 4));
				
				// col = lerp(col, Col3::B, Draw::CircleS(uv, float2(x2*4.0-2.0, 0.0), 2));
				
				// hit.Mask = Draw::Curve(uv, poly);
				// hit.Mask = Draw::CoordSys(uv);
				// hit.Col = Col3::R;
				
				hit.Col = col;
				
				return hit;
			}
			
			return hit;
		}
		
		//linear roots
		{
			Hit hit = Hit::New();
			
			float count = 50.0;
			float lastS = fooPoly.Eval(0.0) > 0.0;
			lastS = 1.0;
			
			for(float i = 0.0; i < count; ++i)
			{
				// if(EvalPoly((i / count)*4.0-2.0, c[0], c[1], c[2], c[3], c[4], c[5], c[6]) < 0.0)
				// if(hitPly.Eval(i / count) < 0.0)
				if(fooPoly.Eval(i / count) < 0.0)
				{
					hit.Mask = 1.0;
					break;
				}
			}
			
			if(hit.Mask)
			{
				hit.Col = 0.5;
				hit.Normal = 1.0;
				hit.Depth = 0.0;
			}
			
			return hit;
		}
		// float3 fCol = lerp(saturate(0), 1.0, Draw::CoordSys(uv));
		
		// hit.Col = fCol;
		// return hit;
	}

	
	return hit;
}

#endif