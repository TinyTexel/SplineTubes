#ifndef RAYCAST_DEFS
	#error [Raycast.vs depends on: Raycast.defs]
#endif

void SplitCurve(uint segment, float v0, float d0, float v1, float d1,
out float v0_out, out float d0_out, out float v1_out, out float d1_out)
{
	float m_v, m_d;
	EvalCSplineMidPoint(v0, d0, v1, d1, m_v, m_d);

	if(segment == 0)
	{
		v0_out = v0;
		d0_out = d0 * 0.5;
		
		v1_out = m_v;
		d1_out = m_d * 0.5;
	}
	else
	{
		v0_out = m_v;
		d0_out = m_d * 0.5;	
		
		v1_out = v1;
		d1_out = d1 * 0.5;
	}
}

void SplitCurve(uint segment, float3 v0, float3 d0, float3 v1, float3 d1,
out float3 v0_out, out float3 d0_out, out float3 v1_out, out float3 d1_out)
{
	float3 m_v, m_d;
	EvalCSplineMidPoint(v0, d0, v1, d1, m_v, m_d);

	if(segment == 0)
	{
		v0_out = v0;
		d0_out = d0 * 0.5;
		
		v1_out = m_v;
		d1_out = m_d * 0.5;
	}
	else
	{
		v0_out = m_v;
		d0_out = m_d * 0.5;	
		
		v1_out = v1;
		d1_out = d1 * 0.5;
	}
}

#define SPLIT_CURVE0(SEG, TUBE, MEM_V, MEM_D, OUT_TUBE) SplitCurve(SEG, TUBE.Start.MEM_V, TUBE.Start.MEM_D, TUBE.End.MEM_V, TUBE.End.MEM_D, OUT_TUBE.Start.MEM_V, OUT_TUBE.Start.MEM_D, OUT_TUBE.End.MEM_V, OUT_TUBE.End.MEM_D)

Tube SplitTube0(uint segment, Tube tube)
{
	Tube subTube;
	
	SPLIT_CURVE0(segment, tube, Pos, PTan, subTube);
	SPLIT_CURVE0(segment, tube, Rad, RTan, subTube);
	SPLIT_CURVE0(segment, tube, Col, CTan, subTube);
	
	return subTube;
}


void SplitCurve(uint segment, float v0, float d0, float v1, float d1,
out float v0_out, out float h_out, out float v1_out)
{
	 v0_out = v0;
	 h_out = v0 + d0 * 0.25;
	float h1 = v1 - d1 * 0.25;
	 v1_out = (h_out + h1) * 0.5;
	
	if(segment == 1)
	{
		v0_out = v1_out;
		h_out = h1;
		v1_out = v1;
	}
}

void SplitCurve(uint segment, float3 v0, float3 d0, float3 v1, float3 d1,
out float3 v0_out, out float3 h_out, out float3 v1_out)
{
	 v0_out = v0;
	 h_out = v0 + d0 * 0.25;
	float3 h1 = v1 - d1 * 0.25;
	 v1_out = (h_out + h1) * 0.5;
	
	if(segment == 1)
	{
		v0_out = v1_out;
		h_out = h1;
		v1_out = v1;
	}
}

#define SPLIT_CURVE(SEG, TUBE, MEM_V, MEM_D, QTUBE) SplitCurve(SEG, TUBE.Start.MEM_V, TUBE.Start.MEM_D, TUBE.End.MEM_V, TUBE.End.MEM_D, QTUBE.Start.MEM_V, QTUBE.H.MEM_V, QTUBE.End.MEM_V)

QTube SplitTube(uint segment, Tube tube)
{
	QTube qTube;
	
	SPLIT_CURVE(segment, tube, Pos, PTan, qTube);
	SPLIT_CURVE(segment, tube, Rad, RTan, qTube);
	SPLIT_CURVE(segment, tube, Col, CTan, qTube);
	
	return qTube;
}

void SplinePointsToPolyCoeffs(float p0, float h, float p1, out float o_c[3])
{
	o_c[0] = p0;
	o_c[1] = -2.0 * p0 + 2.0 * h;
	o_c[2] =   p0 + p1 - 2.0 * h;
}

float EvalPoly(float x, float c0, float c1 = 0, float c2 = 0, float c3 = 0, float c4 = 0, float c5 = 0, float c6 = 0)
{
	return x * (x * (x * (x * (x * (x * c6 + c5) + c4) + c3) + c2) + c1) + c0;
}

float EvalPolyD0(float x, float c[3]) { return EvalPoly(x, c[0], c[1], c[2]); }

/* Always works if the input is non-zero.
 * Doesn’t require the input to be normalised.
 * Doesn’t normalise the output. 
 http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts */
float3 GetOrthoVec(float3 v)
{
    return abs(v.x) > abs(v.z) ? float3(-v.y, v.x, 0.0f) : float3(0.0f, -v.z, v.y);
}

layout (location = 0) 
in vec3 Pos;

in int gl_InstanceID;
#define IId gl_InstanceID
// #define VId gl_VertexID
 
out VS_OUT
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
	
} Out;


void main()
{
#ifdef SPLIT_IN_4_SEGMENTS
	uint segId = IId / 4;
	uint subSegId0 = (IId % 4) / 2;
	uint subSegId = IId % 2;
	
	#if SHADING_STYLE == COLOR_HQ || SHADING_STYLE == LIT_COLOR_HQ
		if(subSegId0 == 0)
		Out.LOff = subSegId == 0 ? 0.0 : 0.25;
		else
		Out.LOff = subSegId == 0 ? 0.5 : 0.75;
	#endif
#else
	uint segId = IId / 2;
	uint subSegId = IId % 2;
	
	#if SHADING_STYLE == COLOR_HQ || SHADING_STYLE == LIT_COLOR_HQ
		Out.LOff = subSegId == 0 ? 0.0 : 0.5;
	#endif
#endif
	
	uint segmentOff = segId * SegmentSize;
	uint tubeOff = uint(TubesBuffer[segmentOff].x);
	
	Tube tube = Tube_Load(tubeOff);
	// Tube tube = CreateTestTube();
#if SHADING_STYLE == COLOR_HQ || SHADING_STYLE == LIT_COLOR_HQ	
	Out.ColP0 = tube.Start.Col;
	Out.ColT0 = tube.Start.CTan;
	Out.ColP1 = tube.End.Col;
	Out.ColT1 = tube.End.CTan;
#endif
	
#ifdef SPLIT_IN_4_SEGMENTS
	tube = SplitTube0(subSegId0, tube);
#endif
	

	QTube qTube = SplitTube(subSegId, tube);
	
	// qTube.Start.Pos = float3(0.0, 0.0, 0.0);
	// qTube.H.Pos = float3(0.0, 1.00, 0.0);
	// qTube.End.Pos = float3(1.0, 0.0, 0.0);
	
	// qTube.Start.Rad = 0.1;
	// qTube.H.Rad = 0.1;
	// qTube.End.Rad = 0.1;
	
	float3 x, y, z;
	float xl, yl;
	bool xq = false;
	bool yq = false;
	{
		x = qTube.End.Pos - qTube.Start.Pos;
		xl = length(x);
		
		if(xl < 0.0001)
		{
			y = qTube.H.Pos - qTube.Start.Pos;
			yl = length(y);
			
			if(yl < 0.0001)
			{
				x = float3(1.0, 0.0, 0.0);
				y = float3(0.0, 1.0, 0.0);
				z = float3(0.0, 0.0, 1.0);
				
				xl = 1.0; xq = true;
				yl = 1.0; yq = true;
			}
			else
			{
				x = normalize(GetOrthoVec(x));
				xl = 1.0; xq = true;
				
				z = cross(x, y);
			}
		}
		else
		{
			y = ProjToPlane(qTube.H.Pos - qTube.Start.Pos, x);
			yl = length(y);
			
			if(yl < 0.0001)
			{
				y = normalize(GetOrthoVec(x));
				yl = 1.0; yq = true;
			}
			
			z = cross(x, y);
		}
	}	

	float3 xd = x / xl;
	float3 yd = y / yl;
	float3 zd = normalize(z);
	
	
	float xm, xp, ym, yp, zm;
	{
		float xyl = dot(qTube.H.Pos - qTube.Start.Pos, xd);
		
		float cx[3];
		SplinePointsToPolyCoeffs(0.0, xyl, xl, cx);
		
		float cy[3];
		SplinePointsToPolyCoeffs(0.0, yl, 0.0, cy);
		
		float rc[3];
		SplinePointsToPolyCoeffs(qTube.Start.Rad, qTube.H.Rad, qTube.End.Rad, rc);
		// SplinePointsToPolyCoeffs(abs(qTube.Start.Rad), abs(qTube.H.Rad), abs(qTube.End.Rad), rc);
		
		
		float c_xm[3];
		c_xm[0] = cx[0] - rc[0]; c_xm[1] = cx[1] - rc[1]; c_xm[2] = cx[2] - rc[2];
		
		float c_xp[3];
		c_xp[0] = cx[0] + rc[0]; c_xp[1] = cx[1] + rc[1]; c_xp[2] = cx[2] + rc[2];
		
		xm = min(-qTube.Start.Rad, min(xl - qTube.End.Rad, EvalPolyD0(saturate(-c_xm[1] / c_xm[2] * 0.5), c_xm)));
		xp = max( qTube.Start.Rad, max(xl + qTube.End.Rad, EvalPolyD0(saturate(-c_xp[1] / c_xp[2] * 0.5), c_xp)));
		
		
		float c_ym[3];
		c_ym[0] = cy[0] - rc[0]; c_ym[1] = cy[1] - rc[1]; c_ym[2] = cy[2] - rc[2];
		
		float c_yp[3];
		c_yp[0] = cy[0] + rc[0]; c_yp[1] = cy[1] + rc[1]; c_yp[2] = cy[2] + rc[2];
		
		ym = min(-qTube.Start.Rad, min(-qTube.End.Rad, EvalPolyD0(saturate(-c_ym[1] / c_ym[2] * 0.5), c_ym)));
		yp = max( qTube.Start.Rad, max( qTube.End.Rad, EvalPolyD0(saturate(-c_yp[1] / c_yp[2] * 0.5), c_yp)));
		
		zm = max( qTube.Start.Rad, max( qTube.End.Rad, EvalPolyD0(saturate(-rc[1] / rc[2] * 0.5), rc)));
		// zm = max(0.0, zm);
		
		
		// xm = min(0.0, xyl) - zm; 
		// xp = max(xl, xyl) + zm; 
		
		// ym = min(0.0, yl) - zm;
		// yp = max(0.0, yl) + zm;
		
		if(xq) { xm = -zm; xp = zm; }
		if(yq) { ym = -zm; yp = zm; }
	}
	
	
	float3x3 mat;
	// mat[0] = xd;
	// mat[1] = yd;
	// mat[2] = zd;
	
	mat[0] = (xp - xm) * xd;
	mat[1] = (yp - ym) * yd;
	mat[2] =  zm * 2.0 * zd;
	
	
	float3 center = qTube.Start.Pos;
	center += xd * (xm + xp) * 0.5;
	center += yd * (ym + yp) * 0.5;
	
#ifdef USE_EARLYZ
	float3 wPos = center + mul(mat, Pos.xyz * 0.5);
#else
	float3 wPos = center + mul(mat, -Pos.xyz * 0.5);
#endif
	
	Out.QTube = qTube;	
	Out.Pos = wPos;
	
	gl_Position = mul(float4(wPos, 1.0), ViewProjMat);
}