#ifndef TUBE_DEFS
	#error [Tube.vs depends on: Tube.defs]
#endif

in uint gl_VertexID;
in uint gl_InstanceID;
out vec3 vPos;

#define VId gl_VertexID
#define IId gl_InstanceID

out VS_OUT
{
	float3 Pos; 
	uint Id; 
	
	float3 LPos; 
	float Cull;
	
	float3 OVec1; 	
	
#ifdef CAPS	
	bool CVFlag;
#endif
	
	float EyeFac;
	float MaxEdgeLen;
} Out;


void main()
{
	uint x = VId % VRingVertexCount;
	uint y = VId / VRingVertexCount;
		
	float lx = float(x) * RcpVRingVertexCount_f * Pi2;
	float ly = float(y);
	
	float3 Pos;
	Pos.y = ly;
	Pos.x = cos(lx);
	Pos.z = sin(lx);
	
#ifdef CAPS
	bool bottomCVFlag = VId == BottomCVId;
	bool topCVFlag = VId == TopCVId;

	if(bottomCVFlag) { Pos = float3(0.0); ly = 0.0; }
	else 
	if(topCVFlag)    { Pos = float3(0.0, TopVRingId_f, 0.0); ly = TopVRingId_f; }
	
	Out.CVFlag = bottomCVFlag || topCVFlag;
#endif
	
	uint segmentOff = IId * SegmentSize;
	uint tubeOff = uint(TubesBuffer[segmentOff].x);
	
	Tube tube = Tube_Load(tubeOff);
	
	Out.Id = tubeOff;
	Out.Pos = Pos.xyz;
	Out.LPos = Pos;


	uint lId = y;
	float l = Pos.y * RcpTopVRingId_f; 
	{
		float lId_f = lId;
		
		
		float3 ovec1;
		float3 ovec2;
		float3 tangent;
		
		// #define EVAL_TFRAMES_MANUALLY
		
		#ifdef EVAL_TFRAMES_MANUALLY
		{
			// float3 c[4];
			// CurveHParasToPolyCoeffs(tube.Start.Pos, tube.Start.PTan, tube.End.Pos, tube.End.PTan, c);
		
			ovec1 = float3(1.0, 0.0, 0.0);
			tangent = float3(0.0, 1.0, 0.0);
			
			for(float i = 0.0; i <= lId_f; ++i)
			{
				float l0 = i * RcpTopVRingId_f;
				
				float3 t = EvalPTanCurve(tube, l0);
				// float3 t = EvalPolyD1(l0, c);
				// t = EvalCSplineTangent(tube.Pos1, tube.PTan1, tube.Pos2, tube.PTan2, l0);
				// float3 t = EvalPolyD1(l0, tube.PosA, tube.PosB, tube.PosC, tube.PosD);
				// t.z = 0.0;
				// lt = float3(0.0, 1.0, 0.0);
				// ovec1 = float3(1.0, 0.0, 0.0);
				
				// Quat q = Quat::New(tangent, t);		
				// ovec1 = Rotate(ovec1, q);
				
				ovec1 = ovec1 - t * dot(ovec1, t) / dot(t, t);
				ovec1 = normalize(ovec1);
				
				tangent = t;
			}
		}
		#else
		{
			// tangent = EvalPolyD1(l, c);
			tangent = EvalPTanCurve(tube, l);

			float4 vec = TubesBuffer[segmentOff + 1 + (lId >> 1)];	
			if((lId & 1) != 0) vec.xy = vec.zw;
		
			ovec1 = CompleteVec(float2(vec.x, abs(vec.y) - 2.0));
			if(vec.y < 0.0) ovec1.z = -ovec1.z;
		}
		#endif
		
		// ovec1 = -normalize(GetOrthoVec(tangent));
		
		ovec2 = cross(tangent, ovec1);
		ovec2 = normalize(ovec2);
		
		float3 ovec = ovec2 * Pos.z - ovec1 * Pos.x;

		float r = EvalRadCurve(tube, l);
		float3 p = EvalPosCurve(tube, l);
			
		float3 normal = CalcSurfaceNormal(tube, ovec, l);
		// Out.Normal = ovec;
	
	#ifdef CAPS
		if(Out.CVFlag) 
		{
			r = 0.0;
			
			// bool tan0ZeroCond = tube.Start.PTan.x == 0.0 && tube.Start.PTan.y == 0.0 && tube.Start.PTan.z == 0.0;
			// bool tan1ZeroCond = tube.End.PTan.x   == 0.0 && tube.End.PTan.y   == 0.0 && tube.End.PTan.z   == 0.0;
			
			// if(tan0ZeroCond && tan1ZeroCond)
			// {
				// normal = tube.End.Pos - tube.Start.Pos;
			// }
			// else if(bottomCVFlag && tan0ZeroCond)
			// {
				// normal = (tube.End.Pos - tube.End.PTan * Rcp3) - tube.Start.Pos;
			// }
			// else if(topCVFlag && tan1ZeroCond)
			// {
				// normal = tube.End.Pos - (tube.Start.Pos + tube.Start.PTan * Rcp3);
			// }
			// else
			normal = bottomCVFlag ? -tube.Start.PTan : tube.End.PTan;
			
			// if(bottomCVFlag) normal = -normal;
			normal =  normalize(normal);
		}
	#endif
		
		Out.Pos = p + ovec * r;

		float3 eyeVec = normalize(CamPos - Out.Pos);
		
		// finer tessellation for silhouettes
		Out.EyeFac = 1.0 - abs(dot(eyeVec, normal));
		Out.EyeFac = saturate(Out.EyeFac * 1.2);
		Out.EyeFac *= Out.EyeFac * Out.EyeFac;

		
		Out.OVec1 = ovec1;
	
	#ifdef USE_BACKPATCHCULLING
		// control point normal facing away from camera -> belonging patch migh get culled
		Out.Cull = dot(eyeVec, normal) < 0.0 ? 1.0 : 0.0;
		// Out.Cull = 0.0;
	#else
		Out.Cull = 0.0;
	#endif
		
	#ifdef CAPS
		// force patch cull for invisible caps
		if(bottomCVFlag && tube.Start.EndNodeFlag == 0.0) Out.Cull = 4.0;
		if(topCVFlag    && tube.End.EndNodeFlag   == 0.0) Out.Cull = 4.0;
		// if(bottomCVFlag) Out.Cull = 4.0;
		// if(topCVFlag   ) Out.Cull = 4.0;
	#endif	
	
	#ifdef USE_FINE_CAPS_TESSELLATION
		// finer tessellation for cap seam
		if((ly == 0.0 && tube.Start.EndNodeFlag != 0.0) || 
		   (ly == TopVRingId_f && tube.End.EndNodeFlag  != 0.0))
		   {
				Out.MaxEdgeLen = 4.0;// seam
		   }
		   else
		   {
				Out.MaxEdgeLen = 32.0;// default
		   }
	#else
		Out.MaxEdgeLen = 32.0;
	#endif
	}	
}