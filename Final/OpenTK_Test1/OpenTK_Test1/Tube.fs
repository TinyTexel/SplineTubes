#ifndef TUBE_DEFS
	#error [Tube.fs depends on: Tube.defs]
#endif

#ifdef USE_GEOSHADER
in GS_OUT
{
	float3 Pos;
	float L;

#ifdef USE_PER_VERTEX_NORMALS	
	float3 Normal;
#endif

	flat uint Id; 

#ifdef CAPS	
	flat bool CapPatchFlag;	
#endif

#ifdef CALC_WIREFRAME
	float2 UV;	
	float3 UVW;
#endif
} In;
#else
in TE_OUT
{
	float3 Pos;
	float L;

#ifdef USE_PER_VERTEX_NORMALS	
	float3 Normal;
#endif

	flat uint Id; 

#ifdef CAPS	
	flat bool CapPatchFlag;	
#endif

#ifdef CALC_WIREFRAME
	float2 UV;	
#endif
} In;
#endif

out vec4 FragColor;
 
void main()
{	
	// float3 wPos = In.Pos;
	
	Tube tube = Tube_Load(In.Id);
	// Tube tube = CreateTestTube();

	float3 col  = EvalColCurve(tube, In.L);
	
	float3 p = EvalPosCurve(tube, In.L);
	float3 t = EvalPTanCurve(tube, In.L);
	float  r = EvalRadCurve(tube, In.L); 

	// float3 ovec = normalize(In.Pos - p);
	float3 ovec = normalize(ProjToPlane(In.Pos - p, t));

	float3 wPos = p + ovec * r;
	float3 normal = CalcSurfaceNormal(tube, ovec, In.L);
		
#ifdef CAPS
	if(In.CapPatchFlag)
	{
		wPos = In.Pos;
		normal = normalize(cross(dFdx(wPos), dFdy(wPos)));

	#ifdef INVERT_COLOR_FOR_CAPS
		// col = float3(1.0) - col;
	#endif
	}
#endif

#ifdef USE_PER_VERTEX_NORMALS
	wPos = In.Pos;
	normal = In.Normal;
#endif
	 // normal = normalize(In.Normal);
	// normal = normalize(cross(dFdx(wPos), dFdy(wPos)));
// wPos = In.Pos;
		// normal = In.Normal;
		
	float3 eyeVec = normalize(CamPos - wPos);
	
	float rim = saturate(dot(eyeVec, normal));
	
	// FragColor.rgb = float3(rim);
	// FragColor.rgb = normal*0.5+0.5;
	// return;
	
	
#ifdef CALC_HARD_NORMALS
wPos = In.Pos;
	float3 hnormal = normalize(cross(dFdx(wPos), dFdy(wPos)));
	float hrim = saturate(dot(eyeVec, hnormal));
#endif

#ifdef CALC_WIREFRAME	
	float quadEdge = min(min(RescaleToSSDist(In.UV.x), RescaleToSSDist(In.UV.y)), min(RescaleToSSDist(1.0 - In.UV.x), RescaleToSSDist(1.0 - In.UV.y)));
	float triEdge  = min(RescaleToSSDist(In.UVW.x), min(RescaleToSSDist(In.UVW.y), RescaleToSSDist(In.UVW.z)));
		
	triEdge  = saturate(saturate(triEdge - 0.125) + 0.4);
#endif

	
#if SHADING_STYLE == COLOR
	FragColor.rgb = col;
	FragColor.rgb *= rim;
	
#elif SHADING_STYLE == LIT_COLOR
	float3 lightDir = normalize(float3(1.0, 1.0, 1.0));
	
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
	// FragColor.rgb = col * (diff + col * back * 0.05) + fresnel;
	FragColor.rgb = col * (diff + col * back * 0.05 + fresnel);
	// FragColor.rgb = float3(1.0 - abs(dot(normalize(t), eyeVec)));
	// FragColor.rgb = col * (diff + back * 0.05);
	
#elif SHADING_STYLE == SMOOTH_NORMAL
	FragColor.rgb = GammaDecode(normal * 0.5 + 0.5);
	FragColor.rgb *= rim;
	
#elif SHADING_STYLE == HARD_NORMAL	
	FragColor.rgb = GammaDecode(hnormal * 0.5 + 0.5);
	FragColor.rgb *= hrim;

#elif SHADING_STYLE == HARD_NORMAL_WF
	quadEdge = saturate(quadEdge - 0.5);

	FragColor.rgb = GammaDecode(hnormal * 0.5 + 0.5);	
	FragColor.rgb *= triEdge;
	FragColor.rgb *= quadEdge;
	FragColor.rgb *= hrim;

#elif SHADING_STYLE == BLUE_WF
	// FragColor.rgb = lerp(lerp(float3(0.0, 0.2, 1.0), float3(triEdge), quadEdge),float3(1.0, 0.2, 0.1) * 0.8,  (1.0 - hrim) * (1.0 - hrim));
	quadEdge = saturate(quadEdge - 0.333);
	
	float3 scol = float3(0.0, 0.2, 1.0);
	// float3 scol = float3(1.0, 0.2, 0.1);
	// triEdge = 1.0;
	FragColor.rgb = lerp(lerp(scol, float3(1.0), (triEdge)), scol * 0.8,  (1.0 - hrim)) * quadEdge;

#elif SHADING_STYLE == WIP
	FragColor.rgb = lerp(float3(1.0), float3(1.0, 0.3, 0.1), rim);
	
#endif
		// FragColor.rgb = lerp(float3(0.3, 0.4, 1.0), float3(1.0), triEdge) * quadEdge; 

	// FragColor = float4(GammaEncode(FragColor.rgb), 1.0);
}