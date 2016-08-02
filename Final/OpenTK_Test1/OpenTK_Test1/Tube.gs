#ifndef TUBE_DEFS
	#error [Tube.gs depends on: Tube.defs]
#endif

layout(triangles) in;
layout (triangle_strip, max_vertices = 3) out;
 
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
} In[3];

out GS_OUT
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
} Out;
 
 
void main()
{
	// float3 normal = normalize(cross(In[1].WPos - In[0].WPos, In[2].WPos - In[0].WPos));
	
	// float2 xy0 = gl_in[0].gl_Position.xy / gl_in[0].gl_Position.w;
	// float2 xy1 = gl_in[1].gl_Position.xy / gl_in[1].gl_Position.w;
	// float2 xy2 = gl_in[2].gl_Position.xy / gl_in[2].gl_Position.w;
	
	// float vn0 = dot(In[0].Normal, CamPos - In[0].Pos);
	// float vn1 = dot(In[1].Normal, CamPos - In[1].Pos);
	// float vn2 = dot(In[2].Normal, CamPos - In[2].Pos);
	
	// bool foo = false;
	// if(SinDot(xy1 - xy0, xy2 - xy0) < 0.0 
	// && (vn0 > 0.0 || vn1 > 0.0 || vn2 > 0.0))
	// {
		// foo = true;
	// // normal = -normal;
	// }
	
	
	// float l0 = length(In[1].Pos - In[2].Pos);
	// float l1 = length(In[0].Pos - In[2].Pos);
	// float l2 = length(In[1].Pos - In[0].Pos);
	
	// float lm = max(l0, max(l1, l2));
	
	// uint lmask = 0;
		// if(lm == l0) {lmask = 0;}
		// if(lm == l1) {lmask = 1;}
		// if(lm == l2) {lmask = 2;}
		
		
	// uint fit = 0;
	// if(abs(lm - l0) < 0.001) {lmask = 0; ++fit;}
	// if(abs(lm - l1) < 0.001) {lmask = 1; ++fit;}
	// if(abs(lm - l2) < 0.001) {lmask = 2; ++fit;}

	
	for(uint i = 0; i < gl_in.length(); i++)
	{
		uint j = i;// < 2 && foo ?  1 - i : i;
		
		gl_Position = gl_in[j].gl_Position;

		Out.Pos    = In[j].Pos;
		
	#ifdef USE_PER_VERTEX_NORMALS
		Out.Normal = In[j].Normal;
	#endif
	
	#ifdef CAPS
		Out.CapPatchFlag = In[j].CapPatchFlag;
	#endif
		// Out.Normal = normal;
		Out.L      = In[j].L;
		Out.Id 	   = In[j].Id;
		
	#ifdef CALC_WIREFRAME
		Out.UV     = In[j].UV;	
		Out.UVW    = float3(ToFloat(i == 0), ToFloat(i == 1), ToFloat(i == 2));
		// Out.UVW  = float3(ToFloat(i == 0 && lmask != 0), ToFloat(i == 1 && lmask != 1), ToFloat(i == 2 && lmask != 2));
	#endif

		EmitVertex();
	}
}