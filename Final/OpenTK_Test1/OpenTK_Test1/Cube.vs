#version 400
 
//layout (location = 0) 
in vec3 Pos;
in int gl_VertexID;

uniform float scale;
uniform mat4 viewProjMat;
 
out vec3 vPos;

void main()
{
    vPos = Pos;

    if(gl_VertexID == 0) vPos.y += 0.1;

    vec4 pPos = vec4(vPos, 1.0) * viewProjMat;

    gl_Position = pPos;
	//gl_Position = vec4(scale * -vPos.yxz, 1.0);
}

            // int ubo;
            // GL.GenBuffers(1, out ubo);
            // GL.BindBuffer(BufferTarget.UniformBuffer, ubo);
            // GL.BufferData<Matrix4>(BufferTarget.UniformBuffer,
            //                        new IntPtr(1 * sizeof(float) * 4 * 4),
            //                        ref viewProjMat, BufferUsageHint.DynamicDraw);
            //// GL.BindBuffer(BufferTarget.UniformBuffer, 0);

            // GL.BindBuffer(BufferTarget.UniformBuffer, ubo);
            // IntPtr ptr = GL.MapBuffer(BufferTarget.UniformBuffer, BufferAccess.WriteOnly);

            // GL.UnmapBuffer(BufferTarget.UniformBuffer);
            // int block_index = GL.GetUniformBlockIndex(shaderProgramHandle, "shader_data");

            // GL.BindBufferBase(BufferRangeTarget.UniformBuffer, 2, ubo);
            // GL.UniformBlockBinding(shaderProgramHandle, block_index, 2);


            // GLuint ubo = 0;
            // glGenBuffers(1, &ubo);
            // glBindBuffer(GL_UNIFORM_BUFFER, ubo);
            // glBufferData(GL_UNIFORM_BUFFER, sizeof(shader_data), &shader_data, GL_DYNAMIC_DRAW);
            // glBindBuffer(GL_UNIFORM_BUFFER, 0);

            // glBindBuffer(GL_UNIFORM_BUFFER, gbo);
            // GLvoid* p = glMapBuffer(GL_UNIFORM_BUFFER, GL_WRITE_ONLY);
            // memcpy(p, &shader_data, sizeof(shader_data))
            // glUnmapBuffer(GL_UNIFORM_BUFFER);

            // unsigned int block_index = glGetUniformBlockIndex(program, "shader_data");

            // GLuint binding_point_index = 2;
            // glBindBufferBase(GL_UNIFORM_BUFFER, binding_point_index, ubo);

            // GLuint binding_point_index = 2;
            // glUniformBlockBinding(program, block_index, binding_point_index);