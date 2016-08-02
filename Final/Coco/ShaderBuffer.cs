using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;

namespace Coco
{
    public enum ShaderBufferType
    {
        UBO,
        SSBO,
    }

    public static class ShaderBufferTypeExtensions
    {
        public static BufferTarget ToBufferTarget(this ShaderBufferType value)
        {
            switch (value)
            {
                case ShaderBufferType.UBO: return BufferTarget.UniformBuffer;

                case ShaderBufferType.SSBO: return BufferTarget.ShaderStorageBuffer;

                default: return BufferTarget.UniformBuffer;
            }
        }

        public static BufferRangeTarget ToBufferRangeTarget(this ShaderBufferType value)
        {
            switch (value)
            {
                case ShaderBufferType.UBO: return BufferRangeTarget.UniformBuffer;

                case ShaderBufferType.SSBO: return BufferRangeTarget.ShaderStorageBuffer;

                default: return BufferRangeTarget.UniformBuffer;
            }
        }

        public static ProgramInterface ToProgramInterface(this ShaderBufferType value)
        {
            switch (value)
            {
                case ShaderBufferType.UBO: return ProgramInterface.UniformBlock;

                case ShaderBufferType.SSBO: return ProgramInterface.ShaderStorageBlock;

                default: return ProgramInterface.UniformBlock;
            }
        }
    }

 
    //public static class BufferTargetExtensions
    //{
    //    public static ShaderBufferType ToShaderBufferType(this BufferTarget value)
    //    {
    //        switch(value)
    //        {
    //            case BufferTarget.UniformBuffer: return ShaderBufferType.UBO;

    //            case BufferTarget.ShaderStorageBuffer: return ShaderBufferType.SSBO;

    //            default: throw new Exception("Can not convert from BufferTarget to ShaderBufferType: " + value.ToString());// return ShaderBufferType.UBO;
    //        }
    //    }
    //}

    public class ShaderBuffer
    {
        protected int handle;
        protected int bufferIndex = -1;
        protected ShaderBufferType type;

        protected IntPtr size;
        protected int float4BlockCount;

        protected BufferTarget bufferTargetType;
        protected BufferRangeTarget bufferRangeTargetType;
        protected ProgramInterface programInterfaceType;


        public int Handle { get { return handle; } }
        public int BufferIndex { get { return bufferIndex; } }
        public ShaderBufferType Type { get { return type; } }
        public int Float4BlockCount { get { return float4BlockCount; } }

        public BufferTarget BufferTargetType { get { return bufferTargetType; } }
        public BufferRangeTarget BufferRangeTargetType { get { return bufferRangeTargetType; } }
        public ProgramInterface ProgramInterfaceType { get { return programInterfaceType; } }


        private ShaderBuffer(ShaderBufferType type, int float4BlockCount)
        {
            this.type = type;
            this.bufferTargetType = type.ToBufferTarget();
            this.bufferRangeTargetType = type.ToBufferRangeTarget();
            this.programInterfaceType = type.ToProgramInterface();

            this.float4BlockCount = float4BlockCount;

            this.size = new IntPtr(4 * sizeof(float) * float4BlockCount);
        }

        public ShaderBuffer(ShaderBufferType type, int float4BlockCount, BufferUsageHint usage = BufferUsageHint.DynamicDraw)
            :this(type, float4BlockCount)
        {
            GL.GenBuffers(1, out handle);

            Bind();

            GL.BufferData(bufferTargetType, size, new float[float4BlockCount * 4], usage);// this one works...
            //GL.BufferData(BufferTarget.UniformBuffer, size, IntPtr.Zero, BufferUsageHint.DynamicDraw); // Request the memory to be allocated

            Unbind();
        }

        public ShaderBuffer(ShaderBufferType type, ShaderBufferContent data, BufferUsageHint usage)
            : this(type, data.Float4BlockCount)
        {
            GL.GenBuffers(1, out handle);

            Bind();

            GL.BufferData(bufferTargetType, size, data.Content, usage);

            Unbind();
        }

        public void BindToIndex(int bufferIndex)
        {
            this.bufferIndex = bufferIndex;

            GL.BindBufferBase(bufferRangeTargetType, bufferIndex, handle);// Bind to Buffer Index
        }

        public int GetProgramResourceIndex(int programHandle, string blockName)
        {
            return GL.GetProgramResourceIndex(programHandle, programInterfaceType, blockName);
        }

        //public void BindToProgram(int programHandle, string uniformBlockName)
        //{
        //    var location = GL.GetUniformBlockIndex(programHandle, uniformBlockName);
            
        //    // GL.BindBufferRange(BufferTarget.UniformBuffer, BufferIndex, ubo, IntPtr.Zero, size);
        //    GL.UniformBlockBinding(programHandle, location, bufferIndex);
        //}

        public void Bind()
        {
            GL.BindBuffer(bufferTargetType, handle);
        }

        public void Unbind()
        {
            GL.BindBuffer(bufferTargetType, 0);
        }

        public void Update<T>(T data) where T : struct
        {
            Bind();

            GL.BufferSubData(bufferTargetType, IntPtr.Zero, size, ref data);

            Unbind();
        }

        public void Update(ShaderBufferContent data)
        {
            Bind();

            GL.BufferSubData(bufferTargetType, IntPtr.Zero, data.Size, data.Content);

            Unbind();
        }

        public void Delete()
        {
            GL.DeleteBuffer(handle);
        }
    }


    public class ShaderBufferContent
    {
        public float[] Content;
        protected int float4BlockCount;
        protected IntPtr size;

        public IntPtr Size
        {
            get { return size; }
        }

        public int Float4BlockCount
        {
            get { return float4BlockCount; }
        }

        public ShaderBufferContent(int float4BlockCount)
        {
            this.float4BlockCount = float4BlockCount;

            this.size = new IntPtr(4 * sizeof(float) * float4BlockCount);

            this.Content = new float[float4BlockCount * sizeof(float)];
        }

        #region Fill

        #region Float

        public void Fill(int offset, int subOffset, float newContent) { Content[offset * 4 + subOffset] = newContent; }


        public void Fill(int offset, int subOffset, Vector2 newContent)
        {
            offset = offset * 4 + subOffset;

            Content[offset] = newContent.X;
            Content[++offset] = newContent.Y;
        }

        public void Fill(int offset, int subOffset, Vector3 newContent)
        {
            offset = offset * 4 + subOffset;

            Content[offset] = newContent.X;
            Content[++offset] = newContent.Y;
            Content[++offset] = newContent.Z;
        }

        public void Fill(int offset, Vector4 newContent)
        {
            offset *= 4;

            Content[offset] = newContent.X;
            Content[++offset] = newContent.Y;
            Content[++offset] = newContent.Z;
            Content[++offset] = newContent.W;
        }

        public void Fill(int offset, Matrix4 newContent)
        {
            offset *= 4;

            Content[offset] = newContent.M11;
            Content[++offset] = newContent.M21;
            Content[++offset] = newContent.M31;
            Content[++offset] = newContent.M41;

            Content[++offset] = newContent.M12;
            Content[++offset] = newContent.M22;
            Content[++offset] = newContent.M32;
            Content[++offset] = newContent.M42;

            Content[++offset] = newContent.M13;
            Content[++offset] = newContent.M23;
            Content[++offset] = newContent.M33;
            Content[++offset] = newContent.M43;

            Content[++offset] = newContent.M14;
            Content[++offset] = newContent.M24;
            Content[++offset] = newContent.M34;
            Content[++offset] = newContent.M44;
        }

        public void FillTransposed(int offset, Matrix4 newContent)
        {
            offset *= 4;

            Content[offset] = newContent.M11;
            Content[++offset] = newContent.M12;
            Content[++offset] = newContent.M13;
            Content[++offset] = newContent.M14;

            Content[++offset] = newContent.M21;
            Content[++offset] = newContent.M22;
            Content[++offset] = newContent.M23;
            Content[++offset] = newContent.M24;

            Content[++offset] = newContent.M31;
            Content[++offset] = newContent.M32;
            Content[++offset] = newContent.M33;
            Content[++offset] = newContent.M34;

            Content[++offset] = newContent.M41;
            Content[++offset] = newContent.M42;
            Content[++offset] = newContent.M43;
            Content[++offset] = newContent.M44;
        }

        #endregion

        #region Int

        public void Fill(int offset, int subOffset, int newContent) { unsafe { Content[offset * 4 + subOffset] = *((float*)(&newContent)); } }

        //public void Fill(int offset, int subOffset, Int2 newContent)
        //{
        //    offset = offset * 4 + subOffset;

        //    unsafe
        //    {
        //        content[offset] = *((float*)(&newContent.X));
        //        content[++offset] = *((float*)(&newContent.Y));
        //    }
        //}

        //public void Fill(int offset, int subOffset, Int3 newContent)
        //{
        //    offset = offset * 4 + subOffset;

        //    unsafe
        //    {
        //        content[offset] = *((float*)(&newContent.X));
        //        content[++offset] = *((float*)(&newContent.Y));
        //        content[++offset] = *((float*)(&newContent.Z));
        //    }
        //}

        //public void Fill(int offset, Int4 newContent)
        //{
        //    offset *= 4;

        //    unsafe
        //    {
        //        content[offset] = *((float*)(&newContent.X));
        //        content[++offset] = *((float*)(&newContent.Y));
        //        content[++offset] = *((float*)(&newContent.Z));
        //        content[++offset] = *((float*)(&newContent.W));
        //    }
        //}

        #endregion

        #endregion
    }
}
