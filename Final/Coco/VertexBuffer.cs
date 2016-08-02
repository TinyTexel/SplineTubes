using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;

namespace Coco
{
    static public class VertexBuffer
    {
        static public void Unbind()
        {
            GL.BindBuffer(BufferTarget.ArrayBuffer, 0);
        }
    }

    public class VertexBuffer<T> where T : struct
    {
        protected int handle;
        protected int length;

        public int Handle { get { return handle; } }
        public int Length { get { return length; } }


        public VertexBuffer(params T[] vertices)
        {
            this.length = vertices.Length;

            GL.GenBuffers(1, out handle);

            Bind();

            GL.BufferData(BufferTarget.ArrayBuffer,
                            new IntPtr(vertices.Length * Vector3.SizeInBytes),
                            vertices, BufferUsageHint.StaticDraw);

            VertexBuffer.Unbind();
        }

        public void Bind()
        {
            GL.BindBuffer(BufferTarget.ArrayBuffer, handle);
        }

        public void Delete()
        {
            GL.DeleteBuffer(handle);
        }


        #region resources
        public static Vector3[] CubeVertices3 = new[]
        {
            new Vector3(-1f, +1f, +1f),
            new Vector3(-1f, -1f, +1f),
            new Vector3(+1f, +1f, +1f),
            new Vector3(+1f, -1f, +1f),
            new Vector3(+1f, -1f, -1f),
            new Vector3(+1f, +1f, -1f),
            new Vector3(-1f, +1f, -1f),
            new Vector3(-1f, -1f, -1f)
        };
        #endregion
    }
}
