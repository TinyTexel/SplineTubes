using System;
using OpenTK.Graphics.OpenGL;

namespace Coco
{
    public struct VertexAttribute
    {
        public int Index;
        public int Size;
        public VertexAttribPointerType Type;
        public bool Normalized;
        public int Stride;
        public int Offset;

        public void Bind()
        {
            GL.EnableVertexAttribArray(Index);
            GL.VertexAttribPointer(Index, Size, Type, Normalized, Stride, Offset);
        }
    }

    public class Mesh<T> where T : struct
    {
        protected VertexAttribute[] vertexAttributes;
        protected VertexBuffer<T> vertexBuffer;
        protected IndexBuffer indexBuffer;

        public Mesh(VertexAttribute[] vertexAttributes, T[] vertices, uint[] indices)
        {
            this.vertexAttributes = (VertexAttribute[])vertexAttributes.Clone();
            
            this.vertexBuffer = new VertexBuffer<T>(vertices);
            this.indexBuffer = new IndexBuffer(indices);
        }

        public void Bind()
        {
            indexBuffer.Bind();
            vertexBuffer.Bind();


            var count = vertexAttributes.Length;

            for (var i = 0; i < count; ++i)
                vertexAttributes[i].Bind();
        }

        public void Draw()
        {
            GL.DrawElements(PrimitiveType.Triangles, indexBuffer.Count, DrawElementsType.UnsignedInt, IntPtr.Zero);
        }

        public void Draw(int indicesCountMultiplier)
        {
            GL.DrawElements(PrimitiveType.Triangles, indexBuffer.Count * indicesCountMultiplier, DrawElementsType.UnsignedInt, IntPtr.Zero);
        }

        public void DrawInstanced(int instanceCount)
        {
            GL.DrawElementsInstanced(PrimitiveType.Triangles, indexBuffer.Count, DrawElementsType.UnsignedInt, IntPtr.Zero, instanceCount);
        }

        public void Delete()
        {
            indexBuffer.Delete();
            vertexBuffer.Delete();
        }
    }
}
