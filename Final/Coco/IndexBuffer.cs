using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;

namespace Coco
{
    public class IndexBuffer
    {
        protected int handle;
        protected int count;

        public int Handle { get { return handle; } }
        public int Count { get { return count; } }


        public IndexBuffer(params uint[] indices)
        {
            this.count = indices.Length;

            GL.GenBuffers(1, out handle);

            Bind();

            GL.BufferData<UInt32>(BufferTarget.ElementArrayBuffer,
                                   new IntPtr(indices.Length * sizeof(UInt32)),
                                   indices, BufferUsageHint.StaticDraw);

            IndexBuffer.Unbind();
        }

        public void Bind()
        {
            GL.BindBuffer(BufferTarget.ElementArrayBuffer, handle);
        }

        static public void Unbind()
        {
            GL.BindBuffer(BufferTarget.ElementArrayBuffer, 0);
        }

        public void Delete()
        {
            GL.DeleteBuffer(handle);
        }

        #region CreateQuadCylinder
        static public IndexBuffer CreateQuadCylinder(int vertexCountX, int vertexCountY)
        {
            var indices = new uint[vertexCountX * vertexCountY * 4];

            var quadCountX = vertexCountX;
            var quadCountY = vertexCountY - 1;

            var quadOff = 0;

            for (var y = 0; y < quadCountY; ++y)
            {
                for (var x = 0; x < quadCountX; ++x)
                {
                    var i0 = x + y * vertexCountX;
                    var i1 = x + (y + 1) * vertexCountX;

                    int i2, i3;
                    if (x + 1 < quadCountX)
                    {
                        i2 = i1 + 1;
                        i3 = i0 + 1;
                    }
                    else
                    {
                        i2 = (y + 1) * vertexCountX;
                        i3 = y * vertexCountX;
                    }

                    indices[quadOff + 0] = (uint)i0;
                    indices[quadOff + 1] = (uint)i1;
                    indices[quadOff + 2] = (uint)i2;
                    indices[quadOff + 3] = (uint)i3;

                    quadOff += 4;
                }
            }

            return new IndexBuffer(indices);
        }
        #endregion

        #region CreateClosedQuadCylinder
        static public IndexBuffer CreateClosedQuadCylinder(int vertexCountX, int vertexCountY)
        {
            var indices = new uint[(vertexCountX * vertexCountY + vertexCountX) * 4];

            var capQuadCount = vertexCountX / 2;

            var quadCountX = vertexCountX;
            var quadCountY = vertexCountY - 1;


            var quadOff = 0;

            uint bottomCenterVertexId = (uint)(vertexCountX * vertexCountY);

            for (var i = 0; i < capQuadCount; ++i)
            {
                uint j = (uint)i * 2u;

                indices[quadOff + 0] = bottomCenterVertexId;
                indices[quadOff + 1] = j;                
                indices[quadOff + 2] = j + 1u;
                indices[quadOff + 3] = i == capQuadCount - 1 ? 0u : j + 2u;

                quadOff += 4;
            }


            for (var y = 0; y < quadCountY; ++y)
            {
                for (var x = 0; x < quadCountX; ++x)
                {
                    var i0 = x + y * vertexCountX;
                    var i1 = x + (y + 1) * vertexCountX;

                    int i2, i3;
                    if (x + 1 < quadCountX)
                    {
                        i2 = i1 + 1;
                        i3 = i0 + 1;
                    }
                    else
                    {
                        i2 = (y + 1) * vertexCountX;
                        i3 = y * vertexCountX;
                    }

                    indices[quadOff + 0] = (uint)i0;
                    indices[quadOff + 1] = (uint)i1;
                    indices[quadOff + 2] = (uint)i2;
                    indices[quadOff + 3] = (uint)i3;

                    quadOff += 4;
                }
            }


            uint topCenterVertexId = bottomCenterVertexId + 1u;
            uint topStartVertexId = bottomCenterVertexId - (uint)vertexCountX;

            for (var i = 0; i < capQuadCount; ++i)
            {
                uint j = topStartVertexId + (uint)i * 2u;

                indices[quadOff + 0] = topCenterVertexId;
                indices[quadOff + 1] = i == capQuadCount - 1 ? topStartVertexId : j + 2u;
                indices[quadOff + 2] = j + 1u;
                indices[quadOff + 3] = j;

                quadOff += 4;
            }


            return new IndexBuffer(indices);
        }
        #endregion


        #region CreateTriDisk
        static public IndexBuffer CreateTriDisk(int vertexCount)
        {
            var indices = new uint[vertexCount * 3];

            for (var i = 0; i < vertexCount; ++i)
            {
                var j = i * 3;

                indices[j + 0] = 0u;
                indices[j + 1] = (uint)(j + 2);
                indices[j + 2] = (uint)(j + 1);
            }

            return new IndexBuffer(indices);
        }
        #endregion

        #region resources
        public uint[] CubeIndices = new uint[]
        {
            0, 1, 2,
            2, 1, 3,
            3, 4, 2,
            2, 4, 5,
            5, 4, 6,
            6, 4, 7,
            7, 1, 0,
            0, 6, 7,
            7, 4, 3,
            3, 1, 7,
            6, 0, 2,
            2, 5, 6,
        };

        public uint[] QuadCubeIndices = new uint[]
        {
            0, 1, 3, 2,
            2, 3, 4, 5,
            5, 4, 7, 6,
            6, 7, 1, 0,
            0, 2, 5, 6,
            1, 7, 4, 3,
        };
        #endregion
    }
}
