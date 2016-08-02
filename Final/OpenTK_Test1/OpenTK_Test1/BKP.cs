using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;
using OpenTK_Test1;
using Coco;

namespace OpenTK_Test1_BKP
{
    class Program : GameWindow
    {
        string projectPath = @"..\..\";

        float variableScale = 0.5f;
        float camPhi, camTheta, camZoom;

        double time;
        Matrix4 projMat, viewProjMat;
        ShaderBuffer ubuffer;
        Shader shader;
        IndexBuffer indexBuffer;
        VertexBuffer<Vector3> vertexBuffer;
        PerfQuery query;


        protected override void OnLoad(EventArgs e)
        {
            GL.ClearColor(Color4.Brown);

            GL.Enable(EnableCap.DepthTest);
            GL.DepthFunc(DepthFunction.Less);

            //GL.Enable(EnableCap.CullFace);
            //GL.CullFace(CullFaceMode.Back);

            vertexBuffer = new VertexBuffer<Vector3>
            (
                new Vector3(-1f, +1f, +1f),
                new Vector3(-1f, -1f, +1f),
                new Vector3(+1f, +1f, +1f),
                new Vector3(+1f, -1f, +1f),
                new Vector3(+1f, -1f, -1f),
                new Vector3(+1f, +1f, -1f),
                new Vector3(-1f, +1f, -1f),
                new Vector3(-1f, -1f, -1f)
            );


            //indexBuffer = new IndexBuffer
            //(
            //    0, 1, 3, 2,
            //    2, 3, 4, 5,
            //    5, 4, 7, 6,
            //    6, 7, 1, 0,
            //    0, 2, 5, 6,
            //    1, 7, 4, 3
            //);

            indexBuffer = IndexBuffer.CreateQuadCylinder(8, 16);


            string shaderPath = projectPath + "TessTriCube";
            string shaderPath2 = projectPath + "TessQuadCube";
            string shaderPath3 = projectPath + "TessCyl";

            shader = new Shader(
            new Shader.Desc(ShaderType.VertexShader, shaderPath3 + ".vs"),
            new Shader.Desc(ShaderType.TessControlShader, shaderPath2 + ".tc"),
                //new Shader.Desc(ShaderType.TessControlShader, shaderPath + ".tc"),
            new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath3 + ".te"),
                //new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath + ".te"),
            new Shader.Desc(ShaderType.FragmentShader, shaderPath + ".ps"));

            shader.Bind();

            //uniformScale = GL.GetUniformLocation(shaderProgramHandle, "scale");
            //uniformMatViewProj = GL.GetUniformLocation(shader.Handle, "viewProjMat");

            ubuffer = new ShaderBuffer(ShaderBufferType.UBO, 4, BufferUsageHint.StreamDraw);
            ubuffer.BindToIndex(0);
            //ubuffer.BindToProgram(shader.Handle, "shader_data");
            shader.BindShaderBuffer(ubuffer, "shader_data");

            query = new PerfQuery();

            projMat = Matrix4.CreatePerspectiveFieldOfView((float)Math.PI * 0.25f, Width / (float)Height, 1f, 5000f);
        }



        protected override void OnUpdateFrame(FrameEventArgs e)
        {
            if (Keyboard[Key.Escape])
                Exit();

            //if (Keyboard[Key.PageUp])
            //    cameraPositionY += 100.0f * (float)e.Time;

            //if (Keyboard[Key.PageDown])
            //    cameraPositionY -= 100.0f * (float)e.Time;

            //if (Keyboard[OpenTK.Input.Key.Left])
            //    cameraPositionX -= 100.0f * (float)e.Time;

            //if (Keyboard[OpenTK.Input.Key.Right])
            //    cameraPositionX += 100.0f * (float)e.Time;

            //if (Keyboard[OpenTK.Input.Key.Up])
            //    cameraPositionZ -= 100.0f * (float)e.Time;

            //if (Keyboard[OpenTK.Input.Key.Down])
            //    cameraPositionZ += 100.0f * (float)e.Time;

            if (Keyboard[OpenTK.Input.Key.Down])
                variableScale += 1.0f * (float)e.Time;

        }

        protected override void OnResize(EventArgs e)
        {
            base.OnResize(e);

            GL.Viewport(0, 0, Width, Height);

            //float aspect_ratio = Width / (float)Height;
            //Matrix4 perpective = Matrix4.CreatePerspectiveFieldOfView(MathHelper.PiOver4, aspect_ratio, 1f, 5000f);
        }

        public static Matrix4 NewTransposed(Vector2 angle)
        {
            float s1 = (float)System.Math.Sin(angle.X);
            float s2 = (float)System.Math.Sin(angle.Y);

            float c1 = (float)System.Math.Cos(angle.X);
            float c2 = (float)System.Math.Cos(angle.Y);

            var matrix = Matrix4.Identity;

            matrix.M11 = c1;
            matrix.M21 = s1 * s2;
            matrix.M31 = c2 * s1;

            matrix.M22 = c2;
            matrix.M32 = -s2;

            matrix.M13 = -s1;
            matrix.M23 = c1 * s2;
            matrix.M33 = c1 * c2;

            return matrix;
        }

        public static Vector3 Multiply(Matrix4 matrix, Vector3 vector)
        {
            var result = Vector3.Zero;

            result.X = vector.X * matrix.M11 + vector.Y * matrix.M12 + vector.Z * matrix.M13;
            result.Y = vector.X * matrix.M21 + vector.Y * matrix.M22 + vector.Z * matrix.M23;
            result.Z = vector.X * matrix.M31 + vector.Y * matrix.M32 + vector.Z * matrix.M33;

            return result;
        }

        public static Matrix4 NewOrbitCamMat(Vector3 center, Vector2 angle, float zoom, out Vector3 position)
        {
            var matrix = NewTransposed(angle);

            position = Vector3.Zero;
            position.X = center.X + zoom * matrix.M31;
            position.Y = center.Y + zoom * matrix.M32;
            position.Z = center.Z + zoom * matrix.M33;

            matrix.Column3 = new Vector4(-Multiply(matrix, position), 1f);

            return matrix;
        }

        ShaderBufferContent content = new ShaderBufferContent(4);

        protected override void OnRenderFrame(FrameEventArgs e)
        {
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);

            //time = (time >= Math.PI) ? 0.0 : time + e.Time;
            time += e.Time;

            Vector3 camPos;
            var viewMat = NewOrbitCamMat(Vector3.Zero, new Vector2(-camPhi, -camTheta), camZoom + 4f, out camPos);

            viewMat = Matrix4.Transpose(viewMat);

            viewProjMat = viewMat * projMat;

            //variableScale = (float)(Math.Sin(time));
            // GL.Uniform1(uniformScale, variableScale);
            //GL.UniformMatrix4(uniformMatViewProj, true, ref viewProjMat);

            query.Start();

            content.Fill(0, viewProjMat);
            ubuffer.Update(content);

            // ubuffer.Update(viewProjMat);
            ubuffer.Bind();

            //GL.EnableVertexAttribArray(0);

            //vertexBuffer.Bind();
            VertexBuffer.Unbind();
            indexBuffer.Bind();

            //GL.VertexAttribPointer(0, 3, VertexAttribPointerType.Float, false, Vector3.SizeInBytes, 0);

            GL.PatchParameter(PatchParameterInt.PatchVertices, 4);

            //GL.DrawArrays(PrimitiveType.Triangles, 0, 9);
            GL.DrawElementsInstanced(PrimitiveType.Patches, indexBuffer.Count, DrawElementsType.UnsignedInt, IntPtr.Zero, 2);
            //GL.DrawElements(PrimitiveType.Patches, 3 * 12, DrawElementsType.UnsignedInt, IntPtr.Zero);
            //GL.DrawElements(PrimitiveType.Triangles, 3 * 12, DrawElementsType.UnsignedInt, IntPtr.Zero);

            //GL.DisableVertexAttribArray(0);

            query.Stop();

            this.Title = "GPU: " + "[" + query.RetrieveResultAsString() + "]";

            SwapBuffers();
        }

        public Program()
            : base(640, 480)
        {
            this.Keyboard.KeyUp += new EventHandler<KeyboardKeyEventArgs>((object sender, KeyboardKeyEventArgs e) =>
            {
                if (e.Key == Key.L)
                {
                    Shader.Unbind();
                    shader.Recreate();
                    shader.Bind();
                }
            });

            var moveCam = false;

            this.MouseDown += new EventHandler<MouseButtonEventArgs>((object sender, MouseButtonEventArgs e) =>
            {
                moveCam = e.Button == MouseButton.Left;
            });

            this.MouseUp += new EventHandler<MouseButtonEventArgs>((object sender, MouseButtonEventArgs e) =>
            {
                moveCam &= e.Button != MouseButton.Left;
            });

            this.MouseMove += new EventHandler<MouseMoveEventArgs>((object sender, MouseMoveEventArgs e) =>
            {
                if (moveCam)
                {
                    camPhi += (float)e.XDelta * 0.01f;
                    camTheta += (float)e.YDelta * 0.01f;
                }
            });

            this.MouseWheel += new EventHandler<MouseWheelEventArgs>((object sender, MouseWheelEventArgs e) =>
            {
                camZoom += e.DeltaPrecise;
            });
        }


        public static void Main_()
        {
            using (Program p = new Program())
            {
                p.Run(60);
            }
        }
    }
}
