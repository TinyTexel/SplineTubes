#define RENDER_TYPE_RAYCAST

using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;
using System.Collections.Generic;
using OpenTK.Audio.OpenAL;
using OpenTK.Audio;
using Coco;

namespace OpenTK_Test1
{
    class Program : GameWindow
    {
        string projectPath = @"..\..\";

#if RENDER_TYPE_RAYCAST
        const int aaSamples = 0;
#else
        const int aaSamples = 8;
#endif

        float camPhi, camTheta, camZoom;

        Matrix4 projMat, viewProjMat;

        IndexBuffer indexBuffer;
        PerfQuery query;

        ShaderBuffer perFrameUBO;
        ShaderBufferContent perFrameData;

        ShaderBuffer staticUBO;
        ShaderBufferContent staticData;

#if RENDER_TYPE_RAYCAST
        SplineTubesRaycast tubes;
#else
        SplineTubesTess tubes;
#endif
        //protected string fileName = "Bulbs";
        //protected string fileName = "Cyl";
        protected string fileName = "Fiber.trk";
        //protected string fileName = "hotRoom3D_steady_lines_2010-08-30_1.vtk";
        //protected string fileName = "Demo";

        protected override void OnLoad(EventArgs e)
        {
            TubesData.ConvertBezdatToTubes(projectPath + fileName + ".bezdat", projectPath + fileName + ".tubes", true);

            GL.Enable(EnableCap.FramebufferSrgb);
            GL.ClearColor(0.45f, 0.04f, 0.04f, 1.0f);
            GL.ClearColor(1.0f, 1.0f, 1.0f, 1.0f);

            GL.Enable(EnableCap.DepthTest);
            GL.DepthFunc(DepthFunction.Less);

            GL.Enable(EnableCap.CullFace);
            GL.CullFace(CullFaceMode.Back);


            indexBuffer = IndexBuffer.CreateQuadCylinder(8, 32);


            //string shaderPath = projectPath + "TessTriCube";
            //string shaderPath2 = projectPath + "TessQuadCube";
            //string shaderPath3 = projectPath + "TessCyl";

            //shader = new Shader(
            //new Shader.Desc(ShaderType.VertexShader, shaderPath3 + ".vs"),
            //new Shader.Desc(ShaderType.TessControlShader, shaderPath2 + ".tc"),
            //    //new Shader.Desc(ShaderType.TessControlShader, shaderPath + ".tc"),
            //new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath3 + ".te"),
            //    //new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath + ".te"),
            //new Shader.Desc(ShaderType.FragmentShader, shaderPath + ".ps"));


            Func<float, Vector3> spiral = (float x) =>
            {
                //float r0 = 1.0f;
                //float r1 = 4.0f;

                //float f = 0.1f;
                //float m = 0.1f;

                ////return new Vector3((float)Math.Cos(x * f) * r1, x * m, (float)Math.Sin(x * f) * r1) + new Vector3(x, (float)Math.Sin(x) * r0, (float)Math.Cos(x) * r0);
                ////return new Vector3(
                ////    (float)Math.Cos(x * f) * r1, 
                ////    (float)Math.Sin(x) * r0 + x * m,
                ////    (float)Math.Cos(x) * r0 + (float)Math.Sin(x * f) * r1);
                //return new Vector3(x, (float)Math.Sin(x) * r0, (float)Math.Cos(x) * r0);
                float f = 2.0f;
                return new Vector3(x, (float)Math.Sin(x * f), (float)Math.Cos(x * f));
            };

            Func<float, Vector3> spiralD = (float x) =>
            {
                var d = 0.001f;

                return (spiral(x + d) - spiral(x - d)) / (d * 2.0f) * 1.0f;
            };

            Func<float, float> wave = (float x) =>
            {
                float w = (float)Math.Sin(x) * 0.5f;

                return w * w;
            };

            Func<float, float> waveD = (float x) =>
            {
                //float w = 2.0f * (float)Math.Sin(x) * (float)Math.Cos(x);

                //return w;
                var d = 0.001f;

                return (wave(x + d) - wave(x - d)) / (d * 2.0f);// *3.0f;
            };

            //// Get your sound buffer
            //byte[] sound = GetMySounds();
            //// Place the data into a stream
            //using (var ms = new MemoryStream(sound))
            //{
            //    // Construct the sound player
            //    var player = new System.Media.SoundPlayer(ms);
            //    player.Play();
            //}

            var tubesCount = 16;
            var totalNodeCount = tubesCount + 1;
            var nodeCount = new int[] { 4, 4, 4, 3, 2 };// tubesCount = 16
            //var nodeCount = new int[] { totalNodeCount };
            //var tubes = new SplineTubes(bufferMaxElementCount: 4 * 5,         

            var filePath = projectPath + fileName + ".tubes";

            var file = new System.IO.BinaryReader(File.Open(filePath, FileMode.Open));

            var version = file.ReadInt32();
            tubesCount = file.ReadInt32();
            nodeCount = new int[tubesCount];
            for (var i = 0; i < tubesCount; ++i)
                nodeCount[i] = file.ReadInt32();

            totalNodeCount = file.ReadInt32();


#if RENDER_TYPE_RAYCAST
            tubes = new SplineTubesRaycast(shaderPath: projectPath + "Raycast",
                                            shadingStyle: RaycastTubeShadingStyle.ColorHQ,
                                            renderSettings: RaycastTubeRenderSetting.Null
                //| RaycastTubeRenderSetting.SplitIn4Segments
                | RaycastTubeRenderSetting.ShaderBufferTypeIsSSBO
                //| RaycastTubeRenderSetting.UseEarlyZ
                                            ,
                                            bufferUsageHintType: BufferUsageHint.StreamDraw,
                                            bufferMaxElementCount: -1,
                                            nodeCount: nodeCount);
#else
            tubes = new SplineTubesTess(shaderPath: projectPath + "Tube",
                                        shadingStyle: TessTubeShadingStyle.Color,
                                        renderSettings: TessTubeRenderSetting.Null
                                            | TessTubeRenderSetting.UseCaps
                //| TessTubeRenderSetting.ShaderBufferTypeIsSSBO
                                            | TessTubeRenderSetting.UseBackPatchCulling
                                            | TessTubeRenderSetting.UseFineCapsTessellation
                                            //| TessTubeRenderSetting.InvertColorsForCaps
                //| TessTubeRenderSetting.UseGeoShader
                                            ,
                                        bufferUsageHintType: BufferUsageHint.StreamDraw,
                                        bufferMaxElementCount: -1,
                                        nodeCount: nodeCount,
                                        vertexRingCount: 16,
                                        vertexRingVertexCount: 8);
#endif

            for (var i = 0; i < totalNodeCount; ++i)
            {
                TubeNode node;

                node.Pos = new Vector3(file.ReadSingle(), file.ReadSingle(), file.ReadSingle());
                node.Rad = file.ReadSingle();
                node.Col = new Vector3(file.ReadSingle(), file.ReadSingle(), file.ReadSingle());

                node.PTan = new Vector3(file.ReadSingle(), file.ReadSingle(), file.ReadSingle());
                node.RTan = file.ReadSingle();
                node.CTan = new Vector3(file.ReadSingle(), file.ReadSingle(), file.ReadSingle());

                tubes.SetNode(i, node);
            }

            if (false)
                for (var i = 0; i < totalNodeCount; ++i)
                {
                    var maxX = 20.0f;
                    var xScale = maxX / (float)tubesCount;
                    //xScale = 1.0f;

                    var x = (float)i;// / (float)tubesCount;
                    x -= maxX * 0.5f;
                    x *= xScale;


                    {
                        TubeNode node;

                        node.Pos = spiral(x);
                        node.PTan = spiralD(x) * xScale;
                        //tube.PTan1 = new Vector3(3.0f, 0.0f, 0.0f);

                        node.Rad = 0.3f + wave(x);
                        node.RTan = waveD(x) * xScale;

                        //node.Rad = 0.3f;
                        //node.RTan = 0.0f;

                        node.Col = new Vector3(0.5f);
                        //node.BaseLayerId = 0.0f;

                        node.CTan = new Vector3(6.0f, 4.0f, 3.0f);
                        //node.EndNodeFlag = 0.0f;
                        //node.EndNodeFlag = i == tubesCount ? 1.0f : 0.0f;

                        tubes.SetNode(i, node);
                    }
                }


            //#if RENDER_TYPE_RAYCAST
            //            tubes = new SplineTubesRaycast( shaderPath: projectPath + "Raycast",
            //                                            bufferMaxElementCount: 4096,
            //                                            bufferElementSize: 4,
            //                                            nodeCount: new[]{2});
            //#else
            //            tubes = new SplineTubesTess(shaderPath: projectPath + "Tube",
            //                                        bufferMaxElementCount: 4096,
            //                                        bufferElementSize: 4,
            //                                        nodeCount: new[] { 2 });
            //#endif
            //            var t = new Vector3(1, 1, 0) * 3.0f;
            //            {
            //                TubeNode node;

            //                node.Pos = new Vector3(-1, -1, 0);
            //                node.PTan = t;

            //                node.Rad = 0.1f;
            //                node.RTan = 3.0f;

            //                node.Col = new Vector3(0.5f);

            //                node.CTan = new Vector3(6.0f, 4.0f, 3.0f);

            //                tubes.SetNode(0, node);
            //            }

            //            {
            //                TubeNode node;

            //                node.Pos = new Vector3(1, -1, 0);
            //                node.PTan = t * new Vector3(1, -1, 1);

            //                node.Rad = 0.1f;
            //                node.RTan = -3.0f;

            //                node.Col = new Vector3(0.5f);

            //                node.CTan = new Vector3(6.0f, 4.0f, 3.0f);

            //                tubes.SetNode(1, node);
            //            }

            tubes.FixZeroTangents();
            tubes.SetCenterToOrigin();

#if !RENDER_TYPE_RAYCAST
            tubes.GenTangentFrames();
#endif

            tubes.BindShaderAndShaderBuffer();
            tubes.BindMeshData();


            //var dim = 2;
            //var tubesCount = dim * dim * dim;
            //var nodeCount = new int[tubesCount];

            //for (var i = 0; i < tubesCount; ++i)
            //    nodeCount[i] = 2;

            //tubes = new SplineTubes(shader: shader,
            //                        bufferMaxElementCount: 4096,
            //                        bufferElementSize: 4,
            //                        nodeCount: nodeCount);


            //for (int i = 0; i < tubesCount; ++i)
            //{
            //    var x = i % dim;
            //    var y = (i / dim) % dim;
            //    var z = (i / dim) / dim;

            //    var ovec = (new Vector3(x, y, z) - new Vector3(dim) * 0.5f) * 2.0f;


            //    {
            //        TubeNode node;

            //        node.Pos = new Vector3(0.0f, 0.0f, 0.0f) + ovec;
            //        node.PTan = new Vector3(4.0f, 4.0f, 4.0f);

            //        node.Rad = 0.3f;
            //        node.RTan = 1.0f;

            //        node.Col = new Vector3(0.5f);
            //        node.CTan = new Vector3(6.0f, 4.0f, 3.0f);

            //        tubes.SetNode(i * 2, node);
            //    }

            //    {
            //        TubeNode node;

            //        node.Pos = new Vector3(4.0f, 0.0f, 0.0f) + ovec;
            //        node.PTan = new Vector3(4.0f, 4.0f, 0.0f);

            //        node.Rad = 0.3f;
            //        node.RTan = -1.0f;

            //        node.Col = new Vector3(0.5f);
            //        node.CTan = new Vector3(6.0f, 4.0f, 3.0f);

            //        tubes.SetNode(i * 2 + 1, node);
            //    }
            //}


            perFrameData = new ShaderBufferContent(5);

            perFrameUBO = new ShaderBuffer(ShaderBufferType.UBO, perFrameData.Float4BlockCount, BufferUsageHint.StreamDraw);
            perFrameUBO.BindToIndex(0);
            tubes.Shader.BindShaderBuffer(perFrameUBO, "PerFrameBuffer");


            staticData = new ShaderBufferContent(1);
            staticData.Fill(0, 0, new Vector2(Width, Height));

            staticUBO = new ShaderBuffer(ShaderBufferType.UBO, staticData, BufferUsageHint.StaticDraw);
            staticUBO.BindToIndex(1);
            tubes.Shader.BindShaderBuffer(staticUBO, "StaticBuffer");


            query = new PerfQuery();

            UpdateProjMat();
        }

        public void UpdateProjMat()
        {
            projMat = Matrix4.CreatePerspectiveFieldOfView((float)Math.PI * 0.25f, Width / (float)Height, 0.1f, 5000f);
        }


        protected override void OnUpdateFrame(FrameEventArgs e)
        {
            if (Keyboard[Key.Escape])
                Exit();
        }

        protected override void OnResize(EventArgs e)
        {
            base.OnResize(e);

            GL.Viewport(0, 0, Width, Height);

            UpdateProjMat();

            staticData.Fill(0, 0, new Vector2(Width, Height));
            staticUBO.Update(staticData);
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

        protected bool startProfiling;
        protected bool doProfiling;
        protected int profilingFrameCounter;
        protected int profilingSetup = 0;

        //protected int width0 = 800;
        //protected int height0 = 600;

        //protected int width1 = 1024;
        //protected int height1 = 768;
        protected int width0 = 1280;
        protected int height0 = 720;

        protected int width1 = 1920;
        protected int height1 = 1080;

        protected float camTheta0 = 0.0f;
        protected float camZoom0 = 75.0f;

        protected float camTheta1 = 0.0f;
        protected float camZoom1 = 10.0f;

        protected void SwitchProfilingSetup()
        {
            switch (profilingSetup)
            {
                case 0:
                    camTheta = camTheta0;
                    camZoom = camZoom0;

                    this.Width = width0;
                    this.Height = height0;
                    break;

                case 1:
                    camTheta = camTheta1;
                    camZoom = camZoom1;

                    this.Width = width0;
                    this.Height = height0;
                    break;

                case 2:
                    camTheta = camTheta0;
                    camZoom = camZoom0;

                    this.Width = width1;
                    this.Height = height1;
                    break;

                case 3:
                    camTheta = camTheta1;
                    camZoom = camZoom1;

                    this.Width = width1;
                    this.Height = height1;
                    break;

                default: break;
            }
        }

        protected void DoProfilingStuffs()
        {
            if (!doProfiling && startProfiling)
            {
                SwitchProfilingSetup();

                camPhi = 0.0f;

                profilingFrameCounter = 0;

                doProfiling = true;
                startProfiling = false;

                ResetMsFiltering();

                return;
            }

            if (doProfiling)
            {
                startProfiling = false;

                if (profilingFrameCounter == 400)
                {
                    using (System.IO.StreamWriter file = new System.IO.StreamWriter(projectPath + "ProfilingData.txt", true))
                    {
                        var ms = Math.Round(avgMs * 10.0) * 0.1;

                        string str;

#if RENDER_TYPE_RAYCAST
                        str = "RC ";
#else
                        str = "TE ";
#endif
                        str += fileName;
                        str += " Res: " + Width.ToString() + " x " + Height.ToString();
                        str += " Zoom: " + camZoom.ToString("N1");
                        //str += " Setup: " + profilingSetup.ToString();
                        str += " MS: " + ms.ToString("N1");

                        str += " AASamples: " + aaSamples.ToString();

                        str += tubes.ShadingStyle.ToProfilingString();
                        str += tubes.RenderSettings.ToProfilingString();

                        file.WriteLine(str);
                    }

                    System.Media.SystemSounds.Beep.Play();

                    doProfiling = false;

                    return;
                }

                camPhi += 0.1f;

                ++profilingFrameCounter;
            }
        }

        protected const int msListLen = 10;
        protected double[] msList = new double[msListLen];
        protected List<double> msSortList = new List<double>(msListLen);
        protected int msListMarker = 0;
        protected double avgMs = 0;
        protected int frameCount = 0;

        public void ResetMsFiltering()
        {
            for (var i = 0; i < msListLen; ++i)
                msList[i] = 0;

            msListMarker = 0;
            avgMs = 0;
            frameCount = 0;
        }

        protected override void OnRenderFrame(FrameEventArgs e)
        {
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);

            DoProfilingStuffs();
            //time += e.Time;

            Vector3 camPos;
            var viewMat = NewOrbitCamMat(Vector3.Zero, new Vector2(-camPhi, -camTheta), camZoom + 4f, out camPos);

            viewProjMat = Matrix4.Transpose(viewMat) * projMat;


            query.Start();

            perFrameData.Fill(0, viewProjMat);
            perFrameData.Fill(4, 0, camPos);

            perFrameUBO.Update(perFrameData);

            //perFrameUBO.Bind();


            tubes.UpdateAndDrawAll();
            //tubes.DrawCaps();
            //indexBuffer.Bind();

            //GL.PatchParameter(PatchParameterInt.PatchVertices, 4);

            //GL.DrawElementsInstanced(PrimitiveType.Patches, indexBuffer.Length, DrawElementsType.UnsignedInt, IntPtr.Zero, 2);


            query.Stop();


            #region filtered perf
            ++frameCount;

            msList[msListMarker] = query.RetrieveResult();

            msListMarker = (++msListMarker) % msListLen;

            msSortList.Clear();

            var count = msListLen < frameCount ? msListLen : frameCount;

            for (var i = 0; i < count; ++i)
                msSortList.Add(msList[i]);

            msSortList.Sort();

            var medianMs = (msSortList[count / 2] + msSortList[(count + 1) / 2 - 1]) * 0.5;

            avgMs -= (avgMs - medianMs) / (double)frameCount;
            //avgMs = ((avgMs * (frameCount - 1.0)) + medianMs) / frameCount;
            #endregion

            this.Title = "GPU: " + "[" + query.RetrieveResultAsString() + "]";
            this.Title += " Filtered: " + "[" + avgMs.ToString("N2") + "]";
            this.Title += " CamPhi: " + "[" + camPhi.ToString("N2") + "]";
            this.Title += " CamTheta: " + "[" + camTheta.ToString("N2") + "]";
            this.Title += " CamZoom: " + "[" + camZoom.ToString("N2") + "]";

            SwapBuffers();
        }

        protected Random rnd = new Random();
        public Program()
            : base(640, 480, new GraphicsMode(32, 24, 0, aaSamples), "Spline Tubes")
        {
            this.KeyUp += new EventHandler<KeyboardKeyEventArgs>((object sender, KeyboardKeyEventArgs e) =>
            {
                if (e.Key == Key.L)
                {
                    tubes.RecreateShader();
                }

                if (e.Key == Key.S)
                {
                    if (GraphicsContext.CurrentContext == null)
                        throw new GraphicsContextMissingException();

                    Bitmap bmp = new Bitmap(this.ClientSize.Width, this.ClientSize.Height);

                    System.Drawing.Imaging.BitmapData data =
                        bmp.LockBits(this.ClientRectangle, System.Drawing.Imaging.ImageLockMode.WriteOnly, System.Drawing.Imaging.PixelFormat.Format24bppRgb);

                    GL.ReadPixels(0, 0, this.ClientSize.Width, this.ClientSize.Height, PixelFormat.Bgr, PixelType.UnsignedByte, data.Scan0);

                    bmp.UnlockBits(data);

                    bmp.RotateFlip(RotateFlipType.RotateNoneFlipY);

                    bmp.Save(projectPath + "SHOT" + rnd.Next(4096) + ".png", System.Drawing.Imaging.ImageFormat.Png);

                    bmp.Dispose();
                }

                if (e.Key == Key.Space)
                {
                    if (doProfiling)
                        doProfiling = false;
                    else
                        startProfiling = true;
                }

                if(e.Key == Key.Number1)
                {
                    profilingSetup = 0;
                    SwitchProfilingSetup();
                }

                if (e.Key == Key.Number2)
                {
                    profilingSetup = 1;
                    SwitchProfilingSetup();
                }

                if (e.Key == Key.Number3)
                {
                    profilingSetup = 2;
                    SwitchProfilingSetup();
                }

                if (e.Key == Key.Number4)
                {
                    profilingSetup = 3;
                    SwitchProfilingSetup();
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

                if (e.Button == MouseButton.Right)
                {
                    tubes.RecreateShader();
                }
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


        public static void Main()
        {
            using (Program p = new Program())
            {
                p.Run(60);
            }
        }
    }
}
