using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;
using System.Collections.Generic;

namespace Coco
{
    public class Shader
    {
        public struct Desc
        {
            public ShaderType Type;
            public string Path;

            public Desc(ShaderType type, string path)
            {
                this.Path = path;
                this.Type = type;
            }
        }

        public class ShaderBufferBinding
        {
            public ShaderBuffer Buffer;
            public string BlockName;
        }

        protected int handle;
        protected Desc[] descs;
        protected string versionStr;
        protected string definesStr;
        protected string commonHeaderPath;
        protected int[] shaderHandles;
        protected int shaderCount;
        protected List<ShaderBufferBinding> shaderBufferBindings = new List<ShaderBufferBinding>();

        public int Handle { get { return handle; } }

        private static string concDefinesStrs(string[] defineStrs)
        {
            if (defineStrs != null && defineStrs.Length > 0)
            {
                var res = "";
                foreach (var str in defineStrs)
                    res += str + Environment.NewLine;

                return res;
            }
            else
            {
                return null;
            }
        }

        public Shader(string versionStr, string definesStr, string commonHeaderPath, params Desc[] descs)
        {
            this.versionStr = versionStr;
            this.definesStr = definesStr;
            this.commonHeaderPath = commonHeaderPath;

            this.shaderCount = descs.Length;
            this.descs = new Desc[shaderCount];
            this.shaderHandles = new int[shaderCount];

            for (var i = 0; i < shaderCount; ++i)
                this.descs[i] = descs[i];

            create();
        }

        public Shader(string versionStr, string[] defineStrs, string commonHeaderPath, params Desc[] descs)
            : this(versionStr, Shader.concDefinesStrs(defineStrs), commonHeaderPath, descs)
        { }

        public Shader(string versionStr, string commonHeaderPath, params Desc[] descs)
            : this(versionStr, (string)null, commonHeaderPath, descs)
        { }

        public Shader(string versionStr, params Desc[] descs)
            : this(versionStr, (string)null, null, descs)
        { }

        public Shader(params Desc[] descs)
            : this(null, (string)null, null, descs)
        { }


        protected void create()
        {
            this.handle = GL.CreateProgram();

            for (var i = 0; i < shaderCount; ++i)
            {
                var shaderHandle = createFromFile(descs[i].Type, descs[i].Path);

                shaderHandles[i] = shaderHandle;

                GL.AttachShader(handle, shaderHandle);
            }

            GL.LinkProgram(handle);
            Console.WriteLine(GL.GetProgramInfoLog(handle));

            foreach (var uboBinding in shaderBufferBindings)
                rebindShaderBufferBinding(uboBinding);
        }

        public void BindShaderBuffer(ShaderBuffer buffer, string uniformBlockName)
        {
            var binding = new ShaderBufferBinding()
            {
                Buffer = buffer,
                BlockName = uniformBlockName,
            };

            shaderBufferBindings.Add(binding);

            rebindShaderBufferBinding(binding);
        }

        protected void rebindShaderBufferBinding(ShaderBufferBinding binding)
        {
            var location = binding.Buffer.GetProgramResourceIndex(handle, binding.BlockName);

            // not really nice but works for now
            if (binding.Buffer.Type == ShaderBufferType.UBO)
                GL.UniformBlockBinding(handle, location, binding.Buffer.BufferIndex);
            else
                GL.ShaderStorageBlockBinding(handle, location, binding.Buffer.BufferIndex);
        }

        public void Bind()
        {
            GL.UseProgram(handle);
        }

        static public void Unbind()
        {
            GL.UseProgram(0);
        }

        static public int Create(ShaderType type, string @source)
        {
            int handle = GL.CreateShader(type);

            GL.ShaderSource(handle, source);

            GL.CompileShader(handle);
            Console.WriteLine(GL.GetShaderInfoLog(handle));

            return handle;
        }

        static public int CreateFromFile(ShaderType type, string path)
        {
            //if (File.Exists(path))
            return Shader.Create(type, File.ReadAllText(path));
            //return Shader.Create(type, File.ReadAllText(path));
            //else
            //    return 0;// meh...
        }

        public int createFromFile(ShaderType type, string path, int oldShaderHandle = 0)
        {
            var str1 = "";
            var str3 = "";

            if (versionStr != null)
                str1 = versionStr + Environment.NewLine;

            if (!File.Exists(path))
            {
                Console.WriteLine("SHADER CREATION ERROR: Can not open " + path);

                return oldShaderHandle;
            }

            if (versionStr != null)
                str3 = File.ReadAllText(commonHeaderPath) + Environment.NewLine;

            string @source;
 
            try
            {
                source = File.ReadAllText(path);
            }
            catch(Exception e)
            {
                Console.WriteLine("SHADER CREATION ERROR: " + e.Message);

                return oldShaderHandle;
            }

            return Shader.Create(type, str1 + definesStr + str3 + source);
        }

        public void Delete()
        {
            for (var i = 0; i < shaderCount; ++i)
                GL.DeleteShader(shaderHandles[i]);

            GL.DeleteShader(handle);
        }

        public void Recreate()
        {
            Delete();
            create();
        }
    }
}
