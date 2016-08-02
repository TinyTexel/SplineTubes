using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Diagnostics;
using Coco;

namespace OpenTK_Test1
{
    static public class Utils
    {
        [MethodImpl(MethodImplOptions.NoInlining)]
        public static string GetCurrentMethod(int index = 1)
        {
            StackTrace st = new StackTrace();
            StackFrame sf = st.GetFrame(index);

            return sf.GetMethod().ToString();
        }

        public static void Assert(bool condition, string conditionStr)
        {
            Debug.Assert(condition, Utils.GetCurrentMethod(2) + ": [" + conditionStr + "]");
        }
    }

    public static class Ex
    {
        public static class Quaternion
        {
            public static OpenTK.Quaternion New(OpenTK.Vector3 from, OpenTK.Vector3 to)
            {
                var quat = OpenTK.Quaternion.Identity;

                var xyz = OpenTK.Vector3.Cross(from, to);
                quat.X = xyz.X;
                quat.Y = xyz.Y;
                quat.Z = xyz.Z;
                quat.W = OpenTK.Vector3.Dot(from, to) + (float)(System.Math.Sqrt((double)(OpenTK.Vector3.Dot(from, from) * OpenTK.Vector3.Dot(to, to))));
                // quat.Val.w = dot(from, to) + 1.f;

                return quat.Normalized();
            }
        }
    }

    public struct TubeNode
    {
        public const int Size = 4 * 4;
        public const float Size_f = (float)TubeNode.Size;

        public Vector3 Pos;
        public float Rad;

        public Vector3 PTan;
        public float RTan;

        public Vector3 Col;
        //public float BaseLayerId;

        public Vector3 CTan;
        //public float EndNodeFlag;
    }

    public class SplineTubes
    {
        public static Vector3 EvalCSpline(Vector3 p1, Vector3 t1, Vector3 p2, Vector3 t2, float l)
        {
            Vector3 h1 = p1 + t1 * 0.33333f;
            Vector3 h2 = p2 - t2 * 0.33333f;

            Vector3 a1 = Vector3.Lerp(p1, h1, l);
            Vector3 a2 = Vector3.Lerp(h1, h2, l);
            Vector3 a3 = Vector3.Lerp(h2, p2, l);

            Vector3 b1 = Vector3.Lerp(a1, a2, l);
            Vector3 b2 = Vector3.Lerp(a2, a3, l);

            return Vector3.Lerp(b1, b2, l);
        }

        public static Vector3 EvalCSplineTangent(Vector3 p1, Vector3 t1, Vector3 p2, Vector3 t2, float l)
        {
            // 
            return +6.0f * p1 * (l - 1.0f) * l
                    - 6.0f * p2 * (l - 1.0f) * l
                    + t1
                    - 4.0f * l * t1
                    + 3.0f * l * l * t1
                    - 2.0f * l * t2
                    + 3.0f * l * l * t2;
        }

        protected int bufferElementSize = 4;
        protected int bufferMaxElementCount = 4096;

        protected List<ShaderBufferContent> tubesBufferContent = new List<ShaderBufferContent>(1);
        protected List<int> segmentCount = new List<int>(1);
        protected ShaderBuffer shaderBuffer;

        protected int nodeElementCount;
        protected int segmentElementCount;

        protected int tubeCount;
        protected int totalNodeCount = 0;
        protected int shaderBufferCount;

        public struct NodeMappingData
        {
            public int BufferId;
            public int Offset;
            public int SegmentOffset;
            public bool Dublicated;
            public bool EndNode;
            public bool StartNode;

            //public void ResetEndNode() { EndNode = false; }
        }

        protected List<NodeMappingData> idToMemMap = new List<NodeMappingData>();
        protected int[] nodeCount;
        protected int[] cummNodeCount;

        public void fillContent(int bufferId, int offset, TubeNode node, bool endNodeFlag, int baseLayerId)
        {
            var content = tubesBufferContent[bufferId];

            content.Fill(offset + 0, 0, node.Pos);
            content.Fill(offset + 0, 3, node.Rad);

            content.Fill(offset + 1, 0, node.PTan);
            content.Fill(offset + 1, 3, node.RTan);

            content.Fill(offset + 2, 0, node.Col);
            content.Fill(offset + 2, 3, baseLayerId);

            content.Fill(offset + 3, 0, node.CTan);
            content.Fill(offset + 3, 3, endNodeFlag ? 1.0f : 0.0f);
        }


        public void SetNode(int nodeId, TubeNode node)
        {
            var mapping = idToMemMap[nodeId];
            //var endNodeFlag = nodeId == cummNodeCount[mapping.]

            fillContent(mapping.BufferId, mapping.Offset, node, mapping.EndNode || mapping.StartNode, 0);

            if (mapping.Dublicated)
            {
                fillContent(mapping.BufferId - 1, idToMemMap[nodeId - 1].Offset + nodeElementCount, node, mapping.EndNode || mapping.StartNode, 0);
            }
        }

        //public int GetNodeId(int bufferId, int localNodeId)
        //{
        //    return cummNodeCount[bufferId] + localNodeId;
        //}


        public SplineTubes(ShaderBufferType shaderBufferType, BufferUsageHint bufferUsageHintType, int bufferMaxElementCount, int[] nodeCount, int segmentElementCount)
        {
            Utils.Assert(bufferElementSize == 1 || bufferElementSize == 2 || bufferElementSize == 4 || bufferElementSize == 8,
                        "bufferElementSize == 1 || bufferElementSize == 2 || bufferElementSize == 4 || bufferElementSize == 8");

            Utils.Assert(nodeCount.Length > 0,
                        "nodeCount.Length > 0");

            if (shaderBufferType == ShaderBufferType.UBO)
            {
                if(bufferMaxElementCount > 0)
                    bufferMaxElementCount = Math.Min(bufferMaxElementCount, GL.GetInteger(GetPName.MaxUniformBlockSize) / sizeof(float) / bufferElementSize);
                else
                    bufferMaxElementCount = GL.GetInteger(GetPName.MaxUniformBlockSize) / sizeof(float) / bufferElementSize;
            }
            else
            {
                if (bufferMaxElementCount > 0)
                    bufferMaxElementCount = Math.Min(bufferMaxElementCount, GL.GetInteger((GetPName)All.MaxShaderStorageBlockSize) / sizeof(float) / bufferElementSize);
                else
                    bufferMaxElementCount = GL.GetInteger((GetPName)All.MaxShaderStorageBlockSize) / sizeof(float) / bufferElementSize;
            }

            //if(bufferMaxElementCount < 0)
            //{
            //    if(shaderBufferType == ShaderBufferType.UBO)
            //        bufferMaxElementCount = GL.GetInteger(GetPName.MaxUniformBlockSize) / sizeof(float) / bufferElementSize;
            //    else
            //        bufferMaxElementCount = GL.GetInteger((GetPName)All.MaxShaderStorageBlockSize) / sizeof(float) / bufferElementSize;
            //}

            this.bufferMaxElementCount = bufferMaxElementCount;
            //this.bufferElementSize = bufferElementSize;

            this.nodeElementCount = TubeNode.Size / bufferElementSize;
            this.segmentElementCount = segmentElementCount;

            Utils.Assert(bufferMaxElementCount >= (nodeElementCount + segmentElementCount) * 2,
                        "bufferMaxElementCount >= (nodeElementCount + segmentElementCount) * 2");


            this.tubeCount = nodeCount.Length;


            var nodeCountCapacity = new int[tubeCount];

            this.cummNodeCount = new int[tubeCount + 1];
            this.nodeCount = new int[tubeCount];

            for (int i = 0; i < tubeCount; ++i)
            {
                Utils.Assert(nodeCount[i] > 1,
                            "nodeCount[i] > 1");

                cummNodeCount[i] = totalNodeCount;
                totalNodeCount += nodeCount[i];

                nodeCountCapacity[i] = nodeCount[i];
                this.nodeCount[i] = nodeCount[i];
            }

            cummNodeCount[tubeCount] = totalNodeCount;

            var totalNodeCountCapacity = totalNodeCount;

            var tubeId = 0;
            var bufferId = 0;

            var bufferCapacity = bufferMaxElementCount;
            var currSegmentCount = 0;
            var currNodeCount = 0;
            var maxBufferSize = 0;
            var dublicatedNode = false;

            while (totalNodeCountCapacity > 0)
            {
                var currNodeCountCapacity = nodeCountCapacity[tubeId];

                var firstNode = currNodeCountCapacity == nodeCount[tubeId];
                var lastNode = currNodeCountCapacity == 1;
                var secondlastNode = currNodeCountCapacity == 2;


                var currElementSize = nodeElementCount;
                var nextElementSize = nodeElementCount;

                if (!lastNode)
                    currElementSize += segmentElementCount;

                if (!secondlastNode)
                    nextElementSize += segmentElementCount;


                var currFits = bufferCapacity >= currElementSize;
                var nextFits = bufferCapacity >= currElementSize + nextElementSize;

                var addNewNode = currFits && (!firstNode || nextFits);


                if (addNewNode)
                {
                    bufferCapacity -= currElementSize;

                    if (nextFits || lastNode)
                    {
                        NodeMappingData mapping;
                        mapping.BufferId = bufferId;
                        mapping.Offset = currNodeCount * nodeElementCount;// not finalized
                        mapping.SegmentOffset = -1;
                        mapping.Dublicated = dublicatedNode;
                        mapping.EndNode = lastNode;
                        mapping.StartNode = firstNode;

                        idToMemMap.Add(mapping);

                        if (!lastNode)
                            ++currSegmentCount;

                        ++currNodeCount;

                        --totalNodeCountCapacity;
                        --nodeCountCapacity[tubeId];

                        if (lastNode) ++tubeId;

                        dublicatedNode = false;
                    }
                    else
                        dublicatedNode = true;
                }


                if (!addNewNode || totalNodeCountCapacity == 0)
                {
                    var bufferSize = bufferMaxElementCount - bufferCapacity;

                    Utils.Assert(bufferSize != 0,
                                "bufferSize != 0");

                    tubesBufferContent.Add(new ShaderBufferContent(bufferSize));
                    segmentCount.Add(currSegmentCount);

                    if (bufferSize > maxBufferSize) maxBufferSize = bufferSize;

                    bufferCapacity = bufferMaxElementCount;
                    currSegmentCount = 0;
                    currNodeCount = 0;

                    ++bufferId;
                }
            }


            shaderBufferCount = tubesBufferContent.Count;

            var bId = -1;
            var bOff = 0;
            var addOffset = 0;
            for (var nId = 0; nId < totalNodeCount; ++nId)
            {
                var mapping = idToMemMap[nId];

                if (mapping.BufferId > bId)
                {
                    bId = mapping.BufferId;
                    bOff = 0;

                    addOffset = segmentCount[bId] * segmentElementCount;
                }

                mapping.Offset += addOffset;
                mapping.SegmentOffset = bOff;

                idToMemMap[nId] = mapping;

                if (!mapping.EndNode)
                {
                    tubesBufferContent[bId].Fill(bOff, 0, (float)mapping.Offset);

                    bOff += segmentElementCount;
                }
            }


            shaderBuffer = new ShaderBuffer(shaderBufferType, maxBufferSize, bufferUsageHintType);
        }

        public void FixZeroTangents(int tubeId)
        {
            Utils.Assert(tubeId >= 0 && tubeId < tubeCount,
                        "tubeId >= 0 && tubeId < tubeCount");

            var baseNodeId = cummNodeCount[tubeId];
            var endNodeId = baseNodeId + nodeCount[tubeId] - 1;

            for (var nodeId = baseNodeId; nodeId < endNodeId; ++nodeId)
            {
                var mapping = idToMemMap[nodeId];

                var content = tubesBufferContent[mapping.BufferId];
                var data = content.Content;


                var off1 = mapping.Offset * 4;
                var off2 = off1 + nodeElementCount * 4;

                var pos1 = new Vector3(data[off1 + 0], data[off1 + 1], data[off1 + 2]);
                var tan1 = new Vector3(data[off1 + 4], data[off1 + 5], data[off1 + 6]);

                var pos2 = new Vector3(data[off2 + 0], data[off2 + 1], data[off2 + 2]);
                var tan2 = new Vector3(data[off2 + 4], data[off2 + 5], data[off2 + 6]);


                var tan1ZeroCond = tan1.X == 0.0 && tan1.Y == 0.0 && tan1.Z == 0.0;
                var tan2ZeroCond = tan2.X == 0.0 && tan2.Y == 0.0 && tan2.Z == 0.0;

                var small = 0.0001f;

                if (tan1ZeroCond && tan2ZeroCond)
                {
                    tan1 = Vector3.Normalize(pos1 - pos2) * small;
                    tan2 = Vector3.Normalize(pos2 - pos1) * small;
                }
                else
                {
                    if (tan1ZeroCond)
                        tan1 = Vector3.Normalize(pos1 - (pos2 - tan2 / 3.0f)) * small;

                    if (tan2ZeroCond)
                        tan2 = Vector3.Normalize(pos2 - (pos1 + tan1 / 3.0f)) * small;
                }

                if (tan1ZeroCond)
                { data[off1 + 4] = tan1.X; data[off1 + 5] = tan1.Y; data[off1 + 6] = tan1.Z; }

                if (tan2ZeroCond)
                { data[off2 + 4] = tan2.X; data[off2 + 5] = tan2.Y; data[off2 + 6] = tan2.Z; }
            }
        }

        public void FixZeroTangents()
        {
            for (var tubeId = 0; tubeId < tubeCount; ++tubeId)
                FixZeroTangents(tubeId);
        }

        public void SetCenter(Vector3 newCenter)
        {
            var center = Vector3.Zero;
            var totalSegmentCount = 0.0f;

            for (var tubeId = 0; tubeId < tubeCount; ++tubeId)
            {
                var baseNodeId = cummNodeCount[tubeId];
                var endNodeId = baseNodeId + nodeCount[tubeId] - 1;

                for (var nodeId = baseNodeId; nodeId < endNodeId; ++nodeId)
                {
                    var mapping = idToMemMap[nodeId];

                    var content = tubesBufferContent[mapping.BufferId];
                    var data = content.Content;


                    var off1 = mapping.Offset * 4;
                    var off2 = off1 + nodeElementCount * 4;

                    var pos1 = new Vector3(data[off1 + 0], data[off1 + 1], data[off1 + 2]);
                    var tan1 = new Vector3(data[off1 + 4], data[off1 + 5], data[off1 + 6]);

                    var pos2 = new Vector3(data[off2 + 0], data[off2 + 1], data[off2 + 2]);
                    var tan2 = new Vector3(data[off2 + 4], data[off2 + 5], data[off2 + 6]);

                    //center += (pos1 + pos2) * 0.5f;
                    center += (pos1 + pos2 + (pos1 + tan1 / 3.0f) + (pos2 - tan2 / 3.0f)) * 0.25f;

                    ++totalSegmentCount;
                }
            }


            center /= totalSegmentCount;

            var centerOff = newCenter - center;


            for (var nodeId = 0; nodeId < totalNodeCount; ++nodeId)
            {
                var mapping = idToMemMap[nodeId];

                var content = tubesBufferContent[mapping.BufferId];
                var data = content.Content;


                var off = mapping.Offset * 4;

                data[off + 0] += centerOff.X;
                data[off + 1] += centerOff.Y;
                data[off + 2] += centerOff.Z;

                if (mapping.Dublicated)
                {
                    content = tubesBufferContent[mapping.BufferId - 1];
                    data = content.Content;

                    off = (idToMemMap[nodeId - 1].Offset + nodeElementCount) * 4;

                    data[off + 0] += centerOff.X;
                    data[off + 1] += centerOff.Y;
                    data[off + 2] += centerOff.Z;
                }
            }
        }

        public void SetCenterToOrigin()
        {
            SetCenter(Vector3.Zero);
        }

        public void BindShaderAndShaderBuffer(Shader shader, int bufferIndex, string uniformBlockName)
        {
            shaderBuffer.BindToIndex(bufferIndex);

            shader.BindShaderBuffer(shaderBuffer, uniformBlockName);

            shader.Bind();
        }
    }

    public enum TessTubeShadingStyle
    {
        Wip = 0,
        Color = 1,
        SmoothNormal = 2,
        HardNormal = 3,
        HardNormalWF = 4,
        BlueWF = 5,
        LitColor = 6,
    }

    public static class TessTubeShadingStyleExtensions
    {
        public static string ToDefineString(this TessTubeShadingStyle value)
        {
            return "#define SHADING_STYLE " + ((int)value).ToString() + Environment.NewLine;
        }

        public static string ToProfilingString(this TessTubeShadingStyle value)
        {
            return " SHADING_STYLE " + ((int)value).ToString();
        }
    }

    [Flags]
    public enum TessTubeRenderSetting
    {
        Null = 0,
        UseCaps = 1,
        ShaderBufferTypeIsSSBO = 2,
        UseGeoShader = 4,
        UseBackPatchCulling = 8,
        UsePerVertexNormals = 16,
        UseFineCapsTessellation = 32,
        InvertColorsForCaps = 64,
    }

    public static class TessTubeRenderSettingExtensions
    {
        public static string ToDefineString(this TessTubeRenderSetting value)
        {
            var str = "";

            if (value.HasFlag(TessTubeRenderSetting.UseCaps))
                str = "#define CAPS" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.ShaderBufferTypeIsSSBO))
                str += "#define TUBES_DATA_TYPE_IS_SSBO" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.UseGeoShader))
                str += "#define USE_GEOSHADER" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.UseBackPatchCulling))
                str += "#define USE_BACKPATCHCULLING" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.UsePerVertexNormals))
                str += "#define USE_PER_VERTEX_NORMALS" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.UseFineCapsTessellation))
                str += "#define USE_FINE_CAPS_TESSELLATION" + Environment.NewLine;

            if (value.HasFlag(TessTubeRenderSetting.InvertColorsForCaps))
                str += "#define INVERT_COLOR_FOR_CAPS" + Environment.NewLine;

            return str;
        }

        public static string ToProfilingString(this TessTubeRenderSetting value)
        {
            var str = "";

            if (value.HasFlag(TessTubeRenderSetting.UseCaps))
                str = " CAPS";

            if (value.HasFlag(TessTubeRenderSetting.ShaderBufferTypeIsSSBO))
                str += " TUBES_DATA_TYPE_IS_SSBO";

            if (value.HasFlag(TessTubeRenderSetting.UseGeoShader))
                str += " USE_GEOSHADER";

            if (value.HasFlag(TessTubeRenderSetting.UseBackPatchCulling))
                str += " USE_BACKPATCHCULLING";

            if (value.HasFlag(TessTubeRenderSetting.UsePerVertexNormals))
                str += " USE_PER_VERTEX_NORMALS";

            if (value.HasFlag(TessTubeRenderSetting.UseFineCapsTessellation))
                str += " USE_FINE_CAPS_TESSELLATION";

            if (value.HasFlag(TessTubeRenderSetting.InvertColorsForCaps))
                str += " INVERT_COLOR_FOR_CAPS";

            return str;
        }
    }


    public class SplineTubesTess : SplineTubes
    {
        protected Shader shader;

        public Shader Shader { get { return shader; } }

        protected IndexBuffer cylIB;
        protected IndexBuffer cylCapIB;

        protected TessTubeShadingStyle shadingStyle;
        protected TessTubeRenderSetting renderSettings;

        public TessTubeShadingStyle ShadingStyle { get { return shadingStyle; } }
        public TessTubeRenderSetting RenderSettings { get { return renderSettings; } }

        protected int vertexRingVertexCount = 8;// power-of-two
        protected int vertexRingCount = 16;// power-of-two

        protected int maxLayerId;


        public SplineTubesTess(string shaderPath, TessTubeShadingStyle shadingStyle, TessTubeRenderSetting renderSettings, BufferUsageHint bufferUsageHintType, int bufferMaxElementCount, int[] nodeCount, int vertexRingCount, int vertexRingVertexCount)
            : base(renderSettings.HasFlag(TessTubeRenderSetting.ShaderBufferTypeIsSSBO) ? ShaderBufferType.SSBO : ShaderBufferType.UBO, bufferUsageHintType, bufferMaxElementCount, nodeCount, 1 + vertexRingCount / 2)
        {
            this.shadingStyle = shadingStyle;
            this.renderSettings = renderSettings;

            this.vertexRingVertexCount = vertexRingVertexCount;
            this.vertexRingCount = vertexRingCount;

            this.maxLayerId = vertexRingCount - 1;

            if (shadingStyle == TessTubeShadingStyle.HardNormalWF || shadingStyle == TessTubeShadingStyle.BlueWF)
                renderSettings |= TessTubeRenderSetting.UseGeoShader;

            if (renderSettings.HasFlag(TessTubeRenderSetting.UseCaps))
                cylIB = IndexBuffer.CreateClosedQuadCylinder(vertexRingVertexCount, vertexRingCount);
            else
                cylIB = IndexBuffer.CreateQuadCylinder(vertexRingVertexCount, vertexRingCount);

            var definesStr = shadingStyle.ToDefineString() + renderSettings.ToDefineString();
            definesStr += "#define TUBES_BUFFER_SIZE " + this.shaderBuffer.Float4BlockCount.ToString();
            definesStr += Environment.NewLine;
            definesStr += "#define VRING_VERTEX_COUNT " + vertexRingVertexCount.ToString();
            definesStr += Environment.NewLine;
            definesStr += "#define VRING_COUNT " + vertexRingCount.ToString();
            definesStr += Environment.NewLine;

            if (renderSettings.HasFlag(TessTubeRenderSetting.UseGeoShader))
            {
                shader = new Shader(
                "#version 440",
                definesStr,
                shaderPath + ".defs",
                new Shader.Desc(ShaderType.VertexShader, shaderPath + ".vs"),
                new Shader.Desc(ShaderType.TessControlShader, shaderPath + ".tc"),
                new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath + ".te"),
                new Shader.Desc(ShaderType.GeometryShader, shaderPath + ".gs"),
                new Shader.Desc(ShaderType.FragmentShader, shaderPath + ".fs"));
            }
            else
            {
                shader = new Shader(
                "#version 440",
                definesStr,
                shaderPath + ".defs",
                new Shader.Desc(ShaderType.VertexShader, shaderPath + ".vs"),
                new Shader.Desc(ShaderType.TessControlShader, shaderPath + ".tc"),
                new Shader.Desc(ShaderType.TessEvaluationShader, shaderPath + ".te"),
                new Shader.Desc(ShaderType.FragmentShader, shaderPath + ".fs"));
            }
        }

        public void GenTangentFrames()
        {
            for (var tubeId = 0; tubeId < tubeCount; ++tubeId)
                GenTangentFrames(tubeId);
        }

        private Vector3 qvec1 = Vector3.UnitX;
        public void GenTangentFrames(int tubeId)
        {
            Utils.Assert(tubeId >= 0 && tubeId < tubeCount,
                        "tubeId >= 0 && tubeId < tubeCount");

            //var qvec1 = new Vector3(1.0f, 0.0f, 0.0f);
            var lt = new Vector3(0.0f, 1.0f, 0.0f);

            var baseNodeId = cummNodeCount[tubeId];
            var endNodeId = baseNodeId + nodeCount[tubeId] - 1;

            for (var nodeId = baseNodeId; nodeId < endNodeId; ++nodeId)
            {
                var mapping = idToMemMap[nodeId];

                var content = tubesBufferContent[mapping.BufferId];
                var data = content.Content;


                var off1 = mapping.Offset * 4;
                var off2 = off1 + nodeElementCount * 4;

                var pos1 = new Vector3(data[off1 + 0], data[off1 + 1], data[off1 + 2]);
                var tan1 = new Vector3(data[off1 + 4], data[off1 + 5], data[off1 + 6]);

                var pos2 = new Vector3(data[off2 + 0], data[off2 + 1], data[off2 + 2]);
                var tan2 = new Vector3(data[off2 + 4], data[off2 + 5], data[off2 + 6]);


                lt = EvalCSplineTangent(pos1, tan1, pos2, tan2, 0.0f);

                //if (Math.Abs(Vector3.Dot(lt, qvec1)) > 0.000f)
                if (Math.Abs(Vector3.Dot(Vector3.Normalize(lt), qvec1)) > 0.0001f)
                {
                    /* Always works if the input is non-zero.
                    * Doesn’t require the input to be normalised.
                    * Doesn’t normalise the output. 
                    http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts */
                    qvec1 = Math.Abs(lt.X) > Math.Abs(lt.Z) ? new Vector3(-lt.Y, lt.X, 0.0f) : new Vector3(0.0f, -lt.Z, lt.Y);
                    qvec1 = Vector3.Normalize(qvec1);
                }

                var off3 = mapping.SegmentOffset + 1;

                for (int i = 0; i <= maxLayerId; ++i)
                {
                    var tan = EvalCSplineTangent(pos1, tan1, pos2, tan2, (float)i / (float)maxLayerId);

                    if (tan.X == 0.0 && tan.Y == 0.0 && tan.Z == 0.0) tan = lt;


                    var q = Ex.Quaternion.New(lt, tan);

                    qvec1 = Vector3.Transform(qvec1, q);
                    qvec1 = Vector3.Normalize(qvec1);

                    if (float.IsNaN(qvec1.X) || float.IsNaN(qvec1.Y) || float.IsNaN(qvec1.Z))
                    {
                        throw new Exception("NAAAAAAAAAN!!!");
                    }

                    var vec = new Vector2(qvec1.X, (qvec1.Y + 2.0f) * (qvec1.Z < 0.0f ? -1.0f : 1.0f));
                    //var vec = new Vector2(qvec1.X, qvec1.Y);

                    content.Fill(off3 + (i / 2), (i % 2) * 2, vec);
                    //content.Fill(off + 4 + ((int)i / 2), ((int)i % 2) * 2, qvec1.XY());

                    lt = tan;
                }

            }
        }

        public void RecreateShader()
        {
            Shader.Unbind();
            shader.Recreate();
            shader.Bind();
        }

        public void BindShaderAndShaderBuffer()
        {
            BindShaderAndShaderBuffer(shader, 2, "TubesData");
        }

        public void BindMeshData()
        {
            cylIB.Bind();

            GL.PatchParameter(PatchParameterInt.PatchVertices, 4);
        }


        public void UpdateAndDrawAll()
        {
            for (int c = 0; c < shaderBufferCount; ++c)
            {
                shaderBuffer.Update(tubesBufferContent[c]);

                GL.DrawElementsInstanced(PrimitiveType.Patches, cylIB.Count, DrawElementsType.UnsignedInt, IntPtr.Zero, segmentCount[c]);
            }
        }
    }


    public enum RaycastTubeShadingStyle
    {
        Wip = 0,
        Color = 1,
        ColorHQ = 2,
        SmoothNormal = 3,
        BBox = 4,
        LitColorHQ = 5,
    }

    public static class RaycastTubeShadingStyleExtensions
    {
        public static string ToDefineString(this RaycastTubeShadingStyle value)
        {
            return "#define SHADING_STYLE " + ((int)value).ToString() + Environment.NewLine;
        }

        public static string ToProfilingString(this RaycastTubeShadingStyle value)
        {
            return " SHADING_STYLE " + ((int)value).ToString();
        }
    }

    [Flags]
    public enum RaycastTubeRenderSetting
    {
        Null = 0,
        SplitIn4Segments = 1,
        ShaderBufferTypeIsSSBO = 2,
        UseEarlyZ = 4,
    }

    public static class RaycastTubeRenderSettingExtensions
    {
        public static string ToDefineString(this RaycastTubeRenderSetting value)
        {
            var str = "";

            if (value.HasFlag(RaycastTubeRenderSetting.SplitIn4Segments))
                str = "#define SPLIT_IN_4_SEGMENTS" + Environment.NewLine;

            if (value.HasFlag(RaycastTubeRenderSetting.ShaderBufferTypeIsSSBO))
                str += "#define TUBES_DATA_TYPE_IS_SSBO" + Environment.NewLine;

            if (value.HasFlag(RaycastTubeRenderSetting.UseEarlyZ))
                str += "#define USE_EARLYZ" + Environment.NewLine;

            return str;
        }

        public static string ToProfilingString(this RaycastTubeRenderSetting value)
        {
            var str = "";

            if (value.HasFlag(RaycastTubeRenderSetting.SplitIn4Segments))
                str = " SPLIT_IN_4_SEGMENTS";

            if (value.HasFlag(RaycastTubeRenderSetting.ShaderBufferTypeIsSSBO))
                str += " TUBES_DATA_TYPE_IS_SSBO";

            if (value.HasFlag(RaycastTubeRenderSetting.UseEarlyZ))
                str += " USE_EARLYZ";

            return str;
        }
    }


    public class SplineTubesRaycast : SplineTubes
    {
        protected Shader shader;

        public Shader Shader { get { return shader; } }

        protected RaycastTubeShadingStyle shadingStyle;
        protected RaycastTubeRenderSetting renderSettings;

        public RaycastTubeShadingStyle ShadingStyle { get { return shadingStyle; } }
        public RaycastTubeRenderSetting RenderSettings { get { return renderSettings; } }

        protected Mesh<Vector3> cube;
        protected int qSegmentCount;


        public SplineTubesRaycast(string shaderPath, RaycastTubeShadingStyle shadingStyle, RaycastTubeRenderSetting renderSettings, BufferUsageHint bufferUsageHintType, int bufferMaxElementCount, int[] nodeCount)
            : base(renderSettings.HasFlag(RaycastTubeRenderSetting.ShaderBufferTypeIsSSBO) ? ShaderBufferType.SSBO : ShaderBufferType.UBO, bufferUsageHintType, bufferMaxElementCount, nodeCount, 1)
        {
            this.shadingStyle = shadingStyle;
            this.renderSettings = renderSettings;
            
            cube = new Mesh<Vector3>(
            new[]
            {
                new VertexAttribute()
                {
                    Index = 0,
                    Size = 3,
                    Type = VertexAttribPointerType.Float,
                    Normalized = false,
                    Stride = Vector3.SizeInBytes,
                    Offset = 0,
                }
            },
            new[]
            {
                new Vector3(-1f, +1f, +1f),
                new Vector3(-1f, -1f, +1f),
                new Vector3(+1f, +1f, +1f),
                new Vector3(+1f, -1f, +1f),
                new Vector3(+1f, -1f, -1f),
                new Vector3(+1f, +1f, -1f),
                new Vector3(-1f, +1f, -1f),
                new Vector3(-1f, -1f, -1f)
            },
            new uint[]
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
            });

            if (renderSettings.HasFlag(RaycastTubeRenderSetting.SplitIn4Segments))
                qSegmentCount = 4;
            else
                qSegmentCount = 2;

            var definesStr = shadingStyle.ToDefineString() + renderSettings.ToDefineString();
            definesStr += "#define TUBES_BUFFER_SIZE " + this.shaderBuffer.Float4BlockCount.ToString();
            definesStr += Environment.NewLine;

            shader = new Shader(
                    "#version 440",
                    definesStr,
                    shaderPath + ".defs",
                    new Shader.Desc(ShaderType.VertexShader, shaderPath + ".vs"),
                    new Shader.Desc(ShaderType.FragmentShader, shaderPath + ".fs"));
        }

        public void RecreateShader()
        {
            Shader.Unbind();
            shader.Recreate();
            shader.Bind();
        }

        public void BindShaderAndShaderBuffer()
        {
            BindShaderAndShaderBuffer(shader, 2, "TubesData");
        }

        public void BindMeshData()
        {
            cube.Bind();
        }

        public void UpdateShaderBuffer(int bufferIndex = 0)
        {
            Utils.Assert(bufferIndex < shaderBufferCount,
                        "bufferIndex < shaderBufferCount");

            shaderBuffer.Update(tubesBufferContent[bufferIndex]);
        }

        public void DrawSegments(int bufferIndex = 0)
        {
            Utils.Assert(bufferIndex < shaderBufferCount,
                        "bufferIndex < shaderBufferCount");

            cube.DrawInstanced(segmentCount[bufferIndex] * qSegmentCount);
        }

        public void UpdateAndDrawAll()
        {
            for (int c = 0; c < shaderBufferCount; ++c)
            {
                shaderBuffer.Update(tubesBufferContent[c]);

                cube.DrawInstanced(segmentCount[c] * qSegmentCount);
            }
        }
    }
}
