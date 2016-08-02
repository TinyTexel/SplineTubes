using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OpenTK_Test1
{
    public static class TubesData
    {
        public class Point
        {
            public float Px, Py, Pz;

            public float R;

            public float Cr, Cg, Cb;
        }

        public class SplinePoint
        {
            public Point V;
            public Point T;
        }

        public static float Pow2(float x) { return x * x; }

        public static void ConvertBezdatToTubes(string bezdatFilePath, string tubesFilePath, bool generatePlainTextVersion = false)
        {
            System.IO.StreamReader file = null;

            try
            {
                file = new System.IO.StreamReader(bezdatFilePath);
            }
            catch (Exception e)
            {
                Console.WriteLine("ERROR: Could not open file '" + bezdatFilePath.ToString() + "'.");
                Console.WriteLine(e.Message);

                return;
            }

            var firstLine = file.ReadLine();
            if (firstLine == null && firstLine != "BezDatA 1.0")
            {
                Console.WriteLine("ERROR: Expected to read 'BezDatA 1.0' on first line, but read '" + firstLine + "' instead.");
                return;
            }

            var points = new System.Collections.Generic.List<Point>();
            var indices = new System.Collections.Generic.List<int>();

            var maxIndex = 0;
            var lineCounter = 1;
            string line;
            while ((line = file.ReadLine()) != null)
            {
                string[] strs = line.Split(null);

                var len = strs.Length;
                if (len > 0)
                {
                    if (strs[0] == "PT")
                    {
                        if (len != 8)
                        {
                            Console.WriteLine("ERROR: Expected 7 values on line " + lineCounter.ToString() + ", but read " + (len - 1).ToString() + " values instead.");
                            file.Close();
                            return;
                        }

                        try
                        {
                            var point = new Point();

                            point.Px = Convert.ToSingle(strs[1], System.Globalization.CultureInfo.InvariantCulture);
                            point.Py = Convert.ToSingle(strs[2], System.Globalization.CultureInfo.InvariantCulture);
                            point.Pz = Convert.ToSingle(strs[3], System.Globalization.CultureInfo.InvariantCulture);

                            point.R = Convert.ToSingle(strs[4], System.Globalization.CultureInfo.InvariantCulture);

                            point.Cr = (float)Math.Pow(Convert.ToDouble(strs[5], System.Globalization.CultureInfo.InvariantCulture) / 255.0, 2.2);
                            point.Cg = (float)Math.Pow(Convert.ToDouble(strs[6], System.Globalization.CultureInfo.InvariantCulture) / 255.0, 2.2);
                            point.Cb = (float)Math.Pow(Convert.ToDouble(strs[7], System.Globalization.CultureInfo.InvariantCulture) / 255.0, 2.2);

                            points.Add(point);
                        }
                        catch (FormatException e)
                        {
                            Console.WriteLine("ERROR: FormatException on line " + lineCounter.ToString() + ".");
                            Console.WriteLine(e.Message);
                            file.Close();
                            return;
                        }
                        catch (OverflowException e)
                        {
                            Console.WriteLine("ERROR: OverflowException on line " + lineCounter.ToString() + ".");
                            Console.WriteLine(e.Message);
                            file.Close();
                            return;
                        }
                    }
                    else
                        if (strs[0] == "BC")
                        {
                            if (len != 5)
                            {
                                Console.WriteLine("ERROR: Expected 4 values on line " + lineCounter.ToString() + ", but read " + (len - 1).ToString() + " values instead.");
                                file.Close();
                                return;
                            }

                            try
                            {
                                for (var i = 1; i <= 4; ++i)
                                {
                                    var index = Convert.ToInt32(strs[i]);

                                    if (index < 0)
                                    {
                                        Console.WriteLine("ERROR: Negative index on line " + lineCounter.ToString() + ".");
                                        file.Close();
                                        return;
                                    }

                                    indices.Add(index);

                                    if (index > maxIndex) maxIndex = index;
                                }
                            }
                            catch (FormatException e)
                            {
                                Console.WriteLine("ERROR: FormatException on line " + lineCounter.ToString() + ".");
                                Console.WriteLine(e.Message);
                                file.Close();
                                return;
                            }
                            catch (OverflowException e)
                            {
                                Console.WriteLine("ERROR: OverflowException on line " + lineCounter.ToString() + ".");
                                Console.WriteLine(e.Message);
                                file.Close();
                                return;
                            }
                        }
                        else
                        {
                            Console.WriteLine("ERROR: Expected 'PT' or 'BC' as first string on line " + lineCounter.ToString() + ", but read '" + strs[0] + "' instead.");
                            file.Close();
                            return;
                        }
                }

                lineCounter++;
            }

            if (maxIndex > points.Count - 1)
            {
                Console.WriteLine("ERROR: Largest Index (" + maxIndex.ToString() + ") does not reference point.");
                file.Close();
                return;
            }

            var indexLen = indices.Count;
            if (indexLen == 0)
            {
                Console.WriteLine("ERROR: File does not contain nodes.");
                file.Close();
                return;
            }

            Console.WriteLine("SUCCESS: File '" + bezdatFilePath.ToString() + "' read.");
            file.Close();



            var sPoints = new System.Collections.Generic.List<SplinePoint>();

            var nodeCountPerTube = new System.Collections.Generic.List<int>();
            nodeCountPerTube.Add(0);

            var currTubeId = 0;
            var tangentDiscCount = 0;

            for (var i = 0; i <= indexLen; i += 4)
            {
                var i0 = i != 0 ? i - 2 : i + 0;
                var i1 = i != 0 ? i - 1 : i + 1;

                var i2 = i != indexLen ? i + 0 : i0;
                var i3 = i != indexLen ? i + 1 : i1;


                var index0 = indices[i0];
                var index1 = indices[i1];
                var index2 = indices[i2];
                var index3 = indices[i3];


                var p0 = points[index0];
                var p1 = points[index1];

                var p2 = points[index2];
                var p3 = points[index3];


                var tan0 = new Point();
                tan0.Px = (p1.Px - p0.Px) * 3.0f;
                tan0.Py = (p1.Py - p0.Py) * 3.0f;
                tan0.Pz = (p1.Pz - p0.Pz) * 3.0f;

                tan0.R = (p1.R - p0.R) * 3.0f;

                tan0.Cr = (p1.Cr - p0.Cr) * 3.0f;
                tan0.Cg = (p1.Cg - p0.Cg) * 3.0f;
                tan0.Cb = (p1.Cb - p0.Cb) * 3.0f;

                var tan1 = new Point();
                tan1.Px = (p3.Px - p2.Px) * 3.0f;
                tan1.Py = (p3.Py - p2.Py) * 3.0f;
                tan1.Pz = (p3.Pz - p2.Pz) * 3.0f;

                tan1.R = (p3.R - p2.R) * 3.0f;

                tan1.Cr = (p3.Cr - p2.Cr) * 3.0f;
                tan1.Cg = (p3.Cg - p2.Cg) * 3.0f;
                tan1.Cb = (p3.Cb - p2.Cb) * 3.0f;


                var pTanDiff = Pow2(tan1.Px - tan0.Px) + Pow2(tan1.Py - tan0.Py) + Pow2(tan1.Pz - tan0.Pz);
                var rTanDiff = Pow2(tan1.R - tan0.R);
                //var cTanDiff = Pow2(tan1.Cr - tan0.Cr) + Pow2(tan1.Cg - tan0.Cg) + Pow2(tan1.Cb - tan0.Cb);

                var tanDiffDelta = 0.001;

                var deltaCond = pTanDiff < tanDiffDelta && rTanDiff < tanDiffDelta;// && cTanDiff < tanDiffDelta;// do not check color delta due to high quantization error? or errors in the data...*shrug*
                var sharedNodeCond = index1 == index2;
                var overwrittenNodesCond = index0 == index2 && index1 == index3;

                if (sharedNodeCond && !deltaCond)
                    ++tangentDiscCount;

                if ((sharedNodeCond || overwrittenNodesCond) && deltaCond)
                {
                    var sPoint = new SplinePoint();

                    sPoint.V = i != indexLen ? p2 : p3;

                    sPoint.T = new Point();
                    sPoint.T.Px = (tan0.Px + tan1.Px) * 0.5f;
                    sPoint.T.Py = (tan0.Py + tan1.Py) * 0.5f;
                    sPoint.T.Pz = (tan0.Pz + tan1.Pz) * 0.5f;

                    sPoint.T.R = (tan0.R + tan1.R) * 0.5f;

                    sPoint.T.Cr = (tan0.Cr + tan1.Cr) * 0.5f;
                    sPoint.T.Cg = (tan0.Cg + tan1.Cg) * 0.5f;
                    sPoint.T.Cb = (tan0.Cb + tan1.Cb) * 0.5f;

                    sPoints.Add(sPoint);

                    ++nodeCountPerTube[currTubeId];
                }
                else
                {
                    var sPoint0 = new SplinePoint();

                    sPoint0.V = p1;
                    sPoint0.T = tan0;

                    sPoints.Add(sPoint0);

                    ++nodeCountPerTube[currTubeId];


                    var sPoint1 = new SplinePoint();

                    sPoint1.V = p2;
                    sPoint1.T = tan1;

                    sPoints.Add(sPoint1);

                    nodeCountPerTube.Add(1);
                    ++currTubeId;
                }
            }

            Console.WriteLine("Number of Tangent Discontinuities: " + tangentDiscCount.ToString());


            // create binary file
            {
                System.IO.BinaryWriter newFile = null;

                try
                {
                    newFile = new System.IO.BinaryWriter(System.IO.File.Open(tubesFilePath, System.IO.FileMode.Create));
                }
                catch (Exception e)
                {
                    Console.WriteLine("ERROR: Could not create/open file '" + tubesFilePath.ToString() + "'.");
                    Console.WriteLine(e.Message);

                    return;
                }

                try
                {
                    var formatVersion = 1;
                    newFile.Write(formatVersion);// 1 x Int32

                    var tubeCount = nodeCountPerTube.Count;
                    newFile.Write(tubeCount);// 1 x Int32

                    foreach (int nodeCount in nodeCountPerTube)// tubeCount x Int32
                        newFile.Write(nodeCount);

                    var totalNodeCount = sPoints.Count;
                    newFile.Write(totalNodeCount);// 1 x Int32

                    foreach (SplinePoint node in sPoints)// totalNodeCount x (3 + 1 + 3) * 2 x Single
                    {
                        newFile.Write(node.V.Px);
                        newFile.Write(node.V.Py);
                        newFile.Write(node.V.Pz);

                        newFile.Write(node.V.R);

                        newFile.Write(node.V.Cr);
                        newFile.Write(node.V.Cg);
                        newFile.Write(node.V.Cb);


                        newFile.Write(node.T.Px);
                        newFile.Write(node.T.Py);
                        newFile.Write(node.T.Pz);

                        newFile.Write(node.T.R);

                        newFile.Write(node.T.Cr);
                        newFile.Write(node.T.Cg);
                        newFile.Write(node.T.Cb);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine("ERROR: Failed to write value to file '" + tubesFilePath.ToString() + "'.");
                    Console.WriteLine(e.Message);

                    newFile.Close();
                    return;
                }

                Console.WriteLine("SUCCESS: File '" + tubesFilePath.ToString() + "' generated.");
                newFile.Close();
            }


            // create plain text file for debugging purposes
            if(generatePlainTextVersion)
            {
                var tubesTxtFilePath = tubesFilePath + ".txt";

                System.IO.StreamWriter newFile = null;

                try
                {
                    newFile = new System.IO.StreamWriter(System.IO.File.Open(tubesTxtFilePath, System.IO.FileMode.Create));
                }
                catch (Exception e)
                {
                    Console.WriteLine("ERROR: Could not create/open file '" + tubesTxtFilePath.ToString() + "'.");
                    Console.WriteLine(e.Message);

                    return;
                }

                try
                {
                    var formatVersion = 1;
                    newFile.Write("Version: ");
                    newFile.Write(formatVersion);// 1 x Int32
                    newFile.Write(Environment.NewLine);

                    var tubeCount = nodeCountPerTube.Count;
                    newFile.Write("Tube Count: ");
                    newFile.Write(tubeCount);// 1 x Int32
                    newFile.Write(Environment.NewLine);

                    newFile.Write("Node Counts per Tube: ");
                    foreach (int nodeCount in nodeCountPerTube)// tubeCount x Int32
                    { newFile.Write(nodeCount); newFile.Write(" "); }

                    newFile.Write(Environment.NewLine);

                    var totalNodeCount = sPoints.Count;
                    newFile.Write("Total Node Count: ");
                    newFile.Write(totalNodeCount);// 1 x Int32

                    var n = 0;
                    newFile.Write(Environment.NewLine);
                    foreach (SplinePoint node in sPoints)// totalNodeCount x (3 + 1 + 3) * 2 x Single
                    {
                        newFile.Write(Environment.NewLine);

                        newFile.Write("Node " + n.ToString() + ":");
                        newFile.Write(Environment.NewLine);
                        ++n;

                        newFile.Write(node.V.Px.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.V.Py.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.V.Pz.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);

                        newFile.Write(node.V.R.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);

                        newFile.Write(node.V.Cr.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.V.Cg.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.V.Cb.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);

                        newFile.Write(node.T.Px.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.T.Py.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.T.Pz.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);

                        newFile.Write(node.T.R.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);

                        newFile.Write(node.T.Cr.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.T.Cg.ToString(System.Globalization.CultureInfo.InvariantCulture)); newFile.Write(" ");
                        newFile.Write(node.T.Cb.ToString(System.Globalization.CultureInfo.InvariantCulture));
                        newFile.Write(Environment.NewLine);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine("ERROR: Failed to write value to file '" + tubesTxtFilePath.ToString() + "'.");
                    Console.WriteLine(e.Message);

                    newFile.Close();
                    return;
                }

                Console.WriteLine("SUCCESS: File '" + tubesTxtFilePath.ToString() + "' generated.");
                newFile.Close();
            }
        }
    }
}
