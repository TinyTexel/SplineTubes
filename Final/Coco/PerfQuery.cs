using System;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using System.IO;

namespace Coco
{
    public class PerfQuery
    {
        protected int startQuery, endQuery;
        protected long timerStart, timerEnd;
        protected int done;

        public PerfQuery()
        {
            startQuery = GL.GenQuery();
            endQuery = GL.GenQuery();

            //var queries = new int[2];
            //GL.GenQueries(1, out startQuery);
            //GL.GenQueries(1, out endQuery);

            //startQuery = 
        }

        public void Start()
        {
            GL.QueryCounter(startQuery, QueryCounterTarget.Timestamp);
        }

        public void Stop()
        {
            GL.QueryCounter(endQuery, QueryCounterTarget.Timestamp);
        }

        public double RetrieveResult()
        {
            done = 0;

            while (done == 0)// proper low budget port approach (ring buffers... pfff)
                GL.GetQueryObject(endQuery, GetQueryObjectParam.QueryResultAvailable, out done);

            GL.GetQueryObject(startQuery, GetQueryObjectParam.QueryResult, out timerStart);
            GL.GetQueryObject(endQuery, GetQueryObjectParam.QueryResult, out timerEnd);

            return (timerEnd - timerStart) / 1000000.0;
        }

        public string RetrieveResultAsString(string format = "N2")
        {
            return RetrieveResult().ToString(format);
        }

        public void Delete()
        {
            GL.DeleteQuery(startQuery);
            GL.DeleteQuery(endQuery);
        }
    }
}
