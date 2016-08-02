Consequently, methods that allow to efficiently render a large number of curved tubes with varying radii are needed.

-----

----- 

The only cross-section frame described that can be computed analytically is the frenet frame \cite{frenet}.
 without pre-processing step
In addition, the used Phong lighting model is adaptively supersampled.
Following this approach instead of generating the geometry completely on the CPU, only a minimal amount of data needs to be transferred to the GPU keeping the bandwidth demands comparatively low for complex scenes. The same strategy is applied by the tessellation based approach in this work.
which prevents them from partially disappearing and re-appearing when the camera moves.

Spline tubes form a handy subset of generalized cylinders \cite{GenCyl}. A generalized cylinders is defined by a space curve forming the axis and an arbitrary cross-section. In the case of a spline tube, the axis follows a spline curve and the cross-section is limited to be circular.
They are defined by bent axes following space curves combined with arbitrary cross-sections.

As already mentioned, the algorithms later presented in this work are specialized in rendering spline tubes, which form a handy subset of generalized cylinders. 

Several works concerned themselves with decomposing scenes into generalized cylinders.

Rendering Generalized Cylinders using the A-Buffer
		~ dynamic adaptive screen-space tessellation of curves with varying radii
		~ using an A-Buffer for AA purposes (also adaptively supersamples its Phong lighting model)
		~ software-based implementation
		~ positions und radii are interpolated using cubic splines (Catmull-Rom)
		~ other properties like color, opacity, reflectance are linearly interpolated 
		~ coarse tessellation followed by length- and breadthwise refinement using different constraints
		~ position constraint I: based on the angle between the screenspace tangent vectors of start/end-point of sub-segment
		~ position constraint II: maintains a desired degree of linearity in the z-component
		  ~> necessary to ensure that a curved segment is adequately subdivided even when viewed from an angle 
			 that makes its screen pro jection close to linear (important for quality of normales and texcoords)
		~ radius contraint: ensure smooth variation in the radius of a segment
		~ inflection point constraint: ensures that inflection point lie between sub-segment (on their vertices)
		~ 3 quality levels for breadthwise refinement, higher levels better approximate tube
		  ~> chosen depending on screenspace size and quality requirements
		~ backface culling for each segment
		~ normal calculation must consider changing radius
		~ high screenspace curvature can lead to ill-formed polygons
		  ~> handled by splitting the offending polygon into an equivalent pair of triangles
		~ variety of effects like length- and breadthwise opacity variation and a global illumination model that approx AO based 
		on the assumption that the paintstroke (tube) is embedded in a layer of a homogeneous medium
		~ benchmarks are pretty deprecated (thesis from 1997), render times given in ticks (s/60, still better than fps...)
		~ overall approach appears annoyingly fiddly



are resolved by splitting the resulting ill-formed polygons into pairs of equival 
		

High Performance generalized cylinders visualization
-adaptive curvature sampling
- no rendetime refindement
- arbitrary corss-section
- skinnning for deformation
- frame time guarantee through dynamic adjustment of used segment number and complexity


The key feature of the presented algorithm is the

However, no dynamic refinement with regard to the is performed at render time.
While this method requires only two tangent frames at the beginning and end of each curve segment, it also limits 

Rendering Tubes from Discrete Curves Using Hardware Tessellation
		~ sequence of points define curve -> interpolated with Catmull-Rom splines
		~ points are sent to GPU, geometry created in hull-shader/tessellator/domain-shader
		~ approx the curve with linear segments, adaptively adds more segments in regions with large derivatives
		~ LOD reduces only cross-section vertex count down to three (triangular cross-section)
		~ radius is scaled up over distance to reduce aliasing 
			~> apparently no MSAA is used in addition but AA is supposedly eliminated (?)
			~> also no opacity reduction applied when reaching sub-pixel diameters
		~ benchmark compares the tessellator-approach to a naive one where the geo is generated on CPU and send to GPU
			~> strongly reduced banwidth requirements
		~ AA approach (scaling radius) is compared to no AA
			~> nearly no fps drop
			
			
		Hyper stream lines stuff
		
--------

The visualization of tube-based structures finds widespread application in various fields.
The visualization of tube-based structures has numerous applications in various fields.
The visualization of tube-based structures finds application in various fields like math, physics, and medicine. In this paper, two approaches to the rendering of tubes, also known as generalized cylinders, are explored. The first approach is tessellation based and refines the tube geometry adaptively at render time using the tessellation unit of modern graphics pipelines. The 

This paper describes the rendering of tubes, also know as generalized cylinders, using a tessellation as well as a raycasting approach. Both approaches are profiled and compared with regard to frame time and memory consumption. 


Advantages and shortcomings of both methods with regard to various parameters are discussed and several potential improvements are proposed.

To reach interactive frame rates when visualizing large data sets, efficient rendering strategies for tubes with varying radii are needed.
To reach interactive frame rates for visualizations of large data sets, strategies that allow to efficiently render tubes with varying radii are needed.
To reach interactive frame rates when visualizing large data sets, strategies that allow to efficiently render tubes with varying radii are needed.
Efficient rendering strategies for tubes with varying radii are important to reach interactive frame rates.  
This paper describes and compares two different rendering approaches for tubes with varying radius. , also know as generalized cylinders.
The first is tessellation based and relies on an adaptive refinement at render time to ensure a smooth appearance for diff viewing conditions, the second is uses raycasting and naturally ensures a smooth depiction of the geometry but suffers shortcomings .

approach, method, strategy, algorithm, technique

- trajectories of any kind

-paricle trajectories /real/virtual
- stream lines	
- white matter tracts
- nerve fibers in general
- vectors fields using many short tubes

- connected tubes for
- graph like structures
- f.i. blood vessels


-------


A spline tube is in the following defined to be a generalized cylinder with circular cross-section of varying radius of which the axis follows a spline curve.
that is referred to as spline tube hereafter.
tubes with a specific parametrization hereafter referred to as spline tubes.

Two such approaches specialized in rendering spline tubes, are presented in this work.

 and are therefore referred to as spline tubes hereafter.
This work presents two such approaches . They are specialized in rendering spline tubes.
This work presents two such approaches specialized in rendering tubes or more specific spline tube.
Two such rendering approaches are presented in this work.

Depending on the complexity, visualizing such structures at interactive frame rates using tubes can require algorithms that are capable to efficiently render those in large numbers. 

Due to the potentially high complexity of these structures, visualizing them at interactive frame rates requires algorithms that are capable to efficiently render tube primitives in large numbers. 

They, therefore, find application in scientific visualizations across various disciplines like math, physics, chemistry, biology, and medicine.MEH
large numbers of tube primitives efficiently. 
However, due to the high complexity these structures can exhibit, visualizing them at interactive frame rates requires efficient algorithms specialized in rendering large numbers of tube primitives.
However, due to the potentially high complexity, visualizing these structures at interactive frame rates requires efficient algorithms specialized in rendering large numbers of tube primitives. 

A generalized cylinder is defined by a curved axis and an arbitrarily formed cross-section. 
In this work only circular cross-sections of varying radii and axes following spline curves are considered. 
This sub-set of generalized cylinder is therefore referred to as spline tubes hereafter.

Being as widely applicable as they are, algorithms that render tubes as efficiently as possible
Therefore, efficient rendering of tube-based structures is of interest to visualization applications for various scientific disciplines.
In the context of computer graphics tubes-like geometry is commonly modeled by generalized cylinders.

A generalize cyl is defined by an axis following an arbitrary crve
by a curved axis and an arbitrarily formed cross-section at any point along a curve axis.

advantageous use cases
Tube-like structures occur naturally in a vairety 
Tubes provide a handy primitive to model a variety of natural structures as well as o visualize scientific data.

Their different
Both approaches are profiled and compared with regard to frame time and memory consumption. 
Generalized cylinders are a well known primitve in comp

Applications range from visus of ...  and ... to ... and ... .
Ranging from applications in biology and medicine where nerve fibers like white matter tracts and blood vessels can be represented to the visu of particle trajectories and abstract concepts like vector/tensor fields.

- what is tube, gen cyl
- what applications, exmaples
- 
Often structures can be visualized with \textit{generalized cylinders} . Generalized Cylinders are cylinders that are generalized.

In the following, first blabla is shortly considered before blablubs is described in the subsequent chapter. 

------

------

First, concepts and implementation aspects both approaches share are outlined.
This subsection outlines concepts and implementation aspect both approaches share. //subsequently described share.

-------

Starting from the tessellation of a straight tube, the vertex shader first transforms this cylinder to a coarse representation of

For both approaches the segment data at least consists of a memory offset to the data of the two nodes defining the segment (see figure \ref{foo}).
The base tessellation provides a cylinder mesh that is modeled as a stack of vertex rings forming a quad-based topology (see figure).
The world space position of the vertex is then computed by a translation from the curve position along an offset vector that is orthogonal to the tangent and of which the length equals the radius (see figure).
For the purpose of transforming this vector to the offset vector, two vectors that together with the curve tangent form an orthonormal basis need to be determined (see figure).
However, none of these strategies is or can be free of discontinuties for continuously varying input vectors \cite{hairyball} and therefore lead to bad topology in some cases due to heavily twisted consecutive vertex rings (see figure).
This vector is then repeatedly projected into the \textit{tangent plane} of the next vertex ring and re-normalized (see figure).
Instead, a surface tangent pointing in curve direction is computed using a central difference (see figure). The expansion vector is then projected onto the plane of this tangent. Normalizing the resulting vector gives the correct normal (see figure). 

figures: memory layout, local tangent frames

model: discs with varying radii along curve
base tessellation, cylinder, vertex rings
quad based, better perf, better tess topology after refinement
local tangent frame required to expand vertex rings from central axis
tangent is used as one axis
several closed form approaches exist but all suffer from degenerated cases/discontinuities
therefore iterative technique where a given orthogonal vector is carry along the curve from one end to the other
for each vertex ring one ortho vec is computed in a pre-processing step and passed to the GPU as per-segment data
the third vec is computed by cross product between tangent and ovec
ovec is compressed to 2 components, possible cause normalized
(instanced) indexed draw call without vertex buffer
vertex positions can be calculated from VertexID since vertex count of the whole tube-segment as well as of each vertex ring is known
vertex shader: 
read per-segment data from buffer: node id ... 
evaluate position, tangent and radius curves
...
use base tessellation to determine refinement parameters for on-the-fly tessellation
tess-control:
two criteria for the amount of tessellation

problems - ugly bulbis!
tangent plane
cross-section stuff

Strategies that specifically attempt to derive a vector orthogonal to the tangent of a cubic spline or Bèzier curve exist \cite{overviewpaper}, but ultimately they too suffer from potential discontinuties.

for continuously varying input vectors \cite{hairyball}
Rendering caps to close up the endings of the tubes requires some additional considerations.
Due to their sharp edges, cap seams need to be finer tessellated to provide a smooth appearance.
Since caps provide distinct silhouettes their seam edges are finer tessellated. 

Their interior geometry, on the other hand, being flat in nature is not further refined. 

, while their interior parts, being flat in nature, are not further refined.

Being flat in nature, the interior part is not further refined.




End caps are explicitly culled for the case that a single tube is constructed from multiple tubes.
In the case that a single tube is constructed from multiple tubes, the end caps at the connection points



To allow for discontinuities along tubes it is possible to construct a single tube from multiple ones. Caps in




The outer and inner tessellation levels of the quad patches are computed based on two criteria.
The tessellation evaluation shader now transforms the refined cylinder geometry again to a tube segment.
This ensures that the tessellation factors of patches close to the silhouette are independent of the view space orientation of their edges.
Each of the two inner tessellation levels for a given quad patch is set to the maximum of the two corresponding outer levels (unclear).
The x and y coordinates are scaled by the radius while the z coordinate is set to zero. This position is then rotated into a origin positioned plane of which the normal is given by the curve tangent, hereafter referred to as tangent plane. After applying a translation by the curve position the resulting coordinates are in world space. The rotation into the tangent plane requires two vectors that together with the normalized tangent form an orthonormal matrix. 

As long as the angle between any two successive tangents is different or not too close to 90 degree, this strategy produces a new non-zero vector orthogonal to the current tangent after each step.
After tessellation refinement by the tessellation control shader, it is the task of the tessellation evaluation shader to perform the transformation from refined 
to again transform the geometry into a tube segment.


 based on the interpolated curve parameter
Eventually the fragment shader outputs a basic shading  

At this point it is the task of the tessellation evaluation shader to again perform the transformation from cylinder to tube segment but this time using the adaptively refined tessellation. There are various approaches to this problem. (meh meh meh)

An alternative, more straightforward, approach w
Alternatively, it is also possible to 

Each ortho vec is therefore calculated by first projecting the closest two pre-computed ortho vecs onto the tangent plane of the considered vertex. The new vector is then determined by linear interpolating these two vectors so that 

the new geometry is again transformed into a tube segment.

One possibility is  to push the vertex out of the patch plane along the vector pointing from the position on the curve determined by evaluating the position spline using the interpolated curve parameter to the interpolated world space position of the vertex.

to calculate an interpolated world space position and curve parameter based on the respective attributes of the four control points as well as the bilinear weights provided to the shader. Directly using the new world space position would not lead to a better approximation of the segment since all of the new vertices of a patch would be coplanar.
By assuming that the vertex
Using the curve parameter a position on the curve is determined.

Using an approach derived from the one already used in the vertex shader prevents this issue. However, for this to work it is necessary to determine orthogonal vectors for the new vertices, which might lie in between two vertex rings and therefore can not simply use one of the pre-computed vectors stored in the tubes buffer.

The new vector is therefore calculated by projecting the two ortho vecs of the vertex rings adjacent to the patch onto the tangent plane of the considered vertex.
These vectors are then linearly interpolated to ensure a continuous transition from one orthogonal vector to the other depending on the distance of the vertex to the vertex rings in parameter space. Finally the resulting vector is normalized and used to construct the second orthogonal vector as described before. To reliably work this approach requires the angle between the tangents of the patch adjacent vertex rings to be smaller than 90°. For the tested data sets this requirement never turned out to be a problem.


This can be achieved by projecting the two ortho vecs of the top and bottom vertex rings into the tangent plane of the considered vertex.
The only difference to t

use the interpolated world space position and curve parameter of each vertex

transform the refined tessellation to a tube segment.

However, as long as only the radius spline is $C1$ discontinuous the resulting artifacts are seldom observed.

----------

Instead, ray intersections with infinitely many scaled geometric primitives positioned along the curves serve as mathematical foundation.
In the following, disks and spheres are considered as such sweep primitives (see figure).

*In the following, disks and spheres are considered as such primitives (see figure). x2
*binary search \cite{bisectionmethod} (see figure)
*iterative approach next less derived version f (see figure).
*Since f' is the derivative of f, of which the minima coincide with the front-side ray-tube intersections, there are up to three real roots (see figure) 
*of which the outer two are of interest. Therefore, another 5 binary searches per interval are necessary to find the relevant roots of f' (see figure).
*precision problem for looong segments (figure)


Finally, based on world position, normal, and color, the shading is computed just like for the tessellation approach. Some care with regard to the color spline is required. When using only two quadratic segments for the approximate of the cubic one, the resulting color can differ substantially from the original one. 

using the assumption that the quadratic sub-segments split the original cubic one equally.
In opposite to the disk-based tube model,
To find the roots of f that relate to actual intersections, the intervals on which f is real need to be determined first.
Differentiation of f with regard to the curve parameter allows results in equation of which the roo
The equation f returns the ray parameter resulting in a front side sphere intersection for a given curve parameter. Differentiation with regard to the curve parameter results in equation of which the roots correspond to tube segment intersections. To reduce the complexity of the formula, the tube segment is first transformed into ray space where the ray itself equals the x-axis.

Depending on the primitive, the complexity of the resulting equations differs as well as the geometrical properties of the rendered tube segment.
Not only does the complexity of the resulting equations depend on the chosen primitive  but also the geometrical characteristics/properties/aesthetics of the rendered tube segment.
characteristics/properties/aesthetics
Two possible primitives are considered in the following: disks and spheres (see figrue).
Instead, ray intersections with geometric primitives, of which the scale and position are parametrized by the radius and position splines, serve as starting points.
As an optimization, the algorithm determines the two roots of polyB'' analytically, effectively trading another two binary searches for a square root.
Due to there being two potential front side intersection, F

, delimited by the roots of polyB.
which delimit the real intervals of f, the roots f 


Instead an iterative, numerical strategy is applied. 
Since these in turn correspond to the intervals on which polyB is positive, the roots of polyB have to be computed.
Computing the roots of polyB provides the intervals on which f is real. 
Since only real roots are relevant, the first step is to determine the intervals on which f is real which correspond to the intervals on which polyB is positive.
for all relevant interavals. 
When trying to determine the intervals for which...it can be observed that all roots lie in between two consecutive extrema of f. 
By observing that the roots of f
The intervals for which binary searches are applied, are 

making the latter more attractive .


those can converge slower than the binary search in practice or even
 
the actual convergence behavior of the binary search guarantees an improvement after each iteration is more robust and 

produces results of reliable precision and is cheap.

Instead the following strategy is applied: first the roots of the second derivative of polyB are computed analytically. 
Using those and the knowledge


Alternatively iterative strategies 
and a xm at which f is negative as well as a xp where f evaluates to a positive value, the 
The segment is then intersected with the x-axis.
 
while not having to hassle with precision problems that can occur for the disk 
 
In both cases finding the roots 

 and differentiating it with respect to the ray parameter
Modeling the tube segments by disks
Depending on the primitive used the resulting geometry differs.

two models: discs and spheres (figures)
either way, no closed solution, so we need numerical root finder
discs are what we did with the tessellation approach, same problems (blulbs)
also precision probs when ray ortho to tangent, because double roots of poly

spheres:
use closed solution for ray-sphere intersection
parametrize position and radius using the respect splines of the segment
transform into ray-space, ray equals x-axis, reduces complexity of resulting terms
 polyA - sqrt(polyB) (figure)
we are interested in the real minima
through differentiation we get
 polyA' - polyB' / (sqrt(polyB) * 2)
we need real roots
idea: first find roots of polyB since those define ranges where f' is real
then search roots of f' in these ranges.
how to find the roots of polyB?
observation roots lie between extrema so knowing two subsequent extrema with one extremum being positive and on negative the root in between can be numerically determined by a binary search, although strategies with better theoretical convergence like Newton-Raphson, which also considers derivatives, exist, the actual convergence of the binary search is more robust and produces results of reliable precision and is cheap.
The extrema themselves correspond to roots of the derivative, so an iterative strategy can be defined bla. This strategy can applied to arbitrary continuous functions.
However the number of necessary binary search grows quadratically with the degree of the polynomial.
Knowing the roots of polyB the same root finding strategy can be applied to the real parts of f'. Depending on the number of potential real roots of f', f' needs to be differentiated a certain number of times until only one root can exist.

For a cubic spline there are up to 3 front side intersections and therefore 3 roots in a single real segment of f'. Therefore f' needs to be differentiated 3 times. 
Also 6 roots for polyB -> therefore many binary searches
Overall too expensive for cubic spline

Idea: approximate cubic spline by two quadratic splines (figure). Since polyB in this case is only of degree 4 it is theoretically to evaluate the roots analytically. However, using the closed form solution for equations of degree 4 is potentially numerically unreliable. Since evaluating a polynomial of degree n can be done with n fused multiply-adds and is therefore fast, using the iterative scheme for poly degree 4 reasonable fast.
Also only 2 potential intersections and roots in each real segment of f' need to be found, meaning f' needs to be differentiated twice to apply the iterative root search.
interval

through root preserving conversion:
 polyA' * sqrt(polyB) * 2 - polyB'
 polyA'^2 * polyB * 4 - polyB'^2

 polyA'^2 - (polyB'^2 / polyB / 4)

for cubic position and radius splines: 
polyA 						degree 3
polyB 						degree 6
polyA' 						degree 2
polyB' 						degree 5
polyA'^2 					degree 4
polyB'^2 					degree 10
polyA'^2 * polyB 	degree 10..

for quadratic position and radius splines: 
polyA 						degree 2
polyB 						degree 4
polyA' 						degree 1
polyB' 						degree 3
polyA'^2 					degree 2
polyB'^2 					degree 6
polyA'^2 * polyB 	degree 6

For each quadratic sub-segment, a cube is rendered using one or multiple instanced draw calls. The number of draw calls depends on how often the tubes buffer need to be refilled, which in turn depends on the maximal allowed size of the buffer type and the amount of available GPU memory. In contrast to the tessellation approach, no orthogonal vectors need to be computed and stored in the tubes buffer. This can result in fewer draw calls. The vertex shader then computes a bounding box for each quadratic tube segment and transforms the cube geometry accordingly (see figure \ref{fig:bb_comp}).


Figure \ref{fig:perf4} compares the frame times of the rendering of the Fibers2 data set to variants with skipped rasterization and no geometry processing at all. To prevent geometry from being rasterized, it is  moved behind the clip plane by negating all z-value of the projected vertex positions. The no-rendering variant additionally only uses dummy vertex- and fragment shaders. 

The camera-to-scene distances are chosen per data so as to achieve a medium and a high amount of screen coverage by the geometry for a far away and a closely positioned camera, respectively.
The camera-to-scene distances are chosen per data set so as to achieve different amounts of screen coverage by the geometry.

- stats fibers:
	- tubes count: 241
	- segment count 3530
	- node count 3771
	- tubes buffer size TE: 0.338 MB
	- tubes buffer size RC: 0.215 MB
	- 0,635
	
- stats fibers2:
	- tubes count: 1701
	- segment count 114489
	- node count 116190
	- tubes buffer size TE: 10.628 MB
	- tubes buffer size RC: 6.621 MB
	- 0,623
	
- stats vectubes:
	- tubes count: 722
	- segment count 30948
	- node count 31670
	- tubes buffer size TE: 2.888 MB
	- tubes buffer size RC: 1.805 MB
	- 0,625
	
	- profile pure bb for rc
- and analytical tangent frames/ no trans	

		\item \emph{far}/\emph{close} (blue/orange) corresponds to the camera-to-scene distances (medium/high screen coverage by the geometry)
		\item \emph{small}/\emph{medium} (light/dark) corresponds to the viewport size (800x600 pixels/1200x768 pixels)
		
		The camera-to-scene distances are chosen per data set so as to achieve a medium amount of screen coverage by the geometry for a large distance as well as high coverage for a small distance.


Considering that both, \cite{foo} and \cite{foo}, apply dynamic tessellation strategies that maintain a topology that allows to identify conic sub-segments, it stands to reason that it might be more efficient to ray cast these instead of tessellating the cross-section.

No mathematical proof has been derived to support this claim. However, comparisons with numerically computed normals provide strong empirical evidence. strongly suggest this to be the case. It is worth noting, that this has been only empirically verified.

The interpolated color, in contrast to position and radius, is not computed from the quadratic approximations but by evaluating the original cubic spline.

- adaptive resampling of tubes, curvature stuff like the others
- ray footprint for tess
- cubic patches
- sphere tessellation
- better root finding
- anlytical normals for disk model
- analyt ortho vecs
- less precomp ortho vecs
- compare to interval/affine arithmetic or hybrid, diff root finder in gneral
- raycasr bounding boxes on CPU side, culling stuffs
- precision problem for looong segments, dynamic iteration count, precision bound
- pathtracing
- soft instancing
- tighter bounding geo, screen space bounding geo
- rendering with dynamic resolution for raycasting + supersampling with fixed cam
- pp aa -> geometrical aa, reprojection aa
- proof for rc normals
- AO volumes

In this paper, two spline tube rendering approaches were presented. A tessellation-based technique as well as a novel ray casting method were described and subsequently evaluated. 
using a numerical root finder
with respect to performance and quality characteristics.

The performance and quality characteristics of both techniques were compared an 

In this work, two spline tube rendering approaches were presented. A tessellation-based technique as well as a novel ray casting method were described and subsequently evaluated. 

provide an attractive approach to the exploration of various kinds of scientific data like vector fields


Every real interval of $f_t'$ corresponds to single ray-segment intersection.

Another problem of the ray casting approach are artifacts that occur when rendering long segments. 

They result from insufficient precision of the intersection computation in these cases due to the fixed iteration counts of the binary searches. Subdividing affected segments in a pre-processing step would allow to easily resolve these, though.

Another kind of artifact occurs when rendering long segments using the ray casting approach.

Furthermore, considering the additional implementation effort, the strategy of separating segment and node data and constructing the spline segments at run-time from the node data to reduce bandwidth demand turned not out to be worthwhile. Pre-computing the polynomial coefficients of the spline segments on the CPU and sending them as segment data to the GPU would be a more convenient and perhaps even more efficient approach.

\teaser{
 \includegraphics[width=\linewidth]{fibers_te.png, fibers_rc.png}
\includegraphics[width=\linewidth]{fibers_rc.png}
\includegraphics[width=\linewidth]{vectubes_rc.png}
 \centering
  \caption{New EG Logo NOOT}
\label{fig:teaser}
}

Tessellation-based approach:\begin{itemize}
	\item level-of-detail scheme applied to the base geometry
	\item curvature adaptive sampling for the base geometry
	\item analytical computation of normal vectors (if possible)
	\item reduced memory footprint due to pre-computed orthogonal vectors (smaller format, less vectors)
	\item removal of bulb artifacts
	\item coupling tessellation quality to target frame time
\end{itemize}

Ray casting-based approach:\begin{itemize}
	\item multi-resolution rendering with geometry-aware upsampling scheme
	\item utilization of conservative depth
	\item sorting the tube segments so as to reduce overdraw
	\item coupling rendering resolution to target frame time
	\item exploration of other root finding strategies
	\item tighter bounding geometry
	\item more efficient rendering of bounding geometry (f.e. by using soft instancing \cite{Soft})
	\item antialiasing (f.e. as a post-processing step using edge information that can be derived from the roots of $polyB_t$)
	\item mathematical proof for the correctness of the normal vector computation
	\item distance field-based deferred ambient occlusion
	\item high quality offline rendering using path tracing
\end{itemize}

	
	
	Results of the ray casting-based tubes rendering technique (first two images from the left) and of the tessellation-based approach (third and fourth image).