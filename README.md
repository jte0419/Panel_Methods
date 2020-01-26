# Panel_Methods

This repository contains all the code related to panel methods.  This includes scripts and functions, in both MATLAB and Python.  It will be updated continuously as I finish the video series with the final Source Panel Method (SPM) and Vortex Panel Method (VPM) code.  Explanations of the different files can be found in the sections below, which are pretty much the same as you'll find on my [website launching page](http://www.joshtheengineer.com/panel-methods/).  Links to my YouTube videos and blog posts on my website are included.  Note that functions are always named in all caps.

## UIUC Airfoil Database File Download

This section is for downloading all the airfoil coordinates from the [UIUC Airfoil Database](https://m-selig.ae.illinois.edu/ads/coord_database.html).  The files are all of type *DAT*, which can be opened in any text editor.  This code is run in Python, since it's free to download.  There are two scripts: one that downloads the files in the default coordinate format (which varies quite a bit), and one that downloads the files in Selig format (which is easier to load into your program of choice).  If I were you, I would download them in Selig format.

* [Get_Save_Airfoils.py](./Get_Save_Airfoils.py)
* [Get_Save_Airfoils_Selig.py](./Get_Save_Airfoils_Selig.py)
* [YouTube Video](https://www.youtube.com/watch?v=nILo18DlqAo)
* [Blog Post](http://www.joshtheengineer.com/2019/01/30/uiuc-airfoil-database-file-download/)

## Load Airfoil Coordinates in MATLAB and Python

With the airfoils download from the UIUC database in the previous section, we now need to read them into our MATLAB/Python code in order to use them in the panel code.  Note that if you're planning on using Python, you need to download the files into Selig format from the previous section, since that's all I coded for the loading function in Python.

* [Test_Load_Airfoil_Function.m](./Test_Load_Airfoil_Function.m)
* [LOAD_AIRFOIL.m](./LOAD_AIRFOIL.m)
* [LOAD_AIRFOIL_SELIG.m](./LOAD_AIRFOIL_SELIG.m)
* [LOAD_AIRFOIL_SELIG.py](./LOAD_AIRFOIL_SELIG.py)
* [YouTube Video](https://www.youtube.com/watch?v=xJYxMfGFrk8)
* [Blog Post](http://www.joshtheengineer.com/2019/03/30/load-airfoil-data-into-matlab/)

## Running XFoil from MATLAB and Python

This section is for running XFoil froma  script in MATLAB/Python.  This is convenient if you want to be able to quickly change a parameter (or multiple parameters) and output the results, but you don't want to go into the command line of XFoil every time.  This is useful for getting airfoil data in a MATLAB GUI for instance.  If you have an edit text box for changing the angle of attack, then every time it's changed by the user, the MATLAB program will call XFoil using the updated AoA value and return the airfoil coordinates and pressure coefficient data (for example).

* [MATLAB_XFOIL.m](./MATLAB_XFOIL.m)
* [Python_XFoil.py](./Python_XFoil.py)
* [YouTube Video](https://www.youtube.com/watch?v=bWjo3N9COz4)
* [Blog Post - MATLAB](http://www.joshtheengineer.com/2019/01/30/running-xfoil-from-matlab/)
* [Blog Post - Python](http://www.joshtheengineer.com/2019/02/06/running-xfoil-from-python/)

## Compute Circulation

Circulation is a fundamental concept in aerodynamics (and more generally, in multivariable calculus).  In order to find the lift of an airfoil using the vortex panel method, the circular around the airfoil needs to be computed.  For pretty much all of the remaining codes, you'll need to have either the MATLAB or Python *COMPUTE_CIRCULATION* function in the run directory.  The script *Circulation* is an example file showing how to use the *COMPUTE_CIRCULATION* function, and is not needed for future codes.

* [Circulation.m](./Circulation.m)
* [COMPUTE_CIRCULATION.m](./COMPUTE_CIRCULATION.m)
* [Circulation.py](./Circulation.py)
* [COMPUTE_CIRCULATION.py](./COMPUTE_CIRCULATION.py)
* [YouTube Video](https://www.youtube.com/watch?v=b8EnhiSjL3o)
* [Blog Post](http://www.joshtheengineer.com/2019/04/01/compute-circulation-of-a-vector-field-in-matlab-and-python/)

## Incompressible Potential Flow

When assuming that the flow is both irrotational and incompressible, we can simplify down some complicated equations into a simpler equation (Laplace's equation), and the resulting flow is called incompressible potential flow.  There is no code here, but links to my explanation YouTube video and blog post can be found below.

* [YouTube Video](https://www.youtube.com/watch?v=Zo5XIcX8s2Q)
* [Blog Post](http://www.joshtheengineer.com/2019/05/05/introduction-to-incompressible-potential-flow/)

## Elementary Potential Flows

In order to build more complex flows, we need to understand the simplest types of potential flows. Here, we go through three elementary flows (uniform, source/sink, and vortex), along with two combinations of these flows (uniform + source/sink, uniform + vortex).

* Uniform Flow
  * [Uniform_Flow.m](./Uniform_Flow.m)
  * [Uniform_Flow.py](./Uniform_Flow.py)
  * [YouTube Video](https://www.youtube.com/watch?v=jCTTDclJZEk)
  * [Blog Post](www.joshtheengineer.com/2019/08/25/elementary-flow-uniform-flow/)
* Source/Sink Flow
  * [Source_Sink_Flow.m](./Source_Sink_Flow.m)
  * [Source_Sink_Flow.py](./Source_Sink_Flow.py)
  * [YouTube Video](https://www.youtube.com/watch?v=eLDI_jV3yo0)
  * [Blog Post](http://www.joshtheengineer.com/2019/08/25/elementary-flow-source-sink-flow/)
* Vortex Flow
  * [Vortex_Flow.m](./Vortex_Flow.m)
  * [Vortex_Flow.py](./Vortex_Flow.py)
  * [YouTube Video](https://www.youtube.com/watch?v=61jvr3rtmLE)
  * [Blog Post](http://www.joshtheengineer.com/2019/08/25/elementary-flow-vortex-flow/)
* Uniform + Source/Sink Flow
  * [Uniform_Source_Sink_Flow.m](./Uniform_Source_Sink_Flow.m)
  * [Uniform_Source_Sink_Flow.py](./Uniform_Source_Sink_Flow.py)
  * [YouTube Video](https://www.youtube.com/watch?v=zIvpN9f9dAA)
  * [Blog Post](http://www.joshtheengineer.com/2019/08/25/combined-flow-uniform-and-source-flow/)
* Uniform + Vortex Flow
  * [Uniform_Vortex_Flow.m](./Uniform_Vortex_Flow.m)
  * [Uniform_Vortex_Flow.py](./Uniform_Votex_Flow.py)
  * [YouTube Video](https://www.youtube.com/watch?v=SoMuRp5v16w)
  * [Blog Post](http://www.joshtheengineer.com/2019/08/25/combined-flow-uniform-and-vortex-flow/)

## Panel Method Geometry

The first step in writing your own panel method code is to understand the geometry and its associated variables.  Here, we go through this in detail for an arbitrary shape (circle approximated by eight panels).  In the MATLAB and Python codes, we also show how this works for an airfoil.  Note that for these codes to work, you will need to have the appropriate *[LOAD_AIRFOIL_SELIG](./LOAD_AIRFOIL_SELIG.m)* function downloaded in the directory (unless you want to comment out the airfoil loading section).

* [Panel_Method_Geometry.m](./Panel_Method_Geometry.m)
* [Panel_Method_Geometry.py](./Panel_Method_Geometry.py)
* [YouTube Video](https://www.youtube.com/watch?v=kIqxbd937PI)
* [Blog Post](http://www.joshtheengineer.com/2019/10/09/panel-method-geometry/)

## Building More Complex Potential Flows

In the *Elementary Potential Flows* section, we went over all the relevant elementary incompressible potential flows.  We can now use these to build up a more complex flow step by step.  The end result of this section is to obtain an expression for the velocity potential induced at an arbitrary point P in the flow due to uniform flow and N source panels (that approximate the airfoil geometry).

* [YouTube Video](https://www.youtube.com/watch?v=EKzbwJvKcmw)

## Flow Around an Airfoil

How do we use the velocity potential equation derived in the previous section to compute the flow around an airfoil?  This section explains the necessary conditions that we can use to compute the unknown source panel strengths.  We use the boundary condition for an impermeable (solid) object to stipulate that the normal velocity at each panel control point should be zero.

* [YouTube Video](https://www.youtube.com/watch?v=cLdv1UfX1g8)

## Source Panel Method: Geometric Integrals

When we take the appropriate partial derivative (normal, tangential, X, or Y) of the velocity potential equation, we end up with a complicated integral expression inside the source panel term.  In order to code up our system of equations in the next section, we need to find an explicit expression for this integral (also called the *geometric integral* since it only depends on the airfoil's geometry).  The videos below have the full derivations for the normal, tangential, X, and Y geometric integrals.  You should at least watch the I(ij) video since it includes the entire derivation, whereas the other videos don't repeat the portions that are the same as the I(ij) derivation video (to keep them shorter).

* [YouTube SPM Normal Velocity Geometric Integral I(ij) Derivation](https://www.youtube.com/watch?v=76vPudNET6U)
* [YouTube SPM Tangential Velocity Geometric Integral J(ij) Derivation](https://www.youtube.com/watch?v=JRHnOsueic8)
* [YouTube SPM Streamline X and Y Geometric Integral (Mx(ij), My(ij)) Derivation](https://www.youtube.com/watch?v=BnPZjGCatcg)

## Source Panel Method System of Equations

Now that we have the normal velocity expression for each panel on the airfoil, we can set up a system of equations that can be easily solved.

* [YouTube Video](https://www.youtube.com/watch?v=ep7vPzGYsbw)


