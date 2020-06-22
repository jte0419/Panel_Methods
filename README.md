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
* [YouTube SPM Streamline X and Y Geometric Integral (Mx(pj), My(pj)) Derivation](https://www.youtube.com/watch?v=BnPZjGCatcg)

## Source Panel Method System of Equations

Now that we have the normal velocity expression for each panel on the airfoil, we can set up a system of equations that can be easily solved.

* [YouTube Video](https://www.youtube.com/watch?v=ep7vPzGYsbw)

## Source Panel Method: Circular Cylinder

We have finally finished the derivations needed to code up a working version of the source panel method.  Recall that this implementation uses constant source panel strengths (which can vary from panel to panel).  This first implementation of the method is for the simplest case: flow over a circular cylinder.  The reason to use this simple case is that it has an analytical solution for the pressure coefficient that we can compare to (and it is also generally used as a source panel method comparison test case).  The code is provided for both MATLAB and Python.  The main code is the *SP_Circle* file, which needs both functions to run properly (*COMPUTE_IJ_SPM* and *STREAMLINE_SPM*).

* MATLAB Code
  * [SP_Circle.m](./SP_Circle.m)
  * [COMPUTE_IJ_SPM.m](./COMPUTE_IJ_SPM.m)
  * [STREAMLINE_SPM.m](./STREAMLINE_SPM.m)
* Python Code
  * [SP_Circle.py](./SP_Circle.py)
  * [COMPUTE_IJ_SPM.py](./COMPUTE_IJ_SPM.py)
  * [STREAMLINE_SPM.py](./STREAMLINE_SPM.py)
* [YouTube Video](https://www.youtube.com/watch?v=zIrDfEz-5mc)

## Source Panel Method: Airfoil

After making sure the simple validation case of the circular cylinder worked properly, we updated the code to be able to run with airfoils.  There are a couple more functions, programs, and directories needed when running this code.  You can run my code in either MATLAB or Python.  Whichever you choose, make sure to download all the files with that extension (*.m* or *.py*) from the list below.  You will also need the directory with the airfoils (*Airfoil_DAT_Selig*, make sure to extract after downloading the *zip* file) and the XFOIL executable (*xfoil.exe*) in the directory with all the code.  To run the code, open the *SP_Airfoil* script and press *Run* or *F5*.

* MATLAB Code
  * [SP_Airfoil.m](./SP_Airfoil.m)
  * [XFOIL.m](./XFOIL.m)
  * [COMPUTE_IJ_SPM.m](./COMPUTE_IJ_SPM.m)
  * [STREAMLINE_SPM.m](./STREAMLINE_SPM.m)
  * [COMPUTE_CIRCULATION.m](./COMPUTE_CIRCULATION.m)
* Python Code
  * [SP_Airfoil.py](./SP_Airfoil.py)
  * [XFOIL.py](./XFOIL.py)
  * [COMPUTE_IJ_SPM.py](./COMPUTE_IJ_SPM.py)
  * [STREAMLINE_SPM.py](./STREAMLINE_SPM.py)
  * [COMPUTE_CIRCULATION.py](./COMPUTE_CIRCULATION.py)
* Other
  * [Airfoil Directory](./Airfoil_DAT_Selig.zip)
  * [XFOIL Executable](https://web.mit.edu/drela/Public/web/xfoil/)
* [YouTube Video](https://www.youtube.com/watch?v=fdNOYdwY9Bw)

## Vortex Panel Methods: Geometric Integrals

For the source panel method, we derived the geometric integral expressions for the normal, tangential, X, and Y directions.  Here, we go through similar derivations, but this time for the vortex panel method.

* [YouTube VPM Normal Velocity Geometric Integral K(ij) Derivation](https://www.youtube.com/watch?v=5lmIv2CUpoc)
* [YouTube VPM Tangential Velocity Geometric Integral L(ij) Derivation](https://www.youtube.com/watch?v=IxWJzwIG_gY)
* [YouTube VPM Streamline X and Y Geometric Integral (Nx(pj), Ny(pj)) Derivation](https://www.youtube.com/watch?v=TBwBnW87hso)

## Vortex Panel Method System of Equations

In the same way that we wrote the system of equations for the source panel method, we can write the system of equations for the vortex panel method.  To include the Kutta condition equation in the matrix, we must remove one of the normal velocity equations for one of the control points.

* [YouTube Video](https://www.youtube.com/watch?v=j3ETHFBiYOg)

## Vortex Panel Method: Airfoil

The source panel method code was updated to be able to solve for the vortex panel strengths instead of the source panel strengths.  You can run my code in either MATLAB or Python.  Whichever you choose, make sure to download all the files with that extension (*.m* or *.py*) from the list below.  You will also need the directory with the airfoils (*Airfoil_DAT_Selig*, make sure to extract after downloading the *zip* file) and the XFOIL executable (*xfoil.exe*) in the directory with all the code.  To run the code, open the *VP_Airfoil* script and press *Run* or *F5*.  The limitations of this VPM implementation are shown in my YouTube video, and motivate the need for the combined source/vortex panel method.

* MATLAB Code
  * [VP_Airfoil.m](./VP_Airfoil.m)
  * [XFOIL.m](./XFOIL.m)
  * [COMPUTE_KL_VPM.m](./COMPUTE_KL_VPM.m)
  * [STREAMLINE_VPM.m](./STREAMLINE_VPM.m)
  * [COMPUTE_CIRCULATION.m](./COMPUTE_CIRCULATION.m)
* Python Code
  * [VP_Airfoil.py](./VP_Airfoil.py)
  * [XFOIL.py](./XFOIL.py)
  * [COMPUTE_KL_VPM.py](./COMPUTE_KL_VPM.py)
  * [STREAMLINE_VPM.py](./STREAMLINE_VPM.py)
  * [COMPUTE_CIRCULATION.py](./COMPUTE_CIRCULATION.py)
* Other
  * [Airfoil Directory](./Airfoil_DAT_Selig.zip)
  * [XFOIL Executable](https://web.mit.edu/drela/Public/web/xfoil/)
* [YouTube Video](https://www.youtube.com/watch?v=JL2fz-xTTT0)

## Source/Vortex Panel Method System of Equations

The previous video on the constant strength vortex panel method showed that there are still some limitations with this simplified method, and to improve it we can combine the source and vortex panel methods into one source panel/vortex panel (SPVP) method.  This video covers how to build the system of equations for the SPVP method and how to derive and apply a different version of the Kutta condition.

* [YouTube Video](https://www.youtube.com/watch?v=bc_pkKGEypU)

## Source/Vortex Panel Method: Airfoil (COMING SOON!)

This video will implement the necessary changes in the code for the SPVP method.  We will look at a bunch of different airfoil cases that happened to fail in the original VPM implementation.  We will compare results to XFOIL.

## Multi-Element Source/Vortex Panel Method (COMING SOON!)

This will be the final video of the series, where I adapt the code from the previous section to be able to load multiple airfoils.  I'll show how to load airfoils, manipulate their locations relative to each other, and how to implement the necessary changes to the code to handle multiple airfoils.  We will go over a couple simple examples to test the validity, and then hopefully I'll be able to show a validation case to an airfoil system from a paper.
