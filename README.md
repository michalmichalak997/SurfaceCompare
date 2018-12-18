# SurfaceCompare

# SurfaceCompare
Comparing the orientation of surfaces
The application SurfaceCompare aims to compare the orientation of geological surfaces. 

Developer: Michał Michalak (mimichalak@us.edu.pl), Department of Applied Geology, Faculty of Earth Sciences, University of Silesia, Poland, Będzińska 60, 41-205 Sosnowiec.

1. System requirements

a) System required: Windows 10 Pro 64-bit, Windows 8.1 64-bit. 

b) Hardware required: As specified in the Windows 10 Pro or Windows 8.1 requirements. 

c) Versions the software has been tested on: Windows 10 Pro 64-bit, Windows 8.1 64-bit.

d) Software required (version numbers can be different, those used in our work are given): 

  -Microsoft Visual Studio 2017, 
  
  -CGAL library (ver. 4.8),
  
  -Boost library (ver. 1.59), 
  
  -Microsoft Visual C++ 2017 Redistributable, 
  
  -libgmp-10.dll, libmpfr-4.dll (attached). 
  
e) Program language: C++. 

f) Program size: 1.1 MB. 

g) License: GNU General Public License v3.0.

2. Installation guide

The below instructions are applicable to run the standalone version.

a) Download SurfaceCompare.exe along with libgmp-10.dll and libmpfr-4.dll. 

b) Put all three above files in the same directory.

c) Run SurfaceCompare.exe

3. Demo and Instructions for use

a) To test the program, please download TestDataBottom.txt and TestDataTerrain.txt. For convenience put them
in the same directory. Having run the SurfaceCompare.exe, you should specify the paths of these input files.
For instance if you put them in C:\CGAL\examples, you should type in: C:\CGAL\examples\TestDataBottom.txt and
C:\CGAL\examples\TestDataTerrain.txt. You should also specify the path of the output file e.g. C:\CGAL\examples\TestResultsBottomTerrain.txt

b) The expected output is as follows:

X_C;Y_C;Distance

920580.573333;254324.183333;1.43308

920834.660000;253472.150000;1.29936

920875.600000;253202.870000;0.856687

920680.523333;253888.936667;1.41487

920979.283333;253888.296667;1.8032

920879.333333;254323.543333;1.37324

921159.243333;253340.333333;1.61983

921231.190000;253582.310000;1.76064

921286.213333;253160.110000;1.61581

921216.716667;252929.303333;1.46042

921340.186667;252801.873333;1.50106


c) The expected run time for demo on a "normal" desktop computer should be no longer than five seconds.
