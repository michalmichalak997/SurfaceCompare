# SurfaceCompare. 
This software offers possibility to conduct measurements of angular distances between geological interfaces.

A freely available C++ code (angular_distances.cpp and orientation_maps.cpp) is required to produce files for the statistical and spatial portion of the research. These C++ programs require a text file as input with the following columns: XY coordinates common for all horizons as two first columns, n columns representing n horizons and the id pointing to the row number as the last column. Then, two R files produce statistical and spatial outputs. Topography of the analyzed interfaces as well as the isopach maps were prepared using akima package (Akima, 1978). The above method fall within the scope of non-geostatistical interpolation methods based on bilinear or bicubic splines (Akima, 1978; Bivand et al., 2013). The geostatistical analysis for angular distance was conducted using autokrige function available in automap package (Hiemstra et al., 2009). This package performs an automatic interpolation by estimating the variogram and then calling gstat package (Gräler et al., 2016; Pebesma, 2004). 

Developer: Michał Michalak (michalmichalak@us.edu.pl)

Affiliations: 
1) Faculty of Natural Sciences, University of Silesia in Katowice, Poland, Będzińska 60, 41-205 Sosnowiec.
2) Faculty of Geology, Geophysics and Environmental Protection, Poland, aleja Adama Mickiewicza 30, 30-059 Kraków.

## Input

![input_file](https://user-images.githubusercontent.com/28152295/161118052-c1a3eaf4-55e8-4f92-858a-fb34fb64dff1.png)

The explanation of the structure of the file is given below:

![input_file_explanation](https://user-images.githubusercontent.com/28152295/161118209-17dc606a-e266-4c82-b4e1-cbeca6f4973c.png)

## Processing

The program looks as follows. You would need to specify the number of horizons and the path.

![processing](https://user-images.githubusercontent.com/28152295/161123010-b3f50a13-f1bb-4d2b-9c96-4c537a456ea6.png)


## Output

If you have 4 input horizons, you will get 6 output files. They will have extensions according to the indices of the taken horizons. In this case, these six files will have the following extensions: _01, _02, _03, _12, _13, and _23.

The file is very wide, so we do not attach a screenshot.

You should delete the first line in your file that includes information about the taken horizons. For example, in the file with extension _01 you should delete the first line "Index of the first surface:0Index of the first surface:1".

## Orientation maps

To create grid maps used in the article you need to use the orientation_maps.cpp file. The structure of the input is the same. The difference is that you need to specify the spacing of the regular grid. To obtain the grid map, you will specifically need the file with the extension _grid_locate.

Please note that in the below example the spacing between points in the regular grid will be 150 m.

![processing_gridmaps](https://user-images.githubusercontent.com/28152295/161153466-fa068793-141a-46d4-939c-a58914b2d853.png)

## System requirements

### System required

Windows 10 Home, 64-bit, processor x64.

### Versions the software has been tested on

Windows 10 Home, 64-bit, processor x64.

### Software required 

Version numbers can be different, those used in our work are given: 

  -Microsoft Visual Studio 2017, 
  
  -CGAL library (ver. 4.8),
  
  -Boost library (ver. 1.59), 
  
  -Microsoft Visual C++ 2017 Redistributable
  
### Program language

C++

### License

GNU General Public License v3.0.

# SurfaceCompare. Move to R - part 1 - Statistical analysis.

Use the Statistical_analysis_angular_distances file.

# SurfaceCompare. Move to R - part 2 - Spatial analysis.

Use the Spatial_analysis_angular_distances file.

# SurfaceCompare. Move to R - part3 - Grid map.

Use the gridmap file.
