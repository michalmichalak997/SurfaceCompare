# SurfaceCompare. Introduction.

Measurements of angular distances between geological interfaces.

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

## System requirements

### System required

Windows 10 Home, 64-bit, processor x64.

### Versions the software has been tested on

Windows 10 Home, 64-bit, processor x64.

### Software required 

Vversion numbers can be different, those used in our work are given: 

  -Microsoft Visual Studio 2017, 
  
  -CGAL library (ver. 4.8),
  
  -Boost library (ver. 1.59), 
  
  -Microsoft Visual C++ 2017 Redistributable
  
### Program language

C++

### License

GNU General Public License v3.0.
