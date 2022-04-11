#include <fstream>
#include <sstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>
#include <random>


using namespace std;


typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel> Gt;
typedef CGAL::Triangulation_vertex_base_with_info_2< std::pair<vector<double>, unsigned int>, Gt > Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Gt, Tds>                    Delaunay;

typedef Kernel::Point_3                                                Point;

const int n = 3; //we work in 3D
const double ex = 2; //we introduce the restriction of collinearity


class plane //a class that stores the crucial figures in terms of computing the orientation
{

private:

	double first_vec[n];            //the first edge of a triangle
	double second_vec[n];			//the second edge of a triangle
	double third_vec[n];			//the third edge of a triangle
	double normal_vec[n];			//normal vector of a triangle
	double directional[n];			//the projection of the normal vector onto the horizontal plane
	double z_axis[n] = { 0,0,1 };   //the definition of the z-axis
	double dip_vec[n];
	double doc;						//a variable that contains the collinearity coefficient
	double area;					//a variable that stores the area of a triangle
	bool lin_dependence;		    //a bool variable to check to answer whether points are (too) collinear
	string dip_degrees;             //a text variable to store the dip angle
	string azimuth_degrees;         //a text variable to store the dip direction
	Point thisp1;       //three points representing a surface
	Point thisp2;
	Point thisp3;


public:
	vector<string> get_Points() {
		//vector<Point> points = { thisp1, thisp2, thisp3 };

		vector<string> coordinates{ std::to_string(thisp1.x()), std::to_string(thisp1.y()), std::to_string(thisp1.z()),
									std::to_string(thisp2.x()), std::to_string(thisp2.y()), std::to_string(thisp2.z()),
									std::to_string(thisp3.x()), std::to_string(thisp3.y()), std::to_string(thisp3.z())
		};
		return coordinates;

	}


	double dip_azimuth(double normal[], int n = 2) //a function that computes the dip azimuth
	{
		double coeff = 180 / M_PI;
		double angle = atan2(normal[1], normal[0]);
		angle = angle * coeff;
		if (angle < 0)
		{
			return angle + 360;
		}
		else
		{
			return angle;
		}
	}

	double dip_angle(double z_axis[], double normal_v[]) //function that computes the dip angle
	{
		double angle;
		double expression;
		double coeff = 180 / M_PI;
		expression = abs(dot_product(z_axis, normal_v)) / (length(z_axis)*length(normal_v));
		angle = acos(expression);
		return angle * coeff;
	}

	double dot_product(double vector_line[], double direction[], int n = 3) //function that computes the dot product of vectors
	{
		double product = 0;
		for (int i = 0; i < n; i++)
		{
			product += direction[i] * vector_line[i];
		}
		return product;
	}

	bool dependence(double v1[], double v2[], double v3[]) //function that checks whether the points are collinear
	{
		double len_v1 = length(v1);
		double len_v2 = length(v2);
		double len_v3 = length(v3);
		double lengths[n] = { len_v1, len_v2, len_v3 };

		sort(lengths, lengths + n);
		this->doc = lengths[2] / (lengths[0] + lengths[1]);
		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (lengths[i] == 0)
			{
				//k += 1;
				throw runtime_error("Points coincidence");

			}
		}
		if (k != 0)
		{
			return true;
		}
		else
		{
			if (doc > ex)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}



	double length(double line_vector[], int n = 3) //function that computes the length of a vector
	{
		double vector_length = sqrt(pow(line_vector[0], 2) + pow(line_vector[1], 2) + pow(line_vector[2], 2));
		return vector_length;
	}


	vector<string> center() //a function that computes the geometrical centre of a Delaunay triangle
	{
		double x = (thisp1.x() + thisp2.x() + thisp3.x()) / (3.0);
		double y = (thisp1.y() + thisp2.y() + thisp3.y()) / (3.0);
		double z = (thisp1.z() + thisp2.z() + thisp3.z()) / (3.0);
		//double z = (thisp1[2] + thisp2[2] + thisp3[2]) / (3.0);
		vector<string> center_triangle{ std::to_string(x), std::to_string(y), std::to_string(z) };
		return center_triangle;

	}

	plane(Point point_1, Point point_2, Point point_3) //the class constructor
	{
		double coeff = 180 / M_PI;


		this->thisp1 = Point(point_1.x(), point_1.y(), point_1.z());
		this->thisp2 = Point(point_2.x(), point_2.y(), point_2.z());
		this->thisp3 = Point(point_3.x(), point_3.y(), point_3.z());

		double first_try[n] = { point_2.x() - point_1.x(), point_2.y() - point_1.y(), point_2.z() - point_1.z() };
		double second_try[n] = { point_3.x() - point_1.x(), point_3.y() - point_1.y(), point_3.z() - point_1.z() };
		double third_try[n] = { point_3.x() - point_2.x(), point_3.y() - point_2.y(), point_3.z() - point_2.z() };

		bool test = dependence(first_try, second_try, third_try);
		if (test == true)
		{
			lin_dependence = true;
		}
		else
		{
			lin_dependence = false;
			for (int i = 0; i < n; i++)
			{
				this->first_vec[i] = first_try[i];
				this->second_vec[i] = second_try[i];
				this->third_vec[i] = third_try[i];
			}
			normal_vec[0] = first_vec[1] * second_vec[2] - second_vec[1] * first_vec[2];
			normal_vec[1] = first_vec[2] * second_vec[0] - second_vec[2] * first_vec[0];
			normal_vec[2] = first_vec[0] * second_vec[1] - second_vec[0] * first_vec[1];

			if (normal_vec[2] < 0) {
				normal_vec[0] *= -1;
				normal_vec[1] *= -1;
				normal_vec[2] *= -1;
			}

			double normal_vector_length = length(normal_vec);

			normal_vec[0] /= normal_vector_length;
			normal_vec[1] /= normal_vector_length;
			normal_vec[2] /= normal_vector_length;

			this->dip_vec[0] = cos(dip_angle(z_axis, normal_vec) / coeff)*cos(dip_azimuth(normal_vec) / coeff);
			this->dip_vec[1] = cos(dip_angle(z_axis, normal_vec) / coeff)*sin(dip_azimuth(normal_vec) / coeff);
			this->dip_vec[2] = -sin(dip_angle(z_axis, normal_vec) / coeff);

			double stala = 0.5; //using Heron formula to calculate area 
			double half = stala * (length(first_vec) + length(second_vec) + length(third_vec)); // Let p=0.5(a+b+c) be the half of the circumference.
			double s = sqrt(half*(half - length(first_vec))*(half - length(second_vec))*(half - length(third_vec))); // Then the area is sqrt(p*(p-a)(p-b)(p-c)))
			this->area = s;
		}
	}

	string measure()//function that supplies orientation results also for singularities
	{
		if (lin_dependence)
		{
			azimuth_degrees = ("LT");
			dip_degrees = ("LT");
			return (dip_degrees + ";" + azimuth_degrees);
		}
		else if (normal_vec[0] == 0 && normal_vec[1] == 0 && normal_vec[2] != 0)
		{
			dip_degrees = "0";
			azimuth_degrees = ("x");
			return (dip_degrees + ";" + azimuth_degrees);
		}
		else if (normal_vec[2] == 0)
		{
			dip_degrees = "90";
			azimuth_degrees = to_string(dip_azimuth(normal_vec));
			return dip_degrees + ";" + azimuth_degrees;
		}
		else
		{
			double dipping_angle = dip_angle(z_axis, normal_vec);
			dip_degrees = to_string(dipping_angle);
			azimuth_degrees = to_string(dip_azimuth(normal_vec));
			return dip_degrees + ";" + azimuth_degrees;
		}
	}

	vector<double> get_normal() //function that -computes- returns the normal vector
	{

		vector<double> normal_vector = { normal_vec[0] ,normal_vec[1], normal_vec[2] };

		return normal_vector;
	}

	double get_area() {

		return(area);
	}

	double get_doc() {

		return(doc);
	}

	vector<double> get_dip_vec() {

		vector<double> dip_vector = { dip_vec[0],dip_vec[1] ,dip_vec[2] };
		return(dip_vector);
	}
};


int main()
{

	int n_surfaces;
	cout << "How many surfaces do you have?" << endl;
	cin >> n_surfaces;

	string path_i, path_o, path_del, path_nor, path_grid, path_gridvis; //text variables for input and output paths, respectively
	double resolution_step; //the density of the grid map

	std::cout << "Type in the path of your input data:" << endl; //the user is required to type in the input path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\JurassicBottomInput.txt" << endl << endl;

	std::cin >> path_i;

	ifstream download(path_i);
	vector<double> surface_elevation;
	vector< std::pair<Point, std::pair<vector<double>, unsigned>> > pts; //a variable storing points representing the first surface

	vector <double> coord_tab;


	if (!download) std::cout << "Error in opening file" << endl; //the case when the file cannot be uploaded

	string tempor;//a temporary variable storing figures while uploading

	while (getline(download, tempor)) // loading points line - by - line associated with the first surface
	{
		istringstream convert(tempor);

		//double x, y, z1, z2;
		unsigned int t;
		double data;
		for (int i = 0; i < n_surfaces + 3; i++) {
			if (!(convert >> data)) { break; }
			coord_tab.push_back(data);
		}
		//if (!(convert >> x >> y)) { break; }

		//if (!(convert >> x >> y >> z1 >> z2 >> t)) { break; }
		for (auto i = 3; i < n_surfaces + 2; i++) {
			surface_elevation.push_back(coord_tab[i]); //elevations after first elevation
		}

		pts.push_back(make_pair(Point(coord_tab[0], coord_tab[1], coord_tab[2]), make_pair(surface_elevation, coord_tab[coord_tab.size() - 1])));

		coord_tab.clear(); //avoiding problems with indexing
		surface_elevation.clear();
	}


	double min_x = pts.begin()->first[0];
	double min_y = pts.begin()->first[1];

	double max_x = pts.begin()->first[0];
	double max_y = pts.begin()->first[1];

	for (auto it = pts.begin(); it != pts.end(); it++) { //calculating boundary coordinates

		if (it->first[0] < min_x) {
			min_x = it->first[0];
		}

		if (it->first[1] < min_y) {
			min_y = it->first[1];
		}

		if (it->first[0] > max_x) {
			max_x = it->first[0];
		}

		if (it->first[1] > max_y) {
			max_y = it->first[1];
		}

	}


	std::cout << "Type in the path of the output:" << endl; //the user is required to type in the output path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\JurassicBottomOutput.txt" << endl << endl;

	for (auto i = 0; i < 4; i++) { //remove the extension
		path_i.pop_back();
	}
	path_o = path_i + "_output";
	std::cout << path_o << endl << endl;

	std::cout << "Type in the path of the Delaunay visualization .vtu file:" << endl; //the user is required to type in the output path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\Delaunay.vtu" << endl << endl;

	path_del = path_i + "_delaunay";
	std::cout << path_del << endl << endl;

	std::cout << "Type in the path of the normals .vtu file:" << endl; //the user is required to type in the output path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\normals.vtu" << endl << endl;

	path_nor = path_i + "_normals";
	std::cout << path_nor << endl << endl;

	std::cout << "Type in the path of gridpath .txt file:" << endl; //the user is required to type in the output path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\gridpath.txt" << endl << endl;

	path_grid = path_i + "_grid_locate.txt";
	std::cout << path_grid << endl << endl;

	std::cout << "Type in the path of grid .vtu file:" << endl; //the user is required to type in the output path
	std::cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\gridvis.vtu" << endl << endl;

	path_gridvis = path_i + "_grid_paraview.vtu";
	std::cout << path_gridvis << endl << endl;

	std::cout << "Type in the grid resolution: (e.g. 100.00)" << endl;
	std::cin >> resolution_step;

	Delaunay dt; //a variable storing the geometrical elements of Delaunay triangulation
	dt.insert(pts.begin(), pts.end());

	cout << "The number of points taken:" << pts.size() << ". The number of vertices in triangulation:" << dt.number_of_vertices() << endl;

	if (pts.size() != dt.number_of_vertices()) { throw(runtime_error("Check for duplicates in data!")); }

	ofstream gridsave(path_grid);

	vector<Point> grid_pts;
	double elevation_grid_pts = 0.00; //the z-coordinate is equal to zero (arbitrary decision)

	for (auto i = min_x; i < max_x; i = i + resolution_step) {
		for (auto j = min_y; j < max_y; j = j + resolution_step) {
			grid_pts.push_back(Point(i, j, elevation_grid_pts));
		}
	}

	gridsave << "px;" << "py;" << "IDT1;" << "IDT2;" << "IDT3" << endl;

	vector<Point> grid_pts_finite;

	for (auto regpoint : grid_pts) {
		Delaunay::Face_handle facereg = dt.locate(regpoint);


		int inf1 = facereg->vertex(0)->info().second;
		int inf2 = facereg->vertex(1)->info().second;
		int inf3 = facereg->vertex(2)->info().second;

		if ((inf1 > 0) && (inf2 > 0) && (inf3 > 0)) //we are interested only in points within the convex hull
		{
			grid_pts_finite.push_back(Point(regpoint.x(), regpoint.y(), elevation_grid_pts));
			gridsave << fixed << regpoint.x() << ";" << regpoint.y() << ";" << inf1 << ";" << inf2 << ";" << inf3 << endl;
		}
	}

	ofstream gridvis(path_gridvis);

	gridvis <<
		"<?xml version=\"1.0\"?>" << "\n" <<
		"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << "\n  " <<
		"<UnstructuredGrid>" << "\n    " <<
		"<Piece NumberOfPoints=\"" << grid_pts_finite.size() << "\" NumberOfCells=\"" << grid_pts_finite.size() << "\">" << "\n      "
		"<PointData Scalars=\"scalars\">" << "\n        "
		"<DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">" << "\n          ";

	for (auto it = grid_pts_finite.begin(); it != grid_pts_finite.end(); it++)
	{
		gridvis << it->z() << "\n          "; //the z-coordinate is equal to zero (zero stored in elevation_grid_pts - arbitrary decision)
	}

	gridvis << "\n        " <<
		"</DataArray>" << "\n      " <<
		"</PointData>" << "\n      " <<
		"<Points>" << "\n        " <<
		"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n           ";

	for (auto it = grid_pts_finite.begin(); it != grid_pts_finite.end(); it++)
	{
		gridvis << fixed << it->y() << " " << it->x() << " " << it->z() << "\n           ";
	}

	gridvis << "\n        " <<
		"</DataArray>" << "\n      " <<
		"</Points>" << "\n      " <<
		"<Cells>" << "\n        " <<
		"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n          ";


	for (auto i = 0; i < grid_pts_finite.size(); i++) //zamienic liczbe 24 na grid_pts size
	{
		gridvis << i << "\n          ";
	}

	gridvis << "\n        " <<
		"</DataArray>" << "\n        " <<
		"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n          ";

	for (auto i = 1; i <= 1 * grid_pts_finite.size(); i = i + 1)
	{
		gridvis << i << " ";
	}

	gridvis << "\n        " <<
		"</DataArray>" << "\n        " <<
		"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<
		"\n          ";

	for (auto i = 1; i <= grid_pts_finite.size(); i++)
	{
		gridvis << 1 << " ";
	}

	gridvis << "\n        " <<
		"</DataArray>" << "\n      " <<
		"</Cells>" << "\n    " <<
		"</Piece>" << "\n  " <<
		"</UnstructuredGrid>" << "\n" <<
		"</VTKFile>";


	for (auto i = 0; i < n_surfaces; i++)
	{
		int surface_index = i;

		string path_output = path_o + "_" + to_string(surface_index) + ".txt";
		ofstream saving(path_output); //a stream variable to save output figures

		saving << "Index of the surface:" << surface_index << endl;

		saving << "X1;" << "Y1;" << "Z1;" << "X2;" << "Y2;" << "Z2;" << "X3;" << "Y3;" << "Z3;" << "X_C;" //column names are saved
			<< "Y_C;" << "Z_C;" << "X_N;" << "Y_N;" << "Z_N;" << "X_D;" << "Y_D;" << "Z_D;" << "Dip_ang;" << "Dip_dir;" << "DOC;" << "Area;" << "IDT1;" << "IDT2;" << "IDT3" << endl;


		vector <plane> planes; //a vector variable storing planes representing the second surface
		std::pair<Point, int> point_1; //we have a pair with point and index of the point
		std::pair<Point, int> point_2;
		std::pair<Point, int> point_3;

		auto plane_index = 0;

		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) //a loop for performing the Delaunay triangulation and save the results

		{

			Delaunay::Face_handle face = fit;
			if (surface_index == 0) {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], dt.triangle(face)[0][2]), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], dt.triangle(face)[1][2]), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], dt.triangle(face)[2][2]), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}
			else {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(surface_index - 1)), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(surface_index - 1)), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(surface_index - 1)), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}

			}

			string result = planes.at(plane_index).measure();							//extracting the dip angle and the dip direction
	
			saving <<
				planes.at(plane_index).get_Points()[0] << ";" <<//coordinates start
				planes.at(plane_index).get_Points()[1] << ";" <<
				planes.at(plane_index).get_Points()[2] << ";" <<
				planes.at(plane_index).get_Points()[3] << ";" <<
				planes.at(plane_index).get_Points()[4] << ";" <<
				planes.at(plane_index).get_Points()[5] << ";" <<
				planes.at(plane_index).get_Points()[6] << ";" <<
				planes.at(plane_index).get_Points()[7] << ";" <<
				planes.at(plane_index).get_Points()[8] << ";" << //coordinates end
				planes.at(plane_index).center()[0] << ";" <<
				planes.at(plane_index).center()[1] << ";" <<
				planes.at(plane_index).center()[2] << ";" <<
				planes.at(plane_index).get_normal()[0] << ";" <<  //normal vector start
				planes.at(plane_index).get_normal()[1] << ";" <<
				planes.at(plane_index).get_normal()[2] << ";" << //normal vector end
				planes.at(plane_index).get_dip_vec()[0] << ";" << //dip vector start
				planes.at(plane_index).get_dip_vec()[1] << ";" <<
				planes.at(plane_index).get_dip_vec()[2] << ";" << //dip vector end
				planes.at(plane_index).measure() << ";" <<
				planes.at(plane_index).get_doc() << ";" << //collinearity
				planes.at(plane_index).get_area() << ";" << //area
				face->vertex(0)->info().second << ";" << //id start
				face->vertex(1)->info().second << ";" <<
				face->vertex(2)->info().second << endl; //id end

			plane_index++; //incrementing index;
		}
		plane_index = 0;

		string path_delaunay = path_del + "_" + to_string(surface_index) + ".vtu";
		ofstream delaunays(path_delaunay);

		delaunays <<
			"<?xml version=\"1.0\"?>" << "\n" <<
			"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << "\n  " <<
			"<UnstructuredGrid>" << "\n    " <<
			"<Piece NumberOfPoints=\"" << dt.number_of_vertices() << "\" NumberOfCells=\"" << dt.number_of_faces() << "\">" << "\n      "
			"<PointData Scalars=\"scalars\">" << "\n        "
			"<DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">" << "\n          ";
		if (surface_index == 0) {
			for (auto it = pts.begin(); it != pts.end(); it++)
			{
				delaunays << fixed << (it->first.z()) << "\n          ";
			}
		}

		else {
			for (auto it = pts.begin(); it != pts.end(); it++)
			{
				delaunays << fixed << it->second.first.at(surface_index - 1) << "\n          ";
			}

		}
		delaunays << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</PointData>" << "\n      " <<
			"<Points>" << "\n        " <<
			"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n           ";

		if (surface_index == 0) {
			for (auto it = pts.begin(); it != pts.end(); it++)
			{
				double x = (it->first.y()), y = (it->first.x()), z = (it->first.z());
				delaunays << fixed << x << " " << y << " " << z << "\n           ";
			}
		}
		else {
			for (auto it = pts.begin(); it != pts.end(); it++)
			{
				double x = (it->first.y()), y = (it->first.x()), z = (it->second.first.at(surface_index - 1));
				delaunays << fixed << x << " " << y << " " << z << "\n           ";
			}
		}

		delaunays << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</Points>" << "\n      " <<
			"<Cells>" << "\n        " <<
			"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n          ";

		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit)
		{
			Delaunay::Face_handle face = fit;
			delaunays << face->vertex(0)->info().second - 1 << " " << face->vertex(1)->info().second - 1 << " " << face->vertex(2)->info().second - 1 << "\n          ";
		}

		delaunays << "\n        " <<
			"</DataArray>" << "\n        " <<
			"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n          ";

		for (auto i = 3; i <= 3 * dt.number_of_faces(); i = i + 3)
		{
			delaunays << i << " ";
		}

		delaunays << "\n        " <<
			"</DataArray>" << "\n        " <<
			"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<
			"\n          ";

		for (auto i = 1; i <= dt.number_of_faces(); i++)
		{
			delaunays << 5 << " ";
		}

		delaunays << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</Cells>" << "\n    " <<
			"</Piece>" << "\n  " <<
			"</UnstructuredGrid>" << "\n" <<
			"</VTKFile>";

		string path_normals = path_nor + "_" + to_string(surface_index) + ".vtu";

		ofstream normalvis(path_normals);

		normalvis <<
			"<?xml version=\"1.0\"?>" << "\n" <<
			"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << "\n  " <<
			"<UnstructuredGrid>" << "\n    " <<
			"<Piece NumberOfPoints=\"" << dt.number_of_faces() << "\" NumberOfCells=\"" << dt.number_of_faces() << "\">" << "\n      "

			"<PointData Scalars=\"scalars\">" << "\n        "
			"<DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">" << "\n          ";


		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) //a loop for performing the Delaunay triangulation and save the results

		{
			Delaunay::Face_handle face = fit;
			if (surface_index == 0) {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], dt.triangle(face)[0][2]), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], dt.triangle(face)[1][2]), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], dt.triangle(face)[2][2]), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}
			else {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(surface_index - 1)), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(surface_index - 1)), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(surface_index - 1)), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}

			}


			normalvis << fixed << planes.at(plane_index).center()[2] << "\n          ";

			plane_index++;
		}

		plane_index = 0;

		normalvis << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</PointData>" << "\n      " <<
			"<Points>" << "\n        " <<
			"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n           ";

		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) //a loop for performing the Delaunay triangulation and save the results

		{
			Delaunay::Face_handle face = fit;

			if (surface_index == 0) {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], dt.triangle(face)[0][2]), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], dt.triangle(face)[1][2]), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], dt.triangle(face)[2][2]), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}
			else {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(surface_index - 1)), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(surface_index - 1)), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(surface_index - 1)), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}

			}

		
			normalvis << fixed << planes.at(plane_index).center()[1] << " " << planes.at(plane_index).center()[0] << " " << planes.at(plane_index).center()[2] << "\n           ";

			plane_index++;

		}

		plane_index = 0;

		normalvis << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</Points>" << "\n      " <<
			"<Cells>" << "\n        " <<
			"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << "\n          ";

		for (auto i = 0; i < dt.number_of_faces(); i++)
		{
			normalvis << i << "\n          ";
		}

		normalvis << "\n        " <<
			"</DataArray>" << "\n        " <<
			"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << "\n          ";

		for (auto i = 1; i <= dt.number_of_faces(); i++)
		{
			normalvis << i << " ";
		}

		normalvis << "\n        " <<
			"</DataArray>" << "\n        " <<
			"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<
			"\n          ";

		for (auto i = 1; i <= dt.number_of_faces(); i++)
		{
			normalvis << 1 << " ";
		}

		normalvis << "\n        " <<
			"</DataArray>" << "\n      " <<
			"</Cells>" << "\n    " <<
			"<CellData Normals=\"cell_normals\">" << "\n      " <<
			"<DataArray type=\"Float32\" Name=\"cell_normals\" NumberOfComponents=\"3\" format=\"ascii\">" << "\n          ";


		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit)
		{
			Delaunay::Face_handle face = fit;

			if (surface_index == 0) {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], dt.triangle(face)[0][2]), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], dt.triangle(face)[1][2]), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], dt.triangle(face)[2][2]), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}
			else {
				point_1 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(surface_index - 1)), face->vertex(0)->info().second);
				point_2 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(surface_index - 1)), face->vertex(1)->info().second);
				point_3 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(surface_index - 1)), face->vertex(2)->info().second);

				try {
					planes.push_back(plane(point_1.first, point_2.first, point_3.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}

			}

			normalvis << fixed << planes.at(plane_index).get_normal()[1] << " " << planes.at(plane_index).get_normal()[0] << " " << planes.at(plane_index).get_normal()[2] << "\n          ";
			plane_index++;
		}

		plane_index = 0;

		normalvis <<
			"</DataArray>" << "\n      " <<
			"</CellData>" << "\n      " <<
			"<Cells>" << "\n        " <<
			"</Cells>" << "\n    " <<
			"</Piece>" << "\n  " <<
			"</UnstructuredGrid>" << "\n" <<
			"</VTKFile>";


	}


	std::system("pause");
	return 0;
}

