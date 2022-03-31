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


using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Triangulation_vertex_base_with_info_2< std::pair<vector<double>, unsigned int>, Gt > Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Gt, Tds> Delaunay;

typedef K::Point_3 Point;

const double coeff = 180 / M_PI; //a constant coefficient for transforming radians into degrees
const int n = 3; //we work in 3D

const double ex = 2; //we introduce the restriction of collinearity

std::vector<std::pair<int, int>> pairs(const int n_surfaces)
{
	std::vector<std::pair<int, int>> indicesf;
	const int g = 2;
	int S[g];
	for (int i = 0; i < g; i++) S[i] = i;

	int p = g - 1;

	while (p >= 0)
	{
		indicesf.push_back(std::make_pair(S[0], S[1]));

		if (S[g - 1] == n_surfaces - 1) { p--; }
		else { p = g - 1; }
		if (p >= 0)
		{
			for (int i = g - 1; i >= p; i--)
			{
				S[i] = S[p] + i - p + 1;
			}
		}

	}
	return indicesf;
}
double dot_product(double normal1[], double normal2[], int n = 3) //a function that computes the dot product of vectors
{
	double product = 0;
	for (int i = 0; i < n; i++)
	{
		product += normal1[i] * normal2[i];
	}
	return product;
}
double length(double line_vector[])//a function that computes the length of a vector
{
	double vector_length = sqrt(pow(line_vector[0], 2) + pow(line_vector[1], 2) + pow(line_vector[2], 2));
	return vector_length;
}
class plane //a class that stores the crucial figures for comparing the orientation
{
private:
	double first_vec[::n];			//the first edge of a triangle
	double second_vec[::n];			//the second edge of a triangle
	double third_vec[n];			//the third edge of a triangle
	double z_axis[::n] = { 0,0,1 };   //the definition of the z-axis
	Point thisp1;       //three points representing a surface
	Point thisp2;
	Point thisp3;
	double normal_vec[::n];  //the normal vector to a plane
	double dip_vec[::n];
	double doc;						//a variable that contains the collinearity coefficient
	bool lin_dependence;		    //a bool variable to check to answer whether points are (too) collinear
	string dip_degrees;             //a text variable to store the dip angle
	string azimuth_degrees;         //a text variable to store the dip direction
	double area;

public:
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

	double dip_angle(double z_axis[], double normal_v[]) //a function that computes the dip angle
	{
		double angle;
		double expression;
		double coeff = 180 / M_PI;
		expression = abs(dot_product(z_axis, normal_v)) / (length(z_axis)*length(normal_v));
		angle = acos(expression);
		return angle * coeff;
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
			azimuth_degrees = std::to_string(dip_azimuth(normal_vec));
			return dip_degrees + ";" + azimuth_degrees;
		}
		else
		{
			double dipping_angle = dip_angle(z_axis, normal_vec);
			dip_degrees = std::to_string(dipping_angle);
			azimuth_degrees = std::to_string(dip_azimuth(normal_vec));
			return dip_degrees + ";" + azimuth_degrees;
		}
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
		}

		double stala = 0.5; //using Heron formula to calculate area 
		double half = stala * (length(first_vec) + length(second_vec) + length(third_vec)); // Let p=0.5(a+b+c) be the half of the circumference.
		double s = sqrt(half*(half - length(first_vec))*(half - length(second_vec))*(half - length(third_vec))); // Then the area is sqrt(p*(p-a)(p-b)(p-c)))
		this->area = s;

	}

	bool dependence(double v1[], double v2[], double v3[]) //a function that checks whether the points are collinear
	{
		double len_v1 = length(v1);
		double len_v2 = length(v2);
		double len_v3 = length(v3);
		double lengths[::n] = { len_v1, len_v2, len_v3 };

		sort(lengths, lengths + ::n);
		this->doc = lengths[2] / (lengths[0] + lengths[1]);
		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (lengths[i] == 0)
			{
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

	vector<string> center() //a function that computes the geometrical centre of a Delaunay triangle
	{
		double x = (thisp1.x() + thisp2.x() + thisp3.x()) / (3.0);
		double y = (thisp1.y() + thisp2.y() + thisp3.y()) / (3.0);
		double z = (thisp1.z() + thisp2.z() + thisp3.z()) / (3.0);
		//double z = (thisp1[2] + thisp2[2] + thisp3[2]) / (3.0);
		vector<string> center_triangle{ std::to_string(x), std::to_string(y), std::to_string(z) };
		return center_triangle;

	}

	vector<double> get_normal() //a function that returns the normal vector
	{

		vector<double> normal_vector = { normal_vec[0] , normal_vec[1], normal_vec[2] };
		return normal_vector;
	}

	vector<double> get_dip_vec() {//a function that -computes- returns the dip vector

		vector<double> dip_vector = { dip_vec[0],dip_vec[1] ,dip_vec[2] };
		return(dip_vector);
	}

	double length(double line_vector[])//a function that computes the length of a vector
	{
		double vector_length = sqrt(pow(line_vector[0], 2) + pow(line_vector[1], 2) + pow(line_vector[2], 2));
		return vector_length;
	}

	double get_doc() {

		return(doc);
	}

	double get_area() {

		return(area);
	}

	vector<string> get_Points() {
		//vector<Point> points = { thisp1, thisp2, thisp3 };

		vector<string> coordinates{ std::to_string(thisp1.x()), std::to_string(thisp1.y()), std::to_string(thisp1.z()),
									std::to_string(thisp2.x()), std::to_string(thisp2.y()), std::to_string(thisp2.z()),
									std::to_string(thisp3.x()), std::to_string(thisp3.y()), std::to_string(thisp3.z())
		};
		return coordinates;

	}
};

double dist_ang_normals(plane plane1, plane plane2) //a function that calculates the angular distance between two vectors
{
	double angle;				 //the angle between vectors
	double expression;			 //the cos of the angle between two vectors

	double normal_1[::n]; //the normal vector of the first plane
	double normal_2[::n]; //the normal vector of the second plane

	for (int k = 0; k < ::n; k++)
	{
		normal_1[k] = plane1.get_normal()[k];
		normal_2[k] = plane2.get_normal()[k];
	}

	expression = abs(dot_product(normal_1, normal_2)) / (length(normal_1)*length(normal_2));
	angle = acos(expression);
	return coeff * angle;
}


double manhattan_normals(plane plane1, plane plane2)

{
	double manh_dist;
	double normal_1[::n]; //the normal vector of the first plane
	double normal_2[::n]; //the normal vector of the second plane

	for (int k = 0; k < ::n; k++)
	{
		normal_1[k] = plane1.get_normal()[k];
		normal_2[k] = plane2.get_normal()[k];
	}
	manh_dist = abs(normal_2[0] - normal_1[0]) + abs(normal_2[1] - normal_1[1]) + abs(normal_2[2] - normal_1[2]);
	return manh_dist;
}


double euclidean_normals(plane plane1, plane plane2)

{
	double euclid_dist;
	double normal_1[::n]; //the normal vector of the first plane
	double normal_2[::n]; //the normal vector of the second plane

	for (int k = 0; k < ::n; k++)
	{
		normal_1[k] = plane1.get_normal()[k];
		normal_2[k] = plane2.get_normal()[k];
	}
	euclid_dist = sqrt(pow((normal_2[0] - normal_1[0]), 2) + pow((normal_2[1] - normal_1[1]), 2) + pow((normal_2[2] - normal_1[2]), 2));
	return euclid_dist;
}

double squaredeuclidean_normals(plane plane1, plane plane2)

{
	double squared_euclid_dist;
	double normal_1[::n]; //the normal vector of the first plane
	double normal_2[::n]; //the normal vector of the second plane

	for (int k = 0; k < ::n; k++)
	{
		normal_1[k] = plane1.get_normal()[k];
		normal_2[k] = plane2.get_normal()[k];
	}
	squared_euclid_dist = pow((normal_2[0] - normal_1[0]), 2) + pow((normal_2[1] - normal_1[1]), 2) + pow((normal_2[2] - normal_1[2]), 2);
	return squared_euclid_dist;
}

double euclidean_dip(plane plane1, plane plane2)

{
	double euclid_dist;
	double dip_1[::n]; //the normal vector of the first plane
	double dip_2[::n]; //the normal vector of the second plane

	for (int k = 0; k < ::n; k++)
	{
		dip_1[k] = plane1.get_dip_vec()[k];
		dip_2[k] = plane2.get_dip_vec()[k];
	}
	euclid_dist = sqrt(pow((dip_2[0] - dip_1[0]), 2) + pow((dip_2[1] - dip_1[1]), 2) + pow((dip_2[2] - dip_1[2]), 2));
	return euclid_dist;
}

int main()
{
	string path_i; // path_i2; //text variables for the input path: the first and the second surface
	string path_o;			//text variable for the output path

	int n_surfaces;
	cout << "How many surfaces do you have?" << endl;
	cin >> n_surfaces;
	int combinations = n_surfaces * (n_surfaces - 1) / 2;
	cout << "You will get: " << combinations << " comparisons !" << endl;

	cout << "Type in the path of the input surfaces:" << endl; //the user is required to type in the input path of the surfaces storing in one file
	cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\Surfaces.txt" << endl << endl; //an example of introducing a path is given

	cin >> path_i;

	ifstream SurfacesDownload(path_i);
	if (!SurfacesDownload) cout << "Error in opening file" << endl; //the case when the file cannot be uploaded

	vector<std::pair<int, int>> indices; //indices for creating pairs
	if (n_surfaces == 2) {
		indices.push_back(make_pair(0, 1));
	}
	else {
		indices = pairs(n_surfaces);
	}
	cout << endl;
	cout << "Your indices: " << endl;

	for (auto k = 0; k < indices.size(); k++) { //displaying indices of surfaces
		cout << indices.at(k).first << "   " << indices.at(k).second << endl;
	}

	string lower; //a temporary variable storing figures of the first surface


	vector<double> surface_elevation; //vector storing elevations for different surfaces
	vector< std::pair<Point, std::pair<vector<double>, unsigned>> > pts; //pts is a vector of pairs: the first element is a point and the second is a pair with the first element being a vector of elevations and the second an index to the point

	vector <double> coord_tab; //a vector for storing all data from the table, its last element will be the index
	//int licznik=1;
	while (getline(SurfacesDownload, lower)) // loading points line - by - line associated with the first surface
	{
		istringstream convert(lower);

		unsigned int t;
		double data;
		for (int i = 0; i < n_surfaces + 3; i++) { //n_surfaces + 3 is the number of columns: we have n_surfaces and xy coordinates and the index
			if (!(convert >> data)) { break; }
			coord_tab.push_back(data);
		}

		for (auto i = 3; i < n_surfaces + 2; i++) { //we start loading the elevations from the fourth column (the columns 0,1,2 are in coord_tab[0],coord_tab[1],coord_tab[2])
			surface_elevation.push_back(coord_tab[i]); //elevations after first elevation
		}
		//	surface_elevation.push_back(coord_tab[i+(n_surfaces+3)*(licznik-1)]);
		pts.push_back(make_pair(Point(coord_tab[0], coord_tab[1], coord_tab[2]), make_pair(surface_elevation, coord_tab[coord_tab.size() - 1])));
		coord_tab.clear(); //avoiding problems with indexing in the next iteration
		surface_elevation.clear();
	}


	Delaunay dt; //a variable storing the geometrical elements of Delaunay triangulation regarding the first surface
	dt.insert(pts.begin(), pts.end());


	if (pts.size() != dt.number_of_vertices()) { throw(runtime_error("Check for duplicates in data!")); }

	cout << "Type in the path of the output:" << endl << endl; //the user is required to type in the path of the output

	cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\Surfaces_output" << endl << endl;

	cout << "DO NOT PROVIDE .TXT EXTENSION. IT WILL BE ADDED ALONG WITH THE INDICES OF THE SURFACES" << endl << endl;

	cin >> path_o;
	//a stream variable to save output figures

	for (auto k = 0; k != indices.size(); ++k) { //size of indices is the number of pairs (comparisons)
		int ind1 = indices.at(k).first; //index of the first surface for the kth iteration
		int ind2 = indices.at(k).second; //index of the second surface for the kth iteration
		string path_output = path_o + "_" + to_string(ind1) + to_string(ind2) + ".txt";
		ofstream save(path_output);
		save << "Index of the first surface:" << ind1 << "Index of the first surface:" << ind2 << endl;
		save << "X11;" << "Y11;" << "Z11;" << "X12;" << "Y12;" << "Z12;" << "X13;" << "Y13;" << "Z13;" << "X1_C;" //column names are saved
			<< "Y1_C;" << "Z1_C;" << "X1_N;" << "Y1_N;" << "Z1_N;" << "X1_D;" << "Y1_D;" << "Z1_D;" << "Dip_ang1;" << "Dip_dir1;" << "DOC1;" << "Area1;" << "IDT11;" << "IDT12;" << "IDT13;" <<
			"X21;" << "Y21;" << "Z21;" << "X22;" << "Y22;" << "Z22;" << "X23;" << "Y23;" << "Z23;" << "X2_C;" //column names are saved
			<< "Y2_C;" << "Z2_C;" << "X2_N;" << "Y2_N;" << "Z2_N;" << "X2_D;" << "Y2_D;" << "Z2_D;" << "Dip_ang2;" << "Dip_dir2;" << "DOC2;" << "Area2;" << "IDT21;" << "IDT22;" << "IDT23;"
			<< "AngularDistance_normals;" << "ManhattanDistance_normals;" << "EuclideanDistance_normals;" << "SquaredEuclidean_normals;" << "EuclideanDistance_dip" << endl; //column names are saved


		vector <plane> planes1; //a vector variable storing planes representing the first surface
		auto plane_index = 0; //for using in "at" function to access elements of a vector containing planes
		vector <plane> planes2; //a vector variable storing planes representing the second surface
		std::pair<Point, int> point_11; //we have a pair with point and index of the point
		std::pair<Point, int> point_12;
		std::pair<Point, int> point_13;
		std::pair<Point, int> point_21;
		std::pair<Point, int> point_22;
		std::pair<Point, int> point_23;
		for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) //a loop for performing the Delaunay triangulation and extracting planes representing the first surface
		{
			Delaunay::Face_handle face = fit;



			if (ind1 == 0) {
				point_11 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], dt.triangle(face)[0][2]), face->vertex(0)->info().second);
				point_12 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], dt.triangle(face)[1][2]), face->vertex(1)->info().second);
				point_13 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], dt.triangle(face)[2][2]), face->vertex(2)->info().second);

				point_21 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(ind2 - 1)), face->vertex(0)->info().second); // why (ind2-1)? For example, if the index of the second surface is 1, then it is stored as the first element of the surface_elevation vector (thus the 0th) element in this vector 
				point_22 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(ind2 - 1)), face->vertex(1)->info().second);
				point_23 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(ind2 - 1)), face->vertex(2)->info().second);
				try {
					planes1.push_back(plane(point_11.first, point_12.first, point_13.first));
					planes2.push_back(plane(point_21.first, point_22.first, point_23.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}

			else {

				point_11 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(ind1 - 1)), face->vertex(0)->info().second);
				point_12 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(ind1 - 1)), face->vertex(1)->info().second);
				point_13 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(ind1 - 1)), face->vertex(2)->info().second);
				point_21 = make_pair(Point(dt.triangle(face)[0][0], dt.triangle(face)[0][1], face->vertex(0)->info().first.at(ind2 - 1)), face->vertex(0)->info().second);
				point_22 = make_pair(Point(dt.triangle(face)[1][0], dt.triangle(face)[1][1], face->vertex(1)->info().first.at(ind2 - 1)), face->vertex(1)->info().second);
				point_23 = make_pair(Point(dt.triangle(face)[2][0], dt.triangle(face)[2][1], face->vertex(2)->info().first.at(ind2 - 1)), face->vertex(2)->info().second);
				try {
					planes1.push_back(plane(point_11.first, point_12.first, point_13.first));
					planes2.push_back(plane(point_21.first, point_22.first, point_23.first));
				}
				catch (exception e) {

					cout << e.what() << endl;
				}
			}

			vector<string> centroid1 = planes1.at(plane_index).center();     //extracting the centre of a Delaunay triangle
			vector<string> centroid2 = planes2.at(plane_index).center();     //extracting the centre of a Delaunay triangle

			string result1 = planes1.at(plane_index).measure();
			string result2 = planes2.at(plane_index).measure();								//extracting the dip angle and the dip direction

			double ang_dist_normals = dist_ang_normals(planes1[plane_index], planes2[plane_index]);
			double manh_dist_normals = manhattan_normals(planes1[plane_index], planes2[plane_index]);
			double euclid_dist_normals = euclidean_normals(planes1[plane_index], planes2[plane_index]);
			double sqeuclid_dist_normals = squaredeuclidean_normals(planes1[plane_index], planes2[plane_index]);
			double euclid_dist_dip = euclidean_dip(planes1[plane_index], planes2[plane_index]);

			save <<
				planes1.at(plane_index).get_Points()[0] << ";" <<//coordinates start
				planes1.at(plane_index).get_Points()[1] << ";" <<
				planes1.at(plane_index).get_Points()[2] << ";" <<
				planes1.at(plane_index).get_Points()[3] << ";" <<
				planes1.at(plane_index).get_Points()[4] << ";" <<
				planes1.at(plane_index).get_Points()[5] << ";" <<
				planes1.at(plane_index).get_Points()[6] << ";" <<
				planes1.at(plane_index).get_Points()[7] << ";" <<
				planes1.at(plane_index).get_Points()[8] << ";" << //coordinates end
				centroid1[0] << ";" << //centre start
				centroid1[1] << ";" <<
				centroid1[2] << ";" << //centre end
				planes1.at(plane_index).get_normal()[0] << ";" <<  //normal vector start
				planes1.at(plane_index).get_normal()[1] << ";" <<
				planes1.at(plane_index).get_normal()[2] << ";" << //normal vector end
				planes1.at(plane_index).get_dip_vec()[0] << ";" << //dip vector start
				planes1.at(plane_index).get_dip_vec()[1] << ";" <<
				planes1.at(plane_index).get_dip_vec()[2] << ";" << //dip vector end
				result1 << ";" << //dip and dip direction
				planes1.at(plane_index).get_doc() << ";" <<  //collinearity
				planes1.at(plane_index).get_area() << ";" << //area
				face->vertex(0)->info().second << ";" <<
				face->vertex(1)->info().second << ";" <<
				face->vertex(2)->info().second << ";" <<
				planes2.at(plane_index).get_Points()[0] << ";" <<//coordinates start
				planes2.at(plane_index).get_Points()[1] << ";" <<
				planes2.at(plane_index).get_Points()[2] << ";" <<
				planes2.at(plane_index).get_Points()[3] << ";" <<
				planes2.at(plane_index).get_Points()[4] << ";" <<
				planes2.at(plane_index).get_Points()[5] << ";" <<
				planes2.at(plane_index).get_Points()[6] << ";" <<
				planes2.at(plane_index).get_Points()[7] << ";" <<
				planes2.at(plane_index).get_Points()[8] << ";" << //coordinates end
				centroid2[0] << ";" << //centre start
				centroid2[1] << ";" <<
				centroid2[2] << ";" << //centre end
				planes2.at(plane_index).get_normal()[0] << ";" <<  //normal vector start
				planes2.at(plane_index).get_normal()[1] << ";" <<
				planes2.at(plane_index).get_normal()[2] << ";" << //normal vector end
				planes2.at(plane_index).get_dip_vec()[0] << ";" << //dip vector start
				planes2.at(plane_index).get_dip_vec()[1] << ";" <<
				planes2.at(plane_index).get_dip_vec()[2] << ";" << //dip vector end
				result2 << ";" << //dip and dip direction
				planes2.at(plane_index).get_doc() << ";" << //collinearity
				planes2.at(plane_index).get_area() << ";" << //area
				face->vertex(0)->info().second << ";" << //id start
				face->vertex(1)->info().second << ";" <<
				face->vertex(2)->info().second << ";" << //id end
				ang_dist_normals << ";" << //distances start
				manh_dist_normals << ";" <<
				euclid_dist_normals << ";" <<
				sqeuclid_dist_normals << ";" <<
				euclid_dist_dip << endl; //distances end

			plane_index++; //incrementing index;

		}

	}
	system("pause");
	return 0;
}
