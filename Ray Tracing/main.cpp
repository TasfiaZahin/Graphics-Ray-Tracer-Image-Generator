#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))

double scale = 2.0;

int drawgrid;
int drawaxes;

int flag=0;

double fov_x;
double fov_y;
double aspectRatio;
double _near;
double _far;
int level_of_recursion;
int screen_size;

int number_of_objects;
int number_of_normal_lights;
int number_of_spot_lights;
int produceimage;

class Vector
{
public:
    double x, y, z;
    Vector()
    {

    }

    // constructs a vector with given components
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};

class Point
{
public:
    double x;
    double y;
    double z;
    Point()
    {

    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point addVector(Vector v)
    {
        double x = this->x + v.x;
        double y = this->y + v.y;
        double z = this->z + v.z;
        Point p(x,y,z);
        return p;
    }

    Point subVector(Vector v)
    {
        double x = this->x - v.x;
        double y = this->y - v.y;
        double z = this->z - v.z;
        Point p(x,y,z);
        return p;
    }

    Point addScalar(double s)
    {
        double x = this->x + s;
        double y = this->y + s;
        double z = this->z + s;
        Point p(x,y,z);
        return p;
    }

    // adds two points
    Point operator+ (const Point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        Point p(x, y, z);
        return p;
    }

    // subtracts one point from another
    Point operator- (const Point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        Point p(x, y, z);
        return p;
    }

    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << endl;
    }

};

class Triangle
{
public:
    Point p1,p2,p3,other_one,other_two;

    Triangle(Point p1, Point p2, Point p3, Point other_one, Point other_two)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        this->other_one = other_one;
        this->other_two = other_two;
    }

    void print()
    {
        p1.print();
        p2.print();
        p3.print();
        cout << "\n";
    }
};

class Color {
public:
    double r, g, b;
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color() {
    }
};

class PolishComponent
{
public:
    double lambert;
    double phong;
    Color ref_color;
    PolishComponent()
    {

    }
    PolishComponent(double lambert, double phong)
    {
        this->lambert = lambert;
        this->phong = phong;
    }
};

struct point
{
	double x,y,z;
};

struct vect
{
	double x,y,z;
};

class CheckerBoard
{
public:
    double width;
    Color color;
    double ambient_coeff;
    double diffuse_coeff;
    double reflection_coeff;
    CheckerBoard()
    {

    }
};

class Pyramid
{
public:
    struct point lowest_pt;
    double width;
    double height;
    Color color;
    double ambient_coeff;
    double diffuse_coeff;
    double specular_coeff;
    double reflection_coeff;
    double shininess;
    Pyramid()
    {

    }
};

class Sphere
{
public:
    struct point center;
    double radius;
    Color color;
    double ambient_coeff;
    double diffuse_coeff;
    double specular_coeff;
    double reflection_coeff;
    double shininess;
    Sphere()
    {

    }

};

class SpherePair
{
public:
    Sphere sphere;
    double t;
    Vector normal;
    SpherePair()
    {
        t = -100.0;
    }
    SpherePair(Sphere s, double t, Vector normal)
    {
        this->sphere = s;
        this->t = t;
        this->normal = normal;
    }
};

class PyramidPair
{
public:
    Pyramid pyramid;
    double t;
    Vector normal;
    PyramidPair()
    {
        t = -100.0;
    }
    PyramidPair(Pyramid p, double t, Vector normal)
    {
        this->pyramid = p;
        this->t = t;
        this->normal = normal;
    }
};

class Index_t_pair
{
public:
    int index;
    double t;
    Vector normal;
    double shininess;
    Index_t_pair()
    {

    }
    Index_t_pair(int index, double t, Vector normal)
    {
        this->index = index;
        this->t = t;
        this->normal = normal;
    }
};

class NormalLightSource
{
public:
    Point position;
    double falloff_parameter;
    NormalLightSource()
    {

    }
};

class SpotLightSource
{
public:
    Point position;
    double falloff_parameter;
    struct point look_point;
    double cutoff_angle;
    SpotLightSource()
    {

    }
};

class PointBufferCell
{
public:
    Point p;
    Color color;
    PointBufferCell()
    {

    }
    PointBufferCell(Point p)
    {
        this->p = p;
    }
};

Color getObjectColor(Point R0, Vector ray);

vector<Sphere> spheres;
vector<Pyramid> pyramids;
CheckerBoard checkerBoard;

vector<NormalLightSource> normal_light_sources;
vector<SpotLightSource> spotlight_sources;

vect U; //up vector
vect R; //right vector
vect L; //look vector
point pos;
double rotate_angle;

bitmap_image b_img ("ece1.bmp");
Color **textureBuffer;
int height, width;
int texture_on;

void getTextureBuffer()
{
    if(produceimage == 1 && texture_on == 1)
    {
        height = b_img.height();
        width = b_img.width();
        //cout << width << " " << height << " ";
        textureBuffer = new Color* [width];
        for (int i = 0; i < width; i++) {
            textureBuffer[i] = new Color [height];
            for (int j = 0; j < height; j++) {
                unsigned char r, g, b;
                b_img.get_pixel(i, j, r, g, b);
                Color c(r/255.0, g/255.0, b/255.0);
                textureBuffer[i][j] = c;
            }
        }
    }
}

void drawAxes()
{
	if(drawaxes==1)
	{
		//glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 300,0,0);
			glVertex3f(-300,0,0);

			glVertex3f(0,-300,0);
			glVertex3f(0, 300,0);

			glVertex3f(0,0, 300);
			glVertex3f(0,0,-300);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-9;i<=9;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -100, 0);
				glVertex3f(i*10,  100, 0);

				//lines parallel to X-axis
				glVertex3f(-100, i*10, 0);
				glVertex3f( 100, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double x,double y, double z, double a) //lower left point and width
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f(x,y,z);
		glVertex3f(x + a,y,z);
		glVertex3f(x + a,y + a,z);
		glVertex3f(x, y + a,z);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments; //what is done here :3
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double center_x, double center_y, double center_z, double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=center_x + r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=center_y + r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{

        double col = 1;

		for(j=0;j<slices;j++)
		{
		    //glColor3f(col,col,col);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,center_z+points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,center_z+points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,center_z+points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,center_z+points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,center_z-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,center_z-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,center_z-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,center_z-points[i+1][j].z);
			}glEnd();

			col = 1 - col;
		}
	}
}

void drawPyramid(double x, double y, double z, double width, double height)
{
    drawSquare(x,y,z,width);
    double mp_x = x + width/2.0;
    double mp_y = y + width/2.0;

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(mp_x,mp_y,z + height);
        glVertex3f(x,y,z);
        glVertex3f(x + width,y,z);
    }
    glEnd();

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(mp_x,mp_y,z + height);
        glVertex3f(x + width,y,z);
        glVertex3f(x + width,y + width,z);
    }
    glEnd();

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(mp_x,mp_y,z + height);
        glVertex3f(x + width,y + width,z);
        glVertex3f(x,y + width,z);
    }
    glEnd();

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(mp_x,mp_y,z + height);
        glVertex3f(x,y + width,z);
        glVertex3f(x,y,z);
    }
    glEnd();

}

double get_magnitude(struct vect p)
{
    double ans = (p.x*p.x)+(p.y*p.y)+(p.z*p.z);
    return sqrt(ans);
}

struct vect get_cross_product(struct vect p, struct vect q)
{
    struct vect ans;
    ans.x = (p.y*q.z) - (p.z*q.y);
    ans.y = -((p.x*q.z) - (p.z*q.x));
    ans.z = (p.x*q.y) - (p.y*q.x);

    return ans;
}

void rotate_UL_anticlockwise()
{

    struct vect new_L, new_U;

    double radian = (rotate_angle/360)*2*pi;
    //printf("radian val: %f",radian);

    struct vect N = get_cross_product(R,L);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    L.x = L.x*cos(radian);
    L.y = L.y*cos(radian);
    L.z = L.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_L.x = L.x + N.x;
    new_L.y = L.y + N.y;
    new_L.z = L.z + N.z;

    N = get_cross_product(R,U);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    U.x = U.x*cos(radian);
    U.y = U.y*cos(radian);
    U.z = U.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_U.x = U.x + N.x;
    new_U.y = U.y + N.y;
    new_U.z = U.z + N.z;

    double mag_L = get_magnitude(L);
    double mag_U = get_magnitude(U);

    new_L.x = new_L.x/mag_L;
    new_L.y = new_L.y/mag_L;
    new_L.z = new_L.z/mag_L;

    new_U.x = new_U.x/mag_U;
    new_U.y = new_U.y/mag_U;
    new_U.z = new_U.z/mag_U;

    L.x = new_L.x;
    L.y = new_L.y;
    L.z = new_L.z;

    U.x = new_U.x;
    U.y = new_U.y;
    U.z = new_U.z;

}

void rotate_UL_clockwise()
{

    struct vect new_L, new_U;

    double radian = (rotate_angle/360)*2*pi;

    struct vect N = get_cross_product(L,R);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    L.x = L.x*cos(radian);
    L.y = L.y*cos(radian);
    L.z = L.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_L.x = L.x + N.x;
    new_L.y = L.y + N.y;
    new_L.z = L.z + N.z;



    N = get_cross_product(U,R);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    U.x = U.x*cos(radian);
    U.y = U.y*cos(radian);
    U.z = U.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_U.x = U.x + N.x;
    new_U.y = U.y + N.y;
    new_U.z = U.z + N.z;

    double mag_L = get_magnitude(L);
    double mag_U = get_magnitude(U);

    new_L.x = new_L.x/mag_L;
    new_L.y = new_L.y/mag_L;
    new_L.z = new_L.z/mag_L;

    new_U.x = new_U.x/mag_U;
    new_U.y = new_U.y/mag_U;
    new_U.z = new_U.z/mag_U;

    L.x = new_L.x;
    L.y = new_L.y;
    L.z = new_L.z;

    U.x = new_U.x;
    U.y = new_U.y;
    U.z = new_U.z;

}

void rotate_UR_anticlockwise()
{

    struct vect new_R, new_U;

    double radian = (rotate_angle/360)*2*pi;

    struct vect N = get_cross_product(U,L);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    U.x = U.x*cos(radian);
    U.y = U.y*cos(radian);
    U.z = U.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_U.x = U.x + N.x;
    new_U.y = U.y + N.y;
    new_U.z = U.z + N.z;


    N = get_cross_product(R,L);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    R.x = R.x*cos(radian);
    R.y = R.y*cos(radian);
    R.z = R.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_R.x = R.x + N.x;
    new_R.y = R.y + N.y;
    new_R.z = R.z + N.z;

    double mag_U = get_magnitude(U);
    double mag_R = get_magnitude(R);

    new_U.x = new_U.x/mag_U;
    new_U.y = new_U.y/mag_U;
    new_U.z = new_U.z/mag_U;

    new_R.x = new_R.x/mag_R;
    new_R.y = new_R.y/mag_R;
    new_R.z = new_R.z/mag_R;

    U.x = new_U.x;
    U.y = new_U.y;
    U.z = new_U.z;

    R.x = new_R.x;
    R.y = new_R.y;
    R.z = new_R.z;

}

void rotate_UR_clockwise()
{

    struct vect new_R, new_U;

    double radian = (rotate_angle/360)*2*pi;

    struct vect N = get_cross_product(L,U);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    U.x = U.x*cos(radian);
    U.y = U.y*cos(radian);
    U.z = U.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_U.x = U.x + N.x;
    new_U.y = U.y + N.y;
    new_U.z = U.z + N.z;


    N = get_cross_product(L,R);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    R.x = R.x*cos(radian);
    R.y = R.y*cos(radian);
    R.z = R.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_R.x = R.x + N.x;
    new_R.y = R.y + N.y;
    new_R.z = R.z + N.z;

    double mag_U = get_magnitude(U);
    double mag_R = get_magnitude(R);

    new_U.x = new_U.x/mag_U;
    new_U.y = new_U.y/mag_U;
    new_U.z = new_U.z/mag_U;

    new_R.x = new_R.x/mag_R;
    new_R.y = new_R.y/mag_R;
    new_R.z = new_R.z/mag_R;

    U.x = new_U.x;
    U.y = new_U.y;
    U.z = new_U.z;

    R.x = new_R.x;
    R.y = new_R.y;
    R.z = new_R.z;

}

void rotate_RL_anticlockwise()
{

    struct vect new_R, new_L;

    double radian = (rotate_angle/360)*2*pi;

    struct vect N = get_cross_product(U,L);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    L.x = L.x*cos(radian);
    L.y = L.y*cos(radian);
    L.z = L.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_L.x = L.x + N.x;
    new_L.y = L.y + N.y;
    new_L.z = L.z + N.z;


    N = get_cross_product(U,R);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    R.x = R.x*cos(radian);
    R.y = R.y*cos(radian);
    R.z = R.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_R.x = R.x + N.x;
    new_R.y = R.y + N.y;
    new_R.z = R.z + N.z;

    double mag_L = get_magnitude(L);
    double mag_R = get_magnitude(R);

    new_L.x = new_L.x/mag_L;
    new_L.y = new_L.y/mag_L;
    new_L.z = new_L.z/mag_L;

    new_R.x = new_R.x/mag_R;
    new_R.y = new_R.y/mag_R;
    new_R.z = new_R.z/mag_R;

    L.x = new_L.x;
    L.y = new_L.y;
    L.z = new_L.z;

    R.x = new_R.x;
    R.y = new_R.y;
    R.z = new_R.z;

}

void rotate_RL_clockwise()
{

    struct vect new_R, new_L;

    double radian = (rotate_angle/360)*2*pi;

    struct vect N = get_cross_product(L,U);

    double mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    L.x = L.x*cos(radian);
    L.y = L.y*cos(radian);
    L.z = L.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_L.x = L.x + N.x;
    new_L.y = L.y + N.y;
    new_L.z = L.z + N.z;

    N = get_cross_product(R,U);

    mag_N = get_magnitude(N);

    N.x = N.x/mag_N;
    N.y = N.y/mag_N;
    N.z = N.z/mag_N;

    R.x = R.x*cos(radian);
    R.y = R.y*cos(radian);
    R.z = R.z*cos(radian);

    N.x = N.x*sin(radian);
    N.y = N.y*sin(radian);
    N.z = N.z*sin(radian);

    new_R.x = R.x + N.x;
    new_R.y = R.y + N.y;
    new_R.z = R.z + N.z;

    double mag_L = get_magnitude(L);
    double mag_R = get_magnitude(R);

    new_L.x = new_L.x/mag_L;
    new_L.y = new_L.y/mag_L;
    new_L.z = new_L.z/mag_L;

    new_R.x = new_R.x/mag_R;
    new_R.y = new_R.y/mag_R;
    new_R.z = new_R.z/mag_R;

    L.x = new_L.x;
    L.y = new_L.y;
    L.z = new_L.z;

    R.x = new_R.x;
    R.y = new_R.y;
    R.z = new_R.z;

}

void drawCheckerBoard()
{
    int x,y,p=0,q=0;
    double width = checkerBoard.width;

    for(x=-2*_far;x<2*_far;x+=width)
    {
         for(y=-2*_far;y<2*_far;y+=width)
         {
             //printf("%d %d",x,y);
             if((p+q)%2==0)
             {
                 glColor3f(1, 1, 1);
             }
             else
             {
                 glColor3f(0, 0, 0);
             }
             drawSquare(x, y, 0, width);
             q++;
         }
         p++;
    }
}

void drawObjects()
{
    if(texture_on == 0)
    {
        drawCheckerBoard();
    }

    for(int i=0;i<spheres.size();i++)
    {
        Sphere sp = spheres[i];
        glColor3f(sp.color.r, sp.color.g, sp.color.b);
        drawSphere(sp.center.x, sp.center.y, sp.center.z, sp.radius, 30, 30);
    }
    for(int i=0;i<pyramids.size();i++)
    {
        Pyramid py = pyramids[i];
        glColor3f(py.color.r, py.color.g, py.color.b);
        drawPyramid(py.lowest_pt.x, py.lowest_pt.y, py.lowest_pt.z, py.width, py.height);
    }
    for(int i=0;i<normal_light_sources.size();i++)
    {
        NormalLightSource src = normal_light_sources[i];
        glColor3f(1,1,1);
        drawSphere(src.position.x, src.position.y, src.position.z, 5, 30, 30);
    }
    for(int i=0;i<spotlight_sources.size();i++)
    {
        SpotLightSource src = spotlight_sources[i];
        glColor3f(1,1,1);
        drawSphere(src.position.x, src.position.y, src.position.z, 5, 30, 30);
    }
}

SpherePair get_closest_sphere(Point buffer_cell_point, Vector ray)
{
    //for spheres array get closes t
    SpherePair minPair;
    vector<SpherePair> t_from_all_sps;

    for(int i=0;i<spheres.size();i++)
    {
        vector<SpherePair> valids;

        Sphere org_sp = spheres[i];
        Sphere sp = spheres[i];

        //Vector new_pos = Vector(org_pos.x - sp.center.x, org_pos.y - sp.center.y, org_pos.z - sp.center.z);
        Point new_cell_pt = Point(buffer_cell_point.x - sp.center.x, buffer_cell_point.y - sp.center.y, buffer_cell_point.z - sp.center.z);
        Vector R0 = Vector(new_cell_pt.x,new_cell_pt.y,new_cell_pt.z);
        Vector Rd = ray;

        //Point temp = new_cell_pt.subVector(new_pos);
        //Vector R0 = Vector(new_cell_pt.x,new_cell_pt.y,new_cell_pt.z); //point buffer point
        //Vector Rd = Vector(temp.x,temp.y,temp.z);
        //Rd.normalize();

        sp.center.x -= sp.center.x;
        sp.center.y -= sp.center.y;
        sp.center.z -= sp.center.z;

        double a = 1.0;
        double b = 2.0*Vector::dot(R0,Rd);
        double c = Vector::dot(R0,R0) - (sp.radius*sp.radius);
        double d = b*b - 4.0*a*c;
        //cout << d << " ";

        if(d < 0)
        {
            continue;
        }
        else if(d > 0)
        {
            //cout << d << " ";
            double t1 = (-b + sqrt(d))/2.0*a;
            double t2 = (-b - sqrt(d))/2.0*a;
            //cout << t1 << " " << t2 << " ";
            if(t1 > 0.0 && t1 < (_far - _near))
            {
                //cout << "yo" << " ";
                Vector P = R0 + Rd*t1; // intersecting point to find normal
                Vector normal = P - Vector(sp.center.x, sp.center.y, sp.center.z);
                normal.normalize();
                valids.push_back(SpherePair(org_sp,t1,normal));
            }
            if(t2 > 0.0 && t2 < (_far - _near))
            {
                //cout << "yo" << " ";
                Vector P = R0 + Rd*t2; // intersecting point to find normal
                Vector normal = P - Vector(sp.center.x, sp.center.y, sp.center.z);
                normal.normalize();
                valids.push_back(SpherePair(org_sp,t2,normal));
            }

            if(valids.size() >= 1)
            {
                //cout << valids.size() << " ";
                for(int m=0;m<valids.size()-1;m++)
                {
                    for(int n=m+1;n<valids.size();n++)
                    {
                        if(valids[n].t < valids[m].t)
                        {
                            SpherePair temp = valids[m];
                            valids[m] = valids[n];
                            valids[n] = temp;
                        }
                    }
                }
                t_from_all_sps.push_back(valids[0]);
            }
        }
        else
        {
            //cout << "hi" << " ";
            double t = -b/2.0*a;
            //cout << t << " ";
            if(t > 0 && t < (_near - _far))
            {
                Vector P = R0 + Rd*t; // intersecting point to find normal
                Vector normal = P - Vector(sp.center.x, sp.center.y, sp.center.z);
                normal.normalize();
                t_from_all_sps.push_back(SpherePair(org_sp,t,normal));
            }
        }
    }

    if(t_from_all_sps.size() >= 1)
    {
        //cout << t_from_all_sps.size() << " ";
        for(int m=0;m<t_from_all_sps.size()-1;m++)
        {
            for(int n=m+1;n<t_from_all_sps.size();n++)
            {
                if(t_from_all_sps[n].t < t_from_all_sps[m].t)
                {
                    SpherePair temp = t_from_all_sps[m];
                    t_from_all_sps[m] = t_from_all_sps[n];
                    t_from_all_sps[n] = temp;
                }
            }
        }
        minPair = t_from_all_sps[0];
    }
    else
    {
        minPair = SpherePair();
    }

    /*if(minPair.t != -100.0)
    {
        cout << minPair.t << " ";
    }*/

    return minPair;
}

double get_determinant(double arr[3][3])
{
    double ans = arr[0][0]*(arr[1][1]*arr[2][2]-arr[1][2]*arr[2][1]) - arr[0][1]*(arr[1][0]*arr[2][2]-arr[1][2]*arr[2][0]) + arr[0][2]*(arr[1][0]*arr[2][1]-arr[1][1]*arr[2][0]);
    //cout << ans << " ";
    return ans;
}

PyramidPair get_closest_pyramid(Vector R0, Vector Rd)
{
    //for pyramids array get closes t

    PyramidPair minPair;
    vector<PyramidPair> t_from_all_pys;

    for(int i=0;i<pyramids.size();i++)
    {
        vector<PyramidPair> valids;

        Pyramid py = pyramids[i];
        double mp_x = py.lowest_pt.x + py.width/2.0;
        double mp_y = py.lowest_pt.y + py.width/2.0;

        Point ll = Point(py.lowest_pt.x,py.lowest_pt.y,py.lowest_pt.z);
        Point lr = Point(py.lowest_pt.x + py.width,py.lowest_pt.y,py.lowest_pt.z);
        Point tr = Point(py.lowest_pt.x + py.width,py.lowest_pt.y + py.width,py.lowest_pt.z);
        Point tl = Point(py.lowest_pt.x,py.lowest_pt.y + py.width,py.lowest_pt.z);
        Point mp = Point(mp_x,mp_y,py.lowest_pt.z + py.height);

        //mp.print();

        vector<Triangle>triangles;
        triangles.push_back(Triangle(mp,ll,lr,tl,tr));
        triangles.push_back(Triangle(mp,lr,tr,ll,tl));
        triangles.push_back(Triangle(mp,tr,tl,ll,lr));
        triangles.push_back(Triangle(mp,tl,ll,lr,tr));


        //cout << triangles.size();
        for(int j=0;j<triangles.size();j++)
        {
            Triangle tr = triangles[j];

            /*Vector a = Vector(tr.p1.x, tr.p1.y, tr.p1.z);
            Vector b = Vector(tr.p2.x, tr.p2.y, tr.p2.z);
            Vector c = Vector(tr.p3.x, tr.p3.y, tr.p3.z);

            Vector n = Vector::cross(b-a,c-a);
            double d1 = -Vector::dot(n,a);
            double d2 = Vector::dot(n,Vector(tr.other_one.x,tr.other_one.y,tr.other_one.z));
            double d = d1 + d2;
            if(d > 0)
            {
                n = n*(-1);
            }
            n.normalize();*/


            double A[3][3] = {
                                    {tr.p1.x-tr.p2.x, tr.p1.x-tr.p3.x, Rd.x},
                                    {tr.p1.y-tr.p2.y, tr.p1.y-tr.p3.y, Rd.y},
                                    {tr.p1.z-tr.p2.z, tr.p1.z-tr.p3.z, Rd.z},
                                };
            double beta_A[3][3] = {
                                    {tr.p1.x-R0.x, tr.p1.x-tr.p3.x, Rd.x},
                                    {tr.p1.y-R0.y, tr.p1.y-tr.p3.y, Rd.y},
                                    {tr.p1.z-R0.z, tr.p1.z-tr.p3.z, Rd.z},
                                };

            double gamma_A[3][3] = {
                                    {tr.p1.x-tr.p2.x, tr.p1.x-R0.x, Rd.x},
                                    {tr.p1.y-tr.p2.y, tr.p1.y-R0.y, Rd.y},
                                    {tr.p1.z-tr.p2.z, tr.p1.z-R0.z, Rd.z},
                                };
            double t_A[3][3] = {
                                    {tr.p1.x-tr.p2.x, tr.p1.x-tr.p3.x, tr.p1.x-R0.x},
                                    {tr.p1.y-tr.p2.y, tr.p1.y-tr.p3.y, tr.p1.y-R0.y},
                                    {tr.p1.z-tr.p2.z, tr.p1.z-tr.p3.z, tr.p1.z-R0.z},
                                };

            double det_A = get_determinant(A);
            double det_beta_A = get_determinant(beta_A);
            double det_gamma_A = get_determinant(gamma_A);
            double det_t_A = get_determinant(t_A);

            double beta = det_beta_A/det_A;
            double gamma = det_gamma_A/det_A;
            double t = det_t_A/det_A;

            /*if(k < 5)
            {
                cout<<det_A<<" "<< det_beta_A<<" "<<det_gamma_A<<" "<<det_t_A<<" ";
                k++;
                cout << "\n";
            }*/

            /*if(k < 5)
            {
                cout<<beta<<" "<<gamma<<" "<<t<<" ";
                k++;
            }*/

            //cout << beta + gamma << " ";
            if((beta+gamma) < 1.0 && beta > 0.0 && gamma > 0.0 && t > 0 && t < (_far - _near))
            {
                Vector a = Vector(tr.p1.x, tr.p1.y, tr.p1.z);
                Vector b = Vector(tr.p2.x, tr.p2.y, tr.p2.z);
                Vector c = Vector(tr.p3.x, tr.p3.y, tr.p3.z);
                Vector normal = Vector::cross(c-b,a-b);
                normal.normalize();

                //cout << "yess" << " ";

                PyramidPair p_pair = PyramidPair(py,t,normal);
                valids.push_back(p_pair);
            }

            /*if(valids.size() != 0)
            {
                //cout << valids.size() << " ";
                //cout << valids[0].t << " ";
            }*/
        }

        // check base intersection with ray
        Vector point_on_plane = Vector(mp_x,mp_y,py.lowest_pt.z);
        //point_on_plane.print();
        Vector plane_outer_normal = Vector(0,0,-1);
        plane_outer_normal.normalize();
        //D = -n.P
        double D = -Vector::dot(plane_outer_normal,point_on_plane);
        //cout << D << " ";
        double intersect_t = -(D + Vector::dot(plane_outer_normal,R0))/Vector::dot(plane_outer_normal,Rd);

        if(intersect_t > 0.0 && intersect_t < (_far - _near))
        {
            //cout << intersect_t << " ";
            Vector intersect_point = R0 + Rd*intersect_t;
            if(intersect_point.x >= ll.x && intersect_point.x <= lr.x && intersect_point.y >= ll.y && intersect_point.y <= tl.y)
            {
                //cout << intersect_t << " ";
                PyramidPair p_pair = PyramidPair(py,intersect_t,Vector(0,0,-1));
                valids.push_back(p_pair);
            }
        }

        //take minimum t from valids and return for this ray
        if(valids.size() >= 1)
        {
            for(int m=0;m<valids.size()-1;m++)
            {
                for(int n=m+1;n<valids.size();n++)
                {
                    if(valids[n].t < valids[m].t)
                    {
                        PyramidPair temp = valids[m];
                        valids[m] = valids[n];
                        valids[n] = temp;
                    }
                }
            }
            t_from_all_pys.push_back(valids[0]);
        }
    }

    if(t_from_all_pys.size() >= 1)
    {
        for(int m=0;m<t_from_all_pys.size()-1;m++)
        {
            for(int n=m+1;n<t_from_all_pys.size();n++)
            {
                if(t_from_all_pys[n].t < t_from_all_pys[m].t)
                {
                    PyramidPair temp = t_from_all_pys[m];
                    t_from_all_pys[m] = t_from_all_pys[n];
                    t_from_all_pys[n] = temp;
                }
            }
        }
        minPair = t_from_all_pys[0];
    }
    else
    {
        minPair = PyramidPair();
    }

    /*if(minPair.t != -100.0)
    {
        cout << minPair.t << " ";
    }*/

    return minPair;
}

double getDistance(Vector v1, Vector v2)
{
    double ans = (v1.x-v2.x)*(v1.x-v2.x) + (v1.y-v2.y)*(v1.y-v2.y) + (v1.z-v2.z)*(v1.z-v2.z);
    return sqrt(ans);
}

int checkIfIntersects(Point R0, Vector ray, double dist)
{
    //check for checker board z=0
    double intersect_t = -R0.z/ray.z;
    //cout << intersect_t << " ";

    if(intersect_t > 0 && intersect_t < (_far - _near) && intersect_t < dist)
    {
        return -1;
    }

    //check intersection with pyramid, no intersection means p_pair.t = -100
    PyramidPair p_pair = get_closest_pyramid(Vector(R0.x,R0.y,R0.z),ray);
    if(p_pair.t != -100.0 && p_pair.t < dist)
    {
        return -1;
    }

    //check intersection with sphere
    SpherePair s_pair = get_closest_sphere(R0,ray);
    if(s_pair.t != -100.0 && s_pair.t < dist)
    {
        return -1;
    }

    return 1;
}

PolishComponent getPolishedComponent(Index_t_pair winner, Vector R0, Vector Rd)
{
    double t = winner.t;
    Vector P = R0 + Rd*t; // P is the intersecting point and S is the light source
    Vector N = winner.normal;

    double lambert = 0;
    double phong = 0;

    for(int i=0;i<normal_light_sources.size();i++)
    {
        NormalLightSource S = normal_light_sources[i];

        Vector S_point = Vector(S.position.x, S.position.y, S.position.z);
        Vector toSource = S_point - P;
        toSource.normalize();
        N.normalize();
        double distance = getDistance(P,S_point);
        //cout << distance << " ";

        //if S does not illumunate P continue
        int test = checkIfIntersects(Point(P.x+0.2,P.y+0.2,P.z+0.2), toSource, distance);
        if(test == -1)
        {
            //cout << "yo" << " ";
            continue;
        }

        double scaling_factor = exp(-distance*distance*S.falloff_parameter);
        lambert += Vector::dot(toSource,N)*scaling_factor;

        //find reflected ray at point P
        Vector sourceToP = P - S_point;
        //sourceToP.normalize();
        Vector R = Rd - N*(2*Vector::dot(Rd,N)); // reflected ray Rd or sourcetoP??
        R.normalize();
        phong += pow(Vector::dot(R,toSource),winner.shininess)*scaling_factor; //dot of R with N or toSource??
    }

    for(int i=0;i<spotlight_sources.size();i++)
    {
        SpotLightSource S = spotlight_sources[i];

        Vector S_point = Vector(S.position.x, S.position.y, S.position.z);
        Vector toSource = S_point - P;
        toSource.normalize();
        N.normalize();
        double distance = getDistance(P,S_point);
        //cout << distance << " ";

        //if S does not illumunate P continue
        int test = checkIfIntersects(Point(P.x+0.2,P.y+0.2,P.z+0.2), toSource, distance);
        if(test == -1)
        {
           // cout << "yo" << " ";
            continue;
        }

        //cout << "yo" << " ";

        //second test
        Vector V1 = P - S_point;
        V1.normalize();
        Vector V2 = Vector(S.look_point.x-S.position.x, S.look_point.y-S.position.y, S.look_point.z-S.position.z);
        V2.normalize();

        double angle_rad = acos(Vector::dot(V1,V2));
        double angle = angle_rad/pi*180.0;
        //V1.print();
        //cout << angle << " ";
        if(angle > S.cutoff_angle)
        {
            continue;
        }

        double scaling_factor = exp(-distance*distance*S.falloff_parameter);
        lambert += Vector::dot(toSource,N)*scaling_factor;

        //find reflected ray at point P
        Vector sourceToP = P - S_point;
        Vector R = Rd - N*(2*Vector::dot(Rd,N));
        R.normalize();
        phong += pow(Vector::dot(R,toSource),winner.shininess)*scaling_factor;
    }

    PolishComponent pc = PolishComponent(lambert,phong);
    return pc;
}

Color getObjectColor(Point R0, Vector ray, int level)
{

    if(level == 0)
    {
        return Color(0,0,0);
    }

     // 0 for checker, 1 for pyramid, 2 for sphere
    Index_t_pair index_t_pairs[3];
    for(int q=0;q<3;q++)
    {
        index_t_pairs[q].index = -1;
        index_t_pairs[q].t = 10000;
    }

    //check for checker board z=0
    double intersect_t = -R0.z/ray.z;
    //cout << intersect_t << " ";

    if(intersect_t > 0 && intersect_t < (_far - _near))
    {
        index_t_pairs[0] = Index_t_pair(0,intersect_t,Vector(0,0,1));
    }

    //check intersection with pyramid
    PyramidPair p_pair = get_closest_pyramid(Vector(R0.x,R0.y,R0.z),ray);
    if(p_pair.t != -100.0)
    {
        //cout << p_pair.t << " ";
        //Color clr = p_pair.pyramid.color;
        //pixels[j][i] = Color(clr.r*255.0,clr.g*255.0,clr.b*255.0);
        index_t_pairs[1] = Index_t_pair(1,p_pair.t,p_pair.normal);
        index_t_pairs[1].shininess = p_pair.pyramid.shininess;
    }

    //check intersection with sphere
    SpherePair s_pair = get_closest_sphere(R0,ray);
    if(s_pair.t != -100.0)
    {
        //cout << p_pair.t << " ";
        //Color clr = s_pair.sphere.color;
        //pixels[j][i] = Color(clr.r*255.0,clr.g*255.0,clr.b*255.0);
        index_t_pairs[2] = Index_t_pair(2,s_pair.t,s_pair.normal);
        index_t_pairs[2].shininess = s_pair.sphere.shininess;
    }

    //sort in ascending order of t
    for(int m=0;m<2;m++)
    {
        for(int n=m+1;n<3;n++)
        {
            if(index_t_pairs[n].t < index_t_pairs[m].t)
            {
                Index_t_pair temp = index_t_pairs[m];
                index_t_pairs[m] = index_t_pairs[n];
                index_t_pairs[n] = temp;
            }
        }
    }

    Index_t_pair winner;
    //if three_t has 10000 in it do not use it
    if(index_t_pairs[0].t == 10000)
    {
        // go to second
        if(index_t_pairs[1].t == 10000)
        {
            //go to third
            if(index_t_pairs[2].t == 10000)
            {
                return Color(0,0,0);
            }
            else
            {
                winner = index_t_pairs[2];
            }
        }
        else
        {
            winner = index_t_pairs[1];
        }
    }
    else
    {
        winner = index_t_pairs[0];
    }

    Color new_clr = Color(0,0,0);
    double final_r_coeff;

    if(winner.index == 0)
    {
        Color clr;

        if(texture_on == 0)
        {
            double check_t = winner.t;
            double intersect_x = R0.x + check_t*ray.x;
            double intersect_y = R0.y + check_t*ray.y;

            int checker_index_x = floor(intersect_x/checkerBoard.width);
            int checker_index_y = floor(intersect_y/checkerBoard.width);

            //cout << checker_index_x << " " << checker_index_y;

            if((checker_index_x+checker_index_y)%2 == 0)
            {
                //pixels[j][i] = Color(255,255,255);
                clr = Color(1,1,1);
            }
            else
            {
                //pixels[j][i] = Color(0,0,0);
                clr = Color(0,0,0);
            }

            winner.shininess = 0;
            PolishComponent pc = getPolishedComponent(winner,Vector(R0.x,R0.y,R0.z),ray);

            new_clr.r = checkerBoard.ambient_coeff*clr.r + checkerBoard.diffuse_coeff*pc.lambert*clr.r;
            new_clr.g = checkerBoard.ambient_coeff*clr.g + checkerBoard.diffuse_coeff*pc.lambert*clr.g;
            new_clr.b = checkerBoard.ambient_coeff*clr.b + checkerBoard.diffuse_coeff*pc.lambert*clr.b;
        }
        else
        {
            double check_t = winner.t;
            double intersect_x = R0.x + check_t*ray.x;
            double intersect_y = R0.y + check_t*ray.y;

            int checker_index_x = floor(intersect_x/checkerBoard.width);
            int checker_index_y = floor(intersect_y/checkerBoard.width);

            double corner_x = checker_index_x*checkerBoard.width;
            double corner_y = checker_index_y*checkerBoard.width;

            double cell_x = intersect_x - corner_x; //wrt to origin to get mapping from texture
            double cell_y = intersect_y - corner_y;

            double cell_h = checkerBoard.width/height;
            double cell_w = checkerBoard.width/width;

            int tb_x = floor(cell_x/cell_w);
            int tb_y = floor(cell_y/cell_h);

            //cout<<tb_x<<" "<<tb_y<<" ";

            winner.shininess = 0;
            PolishComponent pc = getPolishedComponent(winner,Vector(R0.x,R0.y,R0.z),ray);

            //cout << "hi" << " ";

            /*if((checker_index_x+checker_index_y)%2 == 0)
            {
                //pixels[j][i] = Color(255,255,255);
                clr = Color(1,1,1);
            }
            else
            {
                //pixels[j][i] = Color(0,0,0);
                clr = Color(0,0,0);
            }*/

            clr = Color(1,1,1);

            new_clr.r = textureBuffer[tb_x][tb_y].r*(checkerBoard.ambient_coeff*clr.r + checkerBoard.diffuse_coeff*pc.lambert*clr.r);
            new_clr.g = textureBuffer[tb_x][tb_y].g*(checkerBoard.ambient_coeff*clr.g + checkerBoard.diffuse_coeff*pc.lambert*clr.g);
            new_clr.b = textureBuffer[tb_x][tb_y].b*(checkerBoard.ambient_coeff*clr.b + checkerBoard.diffuse_coeff*pc.lambert*clr.b);
        }

        final_r_coeff = checkerBoard.reflection_coeff;
    }
    else if(winner.index == 1)
    {
        Color clr;
        clr = p_pair.pyramid.color;
        winner.shininess = p_pair.pyramid.shininess;
        PolishComponent pc = getPolishedComponent(winner,Vector(R0.x,R0.y,R0.z),ray);
        //cout << pc.lambert << " " << pc.phong << " ";
        Pyramid py = p_pair.pyramid;
        new_clr.r = py.ambient_coeff*clr.r + py.diffuse_coeff*pc.lambert*clr.r + py.specular_coeff*pc.phong*clr.r;
        new_clr.g = py.ambient_coeff*clr.g + py.diffuse_coeff*pc.lambert*clr.g + py.specular_coeff*pc.phong*clr.g;
        new_clr.b = py.ambient_coeff*clr.b + py.diffuse_coeff*pc.lambert*clr.b + py.specular_coeff*pc.phong*clr.b;
        final_r_coeff = py.reflection_coeff;
    }
    else if(winner.index == 2)
    {
        Color clr;
        clr = s_pair.sphere.color;
        //cout << clr.r << " " << clr.g << " " <<clr.b << " ";
        winner.shininess = s_pair.sphere.shininess;
        PolishComponent pc = getPolishedComponent(winner,Vector(R0.x,R0.y,R0.z),ray);

        Sphere sp = s_pair.sphere;
        new_clr.r = sp.ambient_coeff*clr.r + sp.diffuse_coeff*pc.lambert*clr.r + sp.specular_coeff*pc.phong*clr.r;
        new_clr.g = sp.ambient_coeff*clr.g + sp.diffuse_coeff*pc.lambert*clr.g + sp.specular_coeff*pc.phong*clr.g;
        new_clr.b = sp.ambient_coeff*clr.b + sp.diffuse_coeff*pc.lambert*clr.b + sp.specular_coeff*pc.phong*clr.b;
        final_r_coeff = sp.reflection_coeff;
    }

    Vector P = Vector(R0.x,R0.y,R0.z) + ray*winner.t; // P is the intersecting point to new object
    Vector N = winner.normal;
    N.normalize();
    Vector R = ray - N*(2*Vector::dot(ray,N));
    R.normalize();

    Color ref_clr = getObjectColor(Point(P.x+0.2, P.y+0.2, P.z+0.2), R, level-1);
    new_clr.r += (final_r_coeff*ref_clr.r);
    new_clr.g += (final_r_coeff*ref_clr.g);
    new_clr.b += (final_r_coeff*ref_clr.b);

    //cout << new_clr.r << " " << new_clr.g << " " <<new_clr.b << " ";
    return new_clr;
}

void populatePixelBuffer(PointBufferCell** pointbuffer,Color** pixels)
{
    /*if(flag == 0)
    {
        int cnt = 1;
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<5;j++)
            {
                pointbuffer[i][j].p.print();
                cnt ++;
            }
        }
        flag=1;
    }*/

    //get ray
    for(int i=0;i<screen_size;i++)
    {
        for(int j=0;j<screen_size;j++)
        {
            // 0 for checker, 1 for pyramid, 2 for sphere

            PointBufferCell buffer_cell = pointbuffer[i][j];
            Point temp = buffer_cell.p.subVector(Vector(pos.x,pos.y,pos.z));
            Vector ray = Vector(temp.x,temp.y,temp.z);
            ray.normalize();

            Color clr = getObjectColor(buffer_cell.p, ray, level_of_recursion);
            //cout << clr.r << " " << clr.g << " " <<clr.b << " ";
            pixels[j][i] = Color(clr.r*255.0,clr.g*255.0,clr.b*255.0);
        }
    }
}

void produceImage()
{
    if(produceimage == 1)
    {
        //glColor3f(0.1,0.7,0.1);
        //drawSphere(50,0,5,20,10);
        //drawPyramid(0,0,20,25);

        Color** pixels = new Color*[screen_size];

        for (int i = 0; i < screen_size; i++)
        {
            pixels[i] = new Color [screen_size];
            for (int j = 0; j < screen_size; j++)
            {
                pixels[i][j] = Color(0,0,0);
            }
        }
        //cout << "rendering\n";

        //we need to get colors for the pixel buffer using point buffer

        PointBufferCell** pointbuffer = new PointBufferCell*[screen_size];
        for (int i = 0; i < screen_size; i++)
        {
            pointbuffer[i] = new PointBufferCell [screen_size];
        }

        //sreen size means nymber of pixels

        double screen_height = 2.0*_near*tan((fov_y/2.0)/360.0*2*pi);
        double screen_width = 2.0*_near*tan((fov_x/2.0)/360.0*2*pi);
        double cell_width = screen_width/screen_size;
        double cell_height = screen_height/screen_size;

        Vector look = Vector(L.x,L.y,L.z);
        Vector up = Vector(U.x,U.y,U.z);
        Vector right = Vector(R.x,R.y,R.z);

        look.normalize();
        //cout << look.x << look.y << look.z;
        Point midpoint = Point();

        Vector temp_v = look*_near;
        midpoint = Point(pos.x,pos.y,pos.z).addVector(temp_v);

        Vector t1 = up*((screen_size-1)/2.0)*cell_height;
        Vector t2 = right*-((screen_size-1)/2.0)*cell_width;

        //t1.print();
        //t2.print();
        Point top_mid = midpoint.addVector(t1);
        Point left_corner = top_mid.addVector(t2);

        pointbuffer[0][0].p = left_corner;

        Point temp_p = pointbuffer[0][0].p;
        for(int i=1;i<screen_size;i++)
        {
            pointbuffer[i][0].p = temp_p.addVector(up*-cell_height);
            temp_p = pointbuffer[i][0].p;
        }

        temp_p = pointbuffer[0][0].p;
        for(int j=1;j<screen_size;j++)
        {
            pointbuffer[0][j].p = temp_p.addVector(right*cell_width);
            temp_p = pointbuffer[0][j].p;
        }

        for(int i=1;i<screen_size;i++)
        {
            for(int j=1;j<screen_size;j++)
            {
                pointbuffer[i][j].p = pointbuffer[i][j-1].p.addVector(right*cell_width);
            }
        }

        /*if(flag == 0)
        {
            int cnt = 1;
            for(int i=0;i<2;i++)
            {
                for(int j=0;j<5;j++)
                {
                    pointbuffer[i][j].p.print();
                    cnt ++;
                }
            }
            flag=1;
        }*/

        //cout << cnt;

        //find intersections and add colour to point buffer

        populatePixelBuffer(pointbuffer,pixels);


        // the following code generates a bmp image. do not change this.
        bitmap_image image(screen_size, screen_size);
        for (int x = 0; x < screen_size; x++) {
            for (int y = 0; y < screen_size; y++) {
                image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
            }
        }
        image.save_image("out.bmp");

        // free the dynamically allocated memory

        for (int i = 0; i < screen_size; i++) {
            delete [] pixels[i];
            delete [] pointbuffer[i];
        }
        delete [] pixels;
        delete [] pointbuffer;

        if(texture_on == 1)
        {
            for (int i = 0; i < width; i++) {
                delete [] textureBuffer[i];
            }
            delete [] textureBuffer;
        }

        if(flag == 0)
        {
            cout << "Rendering finished\n";
            flag = 1;
        }
        produceimage = 0;
    }
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '0':
            produceimage = 1;
            break;

		case '1':
			//drawgrid=1-drawgrid;
            rotate_RL_anticlockwise();
			break;

        case '2':
			//drawgrid=1-drawgrid;
            rotate_RL_clockwise();
            break;

        case '3':
            //rotate_angle += 1;
            rotate_UL_anticlockwise();
            break;

        case '4':
            //rotate_angle += 1;
            rotate_UL_clockwise();
            break;

        case '5':
            rotate_UR_clockwise();
            break;

        case '6':
            rotate_UR_anticlockwise();
            break;

        case ' ':
            texture_on = 1 - texture_on;
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){

	    case GLUT_KEY_UP:		// up arrow key
			pos.x += scale*L.x;
			pos.y += scale*L.y;
			pos.z += scale*L.z;
			break;

		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= scale*L.x;
			pos.y -= scale*L.y;
			pos.z -= scale*L.z;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			pos.x += scale*R.x;
			pos.y += scale*R.y;
			pos.z += scale*R.z;
			break;

		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
			pos.x -= scale*R.x;
			pos.y -= scale*R.y;
			pos.z -= scale*R.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += scale*U.x;
			pos.y += scale*U.y;
			pos.z += scale*U.z;
			break;

		case GLUT_KEY_PAGE_DOWN:
			pos.x -= scale*U.x;
			pos.y -= scale*U.y;
			pos.z -= scale*U.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;

		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				//drawaxes=1-drawaxes;

				//glTranslatef(-100,0,0);
				//glColor3f(1,0,0);
				//drawSphere(20,30,24);
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void input()
{
    ifstream scene;
    //ofstream stage1;
    scene.open ("description1.txt");
    //stage1.open ("stage1.txt");
    //stage1 << std::fixed;
    //stage1 << std::setprecision(7);

    string command;

    scene >> _near >> _far >> fov_y >> aspectRatio;
    fov_x = fov_y * aspectRatio;
    //cout << fov_y <<" "<< fov_x <<" "<< aspectRatio <<" "<< _near <<" "<< _far <<" ";

    scene >> level_of_recursion;
    scene >> screen_size;
    scene >> checkerBoard.width;
    scene >> checkerBoard.ambient_coeff >> checkerBoard.diffuse_coeff >> checkerBoard.reflection_coeff;
    scene >> number_of_objects;

    //cout << level_of_recursion <<" "<< screen_size <<" "<< checker_width <<" "<< checker_ambient_coeff <<" "<< checker_diffuse_coeff <<" "<< checker_reflection_coeff<<" "<<number_of_objects;
    for(int i=0;i<number_of_objects;i++)
    {
        scene >> command;
        //cout << command << "\n";
        if(command == "sphere")
        {
            Sphere sphere = Sphere();
            struct point cent;
            double r,g,b;

            scene >> cent.x >> cent.y >> cent.z;
            sphere.center = cent;
            scene >> sphere.radius;

            scene >> r >> g >> b;
            sphere.color = Color(r,g,b);

            scene >> sphere.ambient_coeff >> sphere.diffuse_coeff >> sphere.specular_coeff >> sphere.reflection_coeff;
            scene >> sphere.shininess;

            spheres.push_back(sphere);
        }
        else if (command == "pyramid")
        {
            Pyramid pyramid = Pyramid();
            struct point lowest;
            double r,g,b;

            scene >> lowest.x >> lowest.y >> lowest.z;
            pyramid.lowest_pt = lowest;
            scene >> pyramid.width >> pyramid.height;

            scene >> r >> g >> b;
            pyramid.color = Color(r,g,b);

            scene >> pyramid.ambient_coeff >> pyramid.diffuse_coeff >> pyramid.specular_coeff >> pyramid.reflection_coeff;
            scene >> pyramid.shininess;

            pyramids.push_back(pyramid);
        }
    }

    scene >> number_of_normal_lights;
    for(int i=0;i<number_of_normal_lights;i++)
    {
        NormalLightSource src = NormalLightSource();
        scene >> src.position.x >> src.position.y >> src.position.z >> src.falloff_parameter;
        normal_light_sources.push_back(src);
    }

    scene >> number_of_spot_lights;
    for(int i=0;i<number_of_spot_lights;i++)
    {
        SpotLightSource src = SpotLightSource();
        scene >> src.position.x >> src.position.y >> src.position.z >> src.falloff_parameter;
        scene >> src.look_point.x >> src.look_point.y >> src.look_point.z;
        scene >> src.cutoff_angle;
        spotlight_sources.push_back(src);
    }

    //cout << normal_light_sources.size();
    //cout << spotlight_sources.size();
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x,pos.y,pos.z, pos.x+L.x,pos.y+L.y,pos.z+L.z, U.x,U.y,U.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

    //glColor3f(0.1,0.7,0.1);
	//drawAxes();
	//drawGrid();


	//glPushMatrix();

	//glColor3f(0.2,0.2,0.2);
	//glTranslatef(-200,0,0);
	//glRotatef(90,0,1,0);
	//drawSquare(70);

	//glPopMatrix();

	//glRotatef(whole_gun_angle_ud,0,1,0);
	//glRotatef(whole_gun_angle_lr,0,0,1);
	//drawCapsule(10,40,24,20); //radius, hieght of cyln slices/segments, stacks

    //glColor3f(0.1,0.7,0.1);
    //drawSphere(0,0,15,5,20,10);
    //drawPyramid(0,0,15,20,25);

    drawObjects();
    getTextureBuffer();
    produceImage();


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){

	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	produceimage = 0;
	texture_on = 0;
	//drawgrid=1;
	//drawaxes=1;

    rotate_angle = 1;

    pos.x = 50;
    pos.y = 50;
    pos.z = 50;

	U.x = 0.0;
	U.y = 0.0;
	U.z = 1.0;

	R.x = -1/sqrt(2);
	R.y = 1/sqrt(2);
	R.z = 0;

	L.x = -1/sqrt(2);
	L.y = -1/sqrt(2);
	L.z = 0;


	//printf("%lf, %lf, %lf",U.x,U.y,U.z);

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	//cout << _near << " " << _far;
	gluPerspective(fov_y, aspectRatio, _near, _far);//obj has to be in btwn near and far to view
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

    input();
	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}

