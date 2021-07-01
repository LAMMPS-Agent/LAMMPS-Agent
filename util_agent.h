/* ----------------------------------------------------------------------
    Author: Changhao Li (czl478@psu.edu, changhaoli1997@gmail.com)
    This header file defines some utility classes for Dan Sunday's algorithm of minimum distance
    between two line segments. The algorithm is implemented in myutil.cpp.

    To activate macro functions above, use "#define MY_UTIL" in your cpp file.
------------------------------------------------------------------------- */

#pragma once
#include <cmath>
#include <cstdlib>
#include <cstdio>

// activate macro functions
#ifdef MY_UTIL

#define SMALL_NUM 0.00000001 // anything that avoids division overflow
// dot product (3D) which allows vector operations in arguments
#define dot(u, v) ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define norm(v) sqrt(dot(v, v))        // norm = length of  vector
#define d(u, v) norm(u - v)            // distance = norm of difference
#define abs(x) ((x) >= 0 ? (x) : -(x)) //  absolute value

#endif

// 3-dimensional vector
class Vector
{
public:
    double x;
    double y;
    double z;
    Vector(double xx, double yy, double zz)
    {
        x = xx;
        y = yy;
        z = zz;
    }
    Vector operator+(const Vector &v) { return Vector(x + v.x, y + v.y, z + v.z); }
    Vector operator-(const Vector &v) { return Vector(x - v.x, y - v.y, z - v.z); }
    Vector operator*(double i) { return Vector(i * x, i * y, i * z); }
    Vector operator*(float i) { return Vector(i * x, i * y, i * z); }
    Vector operator*(int i) { return Vector(i * x, i * y, i * z); }

    friend Vector operator*(double, const Vector &);
    friend Vector operator*(float, const Vector &);
    friend Vector operator*(int, const Vector &);

    void print() { printf("[%f %f %f]", x, y, z); }
};

// 3-dimensional point
class Point
{
public:
    double x;
    double y;
    double z;
    Point()
    {
        x = (rand() % 200 - 100) / 10.0;
        y = (rand() % 200 - 100) / 10.0;
        z = (rand() % 200 - 100) / 10.0;
    }
    Point(double xx, double yy, double zz)
    {
        x = xx;
        y = yy;
        z = zz;
    }
    Point operator+(const Vector &v) { return Point(x + v.x, y + v.y, z + v.z); }
    Vector operator-(const Point &p) { return Vector(x - p.x, y - p.y, z - p.z); }

    friend Point operator+(const Vector &, const Point &);
    void print() { printf("[%f %f %f]", x, y, z); }
};

// infinite line, defined by two points
class Line
{
public:
    Point P0;
    Point P1;
};

// line segment, defined by two points
class Segment
{
public:
    Point P0;
    Point P1;
    Segment(Point iP0, Point iP1)
    {
        P0 = iP0;
        P1 = iP1;
    }
    void print()
    {
        P0.print();
        printf(" ");
        P1.print();
        printf("\n");
    }
};

// point with 3-dimensional velocity
class Track
{
public:
    Point P0;
    Vector v;
};

// functions to calculate distances between two lines/segments
float dist3D_Line_to_Line(Line L1, Line L2);
float dist3D_Segment_to_Segment(Segment S1, Segment S2, double &t1, double &t2);
float cpa_time(Track Tr1, Track Tr2);
float cpa_distance(Track Tr1, Track Tr2);

// overloaded operators for Vector and Point classes
inline Point operator+(const Vector &v, const Point &p) { return Point(v.x + p.x, v.y + p.y, v.z + p.z); }
inline Vector operator*(double i, const Vector &v) { return Vector(i * v.x, i * v.y, i * v.z); }
inline Vector operator*(float i, const Vector &v) { return Vector(i * v.x, i * v.y, i * v.z); }
inline Vector operator*(int i, const Vector &v) { return Vector(i * v.x, i * v.y, i * v.z); }

// LAMMPS workhorse function, to calculate distance between two line segments
// used in pair_body_rounded_polyhedron_agent.cpp
void distance_bt_edges_new(const double *x1, const double *x2,
                           const double *x3, const double *x4,
                           double *h1, double *h2, double &t1, double &t2, double &r);

// LAMMPS workhorse function, to calculate eqivalent distance between a line segments and a fixed wall
// used in fix_wall_body_polyhedron_agent.cpp
