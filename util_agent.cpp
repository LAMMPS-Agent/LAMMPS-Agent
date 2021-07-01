// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// Assume that classes are already given for the objects:
//    Point and Vector with
//        coordinates {float x, y, z;}
//        operators for:
//            Point   = Point ± Vector
//            Vector =  Point - Point
//            Vector =  Vector ± Vector
//            Vector =  Scalar * Vector
//    Line and Segment with defining points {Point  P0, P1;}
//    Track with initial position and velocity vector
//            {Point P0;  Vector v;}

/* ----------------------------------------------------------------------
   Modified by Changhao Li for modeling of bacteria film
   Last modified: 07/13/2020
------------------------------------------------------------------------- */

#define MY_UTIL

#include "util_agent.h"

// dist3D_Line_to_Line(): get the 3D minimum distance between 2 lines
//    Input:  two 3D lines L1 and L2
//    Return: the shortest distance between L1 and L2
float dist3D_Line_to_Line(Line L1, Line L2)
{
    Vector u = L1.P1 - L1.P0;
    Vector v = L2.P1 - L2.P0;
    Vector w = L1.P0 - L2.P0;
    float a = dot(u, u); // always >= 0
    float b = dot(u, v);
    float c = dot(v, v); // always >= 0
    float d = dot(u, w);
    float e = dot(v, w);
    float D = a * c - b * b; // always >= 0
    float sc, tc;

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM)
    { // the lines are almost parallel
        sc = 0.0;
        tc = (b > c ? d / b : e / c); // use the largest denominator
    }
    else
    {
        sc = (b * e - c * d) / D;
        tc = (a * e - b * d) / D;
    }

    // get the difference of the two closest points
    Vector dP = w + (sc * u) - (tc * v); // =  L1(sc) - L2(tc)

    return norm(dP); // return the closest distance
}
//===================================================================

// dist3D_Segment_to_Segment(): get the 3D minimum distance between 2 segments
//    Input:  two 3D line segments S1 and S2, and two float pointers of segment ratio (0 < t < 1).
//    Return: the shortest distance between S1 and S2
float dist3D_Segment_to_Segment(Segment S1, Segment S2, double &t1, double &t2)
{
    Vector u = S1.P1 - S1.P0;
    Vector v = S2.P1 - S2.P0;
    Vector w = S1.P0 - S2.P0;
    float a = dot(u, u); // always >= 0
    float b = dot(u, v);
    float c = dot(v, v); // always >= 0
    float d = dot(u, w);
    float e = dot(v, w);
    float D = a * c - b * b; // always >= 0
    float sc, sN, sD = D;    // sc = sN / sD, default sD = D >= 0
    float tc, tN, tD = D;    // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM)
    {             // the lines are almost parallel
        sN = 0.0; // force using point P0 on segment S1
        sD = 1.0; // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else
    { // get the closest points on the infinite lines
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN < 0.0)
        { // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD)
        { // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0)
    { // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
    { // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else
        {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (abs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    t1 = sc;
    t2 = tc;

    // get the difference of the two closest points
    Vector dP = w + (sc * u) - (tc * v); // =  S1(sc) - S2(tc)

    return norm(dP); // return the closest distance
}
//===================================================================

// cpa_time(): compute the time of CPA for two tracks
//    Input:  two tracks Tr1 and Tr2
//    Return: the time at which the two tracks are closest
float cpa_time(Track Tr1, Track Tr2)
{
    Vector dv = Tr1.v - Tr2.v;

    float dv2 = dot(dv, dv);
    if (dv2 < SMALL_NUM) // the  tracks are almost parallel
        return 0.0;      // any time is ok.  Use time 0.

    Vector w0 = Tr1.P0 - Tr2.P0;
    float cpatime = -dot(w0, dv) / dv2;

    return cpatime; // time of CPA
}
//===================================================================

// cpa_distance(): compute the distance at CPA for two tracks
//    Input:  two tracks Tr1 and Tr2
//    Return: the distance for which the two tracks are closest
float cpa_distance(Track Tr1, Track Tr2)
{
    float ctime = cpa_time(Tr1, Tr2);
    Point P1 = Tr1.P0 + (ctime * Tr1.v);
    Point P2 = Tr2.P0 + (ctime * Tr2.v);

    return d(P1, P2); // distance at CPA
}
//===================================================================

// utility function: to substitute PairBodyRoundedPolyhedron::distance_bt_edges(),
// which has bugs and cannot give correct minimum distance in some cases

void distance_bt_edges_new(const double *x1, const double *x2,
                           const double *x3, const double *x4,
                           double *h1, double *h2, double &t1, double &t2, double &r)
{

    Point P1(x1[0], x1[1], x1[2]), P2(x2[0], x2[1], x2[2]);
    Point P3(x3[0], x3[1], x3[2]), P4(x4[0], x4[1], x4[2]);
    Segment S1(P1, P2), S2(P3, P4);
    r = dist3D_Segment_to_Segment(S1, S2, t1, t2); // kernel function
    Point D1 = P1 + t1 * (P2 - P1);
    Point D2 = P3 + t2 * (P4 - P3);
    h1[0] = D1.x;
    h2[0] = D2.x;
    h1[1] = D1.y;
    h2[1] = D2.y;
    h1[2] = D1.z;
    h2[2] = D2.z;
}
