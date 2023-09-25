#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

// Point struct for 3D points
struct Point3D
{
    double x, y, z;
    Point3D(double a, double b, double c) : x(a), y(b), z(c) {}
};

// Function to calculate the distance between two 3D points
double distance(const Point3D &p1, const Point3D &p2) { return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) + pow(p2.z - p1.z, 2)); }

// Function to fit a helix curve to three space points
void fitHelixToPoints(const Point3D &p1, const Point3D &p2, const Point3D &p3, double &R, double &pitch, double &p0, double &z0)
{
    // Calculate the distances between the three points
    double d12 = distance(p1, p2);
    double d23 = distance(p2, p3);
    double d31 = distance(p3, p1);

    // Calculate the angles between the three points
    double cos_alpha = (d23 * d23 + d31 * d31 - d12 * d12) / (2 * d23 * d31);
    double sin_alpha = sqrt(1 - cos_alpha * cos_alpha);
    double cos_beta = (d12 * d12 + d23 * d23 - d31 * d31) / (2 * d12 * d23);
    double sin_beta = sqrt(1 - cos_beta * cos_beta);

    // Calculate the pitch and radius of the helix
    pitch = atan2(sin_alpha, cos_alpha * sin_beta);
    R = d12 / (2 * sin_beta);

    // Calculate the center of the helix
    double xc = (p1.x + p2.x) / 2.0;
    double yc = (p1.y + p2.y) / 2.0;
    double zc = (p1.z + p2.z) / 2.0;

    // Calculate the azimuthal angle of the helix
    double phi0 = atan2(p2.y - yc, p2.x - xc);

    // Calculate the z position of the helix
    double z1 = p1.z - zc;
    double z2 = p2.z - zc;
    double z3 = p3.z - zc;
    double dz1 = z2 - z1;
    double dz2 = z3 - z2;
    double a = dz1 / (dz2 - dz1);
    double z0_est = z2 - a * (z3 - z2);

    // Calculate the azimuthal angle of the first point
    double phi1 = atan2(p1.y - yc, p1.x - xc);

    // Calculate the phase of the helix
    p0 = phi1 - pitch * z1 - phi0;

    // Set the helix center position
    z0 = z0_est + zc;
}

// Function to calculate the x and y position of the helix curve at a given z position
Point3D calculateHelixFunction(double R, double pitch, double p0, double z0, double z)
{
    double x = R * cos(p0 + pitch * (z - z0));
    double y = R * sin(p0 + pitch * (z - z0));
    return Point3D(x, y, z);
}

// Function to extrapolate the x and y position of a point given a z position
Point3D extrapolatePoint(const Point3D &p1, const Point3D &p2, const Point3D &p3, double z)
{
    double R, pitch, p0, z0;
    fitHelixToPoints(p1, p2, p3, R, pitch, p0, z0);
    Point3D closest_point = calculateHelixFunction(R, pitch, p0, z0, z);

    // Calculate the tangent vector at the closest point on the helix curve
    double dx = -R * pitch * sin(p0 + pitch * (closest_point.z - z0));
    double dy = R * pitch * cos(p0 + pitch * (closest_point.z - z0));
    double dz = 1;

    // Calculate the x and y position of the extrapolated point
    double dz1 = closest_point.z - p1.z;
    double dx1 = (p1.x - closest_point.x) / dz1;
    double dy1 = (p1.y - closest_point.y) / dz1;
    double dz2 = z - closest_point.z;
    double dx2 = dx * dz2;
    double dy2 = dy * dz2;

    double x = closest_point.x + dx1 * dz2 + dx2;
    double y = closest_point.y + dy1 * dz2 + dy2;

    return Point3D(x, y, z);
}

int main()
{
    Point3D p1(1, 0, 0);
    Point3D p2(0, 1, 0);
    Point3D p3(0, 0, 1);

    Point3D extrapolated_point = extrapolatePoint(p1, p2, p3, 2);

    cout << "Extrapolated point: (" << extrapolated_point.x << ", " << extrapolated_point.y << ", " << extrapolated_point.z << ")" << endl;

    return 0;
}
