#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include <Windows.h>
#include <algorithm>
#include <array>

extern "C"
{
#include "SpiceUsr.h"
}

// kilometers per astronomic unit
double AU = 149597870.7;

std::vector<std::array<int, 3>> major_body_colors = {
	{255, 245, 200},
	{169, 169, 169},
	{255, 238, 219},
	{100, 149, 237},
	{188, 39, 50},
	{218, 165, 105},
	{210, 180, 140},
	{173, 216, 230},
	{72, 61, 139}
};

std::vector<double> major_body_radii = {
	695508,
	4879 / 2,
	12104 / 2,
	12756 / 2,
	6792 / 2,
	142984 / 2,
	120536 / 2,
	51118 / 2,
	49528 / 2
};

double rad2deg(double x)
{
	return x * 180 / pi_c();
}

double deg2rad(double x)
{
	return x * pi_c() / 180;
}

struct Vec3 // extremely self-explanatory
{
	double x, y, z;

	Vec3()
	{
		x = 0;
		y = 0;
		z = 0;
	}

	Vec3(double xp, double yp, double zp)
	{
		x = xp;
		y = yp;
		z = zp;
	}

	Vec3(std::vector<double> vec)
	{
		x = vec[0];
		y = vec[1];
		z = vec[2];
	}

	Vec3(std::array<double, 3Ui64> vec)
	{
		x = vec[0];
		y = vec[1];
		z = vec[2];
	}

	Vec3(double RA, double DEC) // this is equivalent to a spherical2cartezian() function
	{
		x = cos(deg2rad(DEC)) * cos(deg2rad(RA));
		y = cos(deg2rad(DEC)) * sin(deg2rad(RA));
		z = sin(deg2rad(DEC));
	}

	Vec3 operator+(const Vec3& other) const
	{
		return { x + other.x, y + other.y, z + other.z };
	}

	Vec3 operator-(const Vec3& other) const
	{
		return { x - other.x, y - other.y, z - other.z };
	}

	Vec3 operator*(double scalar) const
	{
		return { x * scalar, y * scalar, z * scalar };
	}

	Vec3 operator/(double scalar) const
	{
		return { x / scalar, y / scalar, z / scalar };
	}

	Vec3& operator+=(const Vec3& other)
	{
		x += other.x; y += other.y; z += other.z;
		return *this;
	}

	Vec3 operator-() const
	{
		return Vec3(-x, -y, -z);
	}

	Vec3 cross(const Vec3& other)
	{
		return Vec3(y * other.z - z * other.y,
			z * other.x - x * other.z,
			x * other.y - y * other.x);
	}

	double dot(const Vec3& other)
	{
		return x * other.x + y * other.y + z * other.z;
	}

	Vec3 normalized()
	{
		return Vec3(x, y, z) / Vec3(x, y, z).mag();
	}

	double mag()
	{
		return sqrt(x * x + y * y + z * z);
	}

	void printout()
	{
		std::cout << "Vec3(" << x << ", " << y << ", " << z << ")\n";
	}
};

class State
{
public:
	std::string desig;
	double JD;
	std::string datetime;
	Vec3 p;
	Vec3 v;
};

using StateMatrix = std::array<std::array<std::array<double, 3>, 2>, 9>;

// Constant: font8x8_basic
// Contains an 8x8 font map for unicode points U+0000 - U+007F (basic latin)
uint8_t font8x8_basic[128][8] = {
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0000 (nul)
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0001
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0002
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0003
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0004
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0005
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0006
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0007
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0008
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0009
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000A
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000B
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000C
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000D
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000E
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+000F
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0010
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0011
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0012
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0013
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0014
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0015
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0016
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0017
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0018
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0019
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001A
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001B
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001C
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001D
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001E
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+001F
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0020 (space)
	{ 0x18, 0x3C, 0x3C, 0x18, 0x18, 0x00, 0x18, 0x00},   // U+0021 (!)
	{ 0x36, 0x36, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0022 (")
	{ 0x36, 0x36, 0x7F, 0x36, 0x7F, 0x36, 0x36, 0x00},   // U+0023 (#)
	{ 0x0C, 0x3E, 0x03, 0x1E, 0x30, 0x1F, 0x0C, 0x00},   // U+0024 ($)
	{ 0x00, 0x63, 0x33, 0x18, 0x0C, 0x66, 0x63, 0x00},   // U+0025 (%)
	{ 0x1C, 0x36, 0x1C, 0x6E, 0x3B, 0x33, 0x6E, 0x00},   // U+0026 (&)
	{ 0x06, 0x06, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0027 (')
	{ 0x18, 0x0C, 0x06, 0x06, 0x06, 0x0C, 0x18, 0x00},   // U+0028 (()
	{ 0x06, 0x0C, 0x18, 0x18, 0x18, 0x0C, 0x06, 0x00},   // U+0029 ())
	{ 0x00, 0x66, 0x3C, 0xFF, 0x3C, 0x66, 0x00, 0x00},   // U+002A (*)
	{ 0x00, 0x0C, 0x0C, 0x3F, 0x0C, 0x0C, 0x00, 0x00},   // U+002B (+)
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C, 0x06},   // U+002C (,)
	{ 0x00, 0x00, 0x00, 0x3F, 0x00, 0x00, 0x00, 0x00},   // U+002D (-)
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C, 0x00},   // U+002E (.)
	{ 0x60, 0x30, 0x18, 0x0C, 0x06, 0x03, 0x01, 0x00},   // U+002F (/)
	{ 0x3E, 0x63, 0x73, 0x7B, 0x6F, 0x67, 0x3E, 0x00},   // U+0030 (0)
	{ 0x0C, 0x0E, 0x0C, 0x0C, 0x0C, 0x0C, 0x3F, 0x00},   // U+0031 (1)
	{ 0x1E, 0x33, 0x30, 0x1C, 0x06, 0x33, 0x3F, 0x00},   // U+0032 (2)
	{ 0x1E, 0x33, 0x30, 0x1C, 0x30, 0x33, 0x1E, 0x00},   // U+0033 (3)
	{ 0x38, 0x3C, 0x36, 0x33, 0x7F, 0x30, 0x78, 0x00},   // U+0034 (4)
	{ 0x3F, 0x03, 0x1F, 0x30, 0x30, 0x33, 0x1E, 0x00},   // U+0035 (5)
	{ 0x1C, 0x06, 0x03, 0x1F, 0x33, 0x33, 0x1E, 0x00},   // U+0036 (6)
	{ 0x3F, 0x33, 0x30, 0x18, 0x0C, 0x0C, 0x0C, 0x00},   // U+0037 (7)
	{ 0x1E, 0x33, 0x33, 0x1E, 0x33, 0x33, 0x1E, 0x00},   // U+0038 (8)
	{ 0x1E, 0x33, 0x33, 0x3E, 0x30, 0x18, 0x0E, 0x00},   // U+0039 (9)
	{ 0x00, 0x0C, 0x0C, 0x00, 0x00, 0x0C, 0x0C, 0x00},   // U+003A (:)
	{ 0x00, 0x0C, 0x0C, 0x00, 0x00, 0x0C, 0x0C, 0x06},   // U+003B (;)
	{ 0x18, 0x0C, 0x06, 0x03, 0x06, 0x0C, 0x18, 0x00},   // U+003C (<)
	{ 0x00, 0x00, 0x3F, 0x00, 0x00, 0x3F, 0x00, 0x00},   // U+003D (=)
	{ 0x06, 0x0C, 0x18, 0x30, 0x18, 0x0C, 0x06, 0x00},   // U+003E (>)
	{ 0x1E, 0x33, 0x30, 0x18, 0x0C, 0x00, 0x0C, 0x00},   // U+003F (?)
	{ 0x3E, 0x63, 0x7B, 0x7B, 0x7B, 0x03, 0x1E, 0x00},   // U+0040 (@)
	{ 0x0C, 0x1E, 0x33, 0x33, 0x3F, 0x33, 0x33, 0x00},   // U+0041 (A)
	{ 0x3F, 0x66, 0x66, 0x3E, 0x66, 0x66, 0x3F, 0x00},   // U+0042 (B)
	{ 0x3C, 0x66, 0x03, 0x03, 0x03, 0x66, 0x3C, 0x00},   // U+0043 (C)
	{ 0x1F, 0x36, 0x66, 0x66, 0x66, 0x36, 0x1F, 0x00},   // U+0044 (D)
	{ 0x7F, 0x46, 0x16, 0x1E, 0x16, 0x46, 0x7F, 0x00},   // U+0045 (E)
	{ 0x7F, 0x46, 0x16, 0x1E, 0x16, 0x06, 0x0F, 0x00},   // U+0046 (F)
	{ 0x3C, 0x66, 0x03, 0x03, 0x73, 0x66, 0x7C, 0x00},   // U+0047 (G)
	{ 0x33, 0x33, 0x33, 0x3F, 0x33, 0x33, 0x33, 0x00},   // U+0048 (H)
	{ 0x1E, 0x0C, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00},   // U+0049 (I)
	{ 0x78, 0x30, 0x30, 0x30, 0x33, 0x33, 0x1E, 0x00},   // U+004A (J)
	{ 0x67, 0x66, 0x36, 0x1E, 0x36, 0x66, 0x67, 0x00},   // U+004B (K)
	{ 0x0F, 0x06, 0x06, 0x06, 0x46, 0x66, 0x7F, 0x00},   // U+004C (L)
	{ 0x63, 0x77, 0x7F, 0x7F, 0x6B, 0x63, 0x63, 0x00},   // U+004D (M)
	{ 0x63, 0x67, 0x6F, 0x7B, 0x73, 0x63, 0x63, 0x00},   // U+004E (N)
	{ 0x1C, 0x36, 0x63, 0x63, 0x63, 0x36, 0x1C, 0x00},   // U+004F (O)
	{ 0x3F, 0x66, 0x66, 0x3E, 0x06, 0x06, 0x0F, 0x00},   // U+0050 (P)
	{ 0x1E, 0x33, 0x33, 0x33, 0x3B, 0x1E, 0x38, 0x00},   // U+0051 (Q)
	{ 0x3F, 0x66, 0x66, 0x3E, 0x36, 0x66, 0x67, 0x00},   // U+0052 (R)
	{ 0x1E, 0x33, 0x07, 0x0E, 0x38, 0x33, 0x1E, 0x00},   // U+0053 (S)
	{ 0x3F, 0x2D, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00},   // U+0054 (T)
	{ 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x3F, 0x00},   // U+0055 (U)
	{ 0x33, 0x33, 0x33, 0x33, 0x33, 0x1E, 0x0C, 0x00},   // U+0056 (V)
	{ 0x63, 0x63, 0x63, 0x6B, 0x7F, 0x77, 0x63, 0x00},   // U+0057 (W)
	{ 0x63, 0x63, 0x36, 0x1C, 0x1C, 0x36, 0x63, 0x00},   // U+0058 (X)
	{ 0x33, 0x33, 0x33, 0x1E, 0x0C, 0x0C, 0x1E, 0x00},   // U+0059 (Y)
	{ 0x7F, 0x63, 0x31, 0x18, 0x4C, 0x66, 0x7F, 0x00},   // U+005A (Z)
	{ 0x1E, 0x06, 0x06, 0x06, 0x06, 0x06, 0x1E, 0x00},   // U+005B ([)
	{ 0x03, 0x06, 0x0C, 0x18, 0x30, 0x60, 0x40, 0x00},   // U+005C (\)
	{ 0x1E, 0x18, 0x18, 0x18, 0x18, 0x18, 0x1E, 0x00},   // U+005D (])
	{ 0x08, 0x1C, 0x36, 0x63, 0x00, 0x00, 0x00, 0x00},   // U+005E (^)
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF},   // U+005F (_)
	{ 0x0C, 0x0C, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+0060 (`)
	{ 0x00, 0x00, 0x1E, 0x30, 0x3E, 0x33, 0x6E, 0x00},   // U+0061 (a)
	{ 0x07, 0x06, 0x06, 0x3E, 0x66, 0x66, 0x3B, 0x00},   // U+0062 (b)
	{ 0x00, 0x00, 0x1E, 0x33, 0x03, 0x33, 0x1E, 0x00},   // U+0063 (c)
	{ 0x38, 0x30, 0x30, 0x3e, 0x33, 0x33, 0x6E, 0x00},   // U+0064 (d)
	{ 0x00, 0x00, 0x1E, 0x33, 0x3f, 0x03, 0x1E, 0x00},   // U+0065 (e)
	{ 0x1C, 0x36, 0x06, 0x0f, 0x06, 0x06, 0x0F, 0x00},   // U+0066 (f)
	{ 0x00, 0x00, 0x6E, 0x33, 0x33, 0x3E, 0x30, 0x1F},   // U+0067 (g)
	{ 0x07, 0x06, 0x36, 0x6E, 0x66, 0x66, 0x67, 0x00},   // U+0068 (h)
	{ 0x0C, 0x00, 0x0E, 0x0C, 0x0C, 0x0C, 0x1E, 0x00},   // U+0069 (i)
	{ 0x30, 0x00, 0x30, 0x30, 0x30, 0x33, 0x33, 0x1E},   // U+006A (j)
	{ 0x07, 0x06, 0x66, 0x36, 0x1E, 0x36, 0x67, 0x00},   // U+006B (k)
	{ 0x0E, 0x0C, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00},   // U+006C (l)
	{ 0x00, 0x00, 0x33, 0x7F, 0x7F, 0x6B, 0x63, 0x00},   // U+006D (m)
	{ 0x00, 0x00, 0x1F, 0x33, 0x33, 0x33, 0x33, 0x00},   // U+006E (n)
	{ 0x00, 0x00, 0x1E, 0x33, 0x33, 0x33, 0x1E, 0x00},   // U+006F (o)
	{ 0x00, 0x00, 0x3B, 0x66, 0x66, 0x3E, 0x06, 0x0F},   // U+0070 (p)
	{ 0x00, 0x00, 0x6E, 0x33, 0x33, 0x3E, 0x30, 0x78},   // U+0071 (q)
	{ 0x00, 0x00, 0x3B, 0x6E, 0x66, 0x06, 0x0F, 0x00},   // U+0072 (r)
	{ 0x00, 0x00, 0x3E, 0x03, 0x1E, 0x30, 0x1F, 0x00},   // U+0073 (s)
	{ 0x08, 0x0C, 0x3E, 0x0C, 0x0C, 0x2C, 0x18, 0x00},   // U+0074 (t)
	{ 0x00, 0x00, 0x33, 0x33, 0x33, 0x33, 0x6E, 0x00},   // U+0075 (u)
	{ 0x00, 0x00, 0x33, 0x33, 0x33, 0x1E, 0x0C, 0x00},   // U+0076 (v)
	{ 0x00, 0x00, 0x63, 0x6B, 0x7F, 0x7F, 0x36, 0x00},   // U+0077 (w)
	{ 0x00, 0x00, 0x63, 0x36, 0x1C, 0x36, 0x63, 0x00},   // U+0078 (x)
	{ 0x00, 0x00, 0x33, 0x33, 0x33, 0x3E, 0x30, 0x1F},   // U+0079 (y)
	{ 0x00, 0x00, 0x3F, 0x19, 0x0C, 0x26, 0x3F, 0x00},   // U+007A (z)
	{ 0x38, 0x0C, 0x0C, 0x07, 0x0C, 0x0C, 0x38, 0x00},   // U+007B ({)
	{ 0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x00},   // U+007C (|)
	{ 0x07, 0x0C, 0x0C, 0x38, 0x0C, 0x0C, 0x07, 0x00},   // U+007D (})
	{ 0x6E, 0x3B, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},   // U+007E (~)
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}    // U+007F
};

bool createDirectoryIfNotExists(const std::string& dir_name)
{
	DWORD attribs = GetFileAttributesA(dir_name.c_str());
	if (attribs != INVALID_FILE_ATTRIBUTES && (attribs & FILE_ATTRIBUTE_DIRECTORY))
	{
		return true;
	}

	if (CreateDirectoryA(dir_name.c_str(), NULL))
	{
		return true;
	}

	return false;
}

void loadAllKernels(const std::string& directory)
{
	std::string search_path = directory + "\\*.*";
	WIN32_FIND_DATAA fd;
	HANDLE hFind = ::FindFirstFileA(search_path.c_str(), &fd);

	if (hFind == INVALID_HANDLE_VALUE) {
		std::cerr << "Unable to open directory: " << directory << '\n';
		return;
	}

	do {
		if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			std::string filepath = directory + "\\" + fd.cFileName;
			furnsh_c(filepath.c_str());
			// std::cout << "Loaded kernel: " << filepath << '\n';
		}
	} while (::FindNextFileA(hFind, &fd));

	::FindClose(hFind);
}

StateMatrix getSolarSystemStates(SpiceDouble et)
{
	const char* bodies[9] = {
		"SUN",                 // index 0
		"MERCURY BARYCENTER",  // index 1
		"VENUS BARYCENTER",    // index 2
		"EARTH BARYCENTER",    // index 3
		"MARS BARYCENTER",     // index 4
		"JUPITER BARYCENTER",  // index 5
		"SATURN BARYCENTER",   // index 6
		"URANUS BARYCENTER",   // index 7
		"NEPTUNE BARYCENTER"   // index 8
	};

	StateMatrix states{};

	for (int i = 0; i < 9; ++i)
	{
		SpiceDouble state[6];
		SpiceDouble lt;

		spkezr_c(bodies[i], et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);

		for (int j = 0; j < 3; ++j)
		{
			states[i][0][j] = state[j];     // Position (km)
			states[i][1][j] = state[j + 3]; // Velocity (km/s)
		}
	}

	return states;
}

// feed ecliptic state vectors to this!!
std::vector<double> eclStateVector2Kepler(Vec3 r, Vec3 v, double mu = 1.3271244004193938E+11)
{
	double r_mag = r.mag();
	double v_mag = v.mag();

	Vec3 h = r.cross(v);
	double h_mag = h.mag();

	double inclination = rad2deg(acos(h.z / h_mag));

	Vec3 k = Vec3(0, 0, 1);
	Vec3 n = k.cross(h);
	double n_mag = n.mag();

	double omega;
	double eccentricity;
	double arg_periapsis;
	double true_anomaly;
	double sma;
	double mean_anomaly;

	if (n_mag != 0)
	{
		omega = rad2deg(acos(n.x / n_mag));
		if (n.y < 0)
		{
			omega = 360 - omega;
		}
	}
	else
	{
		omega = 0;
	}

	Vec3 e_vec = (v.cross(h) - r * mu / r_mag) * (1 / mu);
	eccentricity = e_vec.mag();

	if (n_mag != 0)
	{
		if (eccentricity != 0)
		{
			arg_periapsis = rad2deg(acos(n.dot(e_vec) / (n_mag * eccentricity)));
			if (e_vec.z < 0)
			{
				arg_periapsis = 360 - arg_periapsis;
			}
		}
		else
		{
			arg_periapsis = 0;
		}
	}
	else
	{
		arg_periapsis = 0;
	}

	if (eccentricity != 0)
	{
		true_anomaly = rad2deg(acos(e_vec.dot(r) / (eccentricity * r_mag)));
		if (r.dot(v) < 0)
		{
			true_anomaly = 360 - true_anomaly;
		}
	}
	else
	{
		true_anomaly = rad2deg(acos(r.normalized().dot(v.normalized())));
	}

	double specific_energy = v_mag * v_mag / 2 - mu / r_mag;
	if (abs(eccentricity - 1) > 1e-8)
	{
		sma = -mu / (2 * specific_energy);
	}
	else
	{
		sma = 999999;
	}

	if (eccentricity < 1)
	{
		double E = 2 * atan(tan(deg2rad(true_anomaly) / 2) * sqrt((1 - eccentricity) / (1 + eccentricity)));
		if (E < 0)
		{
			E = E + 2 * pi_c();
		}

		mean_anomaly = rad2deg(E - eccentricity * sin(E));
	}
	else if (eccentricity > 1)
	{
		double F = 2 * atanh(tan(deg2rad(true_anomaly) / 2) * sqrt((eccentricity - 1) / (eccentricity + 1)));
		mean_anomaly = rad2deg(eccentricity * sinh(F) - F);
	}
	else
	{
		mean_anomaly = -1.0; // random val.
	}

	return std::vector<double> {sma, eccentricity, inclination, omega, arg_periapsis, true_anomaly, mean_anomaly};

}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> readTycho2(const std::string& filename = "data/Tycho2.csv")
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		throw std::runtime_error("Cannot open star catalog: " + filename);
	}

	std::string line;

	// skip header
	if (!std::getline(file, line))
	{
		throw std::runtime_error("File is empty or cannot read header!");
	}

	std::vector<double> mags, RAs, DECs;

	while (std::getline(file, line))
	{
		std::stringstream ss(line);
		std::string cell;
		int col_index = 0;
		double magval = 0, RAval = 0, DECval = 0;
		bool gotmag = false, gotRA = false, gotDEC = false;

		while (std::getline(ss, cell, ','))
		{
			++col_index;
			try
			{
				if (col_index == 3)
				{
					magval = std::stod(cell);
					gotmag = true;
				}
				else if (col_index == 4)
				{
					RAval = std::stod(cell);
					gotRA = true;
				}
				else if (col_index == 5)
				{
					DECval = std::stod(cell);
					gotDEC = true;
				}
			}
			catch (const std::invalid_argument&)
			{
				throw std::runtime_error("Invalid double value in column " + std::to_string(col_index) + ": " + cell);
			}
		}

		if (gotmag && gotRA && gotDEC) {
			mags.push_back(magval);
			RAs.push_back(RAval);
			DECs.push_back(DECval);
		}
		else {
			throw std::runtime_error("Missing expected columns in line: " + line);
		}
	}

	return std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>{ mags, RAs, DECs };
}

std::vector<State> readStateVectorFile(const std::string& filename)
{
	std::ifstream infile(filename);
	std::string line;

	std::vector<State> states;

	if (!infile) {
		std::cerr << "Failed to open file: " << filename << "\n";
		{
			return {};
		}
	}

	// skip header
	while (std::getline(infile, line)) {
		if (line.find('*') != std::string::npos)
			break;
	}

	while (std::getline(infile, line)) {
		if (line.find('*') != std::string::npos)
			break;

		std::istringstream iss(line);
		double jd, x, y, z, vx, vy, vz;
		std::string utc;

		if (!(iss >> jd >> utc >> x >> y >> z >> vx >> vy >> vz)) {
			std::cerr << "Skipping bad line:\n" << line << '\n';
			continue;
		}

		State new_state;
		new_state.JD = jd;
		new_state.datetime = utc;
		new_state.p = Vec3(x, y, z);
		new_state.v = Vec3(vx, vy, vz);

		states.push_back(new_state);
	}

	return states;
}

double sexRAToDeg(const std::string& RA_str)
{
	int hours, minutes;
	double seconds;
	char sep1, sep2;
	std::istringstream iss(RA_str);
	iss >> hours >> sep1 >> minutes >> sep2 >> seconds;

	if (sep1 != ':' || sep2 != ':') {
		throw std::runtime_error("Invalid RA format: " + RA_str);
	}

	return 15.0 * (hours + minutes / 60.0 + seconds / 3600.0);
}

double sexDECToDeg(const std::string& DEC_str)
{
	int degrees, arcmin;
	double arcsec;
	char sep1, sep2;
	char sign = DEC_str[0];
	std::istringstream iss(DEC_str);
	iss >> degrees >> sep1 >> arcmin >> sep2 >> arcsec;

	if (sep1 != ':' || sep2 != ':') {
		throw std::runtime_error("Invalid DEC format: " + DEC_str);
	}

	double abs_deg = std::abs(degrees) + arcmin / 60.0 + arcsec / 3600.0;
	return degrees < 0 ? -abs_deg : abs_deg;
}

void drawRedCrosshair(std::vector<std::array<int, 3>>& img, int screen_size)
{
	int cx = screen_size / 2;
	int cy = screen_size / 2;
	const std::array<int, 3> red = { 255, 0, 0 };

	// Draw upward arm
	for (int dy = -3; dy > -3 - 5; --dy)
	{
		int y = cy + dy;
		if (y >= 0 && y < screen_size)
			img[y * screen_size + cx] = red;
	}

	// Draw downward arm
	for (int dy = 3; dy < 3 + 5; ++dy)
	{
		int y = cy + dy;
		if (y >= 0 && y < screen_size)
			img[y * screen_size + cx] = red;
	}

	// Draw left arm
	for (int dx = -3; dx > -3 - 5; --dx)
	{
		int x = cx + dx;
		if (x >= 0 && x < screen_size)
			img[cy * screen_size + x] = red;
	}

	// Draw right arm
	for (int dx = 3; dx < 3 + 5; ++dx)
	{
		int x = cx + dx;
		if (x >= 0 && x < screen_size)
			img[cy * screen_size + x] = red;
	}
}

void drawChar(std::vector<std::array<int, 3>>& img, int screen_x, int screen_y,
	int x0, int y0, const uint8_t bitmap[8], std::array<int, 3> color)
{
	for (int row = 0; row < 8; ++row)
	{
		for (int col = 0; col < 8; ++col)
		{
			if (bitmap[row] & (1 << col)) // LSB-left: leftmost pixel = bit 0
			{
				int x = x0 + col;
				int y = y0 + row;
				if (x >= 0 && x < screen_x && y >= 0 && y < screen_y)
				{
					img[y * screen_x + x] = color;
				}
			}
		}
	}
}

void drawText(std::vector<std::array<int, 3>>& img, int screen_x, int screen_y,
	int x, int y, const std::string& text, std::array<int, 3> color)
{
	for (char c : text)
	{
		const uint8_t* bitmap = font8x8_basic[(unsigned char)c]; // assuming it's defined
		drawChar(img, screen_x, screen_y, x, y, bitmap, color);
		x += 8; // fixed spacing
	}
}

void drawCircle(std::vector<std::array<int, 3>>& img, int screen_x, int screen_y, int cx, int cy, int radius, std::array<int, 3> color = {255, 255, 255})
{
	for (int dy = -radius; dy <= radius; ++dy)
	{
		for (int dx = -radius; dx <= radius; ++dx)
		{
			if (dx * dx + dy * dy <= radius * radius)
			{
				int x = cx + dx;
				int y = cy + dy;
				if (x >= 0 && x < screen_x && y >= 0 && y < screen_y)
				{
					int idx = y * screen_x + x;
					img[idx] = color;
				}
			}
		}
	}
}

void drawLine(std::vector<std::array<int, 3>>& img, int screen_x, int screen_y,
	int x0, int y0, int x1, int y1, std::array<int, 3> color)
{
	if (x0 < 0 || y0 < 0 || x1 < 0 || y1 < 0)
	{
		return;
	}

	int dx = std::abs(x1 - x0);
	int dy = -std::abs(y1 - y0);
	int sx = (x0 < x1) ? 1 : -1;
	int sy = (y0 < y1) ? 1 : -1;
	int err = dx + dy; // error value

	int iters = 0;
	const int MAXITERS = max(screen_x, screen_y) * 2;

	while (true)
	{
		if (x0 >= 0 && x0 < screen_x && y0 >= 0 && y0 < screen_y)
		{
			int idx = y0 * screen_x + x0;
			img[idx] = color;
		}

		if (x0 == x1 && y0 == y1) 
		{
			break;
		}

		int e2 = 2 * err;

		if (e2 >= dy)
		{
			err += dy;
			x0 += sx;
		}
		if (e2 <= dx)
		{
			err += dx;
			y0 += sy;
		}

		iters++;

		if (iters > MAXITERS) // failsafe
		{
			break;
		}
	}
}

std::vector<Vec3> getKeplerOrbitPoints(Vec3 p, Vec3 v, int N_points = 720)
{
	// returns {sma, eccentricity, inclination, omega, arg_periapsis, true_anomaly, mean_anomaly}
	std::vector<double> orbital_elems = eclStateVector2Kepler(p, v);

	double a = orbital_elems[0];
	double e = orbital_elems[1];
	double i = deg2rad(orbital_elems[2]);
	double omega = deg2rad(orbital_elems[3]);
	double arg_periapsis = deg2rad(orbital_elems[4]);

	std::vector<Vec3> orbit_points;

	// use (N_points + 1) points to close the curve
	for (int k = 0; k < N_points + 1; ++k)
	{
		double nu = 2.0 * pi_c() * k / N_points; // true anomaly

		// Compute radius in orbital plane
		double r = a * (1 - e * e) / (1 + e * cos(nu));
		double x_orb = r * cos(nu);
		double y_orb = r * sin(nu);
		double z_orb = 0;

		// orbital transformation
		double cos_o = cos(omega);
		double sin_o = sin(omega);
		double cos_i = cos(i); 
		double sin_i = sin(i);
		double cos_w = cos(arg_periapsis);
		double sin_w = sin(arg_periapsis);

		// position in orbital plane
		double x1 = cos_w * x_orb - sin_w * y_orb;
		double y1 = sin_w * x_orb + cos_w * y_orb;
		double z1 = 0;

		// rotate by inclination
		double x2 = x1;
		double y2 = cos_i * y1;
		double z2 = sin_i * y1;

		// rotate by longitude of ascending node
		double x = cos_o * x2 - sin_o * y2;
		double y = sin_o * x2 + cos_o * y2;
		double z = z2;

		orbit_points.push_back(Vec3(x, y, z));
	}

	return orbit_points;
}

std::vector<int> space2screen(Vec3 pos, Vec3 cam_pos, std::vector<Vec3> cam_orient, double f, int screen_x, int screen_y)
{
	// OpenGL-esque
	Vec3 cam_right = cam_orient[0];
	Vec3 cam_up = cam_orient[1];
	Vec3 cam_forward = -cam_orient[2];

	Vec3 rel_pos = pos - cam_pos;

	if (rel_pos.dot(cam_forward) < 0)
	{
		return { -1, -1 };
	}

	double px = f * rel_pos.dot(cam_right) / rel_pos.dot(cam_forward);
	double py = f * rel_pos.dot(cam_up) / rel_pos.dot(cam_forward);

	int pix_x = screen_x / 2 + px + 0.5;
	int pix_y = screen_y / 2 - py + 0.5;

	return std::vector<int> {pix_x, pix_y};
}

// render a single individual image
void renderSolarSystem(State st, Vec3 mp_pos, std::vector<Vec3> mp_orbit,
	std::vector<Vec3> major_pos, std::vector<std::vector<Vec3>> major_orbits,
	std::string cam_mode, double fov, int screen_x, int screen_y,
	Vec3 cam_pos, std::vector<Vec3> cam_orient,
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> starfield,
	std::string save_name)
{
	std::vector<std::array<int, 3>> img(screen_x * screen_y, { 0, 0, 0 });
	int screen_short = min(screen_x, screen_y);

	std::string screen_short_dir = "y";
	if (screen_short == screen_x)
	{
		screen_short_dir = "x";
	}

	// perspective offset factor
	double f = screen_y / (2 * tan(fov / 2));
	if (!strcmp(screen_short_dir.c_str(), "x"))
	{
		f = screen_x / (2 * tan(fov / 2));
	}

	// OpenGL-esque
	Vec3 cam_right = cam_orient[0];
	Vec3 cam_up = cam_orient[1];
	Vec3 cam_forward = -cam_orient[2];

	// now, we render things from back to front as basic renderers do
	// so...
	// draw starfield first
	SpiceDouble et;
	utc2et_c(st.datetime.c_str(), &et);

	// gotta get star coordinates in ecliptic now
	double equ_ecl_rot[3][3];
	pxform_c("J2000", "ECLIPJ2000", et, equ_ecl_rot);

	for (int idx_star = 0; idx_star < std::get<0>(starfield).size(); idx_star++)
	{
		if (std::get<0>(starfield)[idx_star] > 6) // don't draw too dim stars, clutters the background
		{
			continue;
		}

		double radius = 1; // normally in EVIS we have a mag2radius but we don't really need that here
		double RA = std::get<1>(starfield)[idx_star];
		double DEC = std::get<2>(starfield)[idx_star];

		// assign a 3D vector to the star at pseudo-infinite distance
		Vec3 star_pos_equ = Vec3(RA, DEC); // in case of doubt, yes, this takes RA and DEC in degrees
		
		SpiceDouble star_pos_equ_dbl[3] = { star_pos_equ.x, star_pos_equ.y, star_pos_equ.z };
		SpiceDouble star_pos_ecl_dbl[3];
		mxv_c(equ_ecl_rot, star_pos_equ_dbl, star_pos_ecl_dbl);

		Vec3 star_pos = Vec3(star_pos_ecl_dbl[0], star_pos_ecl_dbl[1], star_pos_ecl_dbl[2]);

		// do not do this! no need! star is at infinite distance anyway!!
		// Vec3 star_rel_pos = star_pos - cam_pos;

		if (star_pos.dot(cam_forward) > 0)
		{
			double px = f * star_pos.dot(cam_right) / star_pos.dot(cam_forward);
			double py = f * star_pos.dot(cam_up) / star_pos.dot(cam_forward);

			int pix_x = screen_x / 2 + px + 0.5;
			int pix_y = screen_y / 2 - py + 0.5;

			drawCircle(img, screen_x, screen_y, pix_x, pix_y, radius, {200, 200, 200});
		}
	}

	// ok, next thing, orbit ellipses!
	// minor planet orbit first
	//                                            vvv -- there is one less line than there are points
	for (int idx_op = 0; idx_op < mp_orbit.size() - 1; idx_op++)
	{
		Vec3 p1 = mp_orbit[idx_op];
		Vec3 p2 = mp_orbit[idx_op + 1];

		std::vector<int> p1_scrpos = space2screen(p1, cam_pos, cam_orient, f, screen_x, screen_y);
		std::vector<int> p2_scrpos = space2screen(p2, cam_pos, cam_orient, f, screen_x, screen_y);

		drawLine(img, screen_x, screen_y, p1_scrpos[0], p1_scrpos[1], p2_scrpos[0], p2_scrpos[1], { 0, 255, 0 });
	}

	// now the orbits of major planets (Sun orbit is not drawn, therefore index starts at 1)
	for (int idx_major = 1; idx_major < major_orbits.size(); idx_major++)
	{
		for (int idx_op = 0; idx_op < major_orbits[idx_major].size() - 1; idx_op++)
		{
			Vec3 p1 = major_orbits[idx_major][idx_op];
			Vec3 p2 = major_orbits[idx_major][idx_op + 1];

			std::vector<int> p1_scrpos = space2screen(p1, cam_pos, cam_orient, f, screen_x, screen_y);
			std::vector<int> p2_scrpos = space2screen(p2, cam_pos, cam_orient, f, screen_x, screen_y);

			drawLine(img, screen_x, screen_y, p1_scrpos[0], p1_scrpos[1], p2_scrpos[0], p2_scrpos[1], major_body_colors[idx_major]);
		}
	}

	// now draw the objects themselves
	// starting with the minor planet...
	std::vector<int> mp_scrpos = space2screen(mp_pos, cam_pos, cam_orient, f, screen_x, screen_y);
	if (!(mp_scrpos[0] == -1 && mp_scrpos[1] == -1))
	{
		drawCircle(img, screen_x, screen_y, mp_scrpos[0], mp_scrpos[1], 3);
	}

	// now the major bodies (this time including the Sun, of course)
	for (int idx_major = 0; idx_major < major_orbits.size(); idx_major++)
	{
		std::vector<int> mp_scrpos = space2screen(major_pos[idx_major], cam_pos, cam_orient, f, screen_x, screen_y);
		if (idx_major == 0)
		{
			// compute real angular size in pixels
			double ang_radius = asin(major_body_radii[idx_major] / (major_pos[idx_major] - cam_pos).mag());
			double pix_radius = f * tan(ang_radius);
			double draw_radius = max(5, pix_radius);
			drawCircle(img, screen_x, screen_y, mp_scrpos[0], mp_scrpos[1], draw_radius, major_body_colors[idx_major]);
		}
		else
		{
			double ang_radius = asin(major_body_radii[idx_major] / (major_pos[idx_major] - cam_pos).mag());
			double pix_radius = f * tan(ang_radius);
			double draw_radius = max(3, pix_radius);
			drawCircle(img, screen_x, screen_y, mp_scrpos[0], mp_scrpos[1], draw_radius, major_body_colors[idx_major]);
		}
	}

	drawText(img, screen_x, screen_y, 10, 10, st.datetime, { 255, 0, 0 });

	// now save it to file
	std::string outfilename = save_name;
	std::ofstream outfile(outfilename);

	outfile << "P3\n" << screen_x << " " << screen_y << "\n255\n";

	for (int y = 0; y < screen_y; ++y)
	{
		for (int x = 0; x < screen_x; ++x)
		{
			std::array<int, 3> color = img[y * screen_x + x];
			outfile << color[0] << " " << color[1] << " " << color[2] << " ";
		}
		outfile << "\n";
	}

	outfile.close();
}

// don't ask
double getNextLargerOrRetain(double value, const std::vector<double>& sorted_arr)
{
	std::vector<double>::const_iterator it = std::upper_bound(sorted_arr.begin(), sorted_arr.end(), value);

	if (it != sorted_arr.end())
	{
		return *it;
	}
	else
	{
		return value;
	}
}

// s, starfield, cam_mode, cam_dist, cam_theta, cam_phi, fov, map_name
void mapSS3D(State st, std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> starfield,
	std::string cam_mode, double cam_dist, double cam_theta, double cam_phi, double fov_deg,
	std::string center_obj, std::string carrier_obj,
	std::string map_name, int screen_x, int screen_y)
{
	double fov = deg2rad(fov_deg);

	SpiceDouble et;
	utc2et_c(st.datetime.c_str(), &et);

	// get planet positions
	StateMatrix SolarSystemState = getSolarSystemStates(et);
	// access is StateMatrix[planet idx][pos/vel idx][vector component (x,y,z) idx]

	// convert major body state vectors to J2000 ecliptic version (rather than standard equatorial J2000)
	SpiceDouble rotate[3][3];
	pxform_c("J2000", "ECLIPJ2000", et, rotate);

	std::vector<Vec3> major_pos_eclip;
	std::vector<Vec3> major_vel_eclip;

	for (int idx_major = 0; idx_major < SolarSystemState.size(); idx_major++)
	{
		SpiceDouble equ_pos[3] = { SolarSystemState[idx_major][0][0], SolarSystemState[idx_major][0][1], SolarSystemState[idx_major][0][2] };
		SpiceDouble ecl_pos[3];
		mxv_c(rotate, equ_pos, ecl_pos);

		SpiceDouble equ_vel[3] = { SolarSystemState[idx_major][1][0], SolarSystemState[idx_major][1][1], SolarSystemState[idx_major][1][2] };
		SpiceDouble ecl_vel[3];
		mxv_c(rotate, equ_vel, ecl_vel);

		Vec3 new_pos = Vec3(ecl_pos[0], ecl_pos[1], ecl_pos[2]);
		Vec3 new_vel = Vec3(ecl_vel[0], ecl_vel[1], ecl_vel[2]);
		major_pos_eclip.push_back(new_pos);
		major_vel_eclip.push_back(new_vel);
	}

	// also convert minor planet state
	SpiceDouble mp_equ_pos[3] = { st.p.x, st.p.y, st.p.z };
	SpiceDouble mp_ecl_pos[3];
	mxv_c(rotate, mp_equ_pos, mp_ecl_pos);

	SpiceDouble mp_equ_vel[3] = { st.v.x, st.v.y, st.v.z };
	SpiceDouble mp_ecl_vel[3];
	mxv_c(rotate, mp_equ_vel, mp_ecl_vel);

	Vec3 mp_pos = Vec3(mp_ecl_pos[0], mp_ecl_pos[1], mp_ecl_pos[2]);
	Vec3 mp_vel = Vec3(mp_ecl_vel[0], mp_ecl_vel[1], mp_ecl_vel[2]);

	// get sampled two-body ellipse for the minor planet
	std::vector<Vec3> mp_orbit = getKeplerOrbitPoints(mp_pos, mp_vel);

	// get them for major bodies too
	std::vector<std::vector<Vec3>> major_orbits;
	for (int idx_major = 0; idx_major < SolarSystemState.size(); idx_major++)
	{
		if (idx_major < 3) // having vectors relative to Sun instead of the barycenter makes some less wobbly
		{
			major_orbits.push_back(getKeplerOrbitPoints(major_pos_eclip[idx_major] - major_pos_eclip[0], major_vel_eclip[idx_major] - major_vel_eclip[0]));
		}
		else
		{
			major_orbits.push_back(getKeplerOrbitPoints(major_pos_eclip[idx_major], major_vel_eclip[idx_major]));
		}
	}

	// now we can draw images
	// ========== TOP-DOWN ==========

	// get the extents of the orbit of the minor planet
	double R_mp_max = 0;
	for (int idx_op = 0; idx_op < mp_orbit.size(); idx_op++)
	{
		double Rsq_current = mp_orbit[idx_op].x * mp_orbit[idx_op].x + mp_orbit[idx_op].y * mp_orbit[idx_op].y;
		if (Rsq_current > R_mp_max * R_mp_max)
		{
			R_mp_max = sqrt(Rsq_current);
		}
	}

	// we will push the camera as far back to include the next planet's orbit (unless the minor planet's orbit 
	// is larger than Neptune's, in which case we will go even farther out)
	std::vector<double> planet_sma = { 69.8e6, 108.9e6, 152.1e6, 249.3e6, 816.4e6, 1506.5e6, 3001.4e6, 4558.9e6 };
	double R_max = getNextLargerOrRetain(R_mp_max, planet_sma) * 1.33;
	double fit_dist = R_max / tan(fov / 2);

	Vec3 cam_pos = Vec3(0, 0, fit_dist);
	std::vector<Vec3> cam_orient = {
		Vec3(1, 0, 0),
		Vec3(0, 1, 0),
		Vec3(0, 0, 1)
	};

	std::string save_name = "map_topdown/" + map_name + "_topdown.ppm";
	renderSolarSystem(st, mp_pos, mp_orbit, major_pos_eclip, major_orbits, cam_mode, fov, screen_x, screen_y, cam_pos, cam_orient, starfield, save_name);

	// ========== EDGE-ON ==========
	cam_pos = Vec3(fit_dist, 0, 0);
	cam_orient = {
		Vec3(0, 1, 0),
		Vec3(0, 0, 1),
		Vec3(1, 0, 0)
	};

	save_name = "map_edgeon/" + map_name + "_edgeon.ppm";
	renderSolarSystem(st, mp_pos, mp_orbit, major_pos_eclip, major_orbits, cam_mode, fov, screen_x, screen_y, cam_pos, cam_orient, starfield, save_name);

	// ========== CUSTOM ==========
	cam_pos = -Vec3(cam_theta, cam_phi) * cam_dist;

	if (strcmp(carrier_obj.c_str(), "None")) // NOT equal to "None"
	{
		if (!strcmp(carrier_obj.c_str(), "SOLAR_SYSTEM_BARYCENTER"))
		{
			cam_pos = Vec3(0, 0, 0);
		}
		else if (!strcmp(carrier_obj.c_str(), "SUN"))
		{
			cam_pos = major_pos_eclip[0];
		}
		else if (!strcmp(carrier_obj.c_str(), "MERCURY_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[1];
		}
		else if (!strcmp(carrier_obj.c_str(), "VENUS_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[2];
		}
		else if (!strcmp(carrier_obj.c_str(), "EARTH_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[3];
		}
		else if (!strcmp(carrier_obj.c_str(), "MARS_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[4];
		}
		else if (!strcmp(carrier_obj.c_str(), "JUPITER_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[5];
		}
		else if (!strcmp(carrier_obj.c_str(), "SATURN_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[6];
		}
		else if (!strcmp(carrier_obj.c_str(), "URANUS_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[7];
		}
		else if (!strcmp(carrier_obj.c_str(), "NEPTUNE_BARYCENTER"))
		{
			cam_pos = major_pos_eclip[8];
		}
		else if (!strcmp(carrier_obj.c_str(), "MP"))
		{
			cam_pos = mp_pos;
		}
	}

	Vec3 target_pos = Vec3(0, 0, 0);
	if (!strcmp(center_obj.c_str(), "SOLAR_SYSTEM_BARYCENTER"))
	{
		target_pos = Vec3(0, 0, 0);
	}
	else if (!strcmp(center_obj.c_str(), "SUN"))
	{
		target_pos = major_pos_eclip[0];
	}
	else if (!strcmp(center_obj.c_str(), "MERCURY_BARYCENTER"))
	{
		target_pos = major_pos_eclip[1];
	}
	else if (!strcmp(center_obj.c_str(), "VENUS_BARYCENTER"))
	{
		target_pos = major_pos_eclip[2];
	}
	else if (!strcmp(center_obj.c_str(), "EARTH_BARYCENTER"))
	{
		target_pos = major_pos_eclip[3];
	}
	else if (!strcmp(center_obj.c_str(), "MARS_BARYCENTER"))
	{
		target_pos = major_pos_eclip[4];
	}
	else if (!strcmp(center_obj.c_str(), "JUPITER_BARYCENTER"))
	{
		target_pos = major_pos_eclip[5];
	}
	else if (!strcmp(center_obj.c_str(), "SATURN_BARYCENTER"))
	{
		target_pos = major_pos_eclip[6];
	}
	else if (!strcmp(center_obj.c_str(), "URANUS_BARYCENTER"))
	{
		target_pos = major_pos_eclip[7];
	}
	else if (!strcmp(center_obj.c_str(), "NEPTUNE_BARYCENTER"))
	{
		target_pos = major_pos_eclip[8];
	}
	else if (!strcmp(center_obj.c_str(), "MP"))
	{
		target_pos = mp_pos;
	}

	// dummy default orientation
	cam_orient = {
		Vec3(1, 0, 0),
		Vec3(0, 1, 0),
		Vec3(0, 0, 1)
	};

	Vec3 forward = (target_pos - cam_pos).normalized();
	Vec3 forward_xy = Vec3(forward.x, forward.y, 0).normalized();
	Vec3 right = Vec3(-forward_xy.y, forward_xy.x, 0);
	Vec3 up = forward.cross(right).normalized();
	right = up.cross(forward).normalized();

	cam_orient[0] = right;
	cam_orient[1] = up;
	cam_orient[2] = -forward;

	save_name = "map_custom/" + map_name + "_custom.ppm";
	renderSolarSystem(st, mp_pos, mp_orbit, major_pos_eclip, major_orbits, cam_mode, fov, screen_x, screen_y, cam_pos, cam_orient, starfield, save_name);
}

void printHelpMsg()
{
	std::cout << "SVIS Help\n\n";

	std::cout << "SVIS is a Heliocentric orbit visualization tool. It generates 3D images using state vectors of a "
		<< "Heliocentric object along with the Sun and 8 planets, optionally with background stars.\n\n";

	std::cout << "It requires SPICE kernels, a state vector file generated by SPRO, and optionally a star catalog.\n\n";

	std::cout << "The minimial invocation is merely 'svis' - it assumes the following defaults:\n";
	std::cout << "    State vector file name: state_vectors.txt\n";
	std::cout << "    SPICE kernels path: data/SPICE/\n";
	std::cout << "    Star catalog path: data/Tycho2.csv\n";
	std::cout << "    Field of view: 60 deg\n";
	std::cout << "    Screen size: 640 x 480\n";
	std::cout << "    Camera target: SOLAR_SYSTEM_BARYCENTER\n";
	std::cout << "    Custom cam. carrier: None (fixed camera)\n";
	std::cout << "    Custom cam. RA: 45 deg\n";
	std::cout << "    Custom cam. DEC: +45 deg\n";
	std::cout << "    Custom cam. distance: 15 AU\n";
	std::cout << "    Image output file prefix: 'map_'\n\n";

	std::cout << "You can adjust each setting by using the following arguments:\n";
	std::cout << "    -sv: state vector file path\n";
	std::cout << "    -spice: SPICE kernels directory path\n";
	std::cout << "    -catalog: Star catalog file path (enter 'None' for no background stars)\n";
	std::cout << "    -fov: Field of view in degrees\n";
	std::cout << "    -center: The object the camera is targeting\n";
	std::cout << "    -carrier: The object that the camera is travelling with (leave blank or enter 'None' for a camera fixed in space)\n";
	std::cout << "    -theta: Target RA in degrees if a carrier doesn't exist\n";
	std::cout << "    -phi: Target DEC in degrees if a carrier doesn't exist\n";
	std::cout << "    -dist: Cam. distance from centered object if a carrier doesn't exist\n";
	std::cout << "    -prefix: Image file output prefix\n\n";

	std::cout << "SVIS always outputs two default maps: the top-down and edge-on view maps of the Solar System. The camera settings only affect a third map called the 'custom' map.\n";
	std::cout << "Output images will be saved on the corresponding directories: map_topdown, map_edgeon, map_custom.\n\n";

	std::cout << "SVIS was developed by H. A. Guler.\n";
	std::cout << "SVIS is licensed under GNU General Public License version 2.0 (GPL-2.0 License)\n\n";
}

int main(int argc, char* argv[])
{
	std::cout << "SVIS v0.2.0\n\n";

	// default parameters
	std::string sv_path = "state_vectors.txt";
	std::string spice_path = "data/SPICE/";

	std::string center_obj = "SOLAR_SYSTEM_BARYCENTER";
	std::string carrier_obj = "None"; // the object which the camera is attached to

	std::string cam_mode = "p"; // p for perspective, o for orthogonal
	double cam_dist = 15 * AU;
	double cam_theta = deg2rad(45); // camera position right ascension
	double cam_phi = deg2rad(45); // camera position declination
	double fov = 60; // field-of-view for perspective projection

	int screen_x = 640;
	int screen_y = 480;
	
	std::string out_prefix = "map_";
	std::string starcatalog_path = "data/Tycho2.csv";

	// handle command line arguments
	// there is a more compact version of doing this but this is easier for my brain
	int argtype = 0;
	for (int idx_cmd = 1; idx_cmd < argc; idx_cmd++)
	{
		if (!strcmp(argv[idx_cmd], "-sv") || !strcmp(argv[idx_cmd], "-state_vectors")) // also the default argument
		{
			argtype = 0;
		}
		else if (!strcmp(argv[idx_cmd], "-spice"))
		{
			argtype = 1;
		}
		else if (!strcmp(argv[idx_cmd], "-center"))
		{
			argtype = 2;
		}
		else if (!strcmp(argv[idx_cmd], "-carrier"))
		{
			argtype = 3;
		}
		else if (!strcmp(argv[idx_cmd], "-mode"))
		{
			argtype = 4;
		}
		else if (!strcmp(argv[idx_cmd], "-dist"))
		{
			argtype = 5;
		}
		else if (!strcmp(argv[idx_cmd], "-theta"))
		{
			argtype = 6;
		}
		else if (!strcmp(argv[idx_cmd], "-phi"))
		{
			argtype = 7;
		}
		else if (!strcmp(argv[idx_cmd], "-fov"))
		{
			argtype = 8;
		}
		else if (!strcmp(argv[idx_cmd], "-prefix"))
		{
			argtype = 9;
		}
		else if (!strcmp(argv[idx_cmd], "-catalog"))
		{
			argtype = 10;
		}
		else if (!strcmp(argv[idx_cmd], "-screen_x"))
		{
			argtype = 11;
		}
		else if (!strcmp(argv[idx_cmd], "-screen_y"))
		{
			argtype = 12;
		}
		else if (!strcmp(argv[idx_cmd], "-h") || !strcmp(argv[idx_cmd], "--help")) // two dashes because people are more used to it
		{
			printHelpMsg();
			return 0;
		}
		else
		{
			switch (argtype)
			{
			case 0:
				sv_path = argv[idx_cmd];
				break;
			case 1:
				spice_path = argv[idx_cmd];
				break;
			case 2:
				center_obj = argv[idx_cmd];
				break;
			case 3:
				carrier_obj = argv[idx_cmd];
				break;
			case 4:
				cam_mode = argv[idx_cmd];
				break;
			case 5:
				cam_dist = strtod(argv[idx_cmd], NULL) * AU; // convert km to AU
				break;
			case 6:
				cam_theta = deg2rad(strtod(argv[idx_cmd], NULL));
				break;
			case 7:
				cam_phi = deg2rad(strtod(argv[idx_cmd], NULL));
				break;
			case 8: 
				fov = strtod(argv[idx_cmd], NULL);
				break;
			case 9:
				out_prefix = argv[idx_cmd];
				break;
			case 10:
				starcatalog_path = argv[idx_cmd];
				break;
			case 11:
				screen_x = atoi(argv[idx_cmd]);;
				break;
			case 12:
				screen_y = atoi(argv[idx_cmd]);;
			}
		}
	}

	std::cout << "Loading SPICE kernels... ";
	loadAllKernels(spice_path);
	std::cout << "Done.\n";

	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> starfield;

	if (strcmp(starcatalog_path.c_str(), "None"))
	{
		std::cout << "Reading star catalogue... ";
		starfield = readTycho2(starcatalog_path);
		std::cout << "Done.\n";
	}
	else
	{
		std::cout << "No star catalog provided, skybox will be empty.\n";
	}

	std::cout << "Reading state vector data... ";
	std::vector<State> states = readStateVectorFile(sv_path);
	std::cout << "Done.\n";

	std::cout << "Mapping the Solar System...\n";
	createDirectoryIfNotExists("map_topdown");
	createDirectoryIfNotExists("map_edgeon");
	createDirectoryIfNotExists("map_custom");
	// sanitize ephemeris point data and generate an image for each ephemeris point
	for (int idx_state = 0; idx_state < states.size(); idx_state++)
	{
		std::cout << "    Map " << idx_state + 1 << " / " << states.size() << "...\n";
		State s = states[idx_state];

		// get ephem time from UTC string
		SpiceDouble et;
		utc2et_c(s.datetime.c_str(), &et);

		std::string suffix = s.datetime;
		std::replace(suffix.begin(), suffix.end(), ':', '_'); // keep the OS happy
		std::string map_name = "map_" + suffix;
		mapSS3D(s, starfield, cam_mode, cam_dist, cam_theta, cam_phi, fov, center_obj, carrier_obj, map_name, screen_x, screen_y);
	}
	std::cout << "Done generating charts.\n";

	std::cout << "Program end.\n\n";
	return 0;
}
