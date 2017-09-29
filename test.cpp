// Filename: test.cpp
#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;

class Camera {
public:
 double pos[3] = {0, 0, 0};
 // The position of the camera.
 double angle[3] = {0, 0, 0};
 // The angle that the camera is pointing in, specified in roll, pitch, and yaw.
 // I would initialize these variables in the constructor, but I've declared them up here for other functions to also use, so I have to define something up here anyway.
 Camera() {

 }
 Camera(double *pos) {
	for (int i = 0; i < sizeof(*pos)/sizeof(double); i++) {
		this->pos[i] = pos[i];
	}
 }
 Camera(double *pos, double *angle, int posElements, int normalElements) {
 	for (int i = 0; i < sizeof(*pos)/sizeof(double); i++) {
		this->pos[i] = pos[i];
	}
	for (int i = 0; i < sizeof(*angle)/sizeof(double); i++) {
		this->angle[i] = angle[i];
	}
 }
};
// There should only ever be one camera, unless you want to use multiple screens.
class Viewport {
public:
 double size = 1;
 Camera *cam;
 double initpos[4][3] = {{-1, 1, 1},
{1, 1, 1},
{-1, -1, 1},
{1, -1, 1}};
 double pos[4][3] = {{-1, 1, 1},
{1, 1, 1},
{-1, -1, 1},
{1, -1, 1}};
 // By default, the Viewport should be a rectangle.
 void setVPosToCam () {
 	// This function should only be called if this class has a cam pointer.
	// Specifying x y and z of the camera: cam->pos[0-2]
	// Specying roll, pitch, and yaw of the camera: cam->angle[0-2]
	// Now, I need to rotate each of the points around each of the angles:
	// roll, around the x pole at the camera's position,
	// yaw, around the y pole at the camera's position,
	// and pitch, around the z pole at the camera's position.
	// Now, let's go through the entire array and set the coordinates of this rectangle.
	// First, scale, then, rotate, then add to pos of cam.
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			pos[i][j] = initpos[i][j] * size;
			cout << "initpos[i][j] * size = " << initpos[i][j] * size << endl;
		}
	}
	double magnitude[sizeof(pos)/sizeof(pos[0])] = {0};
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			magnitude[i] += pos[i][j]*pos[i][j];
		}
		magnitude[i] = sqrt(magnitude[i]);
	}
	// To rotate around just one plane, we would need a flat magnitude every time.
	// And unfortunately, since I decided to do it in degrees instead of direction vectors (omg how do you do it in direction vectors), I need to rotate on flat planes, one for each pair of dimensions.
	// This is convenient in 2D and 3D; you have one pair of dimensions, and so need one plane to rotate on, and one measure of degrees. In 3D, you have 3 pairs of dimensions (xy, yz, zx), and three measures of degrees. In 4D however, there are 6 degrees. xy, xz, xm, yz, ym, zm.
	// So, this means: I need to calculate all of the pairs of variables in cam->pos , I need to make sure that there are that many variables in cam->degrees, then I need to calculate a pair's flat magnitude, then I need to set those two coordinates to magnitude[i] * sin(degrees[i]) and magnitude[i] * cos(degrees[i]) respectively. Uh, how do I calculate n choose two again?
	// To iterate through all of the pairs: iterate through two for loops of i and j < the number of dimensions in pos, and don't let j get bigger than or equal to i.
	int coordinatePair = 0;
	for (int i = 0; i < sizeof(pos[0])/sizeof(double); i++) {
		for (int j = 0; j < i; j++) {
			coordinatePair++;
		}
	}
	// The above will iterate over all of the unique coordinates in a ordered double/triple/etc, and make two arrays with that size.
	double flatMagnitude[coordinatePair];
	double flatAngle[coordinatePair];
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			cout << ", " << pos[i][j] << flush;
		}
		cout << endl;
	}
	// The above will then go through those pairs and print them.
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		coordinatePair = 0;	
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			for (int k = 0; k < j; k++) {
				cout << "i is: " << i << ", j is: " << j << ", k is: " << k << ", coordinatePair is: " << coordinatePair << ", cam->angle[coordinatePair] is: " << cam->angle[coordinatePair] << flush;
				flatMagnitude[coordinatePair] = sqrt((pos[i][j]*pos[i][j])+(pos[i][k]*pos[i][k]));
				// The flat magnitude is the magnitude of a vector based on one coordinate pair.
				flatAngle[coordinatePair] = atan2(pos[i][j],pos[i][k]);
				// We need to find the flat angle of the coordinate pairs and then add them to cam->angle, so that all of the points don't end up pointing in the same direction.
				pos[i][k] = cos(cam->angle[coordinatePair]+flatAngle[coordinatePair])*flatMagnitude[coordinatePair];
				pos[i][j] = sin(cam->angle[coordinatePair]+flatAngle[coordinatePair])*flatMagnitude[coordinatePair];
				// If you multiply the flat magnitude by the sin of the degree of rotation, you can get the y coordinate relative to those two coordinates, and if you multiply by the cos you get the x.
				cout << ", pos[" << i << "][" << k << "] is: " << pos[i][k] << ", pos[" << i << "][" << j << "] is: " << pos[i][j] << ", flatMagnitude[" << coordinatePair << "] is: " << flatMagnitude[coordinatePair] << ", flatAngle[" << coordinatePair << "] is: " << flatAngle[coordinatePair] << endl;
				coordinatePair++;
				// The above will rotate each of the 4 points in the rectangle, based on each unique coordinate pair and their given angle (roll, pitch, yaw, etc.) sequentially.
	
			}
		}
	}
	// The result of all of this coordinate pairing is that in 3D, there will be 3 different planes to rotate on: the xy plane, the xz plane, and the yz plane; we will take the magnitude of the 4 points of the rectangle, and then rotate them 3 times according to each of those angles.
	cout << "sizeof(pos)/sizeof(pos[0]): " << sizeof(pos)/sizeof(pos[0]) << endl;
	cout << "sizeof(pos[0])/sizeof(double): " << sizeof(pos[0])/sizeof(double) << endl;
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			pos[i][j] += cam->pos[j];
			// This will move the scaled and rotated viewport to where the camera is.
		}
	}
	
 }
 Viewport () {
	size = 1;
	// If a camera is not specified, nor the size, the viewport should have a size of one, and be pointing toward the x axis.
	// This also means that the default values of pos are okay.
 }
 Viewport (Camera *cam) {
 	this->cam = cam;
	// If a camera is specified, then set the cam pointer to it.
 	size = 1;
	// If the size is not specified, then the viewport should have a size of one.
	setVPosToCam();
	// I decided to move the code that changes the viewport's position to its own function, because I will want to be calling it from two of the constructors, and possibly every single game frame, depending on the viewport.
 }
 Viewport (Camera *cam, double size) {
 	this->cam = cam;
	// If a camera is specified, then set the cam pointer to it.
	this->size = size;
	// If a size is specified, then set the size of the class to it.
	setVPosToCam();
 }
};
double *findIntersect(Viewport *view, double vertice[], double cameraPos[], int arrayIndex) {
	cout << "Input vertice[]: " << endl;
	for (int i = 0; i < arrayIndex; i++) {
		cout << "vertice[" << i << "] is: " << vertice[i] << endl;
	}
	// To find the equation of a line, you need two points.
	// First, take p1p2->
	double v[arrayIndex];
	for (int i = 0; i < arrayIndex; i++) {
		v[i] = vertice[i] - cameraPos[i];
	}
	// Then, the parametric equations are:
	// x = v[0]*t + cameraPos[0]
	// y = v[1]*t + cameraPos[1]
	// z = v[2]*t + cameraPos[2]
	// And so on. I don't make variables x, y, and z in this function because all I need the parametric equations for are to solve for t.
	
	// To find the equation of a plane, you need three points.
	// We have three points in *view.
	// First, we need to make two vectors out of the differences of three points (the two vectors being perpendicular to the plane with the three points).
	double v2[arrayIndex];
	double v3[arrayIndex];
	for (int i = 0; i < arrayIndex; i++) {
		v2[i] = view->pos[0][i] - view->pos[1][i];
	}
	for (int i = 0; i < arrayIndex; i++) {
		v3[i] = view->pos[1][i] - view->pos[2][i];
	}
	// Now, we have to find the coefficients of a, b, and c.
	
	/*
	 *|  a     b     c   |
	 *|v2[0] v2[1] v2[2] |
	 *|v3[0] v3[1] v3[2] |
	*/
	// That will return a vector that is perpendicular to both v2 and v3.
	// Not sure how to find cross product of two 3D vectors or up.
	double a = (v2[1]*v3[2]) - (v2[2]*v3[1]);
	double b = (v2[0]*v3[2]) - (v2[2]*v3[0]);
	double c = (v2[0]*v3[1]) - (v2[1]*v3[0]);

	// Then, we use another point, which will also be on the plane, to find D.
	double d = (c*view->pos[0][0]) + (b*view->pos[0][1]) + (c*view->pos[0][2]);

	// The equation of the plane will come out to be ax + by + cz = d
	// Let's plug in the parametric equations for the line and solve for t.
	// a(v[0]*t + cameraPos[0]) + b(v[1]*t + cameraPos[1]) + c(v[2]*t + cameraPos[2]) = d
	// t = (d - (a * cameraPos[0] + b * cameraPos[1] + c * cameraPos[2]))/(a*v[0] + b*v[1] + c*v[2])
	double t;
	t = (d - (a * cameraPos[0] + b * cameraPos[1] + c * cameraPos[2]))/(a*v[0] + b*v[1] + c*v[2]);
	double intersect[arrayIndex];
	for (int i = 0; i < arrayIndex; i++) {
		intersect[i] = (v[i]*t) + cameraPos[i];
		cout << "i in intersect[i] is: " << i << ", v[i] is: " << v[i] << ", t is: " << t << ", and intersect[" << i << "] is: " << intersect[i] << endl;
	}
	// So, for some reason, this is returning ludicrous values, and only one coordinate of them. Trying to manually set the counter to 3 will return three coordinates of inf, inf, inf.
	return intersect;
}
double* pointsToMatrix(double v[4][3]) {
	// If we get the equation for the matrix, we need to solve for 4 points on the matrix, so we might as well just start with the 4 points from viewport.pos[][]
	// pos[0][] will map to 0,0
	// pos[1][] will map to W_WIDTH,0
	// pos[2][] will map to W_WIDTH,W_HEIGHT
	// pos[3][] will map to 0,W_HEIGHT
	// Either way, we need to calculate once, and only once (per change of the viewport's angle), the matrix that will map points from a plane to a cartesian grid.
	// With the matrix, we will only need to apply the matrix to any position vector in space to get where it would be "on the screen" that is represented by the viewport.
	// The use of the matrix will look like this:
	/*
	 *					|0,0	1,0|
	 *	|0,0	1,0	2,0|	  X	|0,1	1,1|	=	|0,0	1,0|
	 *					|0,2	1,2|
	 *
	*/
	// With the first matrix being the input position, the second matrix being the transformation matrix, and the third matrix being the new coordinates.
	// So let's make a matrix like the one in the middle.
	transformationMatrix = new transformationMatrix[2][3];
	// Now, we need to figure out how to calculate this.
	// Done?:
	// v1x*a + v1y*b + v1c*z = 0
	// v2x*a + v2y*b + v2c*z = W_WIDTH
	// v3x*a + v3y*b + v3c*z = W_WIDTH
	// and
	// v1x*d + v1y*e + v1c*f = 0
	// v2x*d + v2y*e + v2c*f = 0
	// v3x*d + v3y*e + v3c*f = W_HEIGHT 
	// , from the matrices above.
	// We will already have: v1 and v2 values. We will already have 100% of the information about rotation and x+y scale just by where v2 and v3 is relative to v1.
	// So now we just solve for the intersection of those 3 planes to get a, b, c, and then we solve the second equation triplet to get d, e, and f.
	return transformationMatrix;
}
bool isIntersectInShape(Viewport *view, double *intersect) {
	// Because I used the definition for a plane when calculating the intersection, checking if the array intersect is on a plane containing 3 points of the viewport would be redundant. This means that when calculating if a point is in a rectangle, I could do a simpler check to see if it is in a rectangular prism, instead.
	// The viewport is always going to be a rectangle for now, since pentagon screens and 3D displays don't exist yet.
	// https://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
	// Figure out that one simple vector answer in the above link and then explain it so that my future self can understand it.
	// If the vector answer doesn't work out, we could simply make a transformation that maps point 1 of the rectangle to 0,0 on a coordinate plane and point 3 to 1,1 on the coordinate plane, and then apply that transformation to the intersect point, and see if 1 > x > 0 and 1 > y > 0
}
int main() {
	cout << "Hello World!" << endl;
	double camPos[] = {3.4, 2.1, 4};
	double camDegree[] = {8, 2, 1};
	Camera camera(camPos, camDegree, 3, 3);
	for (int i = 0; i < sizeof(camera.pos)/sizeof(double); i++) {
		cout << camera.pos[i] << endl;
	}
	for (int i = 0; i < sizeof(camera.angle)/sizeof(double); i++) {
		cout << camera.angle[i] << endl;
	}
	Viewport view(&camera);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			cout << view.pos[i][j] << "," << flush;
		}
		cout << endl;
	}
	camera.angle[0] = 14;
	camera.angle[1] = -23;
	camera.angle[2] = 33;
	view.size = 10;
	view.setVPosToCam();
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			cout << view.pos[i][j] << "," << flush;
		}
		cout << endl;
	}
	Viewport *pView = &view;
	double vertice[3] = {34, 1, -4};
	double *intersect = findIntersect(pView, vertice, camera.pos, sizeof(camPos)/sizeof(double));
	cout << "done" << endl;
	return 0;
}
