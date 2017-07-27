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
 Camera(double *pos, int posElements) {
	for (int i = 0; i < posElements; i++) {
		this->pos[i] = pos[i];
	}
 }
 Camera(double *pos, double *angle, int posElements, int normalElements) {
 	for (int i = 0; i < posElements; i++) {
		this->pos[i] = pos[i];
	}
	for (int i = 0; i < normalElements; i++) {
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
	double flatMagnitude[coordinatePair+1];
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++) {
		coordinatePair = 0;	
		for (int j = 0; j < sizeof(pos[0])/sizeof(double); j++) {
			for (int k = 0; k < j; k++) {
				flatMagnitude[coordinatePair] = sqrt((pos[i][j]*pos[i][j])+(pos[i][k]*pos[i][k]));
				pos[i][j] = cos(cam->angle[coordinatePair])*flatMagnitude[coordinatePair];
				pos[i][j] = sin(cam->angle[coordinatePair])*flatMagnitude[coordinatePair];
				coordinatePair++;
				// The above will apply the rotation transformations for all of the coordinate pair rotations. It's both disgusting and beautiful, isn't it?
			}
		}
	}
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
int *findIntersect(Viewport *view, double vertice[], double cameraPos[]) {
	// Maths to find where a line from a vertice of a polygon and the camera meets a plane and return a pointer to an array with the position of that intersect. I figured it out once on paper and I just need to nag someone at home to screenshot that paper and send it to me.
	return 0;
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
	Camera *cam = &camera;
	Viewport view(cam);
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
	return 0;
}
