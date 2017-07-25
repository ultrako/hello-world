// Filename: test.cpp
#include <iostream>
#include <stdio.h>
using namespace std;

class Camera {
public:
 double pos[3] = {0, 0, 0};
 // The position of the camera.
 double normal[3] = {1, 0, 0};
 // The direction that the camera is pointing in.
 // It cannot be 0, 0, 0, because the camera has to point _somewhere_.
 // I would initialize these variables in the constructor, but I've declared them up here for other functions to also use, so I have to define something up here anyway.
 Camera() {

 }
 Camera(double *pos, int posElements) {
	for (int i = 0; i < posElements; i++) {
		this->pos[i] = pos[i];
	}
 }
 Camera(double *pos, double *normal, int posElements, int normalElements) {
 	for (int i = 0; i < posElements; i++) {
		this->pos[i] = pos[i];
	}
	for (int i = 0; i < normalElements; i++) {
		this->normal[i] = normal[i];
	}
 }
};
// There should only ever be one camera, unless you want to use multiple screens.
class Viewport {
public:
 double size = 1;
 Camera *cam;
 double pos[4][3] = {{-1, 1, 1},
{1, 1, 1},
{-1, -1, 1},
{1, -1, 1}};
 // By default, the Viewport should be a rectangle.
 void setVPosToCam () {
 	// This function should only be called if this class has a cam pointer.
	// Specifying x y and z of the camera: cam->pos[0-3]
	// The coordinates of the center of the rectangle should be: camera pos + (vector that is a multiple of cam->normal but has magnitude of size).
	// The coordinates of the rectangle's points should be the center of the rectangle plus or minus size in x and y relative to a xy plane that has a normal of cam->normal and a point on it that is the center of the rectangle.
	// I'm sorry, myself. If only I could word the above better, then I could just code it.
	// Now, let's iterate through the entire array and set the coordinates of this rectangle.
	// Err, I forgot how to iterate through multi-dimensioned arrays. My b, will have to look through my example codes that are back home.
 }
 Viewport () {
	size = 1;
	// If a camera is not specified, nor the size, the viewport should have a size of one, and be pointing toward the x axis.
	// This also means that the default values of pos are okay.
 }
 Viewport (&Camera) {
 	cam = &Camera;
	// If a camera is specified, then set the cam pointer to it.
 	size = 1;
	// If the size is not specified, then the viewport should have a size of one.
	initViewPort();
	// I decided to move the code that changes the viewport's position to its own function, because I will want to be calling it from two of the constructors, and possibly every single game frame, depending on the viewport.
 }
 Viewport (&Camera, double size) {
 	Camera *cam = &Camera;
	// If a camera is specified, then set the cam pointer to it.
	this->size = size;
	// If a size is specified, then set the size of the class to it.
	initViewPort();
 }
};
int* findIntersect(&Viewport, double vertice[], double cameraPos[]) {
	// Maths to find where a line from a vertice of a polygon and the camera meets a plane and return a pointer to an array with the position of that intersect. I figured it out once on paper and I just need to nag someone at home to screenshot that paper and send it to me.
	return 0;
}
int main() {
	cout << "Hello World!" << endl;
	double camPos[] = {3.4, 2.1, 4};
	double camNormal[] = {8, 2, 1};
	Camera camera(camPos, camNormal, 3, 3);
	for (int i = 0; i < sizeof(camera.pos)/sizeof(double); i++) {
		cout << camera.pos[i] << endl;
	}
	for (int i = 0; i < sizeof(camera.normal)/sizeof(double); i++) {
		cout << camera.normal[i] << endl;
	}
	return 0;
}
