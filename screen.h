#ifndef SCREEN_H_
#define SCREEN_H_

#include <iostream>
#include <SDL2/SDL.h>
using namespace std;

class Screen {
public:
 int W_WIDTH = 600;
 int W_HEIGHT = 400;
private:
 SDL_Window *window;
 SDL_Renderer *renderer;
 SDL_Texture *texture;
 Uint32 *buffer;
public:
 Screen(): window(NULL), renderer(NULL), texture(NULL), buffer(NULL) {

 }
 bool init() {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		// This will run SDL_Init, and if it fails, it will return a value < 0 , so the program ends.
		cout << "SDL init failed." << endl;
		return 1;
	}
	cout << "SDL init succeeded." << endl;
	window = SDL_CreateWindow("WindowTitle", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, W_WIDTH, W_HEIGHT, SDL_WINDOW_SHOWN);
	// We now create a pointer to a window that is equal to the pointer returned by the function SDL_CreateWindow, which makes a window and returns a pointer to that window.
	if (window == NULL) {
		// If creating the window fails, however, we need to crash our program.
		cout << "Could not create window." << endl;
		SDL_Quit();
		return 2;
	}
	cout << "Window successfully created." << endl;
	buffer = new Uint32[W_WIDTH*W_HEIGHT];
	// We are going to make a buffer (area in memory) that can hold all of the pixels on the screen.
	// A Uint32 is a standardized size for an int by SDL, that is the same size as a pixel.
	// We are going to make an array of ints that can store every single pixel on our screen.
	renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_PRESENTVSYNC);
	// Now, we make a renderer, which can draw textures to the window.
	texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STATIC, W_WIDTH, W_HEIGHT);
	// Now, we make a texture, so that the renderer can render it to the window.
	
	if (renderer == NULL) { cout << "Render init failed!" << " SDL error: " << SDL_GetError() << endl; SDL_DestroyWindow(window); SDL_Quit(); return 3; } 
	// If we failed to create a renderer, we should return that the initialization of the render failed.
	if (texture == NULL) { cout << "Texture init failed!" << endl; SDL_DestroyWindow(window); SDL_DestroyRenderer(renderer); SDL_Quit(); return 4; }
	// Same goes for the texture.
	else { cout << "Render and texture init succeeded." << endl; }
	// The object event is NULL until an event happens, in which case it will be equal to that event.
	
	memset(buffer, 0, W_WIDTH*W_HEIGHT*sizeof(Uint32));
	// Now, we write 255 to every byte of the buffer array.

	return false;
 }
 bool init(int W_WIDTH, int W_HEIGHT) {
 	this->W_WIDTH = W_WIDTH;
	this->W_HEIGHT = W_HEIGHT;
 	// If the init function is called with width and height parameters, the default width and height are overloaded.
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		// This will run SDL_Init, and if it fails, it will return a value < 0 , so the program ends.
		cout << "SDL init failed." << endl;
		return 1;
	}
	cout << "SDL init succeeded." << endl;
	window = SDL_CreateWindow("WindowTitle", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, W_WIDTH, W_HEIGHT, SDL_WINDOW_SHOWN);
	// We now create a pointer to a window that is equal to the pointer returned by the function SDL_CreateWindow, which makes a window and returns a pointer to that window.
	cout << "W_WIDTH is: " << W_WIDTH << " and W_HEIGHT is: " << W_HEIGHT << endl;
	if (window == NULL) {
		// If creating the window fails, however, we need to crash our program.
		cout << "Could not create window." << endl;
		SDL_Quit();
		return 2;
	}
	cout << "Window successfully created." << endl;
	buffer = new Uint32[W_WIDTH*W_HEIGHT];
	// We are going to make a buffer (area in memory) that can hold all of the pixels on the screen.
	// A Uint32 is a standardized size for an int by SDL, that is the same size as a pixel.
	// We are going to make an array of ints that can store every single pixel on our screen.
	renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_PRESENTVSYNC);
	// Now, we make a renderer, which can draw textures to the window.
	texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STATIC, W_WIDTH, W_HEIGHT);
	// Now, we make a texture, so that the renderer can render it to the window.
	
	if (renderer == NULL) { cout << "Render init failed!" << " SDL error: " << SDL_GetError() << endl; SDL_DestroyWindow(window); SDL_Quit(); return 3; } 
	// If we failed to create a renderer, we should return that the initialization of the render failed.
	if (texture == NULL) { cout << "Texture init failed!" << endl; SDL_DestroyWindow(window); SDL_DestroyRenderer(renderer); SDL_Quit(); return 4; }
	// Same goes for the texture.
	else { cout << "Render and texture init succeeded." << endl; }
	// The object event is NULL until an event happens, in which case it will be equal to that event.
	
	memset(buffer, 0, W_WIDTH*W_HEIGHT*sizeof(Uint32));
	// Now, we write 255 to every byte of the buffer array.

	return false;
 }
 void update() {
	SDL_UpdateTexture(texture, NULL, buffer, W_WIDTH*sizeof(Uint32));
	// We will update the texture so that it is equal to the value returned by the buffer. Before we set the buffer, this would return random garbage.
	SDL_RenderClear(renderer);
	// Now, we clear what is in the renderer.
	SDL_RenderCopy(renderer,texture,NULL,NULL);
	// Then, we copy whatever is in the texture to the renderer.
	SDL_RenderPresent(renderer);
	// Then, we present the renderer to SDL. This draws on the screen.
 }
 void setPixel(int x, int y, Uint32 color) {
 	//if (x < 0 || y < 0 || x >= W_WIDTH || y >= W_HEIGHT) { cout << "tried to draw off of the screen!" << endl; return; }
 	if (x < 0 || y < 0 || x >= W_WIDTH || y >= W_HEIGHT) { return; }
	// The above will end the function if it is requested to draw off of the screen.
	buffer[y*W_WIDTH + x] = color;
	// The above will move y * how wide the screen is across the array, to get to the correct y value, and then it will move over x more pixels to get to the correct x value, before we set that pixel to the color specified by the input to the function.
 }
private:
 void drawStraightLineX(int x, int y1, int y2, Uint32 color) {
	for (int y = y1; y < y2; y++) {
		setPixel(x, y, color);
	}
 }
public:
 void drawLine(int x1, int y1, int x2, int y2, Uint32 color) {
 	// This is Bresenham's line algorithm.
 	int deltaX = x2 - x1;
	int deltaY = y2 - y1;
	if (deltaX == 0) {
		drawStraightLineX(x1, y1, y2, color);
	}
	double deltaError = abs((double)deltaY/(double)deltaX);
	double error = 0.0;
	int y = y1;
	for (int x = x1; x < x2; x++) {
		setPixel(x, y, color);
		error += deltaError;
		if (error >= 0.5) {
			y++;
			error--;
		}
	}

 }
 void close() {
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	delete [] buffer;
 }
};

#endif
