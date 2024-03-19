/*
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <unistd.h>
#include <math.h>
*/

#include "graphics.h"
#define PI (3.141592653589793)
#define MAX_X 2.0
#define MAX_Y 2.0



void draw_graph(Display *display, Window win, int n, double *X, double *u){
	Window root_return;
	int x_return, y_return;
	unsigned int width_return, hight_return;
	unsigned int border_width_return, depth_return;
	GC gc;
	XGCValues gcv;
	int blackpixel = BlackPixel(display, DefaultScreen(display));
	int whitepixel = WhitePixel(display, DefaultScreen(display));

	

	double x1;
	double x2;

	double h1;
	double h2;

	

	gcv.background = whitepixel;
	gcv.foreground = blackpixel;
	gc = XCreateGC(display, DefaultRootWindow(display),
	   GCForeground |GCBackground, &gcv);
	XGetGeometry(display, win, &root_return, &x_return, &y_return,&width_return, &hight_return, &border_width_return, &depth_return);
	XDrawLine(display, win, gc, 0, hight_return/2, width_return, hight_return/2);
	XDrawLine(display, win, gc, width_return/4, 0, width_return/4, hight_return);

	// This is for the plot
	// draw a line between two points:
        	
	int j;
	for (j=0; j<n-2; j++){
		x1 =( X[j]/MAX_X + 1.0)*width_return/4.0;
	       
		h1 = (1-(u[j]/MAX_Y))*hight_return/2.0;
	  					      
	 	x2 = (X[j+1]/MAX_X+1.0)*width_return/4.0;
            
		h2 = (1-(u[j+1]/MAX_Y))*hight_return/2.0;
  		XDrawLine(display,win,gc,x1,h1,x2,h2);

	}
	
}


void Draw(int n, double *X, double *u){
	Display *dpy;
	int screen;
	Window win, root_win;
	XEvent event;
	unsigned int depth;
	XSetWindowAttributes attrs;
	dpy = XOpenDisplay(NULL);

	if (dpy == NULL){
		fprintf(stderr, "Cannot display\n");
		exit(1);
	}
	screen = DefaultScreen(dpy);
	depth = DefaultDepth(dpy,screen);
	root_win = RootWindow(dpy, screen);

	attrs.border_pixel = BlackPixel(dpy, screen);
	attrs.background_pixel = WhitePixel(dpy,screen);
	attrs.override_redirect = True;
	attrs.colormap = CopyFromParent;
	attrs.event_mask = ExposureMask | KeyPressMask;
        
	win = XCreateWindow(dpy, root_win, 700,700,800,500,0,depth,InputOutput,
			CopyFromParent, CWBackPixel | CWColormap | CWBorderPixel | CWEventMask, &attrs);

	XMapWindow(dpy, win);
	while(1){
		XNextEvent(dpy, &event);
		if(event.type == Expose){
			draw_graph(dpy,win,n,X,u); // <<------
		}
	
	if (event.type == KeyPress){
		XDestroyWindow(dpy, win);
		XCloseDisplay(dpy);
		break;
	}
	}
}
/*
int main(){
	int n = 11;
	double *x = (double*)malloc(n*sizeof(double));
	double *ID = (double*)malloc(n*sizeof(double));
	double b = 2*PI;
	double a = 0;
	double h = (b-a)/(n-1);
	int i;
	for ( i=0; i<n; i++ ){
		x[i] = a+i*h;             // independent loops
		ID[i] = sin(x[i]);
//		printf("\n%f\n", x[i]);
	}
	Draw(n,x,ID);
	return 0;

}
*/

