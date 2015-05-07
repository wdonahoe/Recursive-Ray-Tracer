#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "udray.h"
#include "glm.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern Camera *ray_cam;       // camera info
extern int image_i, image_j;  // current pixel being shaded
extern bool wrote_image;      // hasf the last pixel been shaded?

// reflection/refraction recursion control

extern int maxlevel;          // maximum depth of ray recursion 
extern double minweight;      // minimum fractional contribution to color

// these describe the scene

extern vector < GLMmodel * > model_list;
extern vector < Sphere * > sphere_list;
extern vector < Light * > light_list;

double near;
double far;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// intersect a ray with the entire scene (.obj models + spheres)

// x, y are in pixel coordinates with (0, 0) the upper-left hand corner of the image.
// color variable is result of this function--it carries back info on how to draw the pixel

void trace_ray(int level, double weight, Ray *ray, Vect color)
{
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	int i;

	// test for intersection with all .obj models

	for (i = 0; i < model_list.size(); i++) {
		inter = intersect_ray_glm_object(ray, model_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// test for intersection with all spheres

	for (i = 0; i < sphere_list.size(); i++) {
		inter = intersect_ray_sphere(ray, sphere_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// "color" the ray according to intersecting surface properties

	// choose one of the simpler options below to debug or preview your scene more quickly.
	// another way to render faster is to decrease the image size.

	if (nearest_inter) {
		//shade_ray_false_color_normal(nearest_inter, color);
		//shade_ray_intersection_mask(color);  
		//shade_ray_diffuse(ray, nearest_inter, color);
		//shade_ray_depth(nearest_inter,color,near,far);
		shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}

	// color the ray using a default

	else
		shade_ray_background(ray, color); 

/*	if (nearest_inter->t == INT_MAX){
		shade_ray_background(ray, color); 
	}
	else if (nearest_inter) {
		shade_ray_false_color_normal(nearest_inter, color);
		shade_ray_intersection_mask(color);  
		shade_ray_diffuse(ray, nearest_inter, color);
		shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}*/
}

void shade_ray_depth(Intersection* inter,Vect color,double clamp_near, double clamp_far){
	// map [near,far] --> [0,1]
	double c = (inter->t - clamp_far) / (clamp_near - clamp_far);

	color[R] += c;
	color[G] += c;
	color[B] += c;

	VectClamp(color,0,1);
}

// ADD TO README
void sphere_intersection_normal(Intersection* inter, Sphere* s){
	Vect N;

	N[X] = (inter->P[X] - s->P[X]) / s->radius;
	N[Y] = (inter->P[Y] - s->P[Y]) / s->radius;
	N[Z] = (inter->P[Z] - s->P[Z]) / s->radius;
	//N[W] = 1.0;
	VectCopy(inter->N,N);

}

//ADD TO README
void sphere_intersection_surface(Intersection* inter, Sphere* s){
	inter->surf = s->surf;
}

//----------------------------------------------------------------------------

// test for ray-sphere intersection; return details of intersection if true

Intersection *intersect_ray_sphere(Ray *ray, Sphere *S)
{

	Intersection* inter;
	Vect ret;
	Vect v;

	VectSub(ray->orig, S->P, v);

	double b = VectDotProd(ray->dir,v);

	double discriminant = pow(b,2) - VectDotProd(v,v) + pow(S->radius,2);
	inter = make_intersection();

	if (discriminant > 0){
		double x1 = -b + sqrt(discriminant);
		double x2 = -b - sqrt(discriminant);
		double t;
		if (x1 < 0){
			if (x2 < 0)
				return NULL;
			t = x2;
		}
		else if (x2 < 0){
			if (x1 < 0)
				return NULL;
			t = x1;
		}
		else {
			t = (x1 > x2) ? x2 : x1;
		}
		inter->t = t;
		ret[X] = ray->orig[X] + t * ray->dir[X];
		ret[Y] = ray->orig[Y] + t * ray->dir[Y];
		ret[Z] = ray->orig[Z] + t * ray->dir[Z];

		VectCopy(inter->P, ret);
		sphere_intersection_normal(inter,S);
		sphere_intersection_surface(inter,S);

	}
	else if (discriminant == 0){
		inter->t = (-1)*b;
		ret[X] = ray->orig[X] + inter->t * ray->dir[X];
		ret[Y] = ray->orig[Y] + inter->t * ray->dir[Y];
		ret[Z] = ray->orig[Z] + inter->t * ray->dir[Z];

		VectCopy(inter->P, ret);
		sphere_intersection_normal(inter,S);
		sphere_intersection_surface(inter,S);
	}
	else if (discriminant < 0){
		return NULL;
	}

	return inter;
}

//----------------------------------------------------------------------------

// only local, ambient + diffuse lighting (no specular, shadows, reflections, or refractions)

//add to readme
void shade_ray_ambient(Ray *ray, Intersection *inter, Light* light, Vect color){

	color[R] += inter->surf->amb[R] * light->amb[R];
	color[G] += inter->surf->amb[G] * light->amb[G];
	color[B] += inter->surf->amb[B] * light->amb[B];

	VectClamp(color, 0, 1);
}

void shade_ray_diffuse(Ray *ray, Intersection *inter, Light* light, Vect color)
{
	double NdotL;
	Vect l;

	VectSub(light->P,inter->P,l);
	VectUnit(l);

	NdotL = VectDotProd(inter->N,l);
	NdotL = NdotL > 0 ? NdotL : 0;

	color[R] += inter->surf->diff[R] * light->diff[R] * NdotL;
	color[G] += inter->surf->diff[G] * light->diff[G] * NdotL;
	color[B] += inter->surf->diff[B] * light->diff[B] * NdotL; 

	VectClamp(color, 0, 1);
}

void shade_ray_specular(Ray *ray, Intersection *inter, Light* light, Vect color){
	Vect H;
	double NdotH;
	Vect eye;
	Vect l;

	VectSub(ray_cam->eye,inter->P,eye);
	VectSub(light->P,inter->P,l);
	VectUnit(eye);
	VectUnit(l);
	VectAddS(1.0,l,eye,H);
	VectUnit(H);

	NdotH = VectDotProd(inter->N,H);
	NdotH = NdotH > 0 ? NdotH : 0;

	color[R] += inter->surf->spec[R] * light->spec[R] * pow(NdotH, inter->surf->spec_exp);
	color[G] += inter->surf->spec[G] * light->spec[G] * pow(NdotH, inter->surf->spec_exp);
	color[B] += inter->surf->spec[B] * light->spec[B] * pow(NdotH, inter->surf->spec_exp);

}

// add to readme
bool trace_ray_shadow(Ray *ray)
{
	int i;
	bool hit = false;

	// test for intersection with all .obj models until hit.

	while (i < model_list.size() && !hit) {
		if (intersect_ray_glm_object(ray, model_list[i]))
			hit = true;
		i++;
	}
	i = 0;

	// test for intersection with all spheres until hit. 
	while(i < sphere_list.size() && !hit) {
		if (intersect_ray_sphere(ray, sphere_list[i]) != NULL){
			hit = true;
		}
		i++;
	}
	return hit;
}

bool in_shadow(Intersection *inter, Light* light){
	double eps = 1.0e-3;
	Vect D;
	Vect orig;
	Ray* shadow_ray;
	
	VectSub(light->P,inter->P,D);
	VectUnit(D);

	orig[X] = inter->P[X] + eps*D[X];
	orig[Y] = inter->P[Y] + eps*D[Y];
	orig[Z] = inter->P[Z] + eps*D[Z];

	shadow_ray = make_ray(orig,D);
	return trace_ray_shadow(shadow_ray);
}


//----------------------------------------------------------------------------

// same as shade_ray_diffuse(), but add specular lighting + shadow rays (i.e., full Phong illumination model)
void shade_ray_local(Ray *ray, Intersection *inter, Vect color)
{

	for(int i = 0; i < light_list.size(); i++){
		if (in_shadow(inter,light_list[i])) {
			shade_ray_ambient(ray,inter,light_list[i],color);
		}
		else {
			shade_ray_ambient(ray,inter,light_list[i],color);
			shade_ray_diffuse(ray,inter,light_list[i],color);
			shade_ray_specular(ray,inter,light_list[i],color);
		}
	} 

	VectClamp(color, 0, 1);
}

//----------------------------------------------------------------------------

// full shading model: ambient/diffuse/specular lighting, shadow rays, recursion for reflection, refraction

// level = recursion level (only used for reflection/refraction)

void shade_ray_recursive(int level, double weight, Ray *ray, Intersection *inter, Vect color)
{
	Vect dir;
	Vect pos;
	//double NdotD;
	Ray* reflection_ray;
	double eps = 1.0e-9;

	// initialize color to Phong reflectance model
	//if (level == 0){
	//}
	shade_ray_local(ray, inter, color);

	// if not too deep, recurse
	if (level + 1 < maxlevel) {

		// add reflection component to color

		if (inter->surf->reflectivity * weight > minweight) {

			reflection_direction(ray->dir,inter->N,dir);
			VectUnit(dir);

			pos[X] += eps * dir[X];
			pos[Y] += eps * dir[Y];
			pos[Z] += eps * dir[Z];

			reflection_ray = make_ray(pos,dir);

			trace_ray(level+1,weight,reflection_ray,color);
				
		}
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// ray trace another pixel if the image isn't finished yet

void idle()
{
	if (image_j < ray_cam->im->h) {
		raytrace_one_pixel(image_i, image_j);

		image_i++;

		if (image_i == ray_cam->im->w) {
			image_i = 0;
			image_j++;
		}    
	}

	// write rendered image to file when done

	else if (!wrote_image) {

		write_PPM("output.ppm", ray_cam->im);

		wrote_image = true;
	}

	glutPostRedisplay();
}

//----------------------------------------------------------------------------

// show the image so far

void display(void)
{
	// draw it!

	glPixelZoom(1, -1);
	glRasterPos2i(0, ray_cam->im->h);

	glDrawPixels(ray_cam->im->w, ray_cam->im->h, GL_RGBA, GL_FLOAT, ray_cam->im->data);

	glFlush ();
}

//----------------------------------------------------------------------------

void init()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ray_cam->im->w, 0.0, ray_cam->im->h);
}

//----------------------------------------------------------------------------

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// initialize scene (must be done before scene file is parsed)

	init_raytracing();

	if (argc == 2){
		parse_scene_file(argv[1], ray_cam);
		near = ray_cam->clip[NEAR];
		far = ray_cam->clip[FAR];
	}
	else {
		printf("missing .scene file\n");
		exit(1);
	}

	// opengl business

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ray_cam->im->w, ray_cam->im->h);
	glutInitWindowPosition(500, 300);
	glutCreateWindow("hw3");
	init();

	glutDisplayFunc(display); 
	glutIdleFunc(idle);

	glutMainLoop();

	return 0; 
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
