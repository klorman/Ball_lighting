#include "TXlib.h"
#include <stdlib.h>
#include <string>
#include <cmath>
#include <math.h>

RGBQUAD* Video_memory = NULL;

const int SSAA = 1;

struct point {
	double x, y, z;
};

class CoordSys {
public:
	point size_;
	void draw_pixel(point coords, point color_of_point);
	point to_pixels(point coords);

	CoordSys(point size) {
		size_ = size;
	}
};

class Vector {
public:
	point coords_;
	double length();
	void normalization();

	Vector(point coords = { 0, 0, 0 }) {
		coords_ = coords;
	}
};


class Object;
class Light_source {
public:
	Vector light_source_;
	Vector light_source_color_;
	void draw_light_source(CoordSys& vector_space, Object* objects, double R);

	Light_source(Vector light_source, Vector light_source_color) {
		light_source_ = light_source;
		light_source_color_ = light_source_color;
		++count();
	}

	~Light_source() {
		--count();
	}

	size_t& count() {
		static size_t c = 0;
		return c;
	}
};

class Object {
public:
	Vector coords_;
	Vector color_;
	double R_;
	void draw_object(CoordSys& vector_space, double rotation, Light_source& purple, Light_source& rainbow, Object* objects, size_t test);

	Object(Vector coords = Vector({ 0, 0, 0 }), Vector color = Vector({ 0, 0, 0 }), double R = 0) {
		coords_ = coords;
		color_ = color;
		R_ = R;
		++count();
	}

	~Object() {
		--count();
	}

	size_t& count() {
		static size_t c = 0;
		return c;
	}
};

Vector operator + (Vector a, Vector b) {
	return Vector({ a.coords_.x + b.coords_.x,
					a.coords_.y + b.coords_.y,
					a.coords_.z + b.coords_.z });
}

Vector operator * (double koef_of_length, Vector a) {
	return Vector({ koef_of_length * a.coords_.x,
					koef_of_length * a.coords_.y,
					koef_of_length * a.coords_.z });
}

Vector operator * (Vector a, Vector b) {
	return Vector({ a.coords_.x * b.coords_.x,
					a.coords_.y * b.coords_.y,
					a.coords_.z * b.coords_.z});
}

double operator ^ (Vector a, Vector b) {
	return (a.coords_.x * b.coords_.x + 
			a.coords_.y * b.coords_.y + 
			a.coords_.z * b.coords_.z);
}

double Vector::length() {
	return sqrt(coords_.x * coords_.x + 
				coords_.y * coords_.y + 
				coords_.z * coords_.z);
}

void Vector::normalization() {
	Vector temp_vector = *this;
	double length = temp_vector.length();

	if (length != 0) temp_vector = 1. / length * temp_vector;

	coords_.x = temp_vector.coords_.x;
	coords_.y = temp_vector.coords_.y;
	coords_.z = temp_vector.coords_.z;
}

point CoordSys::to_pixels(point coords) {
	point start_of_coord = { size_.x / 2, size_.y / 2 };

	return point({ coords.x + start_of_coord.x, start_of_coord.y - coords.y });
}

void CoordSys::draw_pixel(point coords, point color_of_point) {
	point rec_coords = to_pixels(coords);

	RGBQUAD* pixel = &Video_memory[((int)size_.x - (int)rec_coords.y) * (int)size_.y + (int)rec_coords.x];
	pixel->rgbRed = color_of_point.x;
	pixel->rgbGreen = color_of_point.y;
	pixel->rgbBlue = color_of_point.z;
}


Vector get_circle(Vector& n, Vector& O1, Vector color_of_ball, Vector new_color, double R1, double z) {
	Vector n1({ n.coords_.x - O1.coords_.x, n.coords_.y - O1.coords_.y, 0 });

	if (n1.length() <= R1) {
		double z1 = sqrt(R1 * R1 - n1.coords_.x * n1.coords_.x - n1.coords_.y * n1.coords_.y);
		if (z1 + O1.coords_.z >= z && -z1 + O1.coords_.z <= z) {
			return new_color;
		}
	}
	return color_of_ball;
}

Vector get_color(CoordSys& vector_space, Vector& O, Vector& n0, Vector& n, Vector& observer, Vector& O1, Light_source& source, Vector& color_of_ball, double R) {

	if (Vector({ n.coords_.x, n.coords_.y, 0 }).length() < R) {
		double z = sqrt(R * R - pow(n.coords_.x, 2) - pow(n.coords_.y, 2));
		n.coords_.z = z;

		//n.coords_.x += sin(10000 * n.coords_.x + 100000000) * 10;
		//n.coords_.y += cos(10 * n.coords_.z) * 500;
		n.coords_.z += sin(10000 * n.coords_.z + 100000000) * 10;

		Vector color_of_point = color_of_ball;

		//------------------------------------------------------------------------------------------------------//
		color_of_point = get_circle(n, O1, color_of_point, Vector({ 255, 255, 255 }), SSAA * 50 * R / 155, z);

		Vector O2 = O1 + (R / 155 * SSAA * Vector({ 0, 11, 0 }));
		color_of_point = get_circle(n, O2, color_of_point, Vector({ 10, 10, 10 }), SSAA * 15 * R / 155, z);
		color_of_point = get_circle(n, O2, color_of_point, Vector({ 250, 250, 250 }), SSAA * 10 * R / 155 , z);

		O2 = O1 + (R / 155 * SSAA * Vector({ 0, -11, 0 }));
		color_of_point = get_circle(n, O2, color_of_point, Vector({ 10, 10, 10 }), SSAA * 15 * R / 155, z);
		color_of_point = get_circle(n, O2, color_of_point, Vector({ 250, 250, 250 }), SSAA * 10 * R / 155, z);
		//------------------------------------------------------------------------------------------------------//

		//n = n + O;

		Vector L = source.light_source_ + (-1 * (n + O));

		double brightness = (L ^ n) / (L.length() * n.length());

		if (brightness < 0) brightness = 0;


		Vector L1 = 2 * L.length() * brightness * (1. / n.length() * n) + ((-1) * L);

		double k = 50;
		double cos_beta = (observer ^ L1) / (observer.length() * L1.length()), sin_beta = 1 * sqrt(1 - pow((cos_beta), 2));

		
		if (cos_beta < 0 || source.light_source_.coords_.z < -150 || source.light_source_.length() <= R) cos_beta = 0;

		double I = pow(cos_beta, k);

		
		Vector ambient_lighting_color({ 255, 255, 255 });

		Vector ambient = 1. / source.count() * 255 * ((0.04 / 255 * ambient_lighting_color) * (1. / 255 * color_of_point));
		Vector lambert = 255 * (brightness * 1. / 255 * color_of_point) * (1. / 255 * source.light_source_color_);

		color_of_point = ambient + lambert + I * source.light_source_color_;
		
		return color_of_point;
	}

	return Vector({ 0, 0, 0 });
}

bool in_object(Object* objects, Vector pixel, size_t test) {
	for (size_t number_of_object = 0; number_of_object < objects[0].count(); number_of_object++) {
		if (number_of_object == test) continue;

		double z = 0;

		Vector temp = pixel + (-1 * objects[number_of_object].coords_);
		temp.coords_.z = 0;

		if (temp.length() < objects[number_of_object].R_) z = sqrt(pow(objects[number_of_object].R_, 2) - pow(temp.coords_.x, 2) - pow(temp.coords_.y, 2)) + objects[number_of_object].coords_.coords_.z;
		else continue;

		if (pixel.coords_.z < z) return TRUE;
	}

	return FALSE;
}

void Light_source::draw_light_source(CoordSys& vector_space, Object* objects, double R) {
	light_source_ = 1. / SSAA * light_source_;

	double r = 8 * (light_source_.coords_.z + 450) / 800;

	for (double x = light_source_.coords_.x - r; x < light_source_.coords_.x + r; x++) {
		for (double y = light_source_.coords_.y - r; y < light_source_.coords_.y + r; y++) {

			if (x > -350 &&
				x < 350 &&
				y > -350 &&
				y < 350 &&
				(-1 * Vector({ x, y, 0 }) + Vector({ light_source_.coords_.x, light_source_.coords_.y, 0 })).length() < r &&
				!in_object(objects, Vector({ x, y, light_source_.coords_.z }), -1))
				vector_space.draw_pixel(point({ x, y, 0 }), light_source_color_.coords_);
		}
	}
}

void Object::draw_object(CoordSys& vector_space, double rotation, Light_source& purple, Light_source& rainbow, Object* objects, size_t test) {
	R_ *= 1 + coords_.coords_.z / 4 / R_;

	Vector observer({ 0, 0, 350 });
	Vector O1({ R_, 0, 0 });

	observer = SSAA * (observer + (-1 * coords_));
	O1 = SSAA * O1;


	O1.coords_.x = R_ * sin(0.5 * rotation);
	//O1.coords_.y = R * sin(rotation);
	O1.coords_.z = R_ * cos(0.5 * rotation);




	for (double y = -160; y < 160; y++)
		for (double x = -160; x < 160; x++) {
			Vector n({ x, y, 0 });
			Vector color_of_point({ 0, 0, 0 });

			for (double x1 = 0; x1 < SSAA; x1++)
				for (double y1 = 0; y1 < SSAA; y1++) {

					Vector n1 = SSAA * n + Vector({ x1, y1, 0 });
					
					color_of_point = color_of_point + 1. / (SSAA * SSAA) * get_color(vector_space, coords_, n, n1, observer, O1, purple, color_, R_ * SSAA);		
					color_of_point = color_of_point + 1. / (SSAA * SSAA) * get_color(vector_space, coords_, n, n1, observer, O1, rainbow, color_, R_ * SSAA);

					color_of_point = Vector({ MIN(color_of_point.coords_.x, 255),
											  MIN(color_of_point.coords_.y, 255),
											  MIN(color_of_point.coords_.z, 255) });
				}


			if (n.length() < R_) n.coords_.z = sqrt(R_ * R_ - pow(n.coords_.x, 2) - pow(n.coords_.y, 2));

			if (color_of_point.length() != 0 &&	!in_object(objects, n + coords_, test)) vector_space.draw_pixel((n + coords_).coords_, color_of_point.coords_);
		}
}

void draw_scene(CoordSys& vector_space, double rotation) {
	txClear();

	Light_source purple(Vector({ -320, 120, 320 }), Vector({ 255, 0, 255 }));
	Light_source rainbow(Vector({ 150, 320, 220 }), Vector({ 255, 255, 255 }));

	purple.light_source_ = SSAA * (purple.light_source_);
	rainbow.light_source_ = SSAA * (rainbow.light_source_);
	purple.light_source_.coords_.x *= cos(rotation);
	purple.light_source_.coords_.z *= sin(rotation);

	rainbow.light_source_color_.coords_.x *= abs(sin(rotation));
	rainbow.light_source_color_.coords_.y *= abs(sin(1.5 * rotation));
	rainbow.light_source_color_.coords_.z *= abs(sin(3 * rotation));

	rainbow.light_source_.coords_.x *= cos(0.5 * rotation) * cos(rotation);
	rainbow.light_source_.coords_.z *= cos(0.5 * rotation);
	rainbow.light_source_.coords_.y *= sin(0.5 * rotation);
	 
	Object objects[] = { Object(Vector({ 0, 0, 0 }), Vector({ 128, 64, 64 }), 155), Object(Vector({ 220 * cos(rotation), 220 * sin(3 * rotation), 155 * sin(rotation) }), Vector({ 64, 128, 64 }), 55) };

	for (size_t number_of_object = 0; number_of_object < objects[0].count(); number_of_object++) objects[number_of_object].draw_object(vector_space, rotation, purple, rainbow, objects, number_of_object);

	purple.draw_light_source(vector_space, objects, 155);
	rainbow.draw_light_source(vector_space, objects, 155);
}

int main() {
	CoordSys vector_space({ 700, 700 });
	txCreateWindow(vector_space.size_.x, vector_space.size_.y);
	txSetFillColor(TX_BLACK);
	
	Video_memory = txVideoMemory();

	txBegin();



	for (double rotation = 0; !txGetAsyncKeyState(VK_ESCAPE);) {
		double fps = txGetFPS();

		if (fps != 0) rotation += 0.5 / fps;

		if (rotation >= 4 * txPI) rotation -= 4 * txPI;

		draw_scene(vector_space, rotation);

		


		std::cout << fps << "\n";
		txSleep();
	}

	txEnd();
}
