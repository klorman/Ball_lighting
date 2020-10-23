#include "TXlib.h"
#include <stdlib.h>
#include <string>
#include <cmath>
#include <math.h>

RGBQUAD* Video_memory = NULL;

const int SSAA = 2;

struct point {
	double x, y, z;
};

class CoordSys {
private:
	point coords0_, coords1_;
	double scaleX_, scaleY_;

public:
	void draw_pixel(point coords, point color_of_point);
	void draw_point(point coords);
	void draw_line(point coords0, point coords1);
	point to_pixels(point coords);

	CoordSys(point coords0, point coords1, double scaleX, double scaleY) {
		coords0_ = coords0;
		coords1_ = coords1;
		scaleX_ = scaleX;
		scaleY_ = scaleY;
	}
};

class Vector {
public:
	point coords_;
	double length();
	void normalization();
	void draw_vector(point start_coords, CoordSys& vector_space);

	Vector(point coords) {
		coords_.x = coords.x;
		coords_.y = coords.y;
		coords_.z = coords.z;
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

void Vector::draw_vector(point start_coords, CoordSys& vector_space) {
	point end_coords = { coords_.x + start_coords.x, coords_.y + start_coords.y };

	vector_space.draw_line(start_coords, end_coords);
	vector_space.draw_point(start_coords);

	Vector vector_for_arrow1({ -coords_.y ,  coords_.x });
	Vector vector_for_arrow2({ coords_.y , -coords_.x });
	Vector reverse = -1 * *this;

	reverse.normalization();
	vector_for_arrow1.normalization();
	vector_for_arrow2.normalization();

	Vector arrow1 = reverse + vector_for_arrow1;
	Vector arrow2 = reverse + vector_for_arrow2;

	arrow1.normalization();
	arrow2.normalization();

	point end_of_arrow1 = { arrow1.coords_.x + end_coords.x, arrow1.coords_.y + end_coords.y };
	vector_space.draw_line(end_coords, end_of_arrow1);

	point end_of_arrow2 = { arrow2.coords_.x + end_coords.x, arrow2.coords_.y + end_coords.y };
	vector_space.draw_line(end_coords, end_of_arrow2);
}

point CoordSys::to_pixels(point coords) {
	point start_of_coord = { coords1_.x / 2, coords1_.y / 2 };
	point rec_coord;

	rec_coord.x = coords.x * scaleX_ + start_of_coord.x;
	rec_coord.y = start_of_coord.y - coords.y * scaleY_;

	return rec_coord;
}

void CoordSys::draw_point(point coords) {
	txSetColor(TX_BLACK);
	txSetFillColor(TX_BLACK);

	point rec_coord = to_pixels(coords);

	if (rec_coord.y >= coords0_.y &&
		rec_coord.y <= coords1_.y &&
		rec_coord.x <= coords1_.x &&
		rec_coord.x >= coords0_.x)
		txCircle(rec_coord.x, rec_coord.y, 2);
}

void CoordSys::draw_line(point coords0, point coords1) {
	txSetColor(TX_BLUE);

	point rec_coord0 = to_pixels(coords0);
	point rec_coord1 = to_pixels(coords1);

	if (rec_coord0.y >= coords0_.y &&
		rec_coord0.y <= coords1_.y &&
		rec_coord0.x <= coords1_.x &&
		rec_coord0.x >= coords0_.x &&
		rec_coord1.y >= coords0_.y &&
		rec_coord1.y <= coords1_.y &&
		rec_coord1.x <= coords1_.x &&
		rec_coord1.x >= coords0_.x)
		txLine(rec_coord0.x, rec_coord0.y, rec_coord1.x, rec_coord1.y);
}

void CoordSys::draw_pixel(point coords, point color_of_point) {
	//txSetFillColor(color_of_point);
	point rec_coords = to_pixels(coords);
	//txSetPixel(rec_coords[0], rec_coords[1], color_of_point);
	RGBQUAD* pixel = &Video_memory[(700 - (int)rec_coords.y) * 700 + (int)rec_coords.x];
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

Vector get_color(CoordSys& vector_space, Vector& n, Vector& light_source, Vector& light_source_color,int R) {

	if (n.length() <= R) {
		double z = sqrt(R * R - n.coords_.x * n.coords_.x - n.coords_.y * n.coords_.y);
		n.coords_.z = z;

		//n.coords_.x += sin(100 * n.coords_.z) * 50;
		//n.coords_.y += cos(10 * n.coords_.z) * 500;
		//n.coords_.z += sin(0.1 * n.coords_.x) * 1000;

		Vector color_of_ball({ 128, 64, 64 });

		Vector O1({ -120, 0, 120 });
		O1 = SSAA * O1;

		color_of_ball = get_circle(n, O1, color_of_ball, Vector({ 255, 255, 255 }), SSAA * 50, z);

		Vector O2 = O1 + (SSAA * Vector({ -1, 11, 1 }));
		color_of_ball = get_circle(n, O2, color_of_ball, Vector({ 0, 0, 0 }), SSAA * 22, z);
		color_of_ball = get_circle(n, O2, color_of_ball, Vector({ 250, 250, 250 }), SSAA * 18, z);

		O2 = O1 + (SSAA * Vector({ -1, -11, 1 }));
		color_of_ball = get_circle(n, O2, color_of_ball, Vector({ 0, 0, 0 }), SSAA * 22, z);
		color_of_ball = get_circle(n, O2, color_of_ball, Vector({ 250, 250, 250 }), SSAA * 18, z);

		//O2 = O1 + (SSAA * Vector({ -10, -11, -3 }));
		//color_of_ball = get_circle(n, O2, color_of_ball, Vector({ 250, 250, 250 }), SSAA * 22, z);


		Vector L = light_source + (-1 * n);

		double brightness = (L ^ n) / (L.length() * n.length());

		if (brightness < 0) brightness = 0;

		Vector ambient_lighting({ 255, 255, 255 });

		Vector observer({ 0, 0, 350 });
		Vector L1 = 2 * L.length() * brightness * (1. / n.length() * n) + ((-1) * L);

		double k = 50;
		double cos_beta = (observer ^ L1) / (observer.length() * L1.length()), sin_beta = 1 * sqrt(1 - pow((cos_beta), 2));
		
		if (cos_beta < 0 || light_source.coords_.z < 0) cos_beta = 0;

		double I = pow(cos_beta, k);

		Vector lambert = 255 * (brightness * 1. / 255 * color_of_ball) * (1. / 255 * light_source_color);
		Vector ambient = 255 * ((0.04 / 255 * ambient_lighting) * (1. / 255 * color_of_ball));

		color_of_ball = lambert + ambient + I * light_source_color;
		
		return Vector({MIN(color_of_ball.coords_.x, 255),
					   MIN(color_of_ball.coords_.y, 255),
					   MIN(color_of_ball.coords_.z, 255)});
	}

	return Vector({ 0, 0, 0 });
}

void draw_light_source(Vector& light_source_color, Vector& light_source) {
	txSetFillColor(RGB(light_source_color.coords_.x,
		light_source_color.coords_.y,
		light_source_color.coords_.z));

	txSetColor(RGB(light_source_color.coords_.x,
		light_source_color.coords_.y,
		light_source_color.coords_.z));

	if (light_source.coords_.x > -350,
		light_source.coords_.x < 350,
		light_source.coords_.y > -350,
		light_source.coords_.y < 350) txCircle(350 + light_source.coords_.x, 350 - light_source.coords_.y, 8 * (light_source.coords_.z + 450.) / 800);

}

void draw_object(double rotation) {
	txSetFillColor(TX_BLACK);
	txClear();

	double R = 155;
	CoordSys vector_space({ 0, 0 }, { 700, 700 }, 1, 1);
	Vector light_source({ -320, 120, 320 });
	Vector light_source_color({ 255, 0, 255 });

	light_source.coords_.x *= cos(rotation);
	light_source.coords_.z *= sin(rotation);

	for (double y = -160; y < 160; y++) {
		for (double x = -160; x < 160; x++) {
			Vector n({ x, y, 0 });
			Vector color_of_ball({ 0, 0, 0 });

			Vector light_source1 = SSAA * light_source;
			n = SSAA * n;

			for (double x1 = 0; x1 < SSAA; x1++)
				for (double y1 = 0; y1 < SSAA; y1++) {
					Vector n1 = n + Vector({ x1, y1, 0 });
					
					color_of_ball = color_of_ball + 1. / (SSAA * SSAA) * get_color(vector_space, n1, light_source1, light_source_color, R * SSAA);
				}

			n = 1. / SSAA * n;

			if (color_of_ball.length() != 0) vector_space.draw_pixel(point({ n.coords_.x, n.coords_.y, n.coords_.z }), point({ color_of_ball.coords_.x, color_of_ball.coords_.y, color_of_ball.coords_.z }));
		}
	}

	if (Vector({ light_source.coords_.x, light_source.coords_.y, 0}).length() >= R + 8 || light_source.coords_.z >= 0) draw_light_source(light_source_color, light_source);
}

int main() {
	txCreateWindow(700, 700);

	Video_memory = txVideoMemory();

	txBegin();

	double rotation = 0;



	for (rotation; !txGetAsyncKeyState(VK_ESCAPE);) {
		double fps = txGetFPS();

		if (fps != 0) rotation += 1 / fps;

		draw_object(rotation);

		std::cout << fps << "\n";
		txSleep();
	}

	txEnd();
}
