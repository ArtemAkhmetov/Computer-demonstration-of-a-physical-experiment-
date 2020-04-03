#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

using namespace sf;
using namespace std;

enum type_potential { NUL, GR, K, LJ };
const float screenWidth = (float)VideoMode::getDesktopMode().width;
const float screenHeight = (float)VideoMode::getDesktopMode().height; //константы расширения экрана
const double max_width = 500 * screenWidth / 1920; //макс. ширина частиц
const int max_num = 29;    //макс. число частиц
const int max_temperature = 50;
const int max_mass = 50;
const int max_coef_C = 5;
const int nBins = 16;
const double max_speed = 0.08 * screenWidth / 1920;
const double max_accel = 0.001 * screenWidth / 1920;
const float plotSize = 846;//870
const float plotPosX = 800;//820
const float plotPosY = 1033;//1040
const float binCoef = 330;
const float cGR = 9.0;
const int C1_0 = 255, C2_0 = 255, C3_0 = 255;
const int C1_1 = 192, C2_1 = 192, C3_1 = 192;
const int C1_2 = 224, C2_2 = 238, C3_2 = 238;
const int num_buttons = 28, num_texts = 62;
int coef_C = 1;
int crashed = false;
int state_GR[100] = { 0 };
int r = 255, g = 0, b = 0, c_buf = 0;
const float sigma = 0.005;
const float delta = 600.0;

class Particle {
	double* x, *speed, *accel, w, h, *amount;
	int num, nbins, *amount_now;
	bool interact;
	RectangleShape* p;
	RectangleShape* bin;
public:
	Particle() {
		w = h = 0;
		accel = speed = x = amount = 0;
		amount_now = 0;
		num = nbins = 0;
		p = 0;
		bin = 0;
		interact = true;
	}

	~Particle() {
		if (p) delete[] p;
		if (x) delete[] x;
		if (bin) delete[] bin;
		if (speed) delete[] speed;
		if (accel) delete[] accel;
		if (amount) delete[] amount;
		if (amount_now) delete[] amount_now;
	}

	double dist(int i, int j) {
		if (x[i] - x[j] >= 0)
			return x[i] - x[j] - w;
		else
			return x[j] - x[i] - w;
	}

	int sign(double a) {
		if (a > 0)
			return 1;
		else if (a < 0)
			return -1;
		else
			return 0;
	}

	void swapSpeed(int i, int j) {
		double c = speed[i];
		speed[i] = speed[j];
		speed[j] = c;
	}

	double sampleMean() {
		double sum = 0;
		for (int i = 0; i < num - 1; ++i)
			sum += dist(i, i + 1) / ((int)(plotSize * screenWidth / 1920 / nBins));
		return round(100 * sum / ((double)num - 1)) / 100;
	}

	double sampleSTD() {
		double sum = 0, cur_dist, sMean;
		for (int i = 0; i < num - 1; ++i) {
			cur_dist = dist(i, i + 1) / ((int)(plotSize * screenWidth / 1920 / nBins));
			sMean = sampleMean();
			sum += (cur_dist - sMean) * (cur_dist - sMean);
		}
		return round(100 * sqrt(sum / ((double)num - 1))) / 100;
	}

	void update(float time) {
		if (interact) {
			double cur_dist;
			speed[0] += accel[0] * time;
			if (x[0] < 6 * screenWidth / 1920 && speed[0] < 0) {
				speed[0] *= -1;
				x[0] = 6 * screenWidth / 1920;
			}

			if (abs(speed[0]) > max_speed)
				speed[0] = sign(speed[0]) * max_speed;

			if (screenWidth / 1920 < x[0] + speed[0] * time && x[0] + speed[0] * time + w < x[1] + ((abs(speed[1] + accel[1] * time) > max_speed) ? sign(speed[1] + accel[1] * time) * max_speed : (speed[1] + accel[1] * time))* time) {
				x[0] += speed[0] * time;
			}
			else {
				x[0] = (screenWidth / 1920 + x[1] + ((abs(speed[1] + accel[1] * time) > max_speed) ? sign(speed[1] + accel[1] * time) * max_speed : (speed[1] + accel[1] * time))* time) / 2;
				//cout << "x[0] = (6 + x[1])/2" << endl;
			}
			p[0].setPosition((float)x[0], 70 * screenHeight / 1080);

			for (int i = 1; i < num - 1; ++i) {
				speed[i] += accel[i] * time;
				if (abs(speed[i]) > max_speed)
					speed[i] = sign(speed[i]) * max_speed;
				if (x[i - 1] + w < x[i] + speed[i] * time && x[i] + speed[i] * time + w < x[i + 1] + ((abs(speed[i + 1] + accel[i + 1] * time) > max_speed) ? sign(speed[i + 1] + accel[i + 1] * time) * max_speed : (speed[i + 1] + accel[i + 1] * time))* time) {
					x[i] += speed[i] * time;
				}
				else {
					x[i] = (x[i - 1] + w + x[i + 1] + ((abs(speed[i + 1] + accel[i + 1] * time) > max_speed) ? sign(speed[i + 1] + accel[i + 1] * time) * max_speed : (speed[i + 1] + accel[i + 1] * time))* time) / 2;
					//cout << "x[i] = (x[i-1] + x[i+1])/2" << endl;
				}
				p[i].setPosition(x[i], 70 * screenHeight / 1080);
			}

			speed[num - 1] += accel[num - 1] * time;
			if (x[num - 1] > (screenWidth - w - 6 * screenWidth / 1920) && speed[num - 1] > 0) {
				speed[num - 1] *= -1;
				x[num - 1] = screenWidth - w - 6 * screenWidth / 1920;
			}

			if (abs(speed[num - 1]) > max_speed)
				speed[num - 1] = sign(speed[num - 1]) * max_speed;
			if (x[num - 2] + w < x[num - 1] + speed[num - 1] * time && x[num - 1] + speed[num - 1] * time + w < (screenWidth - screenWidth / 1920)) {
				x[num - 1] += speed[num - 1] * time;
			}
			else {
				x[num - 1] = (x[num - 2] + w + (screenWidth - screenWidth / 1920)) / 2;
				//cout << "x[num-1] = (x[num-2] + 1920)/2" << endl;
			}
			p[num - 1].setPosition(x[num - 1], 70 * screenHeight / 1080);

			for (int i = 0; i < nbins; ++i)
				amount_now[i] = 0;

			for (int i = 0; i < num - 1; ++i) {
				cur_dist = dist(i, i + 1);
				if (cur_dist > plotSize* screenWidth / 1920) {
					++amount[nbins - 1];
					++amount_now[nbins - 1];
				}
				else {
					for (int k = 0; k < nbins - 1; ++k) {
						if (cur_dist >= k * plotSize * screenWidth / 1920 / (nbins - 1) && cur_dist <= (k + 1) * plotSize * screenWidth / 1920 / (nbins - 1)) {
							++amount[k];
							++amount_now[k];
							break;
						}
					}
				}
			}
		}
	}

	void interaction() {
		if (interact) {
			for (int i = 0; i < num - 1; ++i) {
				for (int j = i + 1; j < num; ++j) {
					if (x[i] + w - 2 > x[j]) {
						//cout << "crashed" << endl;
						crashed = true;
					}
				}
			}
			for (int i = 0; i < num - 1; ++i) {
				if (dist(i, i + 1) < 7) {
					if (speed[i] > 0 && speed[i + 1] < 0 || sign(speed[i]) == sign(speed[i + 1])) {
						swapSpeed(i, i + 1);
					}
				}
			}
		}
	}

	void potential_interaction(type_potential p, int temperature, int mass, int* state_GR) { //функция взаимодействия
		if (interact) {
			double coefK = 10.0 * coef_C / mass; //1
			double coefGR = 0.1 * coef_C / mass; //0.1
			double coefLJ = 0.000001 * coef_C / mass / num / num; //0.0000001
			double coefGamma_GR = 0.02;
			double coefGamma_K = 0.002;
			double coefGamma_LJ = 0.001;
			double state_1_temp_coef = 0.004;
			double sigmaLJ = 0.01;
			srand((unsigned int)time(NULL));

			if (p == K) {
				accel[0] = -coefK / pow(dist(0, 1), 2) -
					(sign(accel[0]) * coefGamma_K * speed[0] * sign(speed[0])) +
					sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				accel[num - 1] = coefK / pow(dist(num - 1, num - 2), 2) -
					(sign(accel[num - 1]) * coefGamma_K * speed[num - 1] * sign(speed[num - 1])) +
					sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				for (int i = 1; i < num - 1; ++i) {
					accel[i] = coefK / pow(dist(i, i - 1), 2) -
						coefK / pow(dist(i, i + 1), 2) -
						(sign(accel[i]) * coefGamma_K * speed[i] * sign(speed[i])) +
						sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				}
			}

			if (p == GR) {
				if (state_GR[0] == 0) {
					if (dist(0, 1) > cGR) {
						accel[0] = coefGR / pow(dist(0, 1), 2) -
							(sign(accel[0]) * coefGamma_GR * speed[0] * sign(speed[0])) +
							sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
					}
					else {
						state_GR[0] = 1;
						state_GR[1] = 1;
						speed[0] = 0;
						speed[1] = 0;
						accel[0] = 0;
						accel[1] = 0;
					}
				}
				else {
					if (x[0] > cGR)
						accel[0] = sqrt(temperature) * state_1_temp_coef * ((double)(rand() % 2) - 1);
					else
						accel[0] = 0;
					if (dist(0, 1) > cGR) {
						state_GR[0] = 0;
						state_GR[1] = 0;
					}
				}
				if (state_GR[num - 1] == 0) {
					if (dist(num - 1, num - 2) > cGR) {
						accel[num - 1] = -coefGR / pow(dist(num - 1, num - 2), 2) -
							(sign(accel[num - 1]) * coefGamma_GR * speed[num - 1] * sign(speed[num - 1])) +
							sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
					}
					else {
						state_GR[num - 1] = 1;
						state_GR[num - 2] = 1;
						speed[num - 1] = 0;
						speed[num - 2] = 0;
						accel[num - 1] = 0;
						accel[num - 2] = 0;
					}
				}
				else {
					if (x[num - 1] < (double)screenWidth - cGR)
						accel[num - 1] = sqrt(temperature) * state_1_temp_coef * ((double)(rand() % 2));
					else
						accel[num - 1] = 0;
					if (dist(num - 1, num - 2) > cGR) {
						state_GR[num - 1] = 0;
						state_GR[num - 2] = 0;
					}
				}
				for (int i = 1; i < num - 1; ++i) {
					if (state_GR[i] == 0) {
						if ((dist(i, i + 1) > cGR) && (dist(i, i - 1) > cGR)) {
							accel[i] = -coefGR / pow(dist(i, i - 1), 2) +
								coefGR / pow(dist(i, i + 1), 2) -
								(sign(accel[i]) * coefGamma_GR * speed[i] * sign(speed[i])) +
								sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
						}
						else if ((dist(i, i + 1) < cGR) && (dist(i, i - 1) > cGR)) {
							state_GR[i] = 1;
							state_GR[i + 1] = 1;
							speed[i] = 0;
							speed[i + 1] = 0;
							accel[i] = 0;
							accel[i + 1] = 0;
						}
						else if ((dist(i, i + 1) > cGR) && (dist(i, i - 1) < cGR)) {
							state_GR[i] = 1;
							state_GR[i - 1] = 1;
							speed[i] = 0;
							speed[i - 1] = 0;
							accel[i] = 0;
							accel[i - 1] = 0;
						}
						else if ((dist(i, i + 1) < cGR) && (dist(i, i - 1) < cGR)) {
							state_GR[i] = 1;
							state_GR[i - 1] = 1;
							state_GR[i + 1] = 1;
							speed[i] = 0;
							speed[i - 1] = 0;
							speed[i + 1] = 0;
							accel[i] = 0;
							accel[i - 1] = 0;
							accel[i + 1] = 0;
						}
					}
					else {
						if ((dist(i, i + 1) < (double)cGR + 1) && (dist(i, i - 1) > cGR))
							accel[i] = sqrt(temperature) * state_1_temp_coef * ((double)(rand() % 2) - 1)
							- coefGR / pow(dist(i, i - 1), 2);
						if ((dist(i, i + 1) > cGR) && (dist(i, i - 1) < (double)cGR + 1))
							accel[i] = sqrt(temperature) * state_1_temp_coef * ((double)(rand() % 2))
							+ coefGR / pow(dist(i, i + 1), 2);
						if ((dist(i, i + 1) < (double)cGR + 1) && (dist(i, i - 1) < (double)cGR + 1))
							accel[i] = 0;
						if ((dist(i, i + 1) > cGR) && (dist(i, i - 1) > cGR)) {
							state_GR[i] = 0;
						}
					}
				}
			}

			if (p == LJ) {
				accel[0] = coefLJ / pow(dist(0, 1) * sigmaLJ, 7) - coefLJ / pow(dist(0, 1) * sigmaLJ, 13) -
					(sign(accel[0]) * coefGamma_LJ * speed[0] * sign(speed[0])) +
					sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				accel[num - 1] = -coefLJ / pow(dist(num - 1, num - 2) * sigmaLJ, 7) + coefLJ / pow(dist(num - 1, num - 2) * sigmaLJ, 13) -
					(sign(accel[num - 1]) * coefGamma_LJ * speed[num - 1] * sign(speed[num - 1])) +
					sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				for (int i = 1; i < num - 1; ++i) {
					accel[i] = -coefLJ / pow(dist(i, i - 1) * sigmaLJ, 7) + coefLJ / pow(dist(i, i - 1) * sigmaLJ, 13) +
						coefLJ / pow(dist(i, i + 1) * sigmaLJ, 7) - coefLJ / pow(dist(i, i + 1) * sigmaLJ, 13) -
						(sign(accel[i]) * coefGamma_LJ * speed[i] * sign(speed[i])) +
						sqrt(temperature) * 0.00001 * ((double)(rand() % 3) - 1);
				}
			}

			if (p == NUL) {
				for (int i = 0; i < num; ++i) {
					accel[i] = -(sign(accel[i]) * coefGamma_LJ * speed[i] * sign(speed[i])) +
						sqrt(temperature) * 0.000003 * ((double)(rand() % 3) - 1); //0.00001
				}
			}

			for (int i = 0; i < num; ++i) {
				if (abs(accel[i]) > max_accel)
					accel[i] = sign(accel[i]) * max_accel;
			}

		}
	}


	double& get_w() { return w; }
	double& get_h() { return h; }
	int get_num() { return num; }
	int get_nbins() { return nbins; }
	RectangleShape& get_p(int i) { return p[i]; }
	RectangleShape& get_bin(int i) { return bin[i]; }

	int get_maxAmount_now() {
		int max = 0;
		for (int i = 0; i < nbins; ++i)
			if (max < amount_now[i])
				max = amount_now[i];
		return max;
	}

	void set_params(double w_, double h_, int num_) {
		w = w_; h = h_; num = num_;
		if (p) delete[] p;
		p = new RectangleShape[num];
		if (x) delete[] x;
		x = new double[num];
		if (speed) delete[] speed;
		speed = new double[num]();
		if (accel) delete[] accel;
		accel = new double[num]();
		nbins = nBins;
		if (bin) delete[] bin;
		bin = new RectangleShape[nbins];
		if (amount) delete[] amount;
		amount = new double[nbins]();
		if (amount_now) delete[] amount_now;
		amount_now = new int[nbins]();
		for (int i = 0; i < num; ++i) {
			p[i].setSize(Vector2f((float)w, (float)h));
			p[i].setOutlineThickness(2);
			p[i].setOutlineColor(Color::Black);
			p[i].setPosition((float)(x[i] = (((double)i + 1) * (double)(screenWidth / (num + 1)))), 70 * screenHeight / 1080);
			p[i].setFillColor(Color::Blue);
		}

		for (int i = 0; i < nbins; ++i) {
			bin[i].setOutlineThickness(1);
			bin[i].setOutlineColor(Color::Black);
			bin[i].setRotation(180);
			bin[i].setPosition((plotPosX + 5) * screenWidth / 1920 + (i + 1) * (plotSize * screenWidth / 1920 / nbins), plotPosY * screenHeight / 1080);
			bin[i].setFillColor(Color::Green);
		}
	}

	void plot_hist() {
		double max = 0;
		for (int i = 0; i < nbins; ++i)
			if (max < amount[i])
				max = amount[i];
		for (int i = 0; i < nbins; ++i)
			bin[i].setSize(Vector2f((plotSize - 20) * screenWidth / 1920 / nbins, (float)(amount[i] / max) * binCoef * screenHeight / 1080));
	}

	void pause() { interact = false; }

	void resume() { interact = true; }

};

int main() {

	ContextSettings settings;
	settings.antialiasingLevel = 8;

	int state = 0, num = 15, temperature = 25, mass = 25, width = 22;
	bool button_pressed[27] = { 0 }, sliders_pressed[6] = { 0 }, show = false, stop = false, init = true;
	type_potential potential = NUL;

	button_pressed[26] = true;

	float dX1 = 0, dX2 = 0, dX3 = 0, dX4 = 0, dX5 = 0, dY = 0;

	Font font;
	font.loadFromFile("fonts/ArialCyr.ttf");

	RenderWindow* window = new RenderWindow;
	window->create(VideoMode((unsigned int)screenWidth, (unsigned int)screenHeight), "statfiz", Style::Fullscreen, settings);

	RectangleShape tube(Vector2f(screenWidth, 140 * screenHeight / 1080));
	tube.setPosition(0, 70 * screenHeight / 1080);
	tube.setOutlineColor(Color(0, 0, 0));
	tube.setOutlineThickness(2);

	RectangleShape lineX_graph_K;
	lineX_graph_K.setSize(Vector2f(450 * screenWidth / 1920, 4 * screenHeight / 1080));
	lineX_graph_K.setPosition(1350 * screenWidth / 1920, 690 * screenHeight / 1080);
	lineX_graph_K.setFillColor(Color(0, 0, 0));

	RectangleShape lineX_graph_GR;
	lineX_graph_GR.setSize(Vector2f(450 * screenWidth / 1920, 4 * screenHeight / 1080));
	lineX_graph_GR.setPosition(1350 * screenWidth / 1920, 340 * screenHeight / 1080);
	lineX_graph_GR.setFillColor(Color(0, 0, 0));

	RectangleShape lineX_graph_LJ;
	lineX_graph_LJ.setSize(Vector2f(450 * screenWidth / 1920, 4 * screenHeight / 1080));
	lineX_graph_LJ.setPosition(1350 * screenWidth / 1920, 530 * screenHeight / 1080);
	lineX_graph_LJ.setFillColor(Color(0, 0, 0));

	RectangleShape lineY_graph;
	lineY_graph.setSize(Vector2f(4 * screenWidth / 1920, 400 * screenHeight / 1080));
	lineY_graph.setPosition(1350 * screenWidth / 1920, 290 * screenHeight / 1080);
	lineY_graph.setFillColor(Color(0, 0, 0));

	double xGR, xK, xLJ;
	VertexArray graph_GR(LinesStrip, 3900);//график Притягивание
	VertexArray graph_GR1(LinesStrip, 3900);//график Притягивание
	VertexArray graph_GR2(LinesStrip, 3900);//график Притягивание
	VertexArray graph_K(LinesStrip, 3900);//график K
	VertexArray graph_K1(LinesStrip, 3900);//график K
	VertexArray graph_K2(LinesStrip, 3900);//график K
	VertexArray graph_LJ(LinesStrip, 3900);//график LJ
	VertexArray graph_LJ1(LinesStrip, 3900);//график LJ
	VertexArray graph_LJ2(LinesStrip, 3900);//график LJ

	std::stringstream buf;

	RectangleShape lineX;
	lineX.setSize(Vector2f((plotSize + 50) * screenWidth / 1920, 6 * screenHeight / 1080));
	lineX.setPosition(plotPosX * screenWidth / 1920, plotPosY * screenHeight / 1080); //plotPosX=900 plotPosY = 1000
	lineX.setFillColor(Color(0, 0, 0));

	RectangleShape lineY;
	lineY.setSize(Vector2f(6 * screenWidth / 1920, (binCoef + 60) * screenHeight / 1080));
	lineY.setPosition(plotPosX * screenWidth / 1920, (plotPosY - (binCoef + 60)) * screenHeight / 1080); //plotPosX=900 plotPosY = 1000
	lineY.setFillColor(Color(0, 0, 0));

	Text textX1, textX2, textX3, textX4;
	textX1.setFont(font);
	textX1.setFillColor(Color(0, 0, 0));
	textX1.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textX1.setString(L"Расстояния");
	textX1.setPosition((plotPosX + (plotSize + 57)) * screenWidth / 1920, (plotPosY - 113) * screenHeight / 1080);
	textX2.setFont(font);
	textX2.setFillColor(Color(0, 0, 0));
	textX2.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textX2.setString(L"между");
	textX2.setPosition((plotPosX + (plotSize + 57)) * screenWidth / 1920, (plotPosY - 77) * screenHeight / 1080);
	textX3.setFont(font);
	textX3.setFillColor(Color(0, 0, 0));
	textX3.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textX3.setString(L"соседними");
	textX3.setPosition((plotPosX + (plotSize + 57)) * screenWidth / 1920, (plotPosY - 41) * screenHeight / 1080);
	textX4.setFont(font);
	textX4.setFillColor(Color(0, 0, 0));
	textX4.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textX4.setString(L"частицами R");
	textX4.setPosition((plotPosX + (plotSize + 57)) * screenWidth / 1920, (plotPosY - 5) * screenHeight / 1080);
	Text textY1, textY2, textY3;
	textY1.setFont(font);
	textY1.setFillColor(Color(0, 0, 0));
	textY1.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textY1.setString(L"Число");
	textY1.setPosition((plotPosX + 36) * screenWidth / 1920, (plotPosY - (binCoef + 140)) * screenHeight / 1080);
	textY2.setFont(font);
	textY2.setFillColor(Color(0, 0, 0));
	textY2.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textY2.setString(L"попаданий R");
	textY2.setPosition((plotPosX + 18) * screenWidth / 1920, (plotPosY - (binCoef + 106)) * screenHeight / 1080);
	textY3.setFont(font);
	textY3.setFillColor(Color(0, 0, 0));
	textY3.setCharacterSize((unsigned int)(34 * screenWidth / 1920));
	textY3.setString(L"в интервалы");
	textY3.setPosition((plotPosX + 18) * screenWidth / 1920, (plotPosY - (binCoef + 72)) * screenHeight / 1080);
	Text* xtextBins = new Text[(nBins + 1) / 4];
	for (int i = 0; i < (nBins + 1) / 4; ++i) {
		xtextBins[i].setFont(font);
		xtextBins[i].setFillColor(Color(0, 0, 0));
		xtextBins[i].setCharacterSize((unsigned int)(36 * screenWidth / 1920));
		buf << 4 * i;
		xtextBins[i].setString(buf.str());
		buf.str(std::string());
		xtextBins[i].setPosition(plotPosX * screenWidth / 1920 + 4 * i * (plotSize * screenWidth / 1920 / nBins), (plotPosY + 8) * screenHeight / 1080);
	}

	Text ytextBin;
	ytextBin.setFont(font);
	ytextBin.setFillColor(Color(0, 0, 0));
	ytextBin.setCharacterSize((unsigned int)(36 * screenWidth / 1920));
	ytextBin.setPosition((plotPosX - 48) * screenWidth / 1920, (plotPosY - 20) * screenHeight / 1080 - (nBins) * (binCoef) * (screenHeight / 1080) / nBins);

	RectangleShape* xticks = new RectangleShape[nBins];
	for (int i = 0; i < nBins; ++i) {
		xticks[i].setFillColor(Color(0, 0, 0));
		xticks[i].setSize(Vector2f(3 * screenWidth / 1920, 20 * screenHeight / 1080));
		xticks[i].setPosition((plotPosX + 5) * screenWidth / 1920 + (i + 1) * (plotSize * screenWidth / 1920 / nBins), (plotPosY - 7) * (screenHeight / 1080));
	}

	RectangleShape* yticks = new RectangleShape[nBins];
	for (int i = 0; i < nBins; ++i) {
		yticks[i].setFillColor(Color(0, 0, 0));
		yticks[i].setSize(Vector2f(20 * screenWidth / 1920, 3 * screenHeight / 1080));
		yticks[i].setPosition((plotPosX - 7) * screenWidth / 1920, plotPosY * screenHeight / 1080 - (i + 1) * (binCoef) * (screenHeight / 1080) / nBins);
	}

	RectangleShape* button = new RectangleShape[num_buttons];
	button[0].setSize(Vector2f(478 * screenWidth / 1920, 64 * screenHeight / 1080));//Начальная страница
	button[1].setSize(Vector2f(478 * screenWidth / 1920, 64 * screenHeight / 1080));//Демонстрация
	button[2].setSize(Vector2f(478 * screenWidth / 1920, 64 * screenHeight / 1080));//Теория
	button[3].setSize(Vector2f(478 * screenWidth / 1920, 64 * screenHeight / 1080));//Авторы
	button[4].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//поле для задания ширины частиц
	button[5].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//поле для задания числа частиц
	button[6].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//старт
	button[7].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//выход
	button[8].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//Притягивание
	button[9].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//Кулон
	button[10].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//пауза
	button[11].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//поле для задания температуры
	button[12].setSize(Vector2f(246 * screenWidth / 1920, 0.1 * screenHeight / 1080));//поле для ползунка для задания температуры
	button[13].setSize(Vector2f(246 * screenWidth / 1920, 0.1 * screenHeight / 1080));//поле для ползунка для задания ширины частиц
	button[14].setSize(Vector2f(246 * screenWidth / 1920, 0.1 * screenHeight / 1080));//поле для ползунка для задания числа частиц
	button[18].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//поле для задания массы
	button[19].setSize(Vector2f(246 * screenWidth / 1920, 0.1 * screenHeight / 1080));//поле для ползунка для задания массы
	button[21].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//поле для задания с
	button[22].setSize(Vector2f(246 * screenWidth / 1920, 0.1 * screenHeight / 1080));//поле для ползунка с
	button[24].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//сброс
	button[25].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//Леннард-Джонс
	button[26].setSize(Vector2f(360 * screenWidth / 1920, 64 * screenHeight / 1080));//NUL
	button[27].setSize(Vector2f(64 * screenWidth / 1920, screenHeight));//поле для ползунка для скроллинга теории

	button[0].setPosition(0, 0);                                              //Начальная страница
	button[1].setPosition(484 * screenWidth / 1920, 0);                       //Демонстрация
	button[2].setPosition(968 * screenWidth / 1920, 0);                       //Теория
	button[3].setPosition(1452 * screenWidth / 1920, 0);                       //Авторы
	button[4].setPosition(20 * screenWidth / 1920, 290 * screenHeight / 1080);//поле для задания ширины частиц
	button[5].setPosition(20 * screenWidth / 1920, 360 * screenHeight / 1080);//поле для задания числа частиц
	button[6].setPosition(20 * screenWidth / 1920, 640 * screenHeight / 1080);//старт
	button[7].setPosition(20 * screenWidth / 1920, 990 * screenHeight / 1080);//выход
	button[8].setPosition(900 * screenWidth / 1920, 290 * screenHeight / 1080);//Притягивание
	button[9].setPosition(900 * screenWidth / 1920, 360 * screenHeight / 1080);//Кулон
	button[10].setPosition(20 * screenWidth / 1920, 710 * screenHeight / 1080);//пауза
	button[11].setPosition(20 * screenWidth / 1920, 430 * screenHeight / 1080);//поле для задания температуры
	button[12].setPosition(492 * screenWidth / 1920, 496 * screenHeight / 1080);//поле для ползунка для задания температуры
	button[13].setPosition(492 * screenWidth / 1920, 426 * screenHeight / 1080);//поле для ползунка для задания ширины частиц
	button[14].setPosition(492 * screenWidth / 1920, 356 * screenHeight / 1080);//поле для ползунка для задания числа частиц
	button[18].setPosition(20 * screenWidth / 1920, 500 * screenHeight / 1080);//поле для задания массы
	button[19].setPosition(492 * screenWidth / 1920, 566 * screenHeight / 1080);//поле для ползунка для задания массы
	button[21].setPosition(20 * screenWidth / 1920, 570 * screenHeight / 1080);//поле для задания с
	button[22].setPosition(492 * screenWidth / 1920, 636 * screenHeight / 1080);//поле для ползунка для задания с
	button[24].setPosition(20 * screenWidth / 1920, 780 * screenHeight / 1080);//сброс
	button[25].setPosition(900 * screenWidth / 1920, 430 * screenHeight / 1080);//Леннард-Джонс
	button[26].setPosition(900 * screenWidth / 1920, 500 * screenHeight / 1080);//NUL//500
	button[27].setPosition(1854 * screenWidth / 1920, 70 * screenHeight / 1080);//поле для ползунка для скроллинга теории

	for (int i = 0; i < num_buttons; ++i) {
		button[i].setOutlineThickness(2);
		button[i].setOutlineColor(Color(0, 0, 0));
		button[i].setFillColor(Color(C1_0, C2_0, C3_0));
	}

	button[12].setOutlineThickness(1);
	button[13].setOutlineThickness(1);
	button[14].setOutlineThickness(1);
	button[19].setOutlineThickness(1);
	button[22].setOutlineThickness(1);

	button[0].setFillColor(Color(C1_1, C2_1, C3_1));

	ConvexShape* Sliders = new ConvexShape[6];
	for (int i = 0; i < 6; i++) {
		Sliders[i].setPointCount(5);
		Sliders[i].setOutlineColor(Color(0, 0, 0));
		Sliders[i].setOutlineThickness(2);
		Sliders[i].setPoint(0, Vector2f(0, 0));
		Sliders[i].setPoint(1, Vector2f(64 * screenWidth / 1920, 0));
		Sliders[i].setPoint(2, Vector2f(64 * screenWidth / 1920, 32 * screenHeight / 1080));
		Sliders[i].setPoint(3, Vector2f(32 * screenWidth / 1920, 64 * screenHeight / 1080));
		Sliders[i].setPoint(4, Vector2f(0, 32 * screenHeight / 1080));
	}
	Sliders[0].setPosition(583 * screenWidth / 1920, 360 * screenHeight / 1080);//583
	Sliders[1].setPosition(583 * screenWidth / 1920, 290 * screenHeight / 1080);
	Sliders[2].setPosition(583 * screenWidth / 1920, 430 * screenHeight / 1080);
	Sliders[3].setPosition(583 * screenWidth / 1920, 500 * screenHeight / 1080);
	Sliders[4].setPosition(460 * screenWidth / 1920, 570 * screenHeight / 1080);
	Sliders[5].setPosition(1854 * screenWidth / 1920, 70 * screenHeight / 1080);

	Text* text = new Text[num_texts];

	for (int i = 0; i < num_texts; ++i) {
		text[i].setFont(font);
		text[i].setFillColor(Color(0, 0, 0));
		text[i].setCharacterSize((unsigned int)(36 * screenWidth / 1920));
	}

	text[53].setCharacterSize((unsigned int)(26 * screenWidth / 1920));
	text[54].setCharacterSize((unsigned int)(26 * screenWidth / 1920));

	text[0].setString(L"Начальная страница");
	text[0].setPosition(40, 16 * screenHeight / 1080);
	text[1].setString(L"Демонстрация");
	text[1].setPosition(594 * screenWidth / 1920, 16 * screenHeight / 1080);
	text[2].setString(L"Теория");
	text[2].setPosition(1130 * screenWidth / 1920, 16 * screenHeight / 1080);
	text[3].setString(L"Авторы");
	text[3].setPosition(1635 * screenWidth / 1920, 16 * screenHeight / 1080);
	text[4].setString(L"2019");
	text[4].setPosition(920 * screenWidth / 1920, 1000 * screenHeight / 1080);
	text[5].setString(L"Флуктуации плотности в одномерном газе");
	text[5].setOutlineColor(Color(0, 0, 0));
	text[5].setOutlineThickness(1);
	text[5].setPosition(625 * screenWidth / 1920, 300 * screenHeight / 1080);
	text[6].setString(L"Ширина частиц:");
	text[6].setPosition(44 * screenWidth / 1920, 368 * screenHeight / 1080);
	text[7].setString(L"Число частиц:");
	text[7].setPosition(60 * screenWidth / 1920, 298 * screenHeight / 1080);
	text[8].setString(L"Старт");
	text[8].setPosition(144 * screenWidth / 1920, 648 * screenHeight / 1080);
	text[9].setString(L"Выход");
	text[9].setPosition(144 * screenWidth / 1920, 998 * screenHeight / 1080);
	text[10].setPosition(319 * screenWidth / 1920, 368 * screenHeight / 1080);//ширина частиц
	text[11].setPosition(304 * screenWidth / 1920, 298 * screenHeight / 1080);//число частиц
	text[12].setString(L"Выберите потенциал:");
	text[12].setPosition(904 * screenWidth / 1920, 228 * screenHeight / 1080);
	text[13].setString(L"U(R) = -C / R");
	text[13].setPosition(1500 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[14].setString(L"U(R) = C / R");
	text[14].setPosition(1500 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[15].setString(L"Задайте параметры:");
	text[15].setPosition(24 * screenWidth / 1920, 228 * screenHeight / 1080);
	text[16].setString(L"Пауза");
	text[16].setPosition(144 * screenWidth / 1920, 718 * screenHeight / 1080);
	text[17].setString(L"Температура:");
	text[17].setPosition(64 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[18].setPosition(300 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[19].setString(L"1");
	text[19].setPosition(424 * screenWidth / 1920, 438 * screenHeight / 1080);
	buf << max_temperature;
	text[20].setString(buf.str());
	buf.str(std::string());
	text[20].setPosition(784 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[21].setString(L"График потенциала U(R)");
	text[21].setPosition(1354 * screenWidth / 1920, 228 * screenHeight / 1080);
	text[22].setString(L"1");
	text[22].setPosition(424 * screenWidth / 1920, 368 * screenHeight / 1080);
	text[23].setPosition(784 * screenWidth / 1920, 368 * screenHeight / 1080);
	text[24].setString(L"2");
	text[24].setPosition(424 * screenWidth / 1920, 298 * screenHeight / 1080);
	buf << (max_num + 1);
	text[25].setString(buf.str());
	buf.str(std::string());
	text[25].setPosition(784 * screenWidth / 1920, 298 * screenHeight / 1080);
	text[26].setString(L"Масса частицы:");
	text[26].setPosition(50 * screenWidth / 1920, 508 * screenHeight / 1080);
	text[27].setString(L"1");
	text[27].setPosition(424 * screenWidth / 1920, 508 * screenHeight / 1080);
	buf << max_mass;
	text[28].setString(buf.str());
	buf.str(std::string());
	text[28].setPosition(784 * screenWidth / 1920, 508 * screenHeight / 1080);
	text[29].setPosition(324 * screenWidth / 1920, 508 * screenHeight / 1080);
	text[30].setString(L"E (среднее)");
	text[30].setPosition(410 * screenWidth / 1920, 700 * screenHeight / 1080);
	text[31].setString(L"(cтандартное");
	text[31].setPosition(400 * screenWidth / 1920, 805 * screenHeight / 1080);
	text[32].setPosition(650 * screenWidth / 1920, 700 * screenHeight / 1080);//average
	text[33].setPosition(650 * screenWidth / 1920, 805 * screenHeight / 1080);//stdevp
	text[34].setPosition(42 * screenWidth / 1920, 578 * screenHeight / 1080);//Коэффициент C=
	text[34].setString(L"Коэффициент C:");
	text[35].setPosition(338 * screenWidth / 1920, 578 * screenHeight / 1080);
	text[36].setString(L"1");
	text[36].setPosition(424 * screenWidth / 1920, 578 * screenHeight / 1080);
	buf << max_coef_C;
	text[37].setString(buf.str());
	buf.str(std::string());
	text[37].setPosition(784 * screenWidth / 1920, 578 * screenHeight / 1080);
	text[38].setString(L"МГУ им. М. В. Ломоносова");
	text[38].setPosition(740 * screenWidth / 1920, 90 * screenHeight / 1080);
	text[39].setString(L"Компьютерные демонстрации по курсу лекций");
	text[39].setPosition(590 * screenWidth / 1920, 150 * screenHeight / 1080);
	text[40].setString(L"Статистическая физика");
	text[40].setPosition(765 * screenWidth / 1920, 210 * screenHeight / 1080);
	text[41].setString(L"Артём Ахметов");
	text[41].setPosition(166 * screenWidth / 1920, 770 * screenHeight / 1080);
	text[42].setString(L"Никита Кирюшкин");
	text[42].setPosition(1456 * screenWidth / 1920, 770 * screenHeight / 1080);
	text[43].setString(L"Научный руководитель:");
	text[43].setPosition(785 * screenWidth / 1920, 865 * screenHeight / 1080);
	text[44].setString(L"Ольга Александровна Чичигина");
	text[44].setPosition(715 * screenWidth / 1920, 915 * screenHeight / 1080);
	text[45].setString(L"Факультет ВМК");
	text[45].setPosition(13 * screenWidth / 1920, 240 * screenHeight / 1080);
	text[46].setString(L"Физический факультет");
	text[46].setPosition(1525 * screenWidth / 1920, 240 * screenHeight / 1080);
	text[47].setString(L"отклонение)");
	text[47].setPosition(415 * screenWidth / 1920, 840 * screenHeight / 1080);
	text[48].setString(L"Сброс");
	text[48].setPosition(144 * screenWidth / 1920, 788 * screenHeight / 1080);
	text[49].setString(L"U(R) = C/R  - C/R ");//4*eps*((sigma / R) ^ 12 - (sigma / R) ^ 6)
	text[49].setPosition(1500 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[50].setString(L"U(R) = 0");
	text[50].setPosition(1500 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[51].setString(L"FPS:");
	text[51].setPosition(1750 * screenWidth / 1920, 80 * screenHeight / 1080);
	text[52].setPosition(1850 * screenWidth / 1920, 80 * screenHeight / 1080);
	text[53].setString(L"12");
	text[53].setPosition(1674 * screenWidth / 1920, 428 * screenHeight / 1080);
	text[54].setString(L"6");
	text[54].setPosition(1780 * screenWidth / 1920, 428 * screenHeight / 1080);
	text[55].setString(L"Притяжение");
	text[55].setPosition(974 * screenWidth / 1920, 298 * screenHeight / 1080);
	text[56].setString(L"Отталкивание");
	text[56].setPosition(964 * screenWidth / 1920, 368 * screenHeight / 1080);
	text[57].setString(L"Леннард-Джонс");
	text[57].setPosition(944 * screenWidth / 1920, 438 * screenHeight / 1080);
	text[58].setString(L"Нулевой");
	text[58].setPosition(1004 * screenWidth / 1920, 508 * screenHeight / 1080);
	text[59].setString(L"sigma / E");
	text[59].setPosition(400 * screenWidth / 1920, 910 * screenHeight / 1080);
	text[60].setPosition(650 * screenWidth / 1920, 910 * screenHeight / 1080);
	text[61].setString(L"sigma");
	text[61].setPosition(458 * screenWidth / 1920, 770 * screenHeight / 1080);

	Image* iCMC = new Image;
	iCMC->loadFromFile("images/CMC.png");
	Texture* tCMC = new Texture;
	tCMC->loadFromImage(*iCMC);
	Sprite* sCMC = new Sprite;
	sCMC->setTexture(*tCMC);
	sCMC->setPosition(70 * screenWidth / 1920, 80 * screenHeight / 1080);
	sCMC->setScale(screenWidth / 1920, screenHeight / 1080);

	Image* iFF = new Image;
	iFF->loadFromFile("images/FF.png");
	Texture* tFF = new Texture;
	tFF->loadFromImage(*iFF);
	Sprite* sFF = new Sprite;
	sFF->setTexture(*tFF);
	sFF->setPosition(1620 * screenWidth / 1920, 80 * screenHeight / 1080);
	sFF->setScale(screenWidth / 1920, screenHeight / 1080);

	Image* iTh = new Image;
	iTh->loadFromFile("images/Theory.png");
	Texture* tTh = new Texture;
	tTh->loadFromImage(*iTh);
	Sprite* sTh = new Sprite;
	sTh->setTexture(*tTh);
	sTh->setPosition(370 * screenWidth / 1920, 70 * screenHeight / 1080);

	Image* iPhoto1 = new Image;
	iPhoto1->loadFromFile("images/photo1.png");
	Texture* tPhoto1 = new Texture;
	tPhoto1->loadFromImage(*iPhoto1);
	Sprite* sPhoto1 = new Sprite;
	sPhoto1->setTexture(*tPhoto1);
	sPhoto1->setPosition(166 * screenWidth / 1920, 420 * screenHeight / 1080);
	sPhoto1->setScale(screenWidth / 1920, screenHeight / 1080);

	Image* iPhoto2 = new Image;
	iPhoto2->loadFromFile("images/photo2.png");
	Texture* tPhoto2 = new Texture;
	tPhoto2->loadFromImage(*iPhoto2);
	Sprite* sPhoto2 = new Sprite;
	sPhoto2->setTexture(*tPhoto2);
	sPhoto2->setPosition(1475 * screenWidth / 1920, 420 * screenHeight / 1080);
	sPhoto2->setScale(screenWidth / 1920, screenHeight / 1080);

	Particle p;
	Clock clock4Int, clock4Stat;
	//Clock clock4FPS;
	//int fps = 0;
	int maxAmount;
	float time = 0;

	clock4Int.restart();
	clock4Stat.restart();
	//clock4FPS.restart();
	p.set_params(width * screenWidth / 1920, 140 * screenHeight / 1080, num);
	show = true;

	while (window->isOpen()) {
		Vector2i pixelPos = Mouse::getPosition(*window);
		Event event;
		if (crashed) {
			window->close();
		}
		while (window->pollEvent(event)) {
			if (button[0].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[0].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[0].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[1].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[1].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[1].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[2].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[2].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[2].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[3].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[3].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[3].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[6].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[6].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[6].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[7].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[7].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[7].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[8].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[8].setFillColor(Color(C1_2, C2_2, C3_2));
			else if (button_pressed[8])
				button[8].setFillColor(Color(C1_1, C2_1, C3_1));
			else
				button[8].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[9].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[9].setFillColor(Color(C1_2, C2_2, C3_2));
			else if (button_pressed[9])
				button[9].setFillColor(Color(C1_1, C2_1, C3_1));
			else
				button[9].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[10].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[10].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[10].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[24].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[24].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				button[24].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[25].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[25].setFillColor(Color(C1_2, C2_2, C3_2));
			else if (button_pressed[25])
				button[25].setFillColor(Color(C1_1, C2_1, C3_1));
			else
				button[25].setFillColor(Color(C1_0, C2_0, C3_0));
			if (button[26].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				button[26].setFillColor(Color(C1_2, C2_2, C3_2));
			else if (button_pressed[26])
				button[26].setFillColor(Color(C1_1, C2_1, C3_1));
			else
				button[26].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[0].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[0].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[0].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[1].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[1].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[1].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[2].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[2].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[2].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[3].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[3].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[3].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[4].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[4].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[4].setFillColor(Color(C1_0, C2_0, C3_0));
			if (Sliders[5].getGlobalBounds().contains(pixelPos.x, pixelPos.y))
				Sliders[5].setFillColor(Color(C1_2, C2_2, C3_2));
			else
				Sliders[5].setFillColor(Color(C1_0, C2_0, C3_0));

			if (event.type == Event::MouseButtonPressed) {
				if (event.mouseButton.button == Mouse::Left) {
					if (button[0].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//начальная страница
						button[state].setFillColor(Color(C1_0, C2_0, C3_0));
						state = 0;
						button[0].setFillColor(Color(C1_1, C2_1, C3_1));
					}
					else if (button[1].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//демонстрация
						button[state].setFillColor(Color(C1_0, C2_0, C3_0));
						state = 1;
						button[1].setFillColor(Color(C1_1, C2_1, C3_1));
					}
					else if (button[2].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//теория
						button[state].setFillColor(Color(C1_0, C2_0, C3_0));
						state = 2;
						button[2].setFillColor(Color(C1_1, C2_1, C3_1));
					}
					else if (button[3].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//авторы
						button[state].setFillColor(Color(C1_0, C2_0, C3_0));
						state = 3;
						button[3].setFillColor(Color(C1_1, C2_1, C3_1));
					}
					else if ((!show || stop) && button[6].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//старт
						button[4].setFillColor(Color(C1_0, C2_0, C3_0));
						button[5].setFillColor(Color(C1_0, C2_0, C3_0));
						button[8].setFillColor(Color(C1_0, C2_0, C3_0));
						button[9].setFillColor(Color(C1_0, C2_0, C3_0));
						button[10].setFillColor(Color(C1_0, C2_0, C3_0));
						button[6].setFillColor(Color(C1_1, C2_1, C3_1));
						clock4Int.restart();
						clock4Stat.restart();
						//clock4FPS.restart();
						if (init) {
							p.set_params(width * screenWidth / 1920, 140 * screenHeight / 1080, num);
							init = false;
						}
						p.resume();
						show = true;
						stop = false;
					}
					if (button[7].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//Выход
						window->close();
					}
					else if ((!show || stop) && button[8].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//Притягивание
						button_pressed[8] = true;
						button[8].setFillColor(Color(C1_1, C2_1, C3_1));
						button[9].setFillColor(Color(C1_0, C2_0, C3_0));
						button[25].setFillColor(Color(C1_0, C2_0, C3_0));
						button[26].setFillColor(Color(C1_0, C2_0, C3_0));
						button[10].setFillColor(Color(C1_0, C2_0, C3_0));
						button_pressed[9] = false;
						button_pressed[25] = false;
						button_pressed[26] = false;
						potential = GR;
						init = true;
					}
					else if ((!show || stop) && button[9].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//Кулон
						button_pressed[9] = true;
						button[9].setFillColor(Color(C1_1, C2_1, C3_1));
						button[8].setFillColor(Color(C1_0, C2_0, C3_0));
						button[25].setFillColor(Color(C1_0, C2_0, C3_0));
						button[26].setFillColor(Color(C1_0, C2_0, C3_0));
						button[10].setFillColor(Color(C1_0, C2_0, C3_0));
						button_pressed[8] = false;
						button_pressed[26] = false;
						button_pressed[25] = false;
						potential = K;
						init = true;
					}
					else if ((!show || stop) && button[25].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//Леннард-Джонс////
						button_pressed[25] = true;
						button[25].setFillColor(Color(C1_1, C2_1, C3_1));
						button[8].setFillColor(Color(C1_0, C2_0, C3_0));
						button[9].setFillColor(Color(C1_0, C2_0, C3_0));
						button[26].setFillColor(Color(C1_0, C2_0, C3_0));
						button[10].setFillColor(Color(C1_0, C2_0, C3_0));
						button_pressed[8] = false;
						button_pressed[9] = false;
						button_pressed[26] = false;
						potential = LJ;
						init = true;
					}
					else if ((!show || stop) && button[26].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//NUL
						button_pressed[26] = true;
						button[26].setFillColor(Color(C1_1, C2_1, C3_1));
						button[8].setFillColor(Color(C1_0, C2_0, C3_0));
						button[9].setFillColor(Color(C1_0, C2_0, C3_0));
						button[25].setFillColor(Color(C1_0, C2_0, C3_0));
						button[10].setFillColor(Color(C1_0, C2_0, C3_0));
						button_pressed[8] = false;
						button_pressed[9] = false;
						button_pressed[25] = false;
						potential = NUL;
						init = true;
					}
					else if (show && button[10].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//пауза
						button[10].setFillColor(Color(C1_1, C2_1, C3_1));
						button[6].setFillColor(Color(C1_0, C2_0, C3_0));
						stop = true;
					}
					else if (show && button[24].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {//сброс
						button[24].setFillColor(Color(C1_1, C2_1, C3_1));
						button[6].setFillColor(Color(C1_0, C2_0, C3_0));
						stop = true; show = false; init = true;
					}
					else if ((!show || stop) && Sliders[0].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[0].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[0] = true;
						dX1 = pixelPos.x - Sliders[0].getPosition().x;
					}
					else if ((!show || stop) && Sliders[1].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[1].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[1] = true;
						dX2 = pixelPos.x - Sliders[1].getPosition().x;
					}
					else if ((!show || stop) && Sliders[2].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[2].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[2] = true;
						dX3 = pixelPos.x - Sliders[2].getPosition().x;
					}
					else if ((!show || stop) && Sliders[3].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[3].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[3] = true;
						dX4 = pixelPos.x - Sliders[3].getPosition().x;
					}
					else if ((!show || stop) && Sliders[4].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[4].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[4] = true;
						dX5 = pixelPos.x - Sliders[4].getPosition().x;
					}
					else if (state == 2 && Sliders[5].getGlobalBounds().contains(pixelPos.x, pixelPos.y)) {
						Sliders[5].setFillColor(Color(C1_1, C2_1, C3_1));
						sliders_pressed[5] = true;
						dY = pixelPos.y - Sliders[5].getPosition().y;
					}
				}
			}
			if (event.type == Event::MouseButtonReleased) {
				if (event.mouseButton.button == Mouse::Left) {
					sliders_pressed[0] = false;
					sliders_pressed[1] = false;
					sliders_pressed[2] = false;
					sliders_pressed[3] = false;
					sliders_pressed[4] = false;
					sliders_pressed[5] = false;
					Sliders[0].setFillColor(Color(C1_0, C2_0, C3_0));
					Sliders[1].setFillColor(Color(C1_0, C2_0, C3_0));
					Sliders[2].setFillColor(Color(C1_0, C2_0, C3_0));
					Sliders[3].setFillColor(Color(C1_0, C2_0, C3_0));
					Sliders[4].setFillColor(Color(C1_0, C2_0, C3_0));
					Sliders[5].setFillColor(Color(C1_0, C2_0, C3_0));
				}
			}
			if (event.type == Event::Closed)
				window->close();
			if (state == 2 && event.type == Event::MouseWheelScrolled) {
				if (event.mouseWheelScroll.delta > 0) {
					if (sTh->getPosition().y + 20 * screenHeight / 1080 <= 70 * screenHeight / 1080) {
						sTh->move(0, 20 * screenHeight / 1080);
						Sliders[5].setPosition(1854 * screenWidth / 1920, (sTh->getPosition().y -
							(70 * screenHeight / 1080 - 70 * screenHeight / 1080 * ((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)))) /
							((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)));
					}
				}
				else if (event.mouseWheelScroll.delta < 0) {
					if (sTh->getPosition().y - 20 * screenHeight / 1080 >= -1060) {
						sTh->move(0, -20 * screenHeight / 1080);
						Sliders[5].setPosition(1854 * screenWidth / 1920, (sTh->getPosition().y -
							(70 * screenHeight / 1080 - 70 * screenHeight / 1080 * ((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)))) /
							((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)));
					}
				}
			}
		}
		window->clear(Color(C1_0, C2_0, C3_0));
		if (show) {
			if (!stop) {
				if ((float)clock4Int.getElapsedTime().asMilliseconds() > 11) {//11
					p.interaction();
					p.potential_interaction(potential, temperature, mass, state_GR);
					time = (float)clock4Int.getElapsedTime().asMilliseconds();
					clock4Int.restart();
					p.update(time);
					p.plot_hist();
					//++fps;
				}
				//if ((float)clock4FPS.getElapsedTime().asMilliseconds() > 1000.0) {
					//buf << fps;
					//text[52].setString(buf.str());
					//c_buf = r;
					//r = g;
					//g = b;
					//b = c_buf;
					//text[52].setFillColor(Color(r, g, b));
					//buf.str(std::string());
					//fps = 0;
					//clock4FPS.restart();
				//}
			}
			else {
				p.pause();
				//fps = 0;
			}
		}
		if (state == 0 || state == 3) {
			window->draw(*sCMC);
			window->draw(*sFF);
			window->draw(text[4]);//2019
			window->draw(text[5]);//Флуктуации плотности в одномерном газе
			window->draw(text[45]);//Факультет ВМК
			window->draw(text[46]);//Физический факультет
			window->draw(text[38]);//МГУ им.М.В.Ломоносова
			window->draw(text[39]);//Компьютерные демонстрации по курсу лекций
			window->draw(text[40]);//Статистическая физика
			if (state == 3) {
				window->draw(*sPhoto1);
				window->draw(*sPhoto2);
				window->draw(text[41]);
				window->draw(text[42]);
				window->draw(text[43]);
				window->draw(text[44]);
			}
		}

		if (state == 1) {
			window->draw(tube);
			for (int i = 4; i < 27; ++i) {
				if (i != 7)
					window->draw(button[i]);
			}
			for (int i = 0; i < 5; ++i)
				window->draw(Sliders[i]);//ползунки 
			for (int i = 6; i < 9; ++i)
				window->draw(text[i]);
			for (int i = 12; i < 29; ++i)
				if (i != 13 && i != 14 && i != 18)
					window->draw(text[i]);
			for (int i = 34; i < 38; i++)
				if (i != 35)
					window->draw(text[i]);
			window->draw(text[48]);
			window->draw(text[55]);
			window->draw(text[56]);
			window->draw(text[57]);
			window->draw(text[58]);
			if (button_pressed[8]) {//U = -C / R
				for (int i = 0; i < 3900; ++i) {
					xGR = 0.1 * i + 20;
					graph_GR[i].position = Vector2f((1350 + xGR) * screenWidth / 1920, (344 + (1200.0 + 1000 * (double)coef_C) / xGR) * screenHeight / 1080);
					graph_GR1[i].position = Vector2f((1350 + xGR + 1) * screenWidth / 1920, (344 + (1200.0 + 1000 * (double)coef_C) / xGR) * screenHeight / 1080);
					graph_GR2[i].position = Vector2f((1350 + xGR) * screenWidth / 1920, (344 + (1200.0 + 1000 * (double)coef_C) / xGR + 1) * screenHeight / 1080);
					graph_GR[i].color = graph_GR1[i].color = graph_GR2[i].color = Color::Black;
				}
				window->draw(lineX_graph_GR);
				window->draw(lineY_graph);
				window->draw(graph_GR1);
				window->draw(graph_GR2);
				window->draw(graph_GR);
				window->draw(text[13]);
			}
			if (button_pressed[9]) {//U = C / R
				for (int i = 0; i < 3900; ++i) {
					xK = 0.1 * i + 20;
					graph_K[i].position = Vector2f((1350 + xK) * screenWidth / 1920, (690 - (1200.0 + 1000 * (double)coef_C) / xK) * screenHeight / 1080);
					graph_K1[i].position = Vector2f((1350 + xK + 1) * screenWidth / 1920, (690 - (1200.0 + 1000 * (double)coef_C) / xK) * screenHeight / 1080);
					graph_K2[i].position = Vector2f((1350 + xK) * screenWidth / 1920, (690 - (1200.0 + 1000 * (double)coef_C) / xK + 1) * screenHeight / 1080);
					graph_K[i].color = graph_K1[i].color = graph_K2[i].color = Color::Black;
				}
				window->draw(lineX_graph_K);
				window->draw(lineY_graph);
				window->draw(graph_K);
				window->draw(graph_K1);
				window->draw(graph_K2);
				window->draw(text[14]);
			}

			if (button_pressed[25]) {//U = LJ
				for (int i = 0; i < 3900; ++i) {
					switch (coef_C) {
					case 1:
						xLJ = xLJ = 0.1 * i + 178;
						break;
					case 2:
						xLJ = xLJ = 0.1 * i + 184.5;
						break;
					case 3:
						xLJ = xLJ = 0.1 * i + 187.9;
						break;
					case 4:
						xLJ = xLJ = 0.1 * i + 190;
						break;
					case 5:
						xLJ = xLJ = 0.1 * i + 191.5;
						break;
					}
					graph_LJ[i].position = Vector2f((1200 + xLJ) * screenWidth / 1920, (530 - (delta * coef_C * 0.2 * (1.0 / pow(xLJ * sigma, 12) - 1.0 / pow(xLJ * sigma, 6)))) * screenHeight / 1080);
					graph_LJ1[i].position = Vector2f((1200 + xLJ + 1) * screenWidth / 1920, (530 - (delta * coef_C * 0.2 * (1.0 / pow(xLJ * sigma, 12) - 1.0 / pow(xLJ * sigma, 6)))) * screenHeight / 1080);
					graph_LJ2[i].position = Vector2f((1200 + xLJ) * screenWidth / 1920, (530 - (delta * coef_C * 0.2 * (1.0 / pow(xLJ * sigma, 12) - 1.0 / pow(xLJ * sigma, 6)) + 1)) * screenHeight / 1080);
					graph_LJ[i].color = graph_LJ1[i].color = graph_LJ2[i].color = Color::Black;
				}
				window->draw(lineX_graph_LJ);
				window->draw(lineY_graph);
				window->draw(graph_LJ);
				window->draw(graph_LJ1);
				window->draw(graph_LJ2);
				window->draw(text[49]);
				window->draw(text[53]);
				window->draw(text[54]);
			}
			if (button_pressed[26]) { //U = 0
				window->draw(lineX_graph_K);
				window->draw(lineY_graph);
				window->draw(text[50]);
			}

			if (show) {
				if (!stop && clock4Stat.getElapsedTime().asMilliseconds() > 100) {
					buf << p.sampleMean();
					text[32].setString(buf.str());
					buf.str(std::string());
					buf << p.sampleSTD();
					text[33].setString(buf.str());
					buf.str(std::string());
					buf << round(100 * p.sampleSTD() / p.sampleMean()) / 100;
					text[60].setString(buf.str());
					buf.str(std::string());
					maxAmount = p.get_maxAmount_now();
					buf << maxAmount;
					ytextBin.setString(buf.str());
					buf.str(std::string());
					clock4Stat.restart();
				}
				for (int i = 30; i < 34; ++i)
					window->draw(text[i]);
				window->draw(text[59]);// sigma/E
				//if (text[60].getString() != "inf")
				window->draw(text[60]);
				window->draw(text[61]);// sigma
				window->draw(text[47]);
				window->draw(lineX);//
				window->draw(lineY);//   система координат
				window->draw(textX1);//
				window->draw(textX2);//
				window->draw(textX3);//
				window->draw(textX4);//
				window->draw(textY1);//
				window->draw(textY2);//
				window->draw(textY3);//
				for (int i = 0; i < p.get_num(); ++i)
					window->draw(p.get_p(i));
				for (int i = 0; i < p.get_nbins(); ++i) {
					window->draw(p.get_bin(i));//гистограмма
					window->draw(xticks[i]);
					window->draw(yticks[i]);
				}
				for (int i = 0; i < (p.get_nbins() + 1) / 4; ++i) {
					window->draw(xtextBins[i]);
				}
				window->draw(ytextBin);
			}

			if (sliders_pressed[0]) {
				if (pixelPos.x - dX1 >= 460 * screenWidth / 1920 && pixelPos.x - dX1 <= (770 - 64) * screenWidth / 1920)
					Sliders[0].setPosition(pixelPos.x - dX1, 360 * screenHeight / 1080);
				if (Sliders[0].getPosition().x < 460 * screenWidth / 1920)
					width = 1;
				else if (Sliders[0].getPosition().x > (770 - 64)* screenWidth / 1920)
					width = max_width / num;
				else
					width = (int)(Sliders[0].getPosition().x * (max_width / num - 1) / (246 * screenWidth / 1920) + (738 - 460 * max_width / num) / 246);
				init = true;
			}
			if (sliders_pressed[1]) {
				if (pixelPos.x - dX2 >= 460 * screenWidth / 1920 && pixelPos.x - dX2 <= (770 - 64) * screenWidth / 1920)
					Sliders[1].setPosition(pixelPos.x - dX2, 290 * screenHeight / 1080);
				if (Sliders[1].getPosition().x < 460 * screenWidth / 1920)
					num = 2;
				else if (Sliders[1].getPosition().x > (770 - 64)* screenWidth / 1920)
					num = max_num;
				else
					num = (int)(Sliders[1].getPosition().x * (max_num - 1) / (246 * screenWidth / 1920) + 1 + (738 - 460 * max_num) / 246);
				width = (int)(Sliders[0].getPosition().x * (max_width / num - 1) / (246 * screenWidth / 1920) + (738 - 460 * max_width / num) / 246);
				init = true;
			}
			if (sliders_pressed[2]) {
				if (pixelPos.x - dX3 >= 460 * screenWidth / 1920 && pixelPos.x - dX3 <= (770 - 64) * screenWidth / 1920)
					Sliders[2].setPosition(pixelPos.x - dX3, 430 * screenHeight / 1080);
				if (Sliders[2].getPosition().x < 460 * screenWidth / 1920)
					temperature = 1;
				else if (Sliders[2].getPosition().x > (770 - 64)* screenWidth / 1920)
					temperature = max_temperature;
				else
					temperature = (int)(Sliders[2].getPosition().x * (max_temperature - 1) / (246 * screenWidth / 1920) + (738 - 460 * max_temperature) / 246);
				init = true;
			}
			if (sliders_pressed[3]) {
				if (pixelPos.x - dX4 >= 460 * screenWidth / 1920 && pixelPos.x - dX4 <= (770 - 64) * screenWidth / 1920)
					Sliders[3].setPosition(pixelPos.x - dX4, 500 * screenHeight / 1080);
				if (Sliders[3].getPosition().x < 460 * screenWidth / 1920)
					mass = 1;
				else if (Sliders[3].getPosition().x > (770 - 64)* screenWidth / 1920)
					mass = max_mass;
				else
					mass = (int)((Sliders[3].getPosition().x * (max_mass - 1) / (246 * screenWidth / 1920) + (738 - 460 * max_mass) / 246));
				init = true;
			}
			if (sliders_pressed[4]) {
				if (pixelPos.x - dX5 >= 460 * screenWidth / 1920 && pixelPos.x - dX5 <= (770 - 64) * screenWidth / 1920)
					Sliders[4].setPosition(pixelPos.x - dX5, 570 * screenHeight / 1080);
				if (Sliders[4].getPosition().x < 460 * screenWidth / 1920)
					coef_C = 1;
				else if (Sliders[4].getPosition().x > (770 - 64)* screenWidth / 1920)
					coef_C = max_coef_C;
				else
					coef_C = (int)((Sliders[4].getPosition().x * (max_coef_C - 1) / (246 * screenWidth / 1920) + (738 - 460 * max_coef_C) / 246));
				init = true;
			}

			buf << width;
			text[10].setString(buf.str());
			buf.str(std::string());
			window->draw(text[10]);

			buf << num;
			text[11].setString(buf.str());
			buf.str(std::string());
			window->draw(text[11]);

			buf << temperature;
			text[18].setString(buf.str());
			buf.str(std::string());
			window->draw(text[18]);

			buf << mass;
			text[29].setString(buf.str());
			buf.str(std::string());
			window->draw(text[29]);

			buf << coef_C;
			text[35].setString(buf.str());
			buf.str(std::string());
			window->draw(text[35]);

			buf << (int)max_width / num;
			text[23].setString(buf.str());
			buf.str(std::string());
		}
		else if (state == 2) {
			window->draw(*sTh);
			if (sliders_pressed[5]) {
				if (pixelPos.y - dY >= 70 * screenHeight / 1080 && pixelPos.y - dY <= (1080 - 64) * screenHeight / 1080)
					Sliders[5].setPosition(1854 * screenWidth / 1920, pixelPos.y - dY);
				if (Sliders[5].getPosition().y < 70 * screenHeight / 1080)
					sTh->setPosition(370 * screenWidth / 1920, 70 * screenHeight / 1080);
				else if (Sliders[5].getPosition().y > (1080 - 64)* screenHeight / 1080)
					sTh->setPosition(370 * screenWidth / 1920, -1060);
				else
					sTh->setPosition(370 * screenWidth / 1920, Sliders[5].getPosition().y *
					((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)) +
						70 * screenHeight / 1080 - 70 * screenHeight / 1080 * ((-1060 - 70 * screenHeight / 1080) / (screenHeight - 64 * screenHeight / 1080 - 70 * screenHeight / 1080)));
			}
			window->draw(button[27]);
			window->draw(Sliders[5]);
		}
		for (int i = 0; i < 4; ++i)
			window->draw(button[i]);
		for (int i = 0; i < 4; ++i)
			window->draw(text[i]);
		window->draw(button[7]);//кнопка выхода
		window->draw(text[9]);  //выход
		//window->draw(text[51]); //FPS
		//window->draw(text[52]); //FPS
		window->display();
	}
	delete window;
	delete iCMC;
	delete tCMC;
	delete sCMC;
	delete iFF;
	delete tFF;
	delete sFF;
	delete iTh;
	delete tTh;
	delete sTh;
	delete iPhoto1;
	delete tPhoto1;
	delete sPhoto1;
	delete iPhoto2;
	delete tPhoto2;
	delete sPhoto2;
	delete[] xtextBins;
	delete[] xticks;
	delete[] yticks;
	delete[] button;
	delete[] text;
	delete[] Sliders;
	return 0;
}
