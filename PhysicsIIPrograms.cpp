#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>
constexpr static auto k = 8990000000;

namespace dt {
	struct C {
	private:
		static double nCtoC() {
			return .000000001;
		}
	public:
		double value;

		static C nC(double n) {
			return C((n * nCtoC()));
		}

		C(double a) : value(a) {};
		C() : value(0) {};

		friend double operator"" nC(long double a);
	};


	double operator"" nC(long double a) {
		return (a * C::nCtoC());
	}
}

namespace EP {
	using namespace std;
	using namespace dt;

	void whatIsThisNamespace() {
		cout << "This namespace is dedicated to Equipotential mapping. Primarily by creating an array of particles that take a x, y in a 2D map and a charge of either C or nC \n"
			<< "Then we use calcEletricField or calcVoltage for the particles' bounding box with a unit \"l\" \n"
			<< "You can then use writeCSV to output a CSV file of the files with the origin at top left";
	}

	struct Point {
		double x, y;
		double value;
		C* charge = nullptr;
		Point(double x, double y, double v, C* c = nullptr) : x(x), y(y), value(v), charge(c) {}
	};

	struct Particle {
		double x, y;
		C charge;
		Particle(double x, double y, C c) : x(x), y(y), charge(c) {}
		Particle() : x(0), y(0), charge(C(0)) {};
		friend istream& operator >> (istream &in, Particle &p);
		friend ostream& operator << (ostream &out, Particle &p);
	};

	istream& operator >> (istream& in, Particle& p) {
		cout << "Particle's X\n";
		cin >> p.x;
		cout << "Particle's Y\n";
		cin >> p.y;
		cout << "Particle's charge\n";
		string a;
		cin >> a;
		auto i = a.find(" ");
		string b = a.substr(0, i);
		double v = stod(b);
		return in;
	}

	ostream& operator << (ostream& in, Particle& p) {
		cout << "Particle's X: " << p.x << " Particle's Y: " << p.y << " Particle's Charge: " << p.charge.value << "\n";
		return in;
	}

	double distance(double x1, double x2, double y1, double y2) {
		return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	}


	double eField(Particle a, double x, double y) {
		auto d = distance(a.x, x, a.y, y);
		auto c = (k * a.charge.value) / (d * d);
		return c;
	}

	double voltage(Particle a, double x, double y) {
		auto d = distance(a.x, x, a.y, y);
		return (k * a.charge.value) / d;
	}

	template<typename T> T min(T a, T b) {
		if (a > b) {
			return b;
		}
		else {
			return a;
		}
	}

	template<typename T> T max(T a, T b) {
		if (min(a, b) == a) {
			return b;
		}
		else {
			return a;
		}
	}

	struct BoundingBox {
		double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
		BoundingBox(vector<Particle> a) {
			for (Particle x : a) {
				xmin = min(x.x, xmin);
				xmax = max(x.x, xmax);
				ymin = min(x.y, ymin);
				ymax = max(x.y, ymax);
			}
		}

		BoundingBox(double xm, double xma, double ym, double yma) : xmin(xm), xmax(xma), ymin(ym), ymax(yma) {};

		int countX(double l) {
			return round(xmax * 10);
		}

		int countY(double l) {
			return round(ymax * 10);
		}
	};

	//Eletricfield for boundingbox of particles
	vector<Point> calcEletricField(vector<Particle> a, double l, BoundingBox* bb = nullptr, bool DEBUG = false) {
		double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
		vector<Point> r;
		if (bb == nullptr) {
			for (Particle x : a) {
				xmin = min(x.x, xmin);
				xmax = max(x.x, xmax);
				ymin = min(x.y, ymin);
				ymax = max(x.y, ymax);
			}
		}
		else {
			xmin = bb->xmin;
			xmax = bb->xmax;
			ymin = bb->ymin;
			ymax = bb->ymax;
		}
		for (double y = ymin; y <= ymax + l; y += l) {
			for (double x = xmin; x <= xmax + l; x += l) {
				auto v1 = accumulate(a.begin(), a.end(), 0.0, [x, y](double total, Particle current) {return total + eField(current, x, y); });
				if (DEBUG)
					cout << "At " << x << ", " << y << " its eletric field is " << v1 << "\n";
				r.push_back(Point(x, y, v1));
			}
		}
		return r;
	}

	vector<double> calcVoltage(vector<Particle> a, double l, BoundingBox* bb = nullptr, bool DEBUG = false) {
		double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
		vector<double> r;
		if (bb == nullptr) {
			for (Particle x : a) {
				xmin = min(x.x, xmin);
				xmax = max(x.x, xmax);
				ymin = min(x.y, ymin);
				ymax = max(x.y, ymax);
			}
		}
		else {
			xmin = bb->xmin;
			xmax = bb->xmax;
			ymin = bb->ymin;
			ymax = bb->ymax;
		}
		for (double y = ymin; y <= xmax + l; y += l) {
			for (double x = xmin; x <= ymax + l; x += l) {
				auto v1 = accumulate(a.begin(), a.end(), 0.0, [x, y](double total, Particle current) {return total + voltage(current, x, y); });
				if (DEBUG)
					cout << "At " << x << ", " << y << " its voltage is " << v1 << "\n";
				r.push_back(v1);
			}
		}
		return r;
	}

	void writeCSV(vector<Point> a) {
		ofstream f;
		f.open("output.csv");
		for (size_t i = 0; i < a.size(); i++) {
			while (i + 1 != a.size() && a[i].y == a[i + 1].y) {
				f << a[i].value << ",";
				i++;
			}
			if (i + 1 != a.size() && a[i].y != a[i + 1].y) {
				f << "\n";
			}
		}
		f.close();
	}
};

using namespace std;

double ezCoulombsLaw(float a = 0, float b = 0, int c = 0, int d = 0, double distance = 0) {
    if (a == 0) {
    cout << "First base to be mutiplied by 10^X \n";
    cin >> a;
    cout << "First power \n";
    cin >> c;
    cout << "Second base to be mutiplied by 10^X \n";
    cin >> b;
    cout << "First power \n";
    cin >> d;
    cout << "Distance from each other \n";
    cin >> distance;
    }
    auto p = ((a * pow(10, c)) * (b * pow(10, d)));
    auto f = (k * p) / (distance * distance);
    cout << f << endl;
    return f;
}

double ezCoulombsLawWithTwoFieldsAtAPoint() {
    double a, b;
    double d1, d2;
    int e, f;
    cout << "First field charge \n";
    cin >> a;
    cout << "First power \n";
    cin >> e;
    cout << "Distance 1 \n";
    cin >> d1;
    cout << "Second field charge \n";
    cin >> b;
    cout << "Second power \n";
    cin >> f;
    cout << "Distance 2 \n";
    cin >> d2;

    auto c1 = ezCoulombsLaw(a, 1, e, 0, d1);
    auto c2 = ezCoulombsLaw(b, 1, f, 0, d2);
    cout << c1 - c2 << endl; 
    return c1 - c2;
}

void equiPotentialMapping() {
	size_t num = 0;
	cout << "Number of particles: ";
	cin >> num;
	if (num > 25) {
		cout << "Too many particles";
		return;
	}
	else {
		vector<EP::Particle> p;
		for (size_t i = 0; i < num; i++) {
			EP::Particle a;
			cin >> a;
			cout << a;
		}
	}
}

void parallelPlateCapacitor() {
	
}

int main()
{
reset:
    int i;
    cout << "List of available functions: \n";
    cout << "1. Coulumbs Law between two Particles.\n2. Coulombs Law at a point between two fields.\n3. Equipotential Mapping with CVS output\n";
    cin >> i;
    switch (i) {
    case 1:
        ezCoulombsLaw();
        break;
    case 2:
        ezCoulombsLawWithTwoFieldsAtAPoint();
        break;
	case 3:
		equiPotentialMapping();
		break;
    }
    cout << "Use again?\n";
    char yes;
    cin >> yes;
    if (yes == 'y'|| yes=='Y') {
        goto reset; 
    }
}
