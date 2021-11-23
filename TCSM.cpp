#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <utility>

using namespace std;

enum class elTypes {
	RTSG,
	AWNG,
	MDL,
	DMDL,
	ERC,
	MPCH,
};

class telComSys {
private:

	double _endTime; //Modelling end time

	double _digTimeSlot; //Digit time slot / Clock interval

	double _sampInterval; //Sample interval

	vector<double> _s; //Main signal

	vector<double> _initS; //A copy of initial signal

	vector<bool> _oscilFlags; //Flags for oscilloscope usage

	class element {
	public:

		virtual void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) = 0;

	};

	class RTSG : public element { //Random telegraph signal generator
	public:

		double _prob1; //Probability of 1

		RTSG(double prob1) {
			if (prob1 > 1 || prob1 < 0) throw "Error: invalid probability value";
			_prob1 = prob1;
			return;
		}

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			srand(time(0));
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			if (_prob1 < 1 && _prob1 > 0) {
				for (size_t i = 0; i < s.size(); i += length) {
					double r = round(((double)rand() / (RAND_MAX)) - 0.5 + _prob1);
					if (!r) r = -1;
					cout << i << '\t' << r << endl;
					for (size_t j = 0; j < length; ++j) {
						s.at(i + j) = r;
					}
				}
			}
			else {
				for (size_t i = 0; i < s.size(); ++i) {
					_prob1 ? s.at(i) = 1 : s.at(i) = -1;
				}
			}
		}
	};

	class AWGNG : public element { //Additive white Gaussian noise generator
	public:

		double _deviation;

		AWGNG(double sigma) : _deviation(sigma) {};

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			srand(time(0));
			unsigned seed = static_cast<unsigned>(rand());
			default_random_engine dre(seed);
			normal_distribution<double> noise(0, _deviation);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) += noise(dre);
			}
		}
	};

	class MDL : public element { //Modulator
	public:

		char _type;

		vector<double> _carrier;

		vector<double> _carrier2;

		MDL(char type) : _type(type) {};

		void carrierInit(double endTime, double sampInterval) {
			_carrier.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 4 * _Pi);
			}
		}

		void carriersInitFM(double endTime, double sampInterval) {
			_carrier.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 5 * _Pi);
			}
			_carrier2.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 3 * _Pi);
			}
		}

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) += 1;
				s.at(i) *= 0.5 * _carrier.at(i);
			}
		}

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carriersInitFM(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				double tempP = s.at(i);
				double tempM = -s.at(i);
				s.at(i) = tempP * _carrier.at(i) + tempM * _carrier2.at(i);
			}
		}

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) *= _carrier.at(i);
			}
		}

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			switch (_type) {
			case'A':
				AM(endTime, digTimeSlot, sampInterval, s);
				break;
			case'F':
				FM(endTime, digTimeSlot, sampInterval, s);
				break;
			case'P':
				PM(endTime, digTimeSlot, sampInterval, s);
				break;
			default:
				throw "Error: invalid modulation type";
				break;
			}
		}
	};

	class DMDL : public element { //Demodulator
	public:

		char _type;

		vector<double> _carrier;

		vector<double> _carrier2;

		DMDL(char type) : _type(type) {};

		void carrierInit(double endTime, double sampInterval) {
			_carrier.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 4 * _Pi);
			}
		}

		void carriersInitFM(double endTime, double sampInterval) {
			_carrier.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 5 * _Pi);
			}
			_carrier2.resize(static_cast<size_t>(endTime / sampInterval));
			for (size_t i = 0; i < _carrier.size(); ++i) {
				_carrier.at(i) = sin(i * sampInterval * 3 * _Pi);
			}
		}

		void output(double thresholdLevel, double digTimeSlot, double sampInterval, vector<double>& s) {
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			for (size_t i = 0; i < s.size(); i += length) {
				double sum = -thresholdLevel; //Midpoint Riemann sum is used for integral approximation
				for (size_t j = 0; j < length; ++j) {
					if (i + j + 1 <= s.size()) sum += 0.5 * (s.at(i + j) + s.at(i + j + 1)) * sampInterval;
					else sum += s.at(i + j) * sampInterval;
				}
				if (sum += 0.5) {
					for (size_t j = 0; j < length; ++j) {
						s.at(i + j) = 1;
					}
				}
				else {
					for (size_t j = 0; j < length; ++j) {
						s.at(i + j) = -1;
					}
				}
			}
			return;
		}

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) *= _carrier.at(i);
			}
			output(0.25, digTimeSlot, sampInterval, s);
			return;
		}

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carriersInitFM(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) *= _carrier.at(i) - _carrier2.at(i);
			}
			output(0., digTimeSlot, sampInterval, s);
			return;
		}

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) *= _carrier.at(i);
			}
			output(0., digTimeSlot, sampInterval, s);
			return;
		}

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			switch (_type) {
			case'A':
				AM(endTime, digTimeSlot, sampInterval, s);
				break;
			case'F':
				FM(endTime, digTimeSlot, sampInterval, s);
				break;
			case'P':
				PM(endTime, digTimeSlot, sampInterval, s);
				break;
			default:
				throw "Error: invalid modulation type";
				break;
			}
		}
	};

	class ERC : public element { //Error counter

	};

	class MPCH : public element { //Multipath channel

	};

	vector<pair<element*, elTypes>> _queue; //Queue of elements in the system

public:

	telComSys(double endTime, double digTimeSlot, double sampInterval) : _endTime(endTime), _digTimeSlot(digTimeSlot), _sampInterval(sampInterval) {
		//must be bigger than zero
		//must be a multiple
		if (digTimeSlot < sampInterval) throw "Error: digit time slot must be bigger than sample interval";
		//if (digTimeSlot / sampInterval) throw "Error: digit time slot must be bigger than sample interval";
		_s.resize(static_cast<size_t>(endTime / sampInterval));
		return;
	};

	void appendToQueue(const char* str) { //через пользовательский ввод и enum
		int elType = -1;
		vector<double> param;
		while (*str != '\0') {
			str++;
		}
		RTSG gen(0.5);
		_queue.push_back(pair<element*, elTypes>(&gen, elTypes::RTSG));
	}

	void run() {
		for (auto i : _queue) {
			i.first->runEl(_endTime, _digTimeSlot, _sampInterval, _s);
			if (i.second == elTypes::RTSG) _initS = _s; //сделать как статический элемент класса счетчика ошибок?
		}
	}

	void printSignal() {
		for (size_t i = 0; i < _s.size(); ++i) {
			cout << "s[" << i << "] = " << _s.at(i) << '\t';
			if (((i + 1) % 5) == 0) cout << endl;
		}
	}
};

int main() {
	telComSys t(30, 1, 0.1);
	t.appendToQueue("11");
	t.printSignal();
	t.run();
	t.printSignal();

	return 0;
}