#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <utility>
#include "io.h"

using namespace std;



enum class elTypes {
	RTSG,
	AWGNG,
	MDL,
	DMDL,
	ERC,
	MPCH,
	CRTR,
};

class telComSys {
private:

	double _endTime; //Modeling end time

	double _digTimeSlot; //Digit time slot / Clock interval

	double _sampInterval; //Sample interval

	vector<double> _s; //Main signal

	vector<double> _initS; //A copy of initial signal

	vector<bool> _oscilFlags; //Flags for oscilloscope usage

	vector<double> _gammas; //Weight coefficients for multipath channel

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
			return;
		}
	};

	class AWGNG : public element { //Additive white Gaussian noise generator
	public:

		double _deviation;

		AWGNG(double sigma) : _deviation(sigma) {};

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			if (_deviation == 0) return;
			srand(time(0));
			unsigned seed = static_cast<unsigned>(rand());
			default_random_engine dre(seed);
			normal_distribution<double> noise(0, _deviation);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) += noise(dre);
			}
			return;
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
			return;
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
			return;
		}

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) += 1;
				s.at(i) *= 0.5 * _carrier.at(i);
			}
			return;
		}

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carriersInitFM(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				double tempP = s.at(i);
				double tempM = -s.at(i);
				s.at(i) = tempP * _carrier.at(i) + tempM * _carrier2.at(i);
			}
			return;
		}

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			carrierInit(endTime, sampInterval);
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) *= _carrier.at(i);
			}
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
			return;
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
			return;
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
			return;
		}

		void output(double thresholdLevel, double digTimeSlot, double sampInterval, vector<double>& s) {
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			vector<double> temp = s;
			for (size_t i = 0; i < s.size() - length; i += length) {
				double sum = -thresholdLevel; //Midpoint Riemann sum is used for integral approximation
				for (size_t j = 0; j < length; ++j) {
					sum += (s.at(i + j) + s.at(i + j + 1)) * sampInterval;
				}
				if (sum >= 0.5) {
					for (size_t j = 0; j < length; ++j) {
						temp.at(i + j + length) = 1;
					}
				}
				else {
					for (size_t j = 0; j < length; ++j) {
						temp.at(i + j + length) = -1;
					}
				}
			}
			s = temp;
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

		void pM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
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
			case'p':
				pM(endTime, digTimeSlot, sampInterval, s);
				break;
			default:
				throw "Error: invalid modulation type";
				break;
			}
			return;
		}
	};

	class ERC : public element { //Error counter
	public:

		unsigned _cnt;

		unsigned _delay; //Delay in digit time slots

		vector<double> _initS;

		ERC(unsigned delay, vector<double> initS): _delay(delay), _cnt(0) {
			_initS = initS;
		}

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			_cnt = 0;
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			for (size_t i = 0; i < s.size() - _delay * length; ++i) {
				if (s.at(i + _delay * length) != _initS.at(i)) _cnt++;
			}
			cout << "Number of errors: " << _cnt / length << endl;
			return;
		}
	};

	class MPCH : public element { //Multipath channel
	public:

		unsigned _num; //Number of paths

		vector<vector<double>> _sig;

		vector<double> _gammas;

		MPCH(unsigned num, vector<double> coeffs): _num(num) {
			if (_num > coeffs.size()) throw "Error: insufficient number of coefficients";
			if (_num < coeffs.size()) throw "Error: excessive number of coefficients";
			_gammas = coeffs;
		};
		
		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			_sig.resize(_num);
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			for (size_t i = 0; i < _num; ++i) {
				_sig.at(i).resize(s.size() + i * length);
				for (auto j : _sig.at(i)) {
					j = 0;
				}
			}
			for (size_t i = 0; i < _num; ++i) {
				for (size_t j = 0; j < s.size(); ++j) {
					_sig.at(i).at(j + i * length) = s.at(j);
				}
			}
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) = 0;
				for (size_t j = 0; j < _sig.size(); ++j) {
					s.at(i) += _sig.at(j).at(i) * _gammas.at(j);
				}
			}
			return;
		}
	};

	class CRTR : public element { //Corrector
	public:
		char _type;

		unsigned _num; //Number of elements for non-recursive corrector

		vector<vector<double>> _sig;

		vector<double> _gammas;

		CRTR(char type, unsigned num, vector<double> coeffs) : _type(type), _num(num) {
			_gammas = coeffs;
		};

		void recCRTR(size_t length, vector<double>& s, double val = 0, unsigned depth = 0) {
			if (depth * length >= s.size()) return;
			double temp = val;
			val += s.at(depth * length);
			val *= - (_gammas.at(1) / _gammas.at(0));
			for (size_t i = 0; i < length; ++i) {
				s.at(depth * length + i) += temp;
				s.at(depth * length + i) /= _gammas.at(0);
			}
			depth++;
			recCRTR(length, s, val, depth);
		}

		void nrCRTR(size_t length, vector<double> &s) {
			_sig.resize(_num + 1);
			for (size_t i = 0; i <= _num; ++i) {
				_sig.at(i).resize(s.size() + i * length);
				for (auto j : _sig.at(i)) {
					j = 0;
				}
			}
			double k = _gammas.at(0) / _gammas.at(1);
			for (size_t i = 0; i <= _num; ++i) {
				for (size_t j = 0; j < s.size(); ++j) {
					_sig.at(i).at(j + i * length) += s.at(j) * pow((-k), _num - i);
				}
			}
			for (size_t i = 0; i < s.size(); ++i) {
				s.at(i) = 0;
				for (size_t j = 0; j < _sig.size(); ++j) {
					s.at(i) += _sig.at(j).at(i) / _gammas.at(1);
				}
			}
			return;
		}

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			_type == 'R' ? recCRTR(length, s) : nrCRTR(length, s);
			return;
		}
	};

	vector<pair<element*, elTypes>> _queue; //Queue of elements in the system

public:

	telComSys(double endTime, double digTimeSlot, double sampInterval) : _endTime(endTime), _digTimeSlot(digTimeSlot), _sampInterval(sampInterval) {
		//must be bigger than zero
		//must be a multiple
		if (digTimeSlot < sampInterval) throw "Error: digit time slot must be bigger than sample interval";
		if (endTime < digTimeSlot) throw "Error: modeling end time must be bigger than digit time slot";
		//if (digTimeSlot / sampInterval) throw "Error: digit time slot must be bigger than sample interval";
		_s.resize(static_cast<size_t>(endTime / sampInterval));
		return;
	}

	void initAWNG() {
		double sigma = read_double("Enter noise deviation: ", 0., 500.);
		AWGNG* ptr = new AWGNG(sigma);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::AWGNG));
		return;
	}

	void initCRTR() {
		int t = read_int("Choose corrector type:\n1. Recursive\n2. Non-recursive\n", 1, 2);
		char c;
		unsigned n = 0;
		if (t == 2) {
			n = static_cast<unsigned>(read_int("Enter the number of elements (for non-recursive corrector): ", 1, 40));
			c = 'N';
		}
		else c = 'R';
		vector<double> coeffs;
		for (auto i = 2; i > 0; --i) {
			cout << "Enter coefficient for path " << 2 - i + 1 << ": ";
			coeffs.push_back(read_double("", -15, +15));
		}
		CRTR* ptr = new CRTR(c, n, coeffs);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::CRTR));
		return;
	}

	void initDMDL() {
		int t = read_int("Choose modulation type (for demodulator):\n1. Amplitude\n2. Phase\n3. Frequency\n4. Phase (low frequency)\n", 1, 4);
		char c = 0;
		switch (t) {
		case 1:
			c = 'A';
			break;
		case 2:
			c = 'P';
			break;
		case 3:
			c = 'F';
			break;
		case 4:
			c = 'p';
			break;
		default:
			throw "Error: invalid modulation type";
		}
		DMDL* ptr = new DMDL(c);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::DMDL));
		return;
	}

	void initERC() {
		_queue.push_back(pair<element*, elTypes>(nullptr, elTypes::ERC));
		return;
	}

	void initMDL() {
		int t = read_int("Choose modulation type:\n1. Amplitude\n2. Phase\n3. Frequency\n", 1, 3);
		char c = 0;
		switch (t) {
		case 1:
			c = 'A';
			break;
		case 2:
			c = 'P';
			break;
		case 3:
			c = 'F';
			break;
		default:
			throw "Error: invalid modulation type";
		}
		MDL* ptr = new MDL(c);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::MDL));
		return;
	}

	void initMPCH() {
		unsigned n = static_cast<unsigned>(read_int("Enter the number of paths: ", 2, 3));
		vector<double> coeffs;
		for (auto i = n; i > 0; --i) {
			cout << "Enter coefficient for path " << n - i + 1 << ": ";
			coeffs.push_back(read_double("", -15, +15));
		}
		MPCH* ptr = new MPCH(n, coeffs);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::MPCH));
		return;
	}

	void initRTSG() {
		double p = read_double("Enter probability of '1': ", 0., 1.);
		RTSG* ptr = new RTSG(p);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::RTSG));
		return;
	}
	
	void appendToQueue(elTypes type) {
		switch (type) {
		case elTypes::AWGNG:
			initAWNG();
			break;
		case elTypes::CRTR:
			initCRTR();
			break;
		case elTypes::DMDL:
			initDMDL();
			break;
		case elTypes::ERC:
			initERC();
			break;
		case elTypes::MDL:
			initMDL();
			break;
		case elTypes::MPCH:
			initMPCH();
			break;
		case elTypes::RTSG:
			initRTSG();
			break;
		default:
			throw "Error: invalid element type";
			break;
		}
		return;
	}

	void run() {
		for (size_t i = 0; i < _queue.size(); ++i) {
			if (_queue.at(i).second == elTypes::ERC) {
				unsigned u = static_cast<unsigned>(read_int("Enter system's overall delay: ", 0, static_cast<int>(_endTime / _digTimeSlot)));
				_queue.at(i).first = new ERC(u, _initS);
			}
			_queue.at(i).first->runEl(_endTime, _digTimeSlot, _sampInterval, _s);
			if (_queue.at(i).second == elTypes::RTSG) _initS = _s; //сделать как статический элемент класса счетчика ошибок?
			if (_queue.at(i).second == elTypes::RTSG) printSignal();
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
	try {
		telComSys t(30, 1, 0.1);
		t.appendToQueue(elTypes::RTSG);
		t.appendToQueue(elTypes::MPCH);
		t.appendToQueue(elTypes::AWGNG);
		t.appendToQueue(elTypes::CRTR);
		//t.appendToQueue(elTypes::MDL);
		t.appendToQueue(elTypes::DMDL);
		t.appendToQueue(elTypes::ERC);
		t.run();
		t.printSignal();
	}
	catch (const char* str) {
		cout << "Runtime error:" << endl;
		cout << str << endl;
	}

	return 0;
}
