#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <utility>
#include "io.h"

using namespace std;

enum class elTypes {
	RTSG,
	AWNG,
	MDL,
	DMDL,
	ERC,
	MPCH,
	CRTR,
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
			return;
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
			for (size_t i = 0; i < s.size() - _delay * length - 1; ++i) {
				if (s.at(i + _delay * length) != _initS.at(i)) _cnt++; //_initS имеет размер 0 почему-то
			}
			cout << "Number of errors: " << _cnt << endl;
			return;
		}
	};

	class MPCH : public element { //Multipath channel
	public:

		unsigned _num; //Number of paths

		vector<vector<double>> _sig;

		MPCH(unsigned num): _num(num) {};
		
		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
			_sig.resize(_num);
			size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
			for (size_t i = 0; i < _num; ++i) {
				_sig.at(i).resize(s.size() + i * length);
				for (size_t j = 0; j < _sig.at(i).size(); ++j) {
					_sig.at(i).at(j) = 0;
				}
				for (size_t j = 0; j < s.size(); ++j) {
					_sig.at(i).at(j + i * length) = s.at(j);
				}
			}
			for (size_t i = 0; i < s.size(); ++i) {
				for (size_t j = 0; j < _sig.size(); ++j) {
					s.at(i) = 0;
					s.at(i) += _sig.at(j).at(i);
				}
			}
			return;
		}
	};

	class CRTR : public element { //Corrector
	public:
		char _type;

		unsigned _num; //Number of elements for non-recursive corrector

		CRTR(char type, unsigned num) : _type(type), _num(num) {};
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

	void initAWNG() {
		return;
	}

	void initCRTR() {
		return;
	}

	void initDMDL() {
		int t = read_int("Choose modulation type (for demodulator):\n1. Amplitude\n2. Phase\n3. Frequency\n", 1, 3);
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
		return;
	}

	void initRTSG() {
		double p = read_double("Enter probability of '1': ", 0., 1.);
		RTSG* ptr = new RTSG(p);
		_queue.push_back(pair<element*, elTypes>(ptr, elTypes::RTSG));
		return;
	}

	void appendToQueue(elTypes type) { //через пользовательский ввод и enum
		switch (type) {
		case elTypes::AWNG:

			break;
		case elTypes::CRTR:

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
				ERC* ptr = new ERC(1, _initS);
				_queue.at(i).first = ptr;
			}
			_queue.at(i).first->runEl(_endTime, _digTimeSlot, _sampInterval, _s);
			if (_queue.at(i).second == elTypes::RTSG) _initS = _s; //сделать как статический элемент класса счетчика ошибок?
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
	t.appendToQueue(elTypes::RTSG);
	t.appendToQueue(elTypes::MDL);
	t.appendToQueue(elTypes::DMDL);
	t.appendToQueue(elTypes::ERC);
	t.printSignal();
	t.run();
	t.printSignal();

	return 0;
}
