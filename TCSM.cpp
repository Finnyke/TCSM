#include "TCSM.h"

using namespace std;

telComSys::RTSG::RTSG(double prob1) {
	if (prob1 > 1 || prob1 < 0) throw "Error: invalid probability value";
	_prob1 = prob1;
	return;
}

void telComSys::RTSG::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
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


telComSys::AWGNG::AWGNG(double sigma) : _deviation(sigma) {};

void telComSys::AWGNG::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	if (_deviation == 0) return;
	unsigned seed = static_cast<unsigned>(rand());
	default_random_engine dre(seed);
	normal_distribution<double> noise(0, _deviation);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) += noise(dre);
	}
	return;
}


telComSys::MDL::MDL(char type) : _type(type) {};

void telComSys::MDL::carrierInit(double endTime, double sampInterval) {
	_carrier.resize(static_cast<size_t>(endTime / sampInterval));
	for (size_t i = 0; i < _carrier.size(); ++i) {
		_carrier.at(i) = sin(i * sampInterval * 4 * _Pi);
	}
	return;
}

void telComSys::MDL::carriersInitFM(double endTime, double sampInterval) {
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

void telComSys::MDL::AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carrierInit(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) += 1;
		s.at(i) *= 0.5 * _carrier.at(i);
	}
	return;
}

void telComSys::MDL::FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carriersInitFM(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		double tempP = s.at(i);
		double tempM = -s.at(i);
		s.at(i) = tempP * _carrier.at(i) + tempM * _carrier2.at(i);
	}
	return;
}

void telComSys::MDL::PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carrierInit(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) *= _carrier.at(i);
	}
	return;
}

void telComSys::MDL::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
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

telComSys::DMDL::DMDL(char type) : _type(type) {};

void telComSys::DMDL::carrierInit(double endTime, double sampInterval) {
	_carrier.resize(static_cast<size_t>(endTime / sampInterval));
	for (size_t i = 0; i < _carrier.size(); ++i) {
		_carrier.at(i) = sin(i * sampInterval * 4 * _Pi);
	}
	return;
}

void telComSys::DMDL::carriersInitFM(double endTime, double sampInterval) {
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

void telComSys::DMDL::output(double thresholdLevel, double digTimeSlot, double sampInterval, vector<double>& s) {
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

void telComSys::DMDL::AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carrierInit(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) *= _carrier.at(i);
	}
	output(0.25, digTimeSlot, sampInterval, s);
	return;
}

void telComSys::DMDL::FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carriersInitFM(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) *= _carrier.at(i) - _carrier2.at(i);
	}
	output(0., digTimeSlot, sampInterval, s);
	return;
}

void telComSys::DMDL::PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	carrierInit(endTime, sampInterval);
	for (size_t i = 0; i < s.size(); ++i) {
		s.at(i) *= _carrier.at(i);
	}
	output(0., digTimeSlot, sampInterval, s);
	return;
}

void telComSys::DMDL::pM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	output(0., digTimeSlot, sampInterval, s);
	return;
}

void telComSys::DMDL::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
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


telComSys::ERC::ERC(unsigned delay, vector<double> initS) : _delay(delay), _cnt(0) {
	_initS = initS;
}

void telComSys::ERC::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	_cnt = 0;
	size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
	for (size_t i = 0; i < s.size() - _delay * length; ++i) {
		if (s.at(i + _delay * length) != _initS.at(i)) _cnt++;
	}
	cout << "Number of errors: " << _cnt / length << endl;
	return;
}

telComSys::MPCH::MPCH(unsigned num, vector<double> coeffs) : _num(num) {
	if (_num > coeffs.size()) throw "Error: insufficient number of coefficients";
	if (_num < coeffs.size()) throw "Error: excessive number of coefficients";
	_gammas = coeffs;
};

void telComSys::MPCH::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
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


telComSys::CRTR::CRTR(char type, unsigned num, vector<double> coeffs) : _type(type), _num(num) {
	_gammas = coeffs;
};

void telComSys::CRTR::recCRTR(size_t length, vector<double>& s, double val, unsigned depth) {
	if (depth * length >= s.size()) return;
	double temp = val;
	val += s.at(depth * length);
	val *= -(_gammas.at(1) / _gammas.at(0));
	for (size_t i = 0; i < length; ++i) {
		s.at(depth * length + i) += temp;
		s.at(depth * length + i) /= _gammas.at(0);
	}
	depth++;
	recCRTR(length, s, val, depth);
}

void telComSys::CRTR::nrCRTR(size_t length, vector<double>& s) {
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

void telComSys::CRTR::runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) {
	size_t length = static_cast<size_t>(digTimeSlot / sampInterval);
	_type == 'R' ? recCRTR(length, s) : nrCRTR(length, s);
	return;
}

bool telComSys::cmpd(double lhs, double rhs) {
	return (abs(lhs - rhs) < 0.000001);
}

bool telComSys::checkForMltpl(double x, double y) {
	if (x < y) return false;
	x -= y;
	if (cmpd(x, y)) return true;
	else checkForMltpl(x, y);
}
	

telComSys::telComSys(double endTime, double digTimeSlot, double sampInterval) : _endTime(endTime), _digTimeSlot(digTimeSlot), _sampInterval(sampInterval) {
	if (endTime <= 0 || digTimeSlot <= 0 || sampInterval <= 0) throw "Error: all parameters must be positive";
	if (!checkForMltpl(endTime, digTimeSlot)) throw "Error: modeling end time must be a multiple of digit time slot";
	if (!checkForMltpl(digTimeSlot, sampInterval)) throw "Error: digit time slot must be a multiple of sample interval";
	srand(static_cast<unsigned>(time(nullptr)));
	_s.resize(static_cast<size_t>(endTime / sampInterval));
	return;
}

void telComSys::initAWNG() {
	double sigma = read_double("Enter noise deviation: ", 0., 500.);
	AWGNG* ptr = new AWGNG(sigma);
	_queue.push_back(pair<element*, elTypes>(ptr, elTypes::AWGNG));
	return;
}

void telComSys::initCRTR() {
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

void telComSys::initDMDL() {
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

void telComSys::initERC() {
	_queue.push_back(pair<element*, elTypes>(nullptr, elTypes::ERC));
	return;
}

void telComSys::initMDL() {
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

void telComSys::initMPCH() {
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

void telComSys::initRTSG() {
	double p = read_double("Enter probability of '1': ", 0., 1.);
	RTSG* ptr = new RTSG(p);
	_queue.push_back(pair<element*, elTypes>(ptr, elTypes::RTSG));
	return;
}

void telComSys::appendToQueue(elTypes type) {
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

void telComSys::run() {
	for (size_t i = 0; i < _queue.size(); ++i) {
		if (_queue.at(i).second == elTypes::ERC) {
			unsigned u = static_cast<unsigned>(read_int("Enter system's overall delay: ", 0, static_cast<int>(_endTime / _digTimeSlot)));
			_queue.at(i).first = new ERC(u, _initS);
		}
		_queue.at(i).first->runEl(_endTime, _digTimeSlot, _sampInterval, _s);
		if (_queue.at(i).second == elTypes::RTSG) _initS = _s;
		if (_queue.at(i).second == elTypes::RTSG) printSignal();
	}
}

void telComSys::printSignal() {
	for (size_t i = 0; i < _s.size(); ++i) {
		cout << "s[" << i << "] = " << _s.at(i) << '\t';
		if (((i + 1) % 5) == 0) cout << endl;
	}
}
