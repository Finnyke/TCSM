#pragma once
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

	double _endTime;

	double _digTimeSlot;

	double _sampInterval;

	vector<double> _s;

	vector<double> _initS;

	vector<bool> _oscilFlags;

	vector<double> _gammas;

	class element {
	public:

		virtual void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) = 0;

	};

	class RTSG : public element {
	public:

		double _prob1;

		RTSG(double prob1);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class AWGNG : public element {
	public:

		double _deviation;

		AWGNG(double sigma);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class MDL : public element {
	public:

		char _type;

		vector<double> _carrier;

		vector<double> _carrier2;

		MDL(char type);

		void carrierInit(double endTime, double sampInterval);

		void carriersInitFM(double endTime, double sampInterval);

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class DMDL : public element {
	public:

		char _type;

		vector<double> _carrier;

		vector<double> _carrier2;

		DMDL(char type);

		void carrierInit(double endTime, double sampInterval);

		void carriersInitFM(double endTime, double sampInterval);

		void output(double thresholdLevel, double digTimeSlot, double sampInterval, vector<double>& s);

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void pM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class ERC : public element {
	public:

		unsigned _cnt;

		unsigned _delay;

		vector<double> _initS;

		ERC(unsigned delay, vector<double> initS);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class MPCH : public element {
	public:

		unsigned _num;

		vector<vector<double>> _sig;

		vector<double> _gammas;

		MPCH(unsigned num, vector<double> coeffs);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class CRTR : public element {
	public:

		char _type;

		unsigned _num;

		vector<vector<double>> _sig;

		vector<double> _gammas;

		CRTR(char type, unsigned num, vector<double> coeffs);

		void recCRTR(size_t length, vector<double>& s, double val = 0, unsigned depth = 0);

		void nrCRTR(size_t length, vector<double>& s);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	vector<pair<element*, elTypes>> _queue;

public:

	telComSys(double endTime, double digTimeSlot, double sampInterval);

	void initAWNG();

	void initCRTR();

	void initDMDL();

	void initERC();

	void initMDL();

	void initMPCH();

	void initRTSG();

	void appendToQueue(elTypes type);

	void run();

	void printSignal();
};