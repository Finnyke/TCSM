#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <utility>
#include "io.h"

using namespace std;

enum class elTypes { //Types of available elements
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

	double _digTimeSlot; //Digit time slot

	double _sampInterval; //Smaple interval

	vector<double> _s; //Main signal

	vector<double> _initS; //Initial signal (after RTSG)

	vector<double> _gammas; //Coefficients for multipath channel

	class element {
	public:

		virtual void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s) = 0; //Function which models change to the main signal

	};

	class RTSG : public element { //Random telegraph signal generator
	public:

		double _prob1; //Probability of 1

		RTSG(double prob1);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class AWGNG : public element { //Additive white gaussian noise generator
	public:

		double _deviation; //Standard deviation of gaussian noise

		AWGNG(double sigma);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class MDL : public element { //Modulator
	public:

		char _type; //Modulation type

		vector<double> _carrier; //Carrier signals for different modulation types

		vector<double> _carrier2;

		MDL(char type);

		void carrierInit(double endTime, double sampInterval); //Initialization of the carrier signal for AM and PH

		void carriersInitFM(double endTime, double sampInterval); //Initialization of the carrier signal for FM

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s); //Modulation functions for corresponding modulation types

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class DMDL : public element { //Demodulator
	public:

		char _type; //Modulation type

		vector<double> _carrier; //Carrier signals for different modulation types

		vector<double> _carrier2;

		DMDL(char type);

		void carrierInit(double endTime, double sampInterval); //Initialization of the carrier signal for AM and PH

		void carriersInitFM(double endTime, double sampInterval); //Initialization of the carrier signal for FM

		void output(double thresholdLevel, double digTimeSlot, double sampInterval, vector<double>& s); //Integrator and decision-making device

		void AM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void FM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void PM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);

		void pM(double endTime, double digTimeSlot, double sampInterval, vector<double>& s); //Low-frequency phase modulation

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class ERC : public element { //Error counter
	public:

		unsigned _cnt; //Counter itself

		unsigned _delay; //Total delay in the system, in digit time slots

		vector<double> _initS; //Initial signal (after RTSG)

		ERC(unsigned delay, vector<double> initS);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class MPCH : public element { //Mutipath channel
	public:

		unsigned _num; //Number of paths

		vector<vector<double>> _sig; //Signals on each path

		vector<double> _gammas; //Coefficients

		MPCH(unsigned num, vector<double> coeffs);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	class CRTR : public element { //Corrector
	public:

		char _type; //Corrector type

		unsigned _num; //Number of elements in non-recursive corrector

		vector<vector<double>> _sig; //Signals on each branch of non-recursive corrector

		vector<double> _gammas; //Coefficients

		CRTR(char type, unsigned num, vector<double> coeffs);

		void recCRTR(size_t length, vector<double>& s, double val = 0, unsigned depth = 0);

		void nrCRTR(size_t length, vector<double>& s);

		void runEl(double endTime, double digTimeSlot, double sampInterval, vector<double>& s);
	};

	vector<pair<element*, elTypes>> _queue; //Queue of the elements in the system

public:

	bool cmpd(double lhs, double rhs); //Floating point values comparison

	bool checkForMltpl(double x, double y); //Checks if x is a multiple of y

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
