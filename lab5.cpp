#include "TCSM.h"

int main() {
	try {
		telComSys t(30, 1, 0.1);
		t.appendToQueue(elTypes::RTSG);
		t.appendToQueue(elTypes::MPCH);
		t.appendToQueue(elTypes::AWGNG);
		t.appendToQueue(elTypes::CRTR);
		t.run();
		t.printSignal();
	}
	catch (const char* str) {
		cout << "Runtime error:" << endl;
		cout << str << endl;
	}

	return 0;
}
