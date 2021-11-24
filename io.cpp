#include <io.h>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>

int read_int(const char* prompt, int min, int max) {
	while (true) {
		printf("%s", prompt);
		char str[16];
		fgets(str, sizeof(str), stdin);
		if (strchr(str, '\n') == NULL) {
			bool f = true;
			int c;
			while ((c = fgetc(stdin)) != '\n' && c != EOF) {
				if (!isspace(c)) {
					f = false;
				}
			}
			if (f == false) {
				printf("The entered string is too large.\n");
				continue;
			}
		}
		char* endptr = NULL;
		errno = 0;
		long n = strtol(str, &endptr, 10);
		if (errno == ERANGE) {
			printf("The entered number is too large\n");
			continue;
		}
		if (endptr == str) {
			printf("The entered string does not contain a number.\n");
			continue;
		}
		bool correct = true;
		for (char* p = endptr; *p != 0; ++p) {
			if (!isspace(*p)) {
				correct = false;
				break;
			}
		}
		if (correct == false) {
			printf("The entered string contains invalid character\n");
			continue;
		}
		if (n < min || n > max) {
			printf("Valid range is from %d to %d\n", min, max);
			continue;
		}
		return (int)n;
	}
}

double read_double(const char* prompt, double min, double max) {
	while (true) {
		printf("%s", prompt);
		char str[16];
		fgets(str, sizeof(str), stdin);
		if (strchr(str, '\n') == NULL) {
			bool f = true;
			int c;
			while ((c = fgetc(stdin)) != '\n' && c != EOF) {
				if (!isspace(c)) {
					f = false;
				}
			}
			if (f == false) {
				printf("The entered string is too large\n");
				continue;
			}
		}
		char* endptr = NULL;
		errno = 0;
		double n = strtod(str, &endptr);
		if (errno == ERANGE) {
			printf("The entered number is too large\n");
			continue;
		}
		if (endptr == str) {
			printf("The entered string does not contain a number.\n");
			continue;
		}
		bool correct = true;
		for (char* p = endptr; *p != 0; ++p) {
			if (!isspace(*p) && *p != '.') {
				correct = false;
				break;
			}
		}
		if (correct == false) {
			printf("The entered string contains invalid character\n");
			continue;
		}
		if (n < min || n > max) {
			printf("Valid range is from %.2lf to %.2lf\n", min, max);
			continue;
		}
		return (double)n;
	}
}