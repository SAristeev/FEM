#pragma once
#ifndef BASE64_H
#define BASE64_H

#include <math.h>
#include <iostream>
#include <string>
#include <exception>


namespace base64
{
	constexpr size_t block_size = 1024;

	union buf24
	{
		unsigned int i;
		char c[4];
	};
}

char const c_base64_chars[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
unsigned int const c_base64_ordrs[128] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 0, 0, 0, 63, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 0, 0, 0, 0, 0, 0,
	0,  0 , 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 0, 0, 0, 0,
	0, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 0, 0, 0, 0, 0
};

size_t b64esize(size_t data_size);
size_t b64dsize_host(size_t data_size, char const* host_data);
void b64encode_host(char const* in, size_t in_size, char* out);

void b64decode_host(char const* in, size_t in_size, char* out);

#endif // BASE64_H