#pragma once
#ifndef BASE64_H
#define BASE64_H

#include <iostream>
#include <string>

namespace base64
{
	size_t esize(size_t data_size);
	size_t dsize(size_t data_size, char const* data);

	void decode(char const* in, size_t in_size, char* out);
	void encode(char const* in, size_t in_size, char* out);

	template <std::ranges::sized_range Range>
	void encode(Range in, std::span<char> out);

	namespace detail
	{
		constexpr char         chars[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
		constexpr unsigned int ordrs[128] = {
				0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
				0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 62, 0 , 0 , 0 , 63, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 0 , 0 , 0 , 0 , 0 , 0 ,
				0 ,  0 , 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0 , 0 , 0 , 0 , 0 ,
				0 , 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 0 , 0 , 0 , 0 , 0
			};

		void encode4bytes(size_t idx, char const* in, size_t in_size, char* out);
		void decode3bytes(size_t idx, char const* in, char* out);

		union buf24
		{
			unsigned int i;
			char c[4];
		};
	}
}


inline size_t base64::esize(size_t data_size)
{
	return 4 * (size_t)ceil((double)data_size / 3.0);
}


inline size_t base64::dsize(size_t data_size, char const* data)
{
	size_t size = 3 * data_size / 4;

	if (data[data_size - 1] == '=') size--;
	if (data[data_size - 2] == '=') size--;

	return size;
}


inline void base64::encode(char const* in, size_t in_size, char* out)
{
	for (int i = 0; 3 * i < in_size; i++)
		detail::encode4bytes(i, in, in_size, out);
}


inline void base64::decode(char const* in, size_t in_size, char* out)
{
	for (int i = 0; 4 * i < in_size; i++)
		detail::decode3bytes(i, in, out);
}


inline void base64::detail::encode4bytes(size_t idx, char const* in, size_t in_size, char* out)
{
	unsigned int a = 0, b = 0;
	base64::detail::buf24 buf;

	buf.c[3] = 0;
	buf.c[2] = in[3 * idx + 0];

	if (3 * idx + 1 < in_size) buf.c[1] = in[3 * idx + 1]; else { buf.c[1] = 0; a = 1; }
	if (3 * idx + 2 < in_size) buf.c[0] = in[3 * idx + 2]; else { buf.c[0] = 0; b = 1; }

	out[4 * idx + 0] = chars[(buf.i >> 18) & 0x3f];
	out[4 * idx + 1] = chars[(buf.i >> 12) & 0x3f];

	if (a) out[4 * idx + 2] = '=';
	else   out[4 * idx + 2] = chars[(buf.i >> 6) & 0x3f];

	if (b) out[4 * idx + 3] = '=';
	else   out[4 * idx + 3] = chars[(buf.i >> 0) & 0x3f];
}


inline void base64::detail::decode3bytes(size_t idx, char const* in, char* out)
{
	unsigned int a = 1, b = 1;
	base64::detail::buf24 buf;

	buf.i = 0;

	buf.i |= ordrs[in[4 * idx + 0]] << 18;
	buf.i |= ordrs[in[4 * idx + 1]] << 12;

	if (in[4 * idx + 2] == '=') a = 0; else buf.i |= ordrs[in[4 * idx + 2]] << 6;
	if (in[4 * idx + 3] == '=') b = 0; else buf.i |= ordrs[in[4 * idx + 3]] << 0;

	out[3 * idx + 0] = buf.c[2];
	if (a) out[3 * idx + 1] = buf.c[1];
	if (b) out[3 * idx + 2] = buf.c[0];
}

#endif // BASE64_H