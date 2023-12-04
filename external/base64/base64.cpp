#include "base64.h"

size_t b64esize(size_t data_size){
	return 4 * (size_t)ceil((double)data_size / 3.0);
}


size_t b64dsize_host(size_t data_size, char const* host_data){
	size_t size = 3 * data_size / 4;

	if (host_data[data_size - 1] == '=') size--;
	if (host_data[data_size - 2] == '=') size--;

	return size;
}

void b64encode4bytes(size_t idx, char const* in, size_t in_size, char* out){
	unsigned int a = 0, b = 0;
	base64::buf24 buf;

	buf.c[3] = 0;
	buf.c[2] = in[3 * idx + 0];

	if (3 * idx + 1 < in_size) buf.c[1] = in[3 * idx + 1]; else { buf.c[1] = 0; a = 1; }
	if (3 * idx + 2 < in_size) buf.c[0] = in[3 * idx + 2]; else { buf.c[0] = 0; b = 1; }

	out[4 * idx + 0] = c_base64_chars[(buf.i >> 18) & 0x3f];
	out[4 * idx + 1] = c_base64_chars[(buf.i >> 12) & 0x3f];

	if (a) out[4 * idx + 2] = '=';
	else   out[4 * idx + 2] = c_base64_chars[(buf.i >> 6) & 0x3f];

	if (b) out[4 * idx + 3] = '=';
	else   out[4 * idx + 3] = c_base64_chars[(buf.i >> 0) & 0x3f];
}


void b64encode_host(char const* in, size_t in_size, char* out)
{
	for (int i = 0; 3 * i < in_size; i++)
	{
		b64encode4bytes(i, in, in_size, out);
	}
}


void b64decode3bytes(size_t idx, char const* in, char* out){
	unsigned int a = 1, b = 1;
	base64::buf24 buf;

	buf.i = 0;

	buf.i |= c_base64_ordrs[in[4 * idx + 0]] << 18;
	buf.i |= c_base64_ordrs[in[4 * idx + 1]] << 12;

	if (in[4 * idx + 2] == '=') a = 0; else buf.i |= c_base64_ordrs[in[4 * idx + 2]] << 6;
	if (in[4 * idx + 3] == '=') b = 0; else buf.i |= c_base64_ordrs[in[4 * idx + 3]] << 0;

	       out[3 * idx + 0] = buf.c[2];
	if (a) out[3 * idx + 1] = buf.c[1];
	if (b) out[3 * idx + 2] = buf.c[0];
}


void b64decode_host(char const* in, size_t in_size, char* out)
{
	for (int i = 0; 4 * i < in_size; i++)
	{
		b64decode3bytes(i, in, out);
	}
}