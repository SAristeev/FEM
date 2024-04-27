#pragma once
#ifndef OFMAPSTREAM_H
#define OFMAPSTREAM_H

#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else //POSIX
#include <unistd.h>
#include <cstring>
#include <cerrno>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h> 
#endif //_WIN32

using bool_type =
#ifdef _WIN32
BOOL;
#else //POSIX
bool;
#endif //_WIN32


class ofmapstream
{
	std::string d_filename;

#ifdef _WIN32
	HANDLE d_file_handle = 0;
	HANDLE d_mapping_handle = 0;
	LPVOID d_base_adr = 0;
#else //POSIX
	int d_file_handle = 0;
	void* d_base_adr = 0;
#endif //_WIN32

	uint64_t d_max_size = 0;
	uint64_t d_put_ptr = 0;
	uint64_t d_file_size = 0;

	bool d_mapped = false;

public:
	~ofmapstream();

	void map(std::string const& filename, uint64_t size);
	uint64_t tellp();
	void seekp(uint64_t s);

	char* ptr();

	bool_type flush();
	bool_type unmap();
	//void close();

	ofmapstream& operator<<(std::string const& data);

	struct raw_data_t
	{
		char const* d_ptr = nullptr;
		size_t d_size     = 0;
	};

	ofmapstream& operator<<(raw_data_t const &data);
};


inline ofmapstream::~ofmapstream()
{
	flush();
	unmap();
}


inline void ofmapstream::map(std::string const& filename, uint64_t size)
{
	if (d_mapped)
	{
		throw std::runtime_error("ofmapstream::map(std::string const&, uint64_t) error: a file has already been mapped");
	}

	d_filename = filename;
	d_max_size = size;
	d_file_size = 0;
	d_put_ptr = 0;

#ifdef _WIN32
	d_file_handle = CreateFileA(d_filename.c_str(), GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (d_file_handle == INVALID_HANDLE_VALUE || d_file_handle == NULL)
	{
		throw std::runtime_error("Win32 API error: cannot open file " + d_filename + " error code: " + std::to_string(GetLastError()));
	}

	std::string handle_name = "ofmapstream_" + std::to_string(reinterpret_cast<uint64_t>(this));

	d_mapping_handle = CreateFileMappingA(d_file_handle, NULL, PAGE_READWRITE, d_max_size >> 32, d_max_size & 0x00000000ffffffff, handle_name.c_str());
	if (d_mapping_handle == INVALID_HANDLE_VALUE || d_mapping_handle == NULL)
	{
		CloseHandle(d_file_handle);

		throw std::runtime_error("Win32 API error: cannot create file mapping " + d_filename + " error code: " + std::to_string(GetLastError()));
	}

	d_base_adr = MapViewOfFile(d_mapping_handle, FILE_MAP_WRITE, 0, 0, 0);
	if (d_base_adr == NULL)
	{
		CloseHandle(d_mapping_handle);
		CloseHandle(d_file_handle);

		throw std::runtime_error("Win32 API error: cannot create file view " + d_filename + " error code: " + std::to_string(GetLastError()));
	}
#else //POSIX
	d_file_handle = open(d_filename.c_str(), O_RDWR | O_CREAT);
	if (d_file_handle < 0)
	{
		throw std::runtime_error("open() error: cannot open file " + d_filename + " error code: " + std::strerror(errno));
	}

	if (ftruncate(d_file_handle, d_max_size) != 0)
	{
		throw std::runtime_error("ftruncate() error: cannot truncate file " + d_filename + " error code: " + std::strerror(errno));
	}

	d_base_adr = mmap(0, d_max_size, PROT_READ | PROT_WRITE, MAP_SHARED, d_file_handle, 0);
	if (d_base_adr == MAP_FAILED)
	{
		throw std::runtime_error("mmap() error: cannot create file mapping " + d_filename + " error code: " + std::strerror(errno));
	}

	if (close(d_file_handle) != 0)
	{
		throw std::runtime_error("close() error: cannot close file " + d_filename + " error code: " + std::strerror(errno));
	}

	memset(d_base_adr, 0, d_max_size);
#endif //_WIN32

	d_mapped = true;
}


inline bool_type ofmapstream::flush()
{
	bool_type res = 0;

	if (!d_mapped)
		return res;

#ifdef _WIN32
	if (d_file_size)
	{
		res |= FlushViewOfFile(d_base_adr, d_file_size);
		if (!res)
		{
			throw std::runtime_error("Win32 API error: cannot flush view of file. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}
	}

	if (d_file_handle)
	{
		res |= FlushFileBuffers(d_file_handle);
		if (!res)
		{
			throw std::runtime_error("Win32 API error: cannot flush file buffers. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}
	}
#endif //_WIN32

	return !res;
}


inline bool_type ofmapstream::unmap()
{
	bool_type res = 0;

	if (!d_mapped)
		return res;

#ifdef _WIN32
	if (d_file_handle && d_base_adr && d_mapping_handle)
	{
		res |= UnmapViewOfFile(d_base_adr);
		if (res)
			d_base_adr = NULL;
		else
		{
			throw std::runtime_error("Win32 API error: cannot close file view. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}

		res |= CloseHandle(d_mapping_handle);
		if (res)
			d_mapping_handle = NULL;
		else
		{
			throw std::runtime_error("Win32 API error: cannot close file mapping. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}

		LONG off_hight = d_file_size >> 32;
		SetFilePointer(d_file_handle, d_file_size & 0x00000000ffffffff, &off_hight, 0);
		res |= SetEndOfFile(d_file_handle);
		if (!res)
		{
			throw std::runtime_error("Win32 API error: cannot set end of file. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}

		res |= CloseHandle(d_file_handle);
		if (res)
			d_file_handle = NULL;
		else
		{
			throw std::runtime_error("Win32 API error: cannot close file. error code: " + d_filename + ". error code: " + std::to_string(GetLastError()));
		}
	}

	d_mapping_handle = 0;
#else //POSIX
	if (munmap(d_base_adr, d_max_size) != 0)
	{
		throw std::runtime_error("munmap() error: cannot unmap file " + d_filename + ". error code: " + std::strerror(errno));
	}

	if (truncate(d_filename.c_str(), d_file_size) != 0)
	{
		throw std::runtime_error("truncate() error: cannot truncate file " + d_filename + ". error code: " + std::strerror(errno));
	}
#endif //_WIN32

	d_file_handle = 0;
	d_base_adr = 0;
	d_max_size = 0;
	d_file_size = 0;
	d_put_ptr = 0;
	d_mapped = false;

	return !res;
}


// inline void ofmapstream::close()
// {
// 	flush();
// 	unmap();
// }


inline uint64_t ofmapstream::tellp()
{
	return d_put_ptr;
}


inline char* ofmapstream::ptr()
{
	return (char*)d_base_adr;
}


inline void ofmapstream::seekp(uint64_t s)
{
	if (s <= d_max_size)
	{
		d_put_ptr = s;
		if (d_put_ptr > d_file_size)
			d_file_size = d_put_ptr;
	}
	else
	{
		throw std::runtime_error("ofmapstream::seekp(uint64_t) error: seek exceeds file size " + d_filename);
	}
}


inline ofmapstream& ofmapstream::operator<<(std::string const& data)
{
	if (d_put_ptr + data.size() <= d_max_size)
	{
		memcpy((char*)d_base_adr + d_put_ptr, data.c_str(), data.size());
		d_put_ptr += data.size();

		if (d_put_ptr > d_file_size)
			d_file_size = d_put_ptr;
	}

	return *this;
}


inline ofmapstream& ofmapstream::operator<<(raw_data_t const& data)
{
	if (d_put_ptr + data.d_size <= d_max_size)
	{
		memcpy((char*)d_base_adr + d_put_ptr, data.d_ptr, data.d_size);
		d_put_ptr += data.d_size;

		if (d_put_ptr > d_file_size)
			d_file_size = d_put_ptr;
	}
	else
	{
		throw std::runtime_error("ofmapstream::operator<<(GPUdata const&) error: data exceeds file size " + d_filename);
	}

	return *this;
}


inline std::ofstream& operator<<(std::ofstream &stream, ofmapstream::raw_data_t const& data)
{
	stream.write(data.d_ptr, data.d_size);
	return stream;
}

#endif //OFMAPSTREAM_H
