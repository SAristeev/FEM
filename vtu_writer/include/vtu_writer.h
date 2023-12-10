#pragma once
#ifndef VTU_WRITER_H
#define VTU_WRITER_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <span>
#include <numeric>
#include <concepts>
#include <ranges>
#include <algorithm>
#include <string_view>

namespace vtu 
{
	template <class T> struct type {};
	template <> struct type<int8_t>   { static constexpr std::string_view name = "Int8";    };
	template <> struct type<uint8_t>  { static constexpr std::string_view name = "UInt8";   };
	template <> struct type<int16_t>  { static constexpr std::string_view name = "Int16";   };
	template <> struct type<uint16_t> { static constexpr std::string_view name = "UInt16";  };
	template <> struct type<int32_t>  { static constexpr std::string_view name = "Int32";   };
	template <> struct type<uint32_t> { static constexpr std::string_view name = "UInt32";  };
	template <> struct type<int64_t>  { static constexpr std::string_view name = "Int64";   };
	template <> struct type<uint64_t> { static constexpr std::string_view name = "UInt64";  };
	template <> struct type<double>   { static constexpr std::string_view name = "Float64"; };
	template <> struct type<float>    { static constexpr std::string_view name = "Float32"; };

	class spatial_data_t;
	class writer_t;
}


class vtu::spatial_data_t
{
	char const*              m_data;
	std::string              m_name;
	size_t                   m_type_size;
	std::string              m_type_name;
	std::vector<std::string> m_comp_names;
	size_t                   m_ncomponents;

public:
	template <class T>
	spatial_data_t(T const* data, std::string const& name, std::vector<std::string> const& component_names = {})
	:
		m_data{ reinterpret_cast<char const*>(data) },
		m_name{ name },
		m_type_size{ sizeof(T) },
		m_type_name{ vtu::type<T>::name },
		m_comp_names{ component_names },
		m_ncomponents{ component_names.size() }
	{}

	template <class T>
	spatial_data_t(T const* data, std::string const& name, size_t ncomponents)
	:
		spatial_data_t(data, name)
	{
		m_ncomponents = ncomponents;
	}

	char const*        data()                   const;
	std::string const& name()                   const;
	std::string const& type_name()              const;
	size_t             n_components()           const;
	bool               has_components()         const;
	bool               has_component_names()    const;
	std::string const& component_name(size_t i) const;
	size_t             item_size()              const;
};


class vtu::writer_t
{
	class append_data_block_t
	{
		static constexpr size_t s_header_size = sizeof(uint64_t);

		char const* m_data = nullptr;
		uint64_t    m_data_size = 0;
		size_t      m_offset = 0;
		bool        m_b64 = true;

	public:
		append_data_block_t() = default;
		append_data_block_t(char const* data, size_t data_size, bool base64 = true);

		size_t size() const;
		size_t offset() const;
		void   offset(size_t val);
		void   write(char* dst) const;
	};

	std::string m_filename;
	bool        m_b64;

	std::vector<append_data_block_t> m_vtu_blocks;
	std::vector<append_data_block_t> m_point_blocks;
	std::vector<append_data_block_t> m_cell_blocks;
	std::string m_float_name, m_size_name, m_char_name;

public:
	writer_t(std::string const& filename, bool b64encode = true);

	template <
		class float_range, 
		class sizet_range, 
		class ctype_range
	>
	void write(
		float_range p_cords,
		sizet_range c_conns,
		ctype_range c_types,
		std::vector<spatial_data_t>& p_data = {},
		std::vector<spatial_data_t>& c_data = {}
	) requires 
		std::ranges::sized_range<float_range> && std::ranges::random_access_range<float_range> &&
		std::ranges::sized_range<sizet_range> && std::ranges::random_access_range<sizet_range> &&
		                                         std::ranges::random_access_range<ctype_range>
	{
		using float_type = std::ranges::range_value_t<float_range>;
		using size_type  = std::ranges::range_value_t<sizet_range>;
		using char_type  = std::ranges::range_value_t<ctype_range>;

		m_float_name = std::string{ vtu::type<float_type>::name };
		m_size_name  = std::string{ vtu::type<size_type>::name }; 
		m_char_name  = std::string{ vtu::type<char_type>::name };

		static_assert(std::same_as<char_type, unsigned char>);
		assert(std::ranges::size(p_cords) % 3 == 0);

		size_t const n_points = std::ranges::size(p_cords) / 3;
		size_t const n_cells  = std::ranges::size(c_types);

		// create cell offsets
		auto c_verts = std::views::transform(c_types, 
			[&](char_type const t) -> size_type {
				switch (t)
				{
					case 1:  return 1; // vertex
					case 3:  return 2; // line
					case 5:  return 3; // triangle
					case 8:  return 4; // pixel
					case 9:  return 4; // quad
					case 10: return 4; // tetra
					case 11: return 8; // voxel
					case 12: return 8; // hexahedron
					case 13: return 6; // wedge
					case 14: return 5; // pyramid
					default: return 0;
				}
			}
		);
		std::span<size_type> c_offsets{ new size_type[n_cells], n_cells };
		std::inclusive_scan(c_verts.begin(), c_verts.end(), c_offsets.begin());
		size_t const c_conns_size = c_offsets.back();
		
		// create append data blocks
		m_point_blocks.clear();
		for (auto& pd : p_data)
			m_point_blocks.emplace_back(pd.data(), n_points * pd.item_size(), m_b64);

		m_cell_blocks.clear();
		for (auto& cd : c_data)
			m_cell_blocks.emplace_back(cd.data(), n_cells * cd.item_size(), m_b64);

		m_vtu_blocks.resize(4);
		m_vtu_blocks[0] = append_data_block_t(reinterpret_cast<char const*>(std::ranges::cdata(p_cords))  , n_points * 3 * sizeof(float_type), m_b64);
		m_vtu_blocks[1] = append_data_block_t(reinterpret_cast<char const*>(std::ranges::cdata(c_conns))  , c_conns_size * sizeof(size_type) , m_b64);
		m_vtu_blocks[2] = append_data_block_t(reinterpret_cast<char const*>(std::ranges::cdata(c_offsets)), n_cells      * sizeof(size_type) , m_b64);
		m_vtu_blocks[3] = append_data_block_t(reinterpret_cast<char const*>(std::ranges::cdata(c_types))  , n_cells      * sizeof(char_type) , m_b64);

		// setting append data block offsets
		auto init_seq_data_block = [&](append_data_block_t& block, append_data_block_t const& from) -> size_t {
			block.offset(from.offset() + from.size());
			return block.offset() + block.size();
		};

		auto init_seq_data_blocks = [&](std::vector<append_data_block_t>& data_blocks, size_t init = 0) -> size_t {
			size_t total = init;

			if (!data_blocks.empty())
			{
				data_blocks[0].offset(init);
				total += data_blocks[0].size();

				for (auto it = data_blocks.begin() + 1; it != data_blocks.end(); it++)
					total = init_seq_data_block(*it, *(it - 1));
			}

			return total;
		};

		size_t const total_offset = init_seq_data_blocks(m_vtu_blocks, init_seq_data_blocks(m_cell_blocks, init_seq_data_blocks(m_point_blocks)));
		write_file(p_data, c_data, total_offset, n_points, n_cells);
	}

private:
	std::string get_endian();
	std::string get_encoding();
	void write_file(
		std::vector<spatial_data_t>& p_data,
		std::vector<spatial_data_t>& c_data,
		size_t      total_offset,
		size_t      n_points,
		size_t      n_cells
	);
};


#endif // VTU_WRITER_H
