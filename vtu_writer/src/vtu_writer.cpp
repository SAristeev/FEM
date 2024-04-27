#include <tuple>

#include "vtu_writer.h"
#include "ofmapstream.h"
#include "base64.h"


char const* vtu::spatial_data_t::data() const
{
	return m_data;
}


std::string const& vtu::spatial_data_t::name() const
{
	return m_name;
}


std::string const& vtu::spatial_data_t::type_name() const
{
	return m_type_name;
}


size_t vtu::spatial_data_t::n_components() const
{
	return m_ncomponents;
}


bool vtu::spatial_data_t::has_components() const
{
	return m_ncomponents > 0;
}


bool vtu::spatial_data_t::has_component_names() const
{
	return m_ncomponents > 0 ? m_ncomponents == m_comp_names.size() : false;
}


std::string const& vtu::spatial_data_t::component_name(size_t i) const
{
	return m_comp_names[i];
}


size_t vtu::spatial_data_t::item_size() const
{
	if (has_components())
		return m_ncomponents * m_type_size;
	return m_type_size;
}


vtu::writer_t::append_data_block_t::append_data_block_t(char const* data, size_t data_size, bool base64)
:
	m_data{ data },
	m_data_size{ data_size },
	m_offset{ 0 },
	m_b64{ base64 }
{}


size_t vtu::writer_t::append_data_block_t::size() const
{
	if (m_b64)
		return base64::esize(s_header_size) + base64::esize(m_data_size);

	return s_header_size + m_data_size;
}


size_t vtu::writer_t::append_data_block_t::offset() const
{
	return m_offset;
}


void vtu::writer_t::append_data_block_t::offset(size_t val)
{
	m_offset = val;
}


void vtu::writer_t::append_data_block_t::write(char* dst) const
{
	if (m_b64)
	{
		base64::encode(reinterpret_cast<char const*>(&m_data_size), s_header_size, dst);
		base64::encode(m_data, m_data_size, dst + base64::esize(s_header_size));
	}
	else
	{
		memcpy(dst, reinterpret_cast<char const*>(&m_data_size), s_header_size);
		memcpy(dst + s_header_size, m_data, m_data_size);
	}
}


vtu::writer_t::writer_t(std::string const& filename, bool b64encode)
:
	m_filename{ filename },
	m_b64{ b64encode }
{}


std::string vtu::writer_t::get_endian()
{
	union
	{
		uint32_t i;
		uint8_t c[4];
	} check_endian = { 0x01000000 };

	if (check_endian.c[0] == 1)
		return "BigEndian";

	return "LittleEndian";
}


std::string vtu::writer_t::get_encoding()
{
	if (m_b64)
		return "base64";
	return "raw";
}

void vtu::writer_t::write_file(
	std::vector<spatial_data_t>& p_data,
	std::vector<spatial_data_t>& c_data,
	size_t      total_offset,
	size_t      n_points,
	size_t      n_cells
)
{
	append_data_block_t& p_cords_block   = m_vtu_blocks[0];
	append_data_block_t& c_conns_block   = m_vtu_blocks[1];
	append_data_block_t& c_offsets_block = m_vtu_blocks[2];
	append_data_block_t& c_types_block   = m_vtu_blocks[3];

	ofmapstream out_file;
	out_file.map(m_filename, total_offset + 10'000 + 1000 * (p_data.size() + c_data.size()));

	auto write_spatial_data = [&](std::vector<spatial_data_t> const& datas, std::vector<append_data_block_t> const& blocks) -> void {
		assert(datas.size() == blocks.size());
		for (size_t i = 0; i < datas.size(); i++)
		{
			spatial_data_t const& data = datas[i];

			out_file << "        <DataArray type=\"" << data.type_name() << "\" Name=\"" << data.name() << "\" ";

			if (data.has_components())
			{
				out_file << "NumberOfComponents=\"" << std::to_string(data.n_components()) << "\" ";

				if (data.has_component_names())
					for (int i = 0; i < data.n_components(); i++)
						out_file << "ComponentName" << std::to_string(i) << "=\"" << data.component_name(i) << "\" ";
			}

			out_file << "format=\"appended\" offset=\"" << std::to_string(blocks[i].offset()) << "\"/>\n";
		}
	};

	out_file << "<?xml version=\"1.0\"?>\n"
			 << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + get_endian() + "\" header_type=\"UInt64\">\n"
			 << "  <UnstructuredGrid>\n"
			 << "    <Piece NumberOfPoints=\"" << std::to_string(n_points) << "\" NumberOfCells=\"" << std::to_string(n_cells) << "\">\n";

	if (!p_data.empty())
	{
		out_file << "      <PointData>\n";
		write_spatial_data(p_data, m_point_blocks);
		out_file << "      </PointData>\n";
	}

	if (!c_data.empty())
	{
		out_file << "      <CellData>\n";
		write_spatial_data(c_data, m_cell_blocks);
		out_file << "      </CellData>\n";
	}
	
	out_file << "      <Points>\n"
		     << "        <DataArray type=\"" << m_float_name << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << std::to_string(p_cords_block.offset()) << "\"/>\n"
		     << "      </Points>\n"
		     << "      <Cells>\n"
		     << "        <DataArray type=\"" << m_size_name << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << std::to_string(c_conns_block.offset()) << "\"/>\n"
		     << "        <DataArray type=\"" << m_size_name << "\" Name=\"offsets\" format=\"appended\" offset=\"" << std::to_string(c_offsets_block.offset()) << "\"/>\n"
		     << "        <DataArray type=\"" << m_char_name << "\" Name=\"types\" format=\"appended\" offset=\"" << std::to_string(c_types_block.offset()) << "\"/>\n"
			 << "      </Cells>\n"
			 << "    </Piece>\n"
			 << "  </UnstructuredGrid>\n"
			 << "  <AppendedData encoding=\"" + get_encoding() + "\">\n   _";

	uint64_t append_data_offset = out_file.tellp();

	out_file.seekp(append_data_offset + total_offset);
	out_file << "\n  </AppendedData>\n</VTKFile>\n";

	// shifting blocks by meta data
	auto shift_data_block = [&](append_data_block_t& block) -> void {
		block.offset(block.offset() + append_data_offset);
	};

	std::ranges::for_each(m_point_blocks, shift_data_block);
	std::ranges::for_each(m_cell_blocks , shift_data_block);
	std::ranges::for_each(m_vtu_blocks  , shift_data_block);
	
	// writting append data
	auto write_block = [&](append_data_block_t const& block) -> void {
		block.write(out_file.ptr() + block.offset());
	};

	std::ranges::for_each(m_point_blocks, write_block);
	std::ranges::for_each(m_cell_blocks , write_block);
	std::ranges::for_each(m_vtu_blocks  , write_block);

	// file will closed in destructor
	//out_file.close();
}
