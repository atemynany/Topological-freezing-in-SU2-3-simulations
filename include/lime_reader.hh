// lime_reader.hh - LIME/ILDG gauge field reader
// Reads CL2QCD LIME format configs and converts to your binary format
#ifndef LIME_READER_HH
#define LIME_READER_HH

#include <string>
#include <cstdint>

// Note: LimeReader is defined in lime.h, no forward declaration needed

namespace lime_io {

/**
 * Read a LIME/ILDG format gauge field and convert to our format.
 * 
 * LIME uses index order (x,y,z,t) with mu as innermost
 * Our format uses (t,x,y,z) with mu as part of link index
 * 
 * @param filename  Path to LIME file
 * @param gf        Output buffer (must be pre-allocated: T*L^3*4*18 doubles)
 * @param T         Expected temporal extent
 * @param L         Expected spatial extent
 * @return          true on success, false on error
 */
bool read_lime_gaugefield(const std::string& filename, double* gf, int T, int L);

/**
 * Convert CL2QCD index order to our index order in-place.
 * CL2QCD: site = x + y*L + z*L^2 + t*L^3, link = mu + 4*site
 * Ours:   site = ((t*L+x)*L+y)*L+z, link = (4*site + mu) * 18
 * 
 * @param out       Output buffer in our format
 * @param in        Input buffer in CL2QCD format
 * @param T         Temporal extent
 * @param L         Spatial extent
 */
void convert_cl2qcd_to_our_format(double* out, const double* in, int T, int L);

/**
 * Write gauge field in our binary format.
 * Format: ASCII header line + raw binary doubles
 * 
 * @param filename  Output file path
 * @param gf        Gauge field buffer
 * @param T         Temporal extent
 * @param L         Spatial extent
 * @param conf_num  Configuration number (for header)
 * @return          true on success
 */
bool write_raw_gaugefield(const std::string& filename, const double* gf, 
                          int T, int L, int conf_num);

/**
 * Get lattice dimensions from LIME file header.
 * 
 * @param filename  Path to LIME file
 * @param T         Output: temporal extent
 * @param L         Output: spatial extent
 * @return          true on success
 */
bool get_lime_dimensions(const std::string& filename, int* T, int* L);

} // namespace lime_io

#endif // LIME_READER_HH
