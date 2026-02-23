// lime_reader.cc - LIME/ILDG gauge field reader
// Reads CL2QCD LIME format configs and converts to our binary format
//
// Compile with: g++ -O3 -o lime_to_raw lime_reader.cc -llime
// Requires: c-lime library (https://github.com/usqcd-software/c-lime)

#include "lime_reader.hh"

extern "C" {
#include <lime.h>
}

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <arpa/inet.h>  // For byte swapping

namespace lime_io {

// Check if system is big-endian
static bool is_big_endian() {
    union { uint32_t i; char c[4]; } test = {0x01020304};
    return test.c[0] == 1;
}

// Swap bytes for double (ILDG is big-endian)
static double swap_double(double d) {
    double result;
    char* src = (char*)&d;
    char* dst = (char*)&result;
    for (int i = 0; i < 8; i++) {
        dst[i] = src[7-i];
    }
    return result;
}

void convert_cl2qcd_to_our_format(double* out, const double* in, int T, int L) {
    // CL2QCD index: site_cl2qcd = x + y*L + z*L^2 + t*L^3
    //               link_cl2qcd = mu + 4*site_cl2qcd
    //               elem_cl2qcd = link_cl2qcd * 18 + matrix_element
    //
    // Our index:    site_ours = ((t*L+x)*L+y)*L+z
    //               link_ours = (4*site_ours + mu) * 18
    
    const int volume = T * L * L * L;
    
    for (int t = 0; t < T; t++) {
        for (int x = 0; x < L; x++) {
            for (int y = 0; y < L; y++) {
                for (int z = 0; z < L; z++) {
                    // CL2QCD site index (x,y,z,t order)
                    int site_cl2qcd = x + y*L + z*L*L + t*L*L*L;
                    
                    // Our site index (t,x,y,z order)
                    int site_ours = ((t*L + x)*L + y)*L + z;
                    
                    for (int mu = 0; mu < 4; mu++) {
                        // CL2QCD link index (mu innermost)
                        int link_cl2qcd = (mu + 4*site_cl2qcd) * 18;
                        
                        // Our link index (site*4 + mu)
                        int link_ours = (4*site_ours + mu) * 18;
                        
                        // Copy 18 doubles (SU(3) matrix)
                        memcpy(out + link_ours, in + link_cl2qcd, 18 * sizeof(double));
                    }
                }
            }
        }
    }
}

bool get_lime_dimensions(const std::string& filename, int* T, int* L) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Cannot open LIME file: " << filename << std::endl;
        return false;
    }
    
    LimeReader* reader = limeCreateReader(fp);
    if (!reader) {
        fclose(fp);
        return false;
    }
    
    int status;
    while ((status = limeReaderNextRecord(reader)) != LIME_EOF) {
        if (status != LIME_SUCCESS) continue;
        
        const char* type = limeReaderType(reader);
        
        // Look for ILDG format info
        if (strcmp(type, "ildg-format") == 0) {
            n_uint64_t nbytes = limeReaderBytes(reader);
            std::vector<char> buf(nbytes + 1);
            limeReaderReadData(buf.data(), &nbytes, reader);
            buf[nbytes] = '\0';
            
            // Parse XML for dimensions
            std::string xml(buf.data());
            
            // Simple parsing - look for <lx>, <ly>, <lz>, <lt>
            auto extract_dim = [&xml](const std::string& tag) -> int {
                std::string open = "<" + tag + ">";
                std::string close = "</" + tag + ">";
                size_t start = xml.find(open);
                size_t end = xml.find(close);
                if (start != std::string::npos && end != std::string::npos) {
                    start += open.length();
                    return std::stoi(xml.substr(start, end - start));
                }
                return -1;
            };
            
            int lx = extract_dim("lx");
            int ly = extract_dim("ly");
            int lz = extract_dim("lz");
            int lt = extract_dim("lt");
            
            if (lx > 0 && ly == lx && lz == lx && lt > 0) {
                *L = lx;
                *T = lt;
                limeDestroyReader(reader);
                fclose(fp);
                return true;
            }
        }
    }
    
    limeDestroyReader(reader);
    fclose(fp);
    return false;
}

bool read_lime_gaugefield(const std::string& filename, double* gf, int T, int L) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Cannot open LIME file: " << filename << std::endl;
        return false;
    }
    
    LimeReader* reader = limeCreateReader(fp);
    if (!reader) {
        std::cerr << "Error: Cannot create LIME reader" << std::endl;
        fclose(fp);
        return false;
    }
    
    const size_t volume = T * L * L * L;
    const size_t expected_bytes = volume * 4 * 18 * sizeof(double);
    bool found_data = false;
    
    int status;
    while ((status = limeReaderNextRecord(reader)) != LIME_EOF) {
        if (status != LIME_SUCCESS) continue;
        
        const char* type = limeReaderType(reader);
        
        // Look for binary gauge data record
        if (strcmp(type, "ildg-binary-data") == 0 || 
            strcmp(type, "scidac-binary-data") == 0) {
            
            n_uint64_t nbytes = limeReaderBytes(reader);
            
            if (nbytes != expected_bytes) {
                std::cerr << "Warning: Expected " << expected_bytes 
                          << " bytes, got " << nbytes << std::endl;
            }
            
            // Read into temporary buffer (ILDG order)
            std::vector<double> temp(volume * 4 * 18);
            n_uint64_t bytes_to_read = nbytes;
            status = limeReaderReadData((char*)temp.data(), &bytes_to_read, reader);
            
            if (status != LIME_SUCCESS) {
                std::cerr << "Error reading gauge data" << std::endl;
                continue;
            }
            
            // Byte swap if needed (ILDG is big-endian)
            if (!is_big_endian()) {
                for (size_t i = 0; i < temp.size(); i++) {
                    temp[i] = swap_double(temp[i]);
                }
            }
            
            // Convert from CL2QCD/ILDG order to our order
            convert_cl2qcd_to_our_format(gf, temp.data(), T, L);
            
            found_data = true;
            break;
        }
    }
    
    limeDestroyReader(reader);
    fclose(fp);
    
    if (!found_data) {
        std::cerr << "Error: No gauge data found in LIME file" << std::endl;
        return false;
    }
    
    return true;
}

bool write_raw_gaugefield(const std::string& filename, const double* gf, 
                          int T, int L, int conf_num) {
    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Cannot create output file: " << filename << std::endl;
        return false;
    }
    
    // Write header line
    fprintf(fp, "# SU3 config T=%d L=%d conf=%d\n", T, L, conf_num);
    
    // Write binary data
    const size_t n_doubles = T * L * L * L * 4 * 18;
    size_t written = fwrite(gf, sizeof(double), n_doubles, fp);
    
    fclose(fp);
    
    if (written != n_doubles) {
        std::cerr << "Error: Only wrote " << written << " of " << n_doubles << " doubles" << std::endl;
        return false;
    }
    
    return true;
}

} // namespace lime_io

// =============================================================================
// Main program: lime_to_raw converter
// =============================================================================
#ifdef LIME_READER_MAIN

#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <regex>

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options] <input> <output_dir>\n"
              << "\n"
              << "Convert LIME/ILDG gauge configs to raw binary format.\n"
              << "\n"
              << "Arguments:\n"
              << "  input       LIME file or directory containing LIME files\n"
              << "  output_dir  Directory for converted configs\n"
              << "\n"
              << "Options:\n"
              << "  -T <int>    Expected temporal extent (auto-detected if not given)\n"
              << "  -L <int>    Expected spatial extent (auto-detected if not given)\n"
              << "  -p <str>    Output prefix (default: 'conf_')\n"
              << "  -h          Show this help\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    int T = 0, L = 0;
    std::string prefix = "conf_";
    
    // Parse arguments
    int arg_idx = 1;
    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        std::string opt = argv[arg_idx];
        if (opt == "-T" && arg_idx + 1 < argc) {
            T = std::stoi(argv[++arg_idx]);
        } else if (opt == "-L" && arg_idx + 1 < argc) {
            L = std::stoi(argv[++arg_idx]);
        } else if (opt == "-p" && arg_idx + 1 < argc) {
            prefix = argv[++arg_idx];
        } else if (opt == "-h") {
            print_usage(argv[0]);
            return 0;
        }
        arg_idx++;
    }
    
    if (argc - arg_idx < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string input_path = argv[arg_idx];
    std::string output_dir = argv[arg_idx + 1];
    
    // Create output directory
    mkdir(output_dir.c_str(), 0755);
    
    // Collect input files
    std::vector<std::string> files;
    struct stat st;
    if (stat(input_path.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            // Directory: find all LIME files
            DIR* dir = opendir(input_path.c_str());
            if (dir) {
                struct dirent* entry;
                while ((entry = readdir(dir)) != nullptr) {
                    std::string name = entry->d_name;
                    if (name.find(".lime") != std::string::npos ||
                        name.find("conf") != std::string::npos) {
                        files.push_back(input_path + "/" + name);
                    }
                }
                closedir(dir);
            }
        } else {
            // Single file
            files.push_back(input_path);
        }
    }
    
    std::sort(files.begin(), files.end());
    
    if (files.empty()) {
        std::cerr << "Error: No LIME files found" << std::endl;
        return 1;
    }
    
    std::cout << "Found " << files.size() << " LIME file(s)" << std::endl;
    
    // Auto-detect dimensions from first file if not specified
    if (T == 0 || L == 0) {
        if (!lime_io::get_lime_dimensions(files[0], &T, &L)) {
            std::cerr << "Error: Cannot auto-detect dimensions. Specify -T and -L" << std::endl;
            return 1;
        }
        std::cout << "Auto-detected: T=" << T << ", L=" << L << std::endl;
    }
    
    // Allocate buffer
    const size_t buffer_size = T * L * L * L * 4 * 18;
    std::vector<double> gf(buffer_size);
    
    // Convert each file
    int conf_num = 0;
    std::regex num_regex("([0-9]+)");
    
    for (const auto& file : files) {
        std::cout << "Converting: " << file << std::endl;
        
        // Try to extract config number from filename
        std::smatch match;
        std::string basename = file.substr(file.rfind('/') + 1);
        if (std::regex_search(basename, match, num_regex)) {
            conf_num = std::stoi(match[1]);
        } else {
            conf_num++;
        }
        
        // Read LIME file
        if (!lime_io::read_lime_gaugefield(file, gf.data(), T, L)) {
            std::cerr << "  Failed to read, skipping" << std::endl;
            continue;
        }
        
        // Write in our format
        char outname[256];
        snprintf(outname, sizeof(outname), "%s/%s%05d.bin", 
                 output_dir.c_str(), prefix.c_str(), conf_num);
        
        if (lime_io::write_raw_gaugefield(outname, gf.data(), T, L, conf_num)) {
            std::cout << "  -> " << outname << std::endl;
        }
    }
    
    std::cout << "Done!" << std::endl;
    return 0;
}

#endif // LIME_READER_MAIN
