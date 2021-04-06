#include <fstream>
#include <tuple>
#include <numeric>

#include <omp.h>

#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

namespace vibron {

// Construct `max_phonons`
void Options::construct_phonons() {
    max_phonons = vib_sets[0]->max_phonons();
    for (size_t i = 1; i < NIrreds; i++) {
        const auto & ith_max = vib_sets[i]->max_phonons();
        if (ith_max.empty()) continue;
        for (size_t irred = 0; irred < NIrreds; irred++)
        for (size_t mode = 0; mode < NModes[irred]; mode++)
        if (ith_max[irred][mode] > max_phonons[irred][mode])
        max_phonons[irred][mode] = ith_max[irred][mode];
    }
}
// Construct `NSegs`, `starts`, `stops`
void Options::construct_segmentation() {
    NSegs = omp_get_max_threads();
    // segment starts
    starts.resize(NSegs);
    for (auto & start : starts) start.resize(NStates);
    // 1st segment
    std::fill(starts[0].begin(), starts[0].end(), 0);
    std::vector<size_t> lengths(NStates);
    for (size_t i = 0; i < NStates; i++) lengths[i] = vib_sets[vib_irreds[i]]->size() / NSegs;
    // 2nd to last segments
    for (size_t j = 1; j < NSegs; j++)
    for (size_t i = 0; i < NStates; i++)
    starts[j][i] = starts[j - 1][i] + lengths[i];
    // segment stops
    stops.resize(NSegs);
    // 1st to (n - 1)-th segments
    for (size_t j = 0; j < NSegs - 1; j++) stops[j] = starts[j + 1];
    // last segment
    stops.back().resize(NStates);
    for (size_t i = 0; i < NStates; i++) stops.back()[i] = vib_sets[vib_irreds[i]]->size();
}

Options::Options() {}
Options::Options(const std::string & wfn_file, const std::vector<std::string> & vib_files) {
    // basic information
    std::ifstream ifs; ifs.open(wfn_file);
    if (! ifs.good()) throw CL::utility::file_error(wfn_file);
    else {
        std::string line;
        std::vector<std::string> strs;
        // Number of irreducible representations
        std::getline(ifs, line);
        std::getline(ifs, line);
        if (! ifs.good()) throw CL::utility::file_error(wfn_file);
        NIrreds = std::stoul(line);
        // Point group product table
        std::getline(ifs, line);
        product_table.resize(NIrreds);
        for (auto & row : product_table) {
            std::getline(ifs, line);
            if (! ifs.good()) throw CL::utility::file_error(wfn_file);
            strs = CL::utility::split(line);
            if (strs.size() != NIrreds) throw CL::utility::file_error(wfn_file);
            for (size_t i = 0; i < NIrreds; i++) row[i] = std::stoul(strs[i]) - 1;
        }
        // Number of normal modes per irreducible
        std::getline(ifs, line);
        NModes.resize(NIrreds);
        std::getline(ifs, line);
        if (! ifs.good()) throw CL::utility::file_error(wfn_file);
        strs = CL::utility::split(line);
        if (strs.size() != NIrreds) throw CL::utility::file_error(wfn_file);
        for (size_t i = 0; i < NIrreds; i++) NModes[i] = std::stoul(strs[i]);
        // Number of electronic states
        std::getline(ifs, line);
        std::getline(ifs, line);
        if (! ifs.good()) throw CL::utility::file_error(wfn_file);
        NStates = std::stoul(line);
        // Vibrational irreducible of each electronic state
        std::getline(ifs, line);
        vib_irreds.resize(NStates);
        std::getline(ifs, line);
        if (! ifs.good()) throw CL::utility::file_error(wfn_file);
        strs = CL::utility::split(line);
        if (strs.size() != NStates) throw CL::utility::file_error(wfn_file);
        for (size_t i = 0; i < NStates; i++) vib_irreds[i] = std::stoul(strs[i]) - 1;
    }
    ifs.close();
    // vibrational basis function details
    if (vib_files.size() != NIrreds) throw std::invalid_argument(
    "vibron::Options: The number of vibrational basis files must equal to the number of irreducible representations");
    vib_sets.resize(NIrreds);
    #pragma omp parallel for
    for (size_t i = 0; i < NIrreds; i++) vib_sets[i] = new VibrationSet(vib_files[i], NIrreds);
    // Infer the rest
    this->construct_phonons();
    this->construct_segmentation();
}
Options::~Options() {}

void Options::pretty_print(std::ostream & stream) const {
    stream << "Number of irreducible representations: " << NIrreds << '\n';
    stream << "Point group product table:\n";
    for (const auto & row : product_table) {
        for (const size_t & el : row) stream << "    " << el + 1;
        stream << '\n';
    }
    stream << "Number of normal modes:\n";
    for (size_t i = 0; i < NModes.size(); i++)
    stream << "    irreducible "<< i << ": " << NModes[i] << '\n';
    stream << "Number of electronic states " << NStates << '\n';
    stream << "Vibrational irreducible:\n";
    for (size_t i = 0; i < vib_irreds.size(); i++)
    stream << "    state " << i + 1 << ": " << vib_irreds[i] + 1 << '\n';
    stream << "Number of segmentations = " << NSegs << '\n';
    for (size_t i = 0; i < NSegs; i++) {
        stream << "Segment " << i << " owns:\n";
        for (size_t j = 0; j < NStates; j++)
        stream << "Electronic state " << j << ": [" << starts[i][j]
               << ", " << stops[i][j] << ")\n";
    }
    for (size_t i = 0; i < vib_sets.size(); i++) {
        stream << "Vibrational basis functions of irreducible " << i << ":\n";
        vib_sets[i]->pretty_print(stream);
    }
}

size_t Options::intdim() const {return std::accumulate(NModes.begin(), NModes.end(), 0);}
// Given the index of a normal mode, return its irreducible and index within the irreducible
std::pair<size_t, size_t> Options::vib_irred_mode(const size_t & index) const {
    size_t irred, mode = index;
    for (irred = 0; irred < NIrreds; irred++)
    if (mode < NModes[irred]) return std::pair<size_t, size_t>(irred, mode);
    else mode -= NModes[irred];
    throw std::invalid_argument("vibron::Options::vib_irred_mode: index out of range");
}
// Given the irreducible and the index within the irreducible of a normal mode, return its index
size_t Options::vib_C1_index(const size_t & irred, const size_t & mode) const {
    size_t index = mode;
    for (size_t i = 0; i < irred; i++) index += NModes[i];
    return index;
}

size_t Options::NVibron() const {
    size_t NVibron = 0;
    for (size_t i = 0; i < NStates; i++) NVibron += vib_sets[vib_irreds[i]]->size();
    return NVibron;
}

// Given a vibrational basis function defined by phonons, return its irreducible
size_t Options::vib_irred(const std::vector<std::vector<size_t>> & phonons) const {
    if (phonons.size() != NIrreds) throw std::invalid_argument(
    "vibron::Options::vib_irred: define phonons for every irreducible");
    size_t irred = 0;
    for (size_t i = 1; i < NIrreds; i++)
    if (phonons[i].size() != NModes[i]) throw std::invalid_argument(
    "vibron::Options::vib_irred: define phonons for every mode");
    else {
        size_t order = 0;
        for (const size_t & phonon : phonons[i]) order += phonon;
        if (order % 2 == 1) irred = product_table[irred][i];
    }
    return irred;
}

// Given segment, electronic state, vibrational index in this segment,
// return the vibrational index in the vibrational basis function set
size_t Options::vib_index_abs(const size_t & seg, const size_t & state, const size_t & seg_vib) const {
    return seg_vib + starts[seg][state];
}
// Given electronic state, vibrational index in the vibrational basis function set
// return segment and vibrational index in this segment
std::pair<size_t, size_t> Options::vib_index(const size_t & state, const size_t & abs_vib) const {
    size_t seg;
    for (seg = 0; seg < NSegs; seg++)
    if (starts[seg][state] <= abs_vib && abs_vib < stops[seg][state])
    return std::pair<size_t, size_t>(seg, abs_vib - starts[seg][state]);
    throw std::invalid_argument("vibron::Options::vib_index: absolute vibrational index out of range");
}

} // namespace vibron