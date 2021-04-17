#include <vibron/wfn.hpp>

void output_wfn(const size_t & NVecs, const CL::utility::matrix<double> & eigvec,
const std::shared_ptr<vibron::Options> & op, const std::string & prefix) {
    std::vector<std::vector<std::ifstream>> ifs(op->NSegs);
    for (size_t i = 0; i < op->NSegs; i++) {
        ifs[i].resize(op->NStates);
        for (size_t j = 0; j < op->NStates; j++) {
            std::string file = prefix + "-" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + ".wfn";
            ifs[i][j].open(file, std::ifstream::binary);
            if (! ifs[i][j].good()) throw CL::utility::file_error(file);
        }
    }

    std::cout << "Vibronic wave functions are written to wfn*.wfn\n"
              << "Contribution of each basis can be found in population.txt\n";
    std::ofstream ofs; ofs.open("population.txt");
    for (size_t i = 0; i < NVecs; i++) {
        vibron::Wfn wfn(op), Krylov(op);
        wfn = 0.0;
        for (auto & seg : ifs) for (auto & state : seg) state.seekg(0);
        for (const double & el : eigvec[i]) {
            Krylov.read(ifs);
            wfn.add_(el, Krylov);
        }
        wfn.write("wfn" + std::to_string(i));

        ofs << "!!!!!!!!!! Vibronic wave function " << i << ": !!!!!!!!!!\n";
        at::Tensor v = wfn.cat();
        at::Tensor sorted_v, indices;
        std::tie(sorted_v, indices) = v.sort(0, true);
        for (size_t j = 0; j < 10; j++) {
            double population = sorted_v[j].item<double>();
            population *= population;
            ofs << "population = " << population << '\n';
            auto op = wfn.options();
            auto seg_state_vib = wfn.seg_state_vib(indices[j].item<int64_t>());
            size_t vib_index_abs = op->vib_index_abs(seg_state_vib);
            ofs << "electronic state = " << seg_state_vib.second << '\n'
                << "phonons:\n";
            (*op->vib_sets[op->vib_irreds[seg_state_vib.second]])[vib_index_abs].pretty_print(ofs);
            ofs << "-----------------------------------------------\n";
        }
        ofs << '\n';
    }
    ofs.close();
}