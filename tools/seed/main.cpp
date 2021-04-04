#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

#include <Lanczos/mv.hpp>
#include <Lanczos/Lanczos.hpp>

#include <seed/final.hpp>

argparse::ArgumentParser parse_args(const size_t & argc, const char ** & argv) {
    CL::utility::echo_command(argc, argv, std::cout);
    std::cout << '\n';
    argparse::ArgumentParser parser("Generate Lanczos seed vector");

    // required arguments
    parser.add_argument("-w","--wfn",                 1, false, "vibronic wave function definition file");
    parser.add_argument("-v","--vibration",         '+', false, "vibrational basis definition files");
    parser.add_argument("-i","--initial_frequency", '+', false, "frequency of each initial-state normal mode");
    parser.add_argument("-f","--final_frequency",   '+', false, "frequency of each final-state normal mode");
    parser.add_argument("-t","--tran_matrix",         1, false, "transformation matrix from final-state to initial-state normal mode");
    parser.add_argument("-s","--shift_vector",        1, false, "shift vector from final-state to initial-state normal mode");

    // optional arguments
    parser.add_argument("-d","--dipole", '+', true, "one component or norm of transition dipole");

    parser.parse_args(argc, argv);
    return parser;
}

int main(size_t argc, const char ** argv) {
    std::cout << "Generate Lanczos seed vector\n"
              << "Yifan Shen 2021\n\n";
    argparse::ArgumentParser args = parse_args(argc, argv);
    CL::utility::show_time(std::cout);
    std::cout << '\n';

    auto wfn_file  = args.retrieve<std::string>("wfn");
    auto vib_files = args.retrieve<std::vector<std::string>>("vibration");
    auto op = std::make_shared<vibron::Options>(wfn_file, vib_files);
    op->pretty_print(std::cout);
    vibron::Wfn wfn(op);

    std::vector<double> data;
    std::vector<at::Tensor> tensors;
    // frequency of each initial-state normal mode
    auto init_freq_files = args.retrieve<std::vector<std::string>>("initial_frequency");
    tensors.resize(init_freq_files.size());
    for (size_t i = 0; i < tensors.size(); i++) {
        data = CL::utility::read_vector(init_freq_files[i]);
        tensors[i] = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64)).clone();
    }
    at::Tensor init_freq = at::cat(tensors).clone();
    // frequency of each final-state normal mode
    auto final_freq_files = args.retrieve<std::vector<std::string>>("final_frequency");
    tensors.resize(final_freq_files.size());
    for (size_t i = 0; i < tensors.size(); i++) {
        data = CL::utility::read_vector(final_freq_files[i]);
        tensors[i] = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64)).clone();
    }
    at::Tensor final_freq = at::cat(tensors).clone();
    // transformation matrix
    auto mat_file = args.retrieve<std::string>("tran_matrix");
    data = CL::utility::read_vector(mat_file);
    at::Tensor tran_matrix = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64)).clone();
    int64_t intdim = init_freq.size(0);
    tran_matrix = tran_matrix.view({intdim, intdim});
    // shift vector
    auto vec_file = args.retrieve<std::string>("shift_vector");
    data = CL::utility::read_vector(vec_file);
    at::Tensor shift_vector = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64)).clone();

    std::vector<double> dipole(op->NStates, 1.0);
    if (args.gotArgument("dipole")) {
        dipole = args.retrieve<std::vector<double>>("dipole");
        if (dipole.size() != op->NStates) throw std::invalid_argument(
        "the number of dipole elements must match the number of electronic states");
    }

    seed::Final final_bias(dipole, op, init_freq, final_freq, tran_matrix, shift_vector);
    double norm = final_bias.generate_seed(wfn);
    std::cout << "\nNorm of the seed vector = " << norm << '\n';
    std::cout << "The seed vector can be found in seed-*-*.wfn\n";
    std::vector<std::vector<std::ofstream>> ofs(op->NSegs);
    for (size_t i = 0; i < op->NSegs; i++) {
        ofs[i].resize(op->NStates);
        for (size_t j = 0; j < op->NStates; j++)
        ofs[i][j].open("seed-" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + ".wfn");
    }
    wfn.write(ofs);

    std::cout << '\n';
    CL::utility::show_time(std::cout);
    std::cout << "Mission success\n";
}