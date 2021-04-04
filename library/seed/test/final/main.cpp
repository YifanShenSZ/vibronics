#include <seed/final.hpp>

int main() {
    std::vector<std::string> vib_files = {"vibration1.in"};
    auto op = std::make_shared<vibron::Options>("symmetry.in", vib_files);
    op->pretty_print(std::cout);

    std::vector<double> data;
    data = CL::utility::read_vector("initial-freq.in");
    at::Tensor init_freq = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64))
                           .clone();
    data = CL::utility::read_vector("final-freq.in");
    at::Tensor final_freq = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64))
                            .clone();
    data = CL::utility::read_vector("transformation-matrix.in");
    at::Tensor T = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64))
                   .clone();
    int64_t dim = sqrt(T.size(0));
    T = T.view({dim, dim});
    T.transpose_(0, 1);
    data = CL::utility::read_vector("shift-vector.in");
    at::Tensor d = at::from_blob(data.data(), data.size(), at::TensorOptions().dtype(torch::kFloat64))
                   .clone();

    vibron::Wfn wfn(op);

    seed::Final final_bias({1.0, 1.0}, op, init_freq, final_freq, T, d);
    double norm = final_bias.generate_seed(wfn);
    std::cout << "Norm of the seed vector = " << norm << '\n';
}