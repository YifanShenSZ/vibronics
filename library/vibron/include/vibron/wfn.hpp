#ifndef vibron_wfn_hpp
#define vibron_wfn_hpp

#include <torch/torch.h>

#include <CppLibrary/utility.hpp>

#include <vibron/options.hpp>

namespace vibron {

// We segment the vibronic wave function only by its vibrational part
// i.e. every segment (wfn[i]) contains all electronic states,
// but only a piece of vibration (wfn[i, j]) on each electronic state
class Wfn {
    private:
        std::shared_ptr<Options> op_;
        // data_ holds the underlying data of wfn
        // data_[i][j][k] = wfn[i, j, k] is the k-th vibrational element on the j-th electronic state in the i-th segment
        std::vector<std::vector<at::Tensor>> data_;
        std::vector<std::vector<double *>> data_ptrs_;
        std::vector<std::vector<size_t>> lengthes_;
    public:
        Wfn();
        Wfn(const std::shared_ptr<Options> & _op);
        ~Wfn();

        const std::shared_ptr<Options> & options() const;

        at::Tensor cat() const;

        const std::vector<at::Tensor> & operator[](const size_t & seg) const;
        at::Tensor & operator[](const std::pair<size_t, size_t> & seg_state);
        const at::Tensor & operator[](const std::pair<size_t, size_t> & seg_state) const;

        double & select(const size_t & seg, const size_t & state, const size_t & vib);
        const double & select(const size_t & seg, const size_t & state, const size_t & vib) const;

        // Read the vibronic wave function from files
        void read(const std::string & prefix);
        void read(std::vector<std::vector<std::ifstream>> & ifs);
        void read(std::vector<std::vector<std:: fstream>> & ifs);
        // Write the vibronic wave function to files
        void write(const std::string & prefix) const;
        void write(std::vector<std::vector<std::ofstream>> & ofs) const;
        void write(std::vector<std::vector<std:: fstream>> & ofs) const;

        double dot(const Wfn & other) const;
        void operator=(const double & scalar);
        void mul(const double & scalar, Wfn & result) const;
        void operator*=(const double & scalar);
        void add_(const double & c, const Wfn & add);
        void sub_(const double & c, const Wfn & sub);
        void sub2_(const double & c1, const Wfn & sub1, const double & c2, const Wfn & sub2);
};

} // namespace vibron

#endif