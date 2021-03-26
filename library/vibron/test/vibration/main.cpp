#include <iostream>
#include <memory>
#include <random>

#include <vibron/vibration.hpp>

int main() {
    std::vector<std::shared_ptr<vibron::VibrationSet>> sets(2);
    sets[0] = std::make_shared<vibron::VibrationSet>("vibration1.in", 2);
    sets[1] = std::make_shared<vibron::VibrationSet>("vibration2.in", 2);

    std::cout << "Index a vibrational basis function:\n";
    srand((unsigned)time(NULL));
    for (size_t i = 0; i < 20; i++) {
        size_t index = rand() % sets[0]->size();
        std::cout << sets[0]->index_vibration((*sets[0])[index]) - index << ' ';
    }
    std::cout << '\n';
}