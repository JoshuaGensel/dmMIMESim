#include <chrono>
#include <random>

class Generator
{

  private:
    static Generator* instance;

    Generator(const unsigned int seed) : engine{std::default_random_engine(seed)} {};

    Generator(const Generator&) = delete;
    Generator(Generator&&) = delete;
    void operator=(const Generator&) = delete;

  public:
    static void create_instance(const unsigned int seed);
    static Generator* get_instance();
    static void release_instance();
    std::default_random_engine engine;
};