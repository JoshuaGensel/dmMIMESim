#include "Generator.hpp"

#include <iostream>

// Initialize static member instance
Generator* Generator::instance = NULL;

void Generator::create_instance(const unsigned int seed)
{
    if (instance != NULL)
    {
        std::cerr << "Generator was already instantiated (create after create)." << std::endl;
        throw std::runtime_error("Generator was already instantiated (create after create).");
    }

    instance = new Generator(seed);
}

Generator* Generator::get_instance()
{
    if (instance == NULL)
    {
        std::cerr << "Generator not correctly instantiated (get before create or after release)." << std::endl;
        throw std::runtime_error("Generator not correctly instantiated (get before create or after release).");
    }
    return instance;
}

void Generator::release_instance()
{
    if (instance == NULL)
    {
        std::cerr << "Deletion of not correctly instantiated Generator (release after release)." << std::endl;
        throw std::runtime_error("Deletion of not correctly instantiated Generator (release after release).");
    }
    delete instance;
    instance = NULL;
}