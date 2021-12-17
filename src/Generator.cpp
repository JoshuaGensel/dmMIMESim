#include "Generator.hpp"

#include <iostream>

// Initialize static member instance
Generator* Generator::instance = NULL;

void Generator::create_instance(const unsigned int seed)
{
    if (instance != NULL)
    {
        std::cerr << "Singleton was already instantiated." << std::endl;
        // TODO anderen exception type
        throw std::invalid_argument("Singleton was already instantiated.");
    }

    instance = new Generator(seed);
}

Generator* Generator::get_instance()
{
    if (instance == NULL)
    {
        std::cerr << "Singleton not correctly instantiated." << std::endl;
        // TODO anderen exception type
        throw std::invalid_argument("Singleton not correctly instantiated.");
    }
    return instance;
}