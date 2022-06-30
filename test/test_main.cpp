//
//  test_main.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 24.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

//##define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_RUNNER // To be able to run the main on your own
#include "Constants.hpp"
#include "catch.hpp"

namespace fs = std::filesystem;
int main(int argc, char* argv[])
{
    std::cout << " blub " << std::endl;
    // TODO verschieder mains mit verschiedenen konstanten erstellen
    // Call with parameters: pre run singleton creation, if no parameters: don't create singleton and run constants test
    // There must be exactly one instance

    //    for(auto i = 1; i<argc; ++i) {
    //        std::cout << argv[i] << std::endl;
    //        //if a parameter file is given with the command -f, open it, read it and create the constants instance
    //        if(std::string(argv[i]) == "-f" && argc>=i+1) {
    //            if(fs::is_regular_file(argv[i+1])) {
    //                fs::path constantFile(argv[i + 1]);
    //                utils::readParameters(constantFile);
    //                //break;
    //            }
    //        }
    //    }

    int result = Catch::Session().run(argc, argv);
    std::cout << " Halloooooo " << result << std::endl;
    return result;
}