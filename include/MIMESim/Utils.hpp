//
//  Utils.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 22.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <filesystem>
#include <map>
#include <vector>

namespace fs = std::filesystem;

namespace utils
{
    typedef std::map<unsigned int, unsigned int> mutatedPositions;

    // TODO: refactoring -> besesr zu DCA und in DCABenchmark benutzen
    unsigned long nChoosek(const unsigned n, const unsigned k);
    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L);
    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q);
    // TODO QUESTION ist in Species vorhanden.... umschiften, ebenso wie mutPosToIndex?
    mutatedPositions specIdxToMutPos(const unsigned long specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int>& nMutRange);

}

#endif /* Utils_hpp */
