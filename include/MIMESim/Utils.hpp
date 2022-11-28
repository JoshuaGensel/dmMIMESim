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

    long double nChoosek(const unsigned n, const unsigned k);
    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L);
    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q);
    mutatedPositions specIdxToMutPos(const unsigned long specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int>& nMutRange);

    enum SampleID : int
    {
        wt_bound,
        wt_unbound,
        mut_bound,
        mut_unbound,
        mut_bound_bound,
        mut_bound_unbound,
        mut_unbound_bound,
        mut_unbound_unbound,
    };

    // return string identifier for sample as defined in SampleID enum
    std::string SampleIDStr(SampleID id);

    /**
     * Split string s into elements separated by delim
     */
    void split_string(const std::string& s, char delim, std::vector<std::string>& elems);
}

#endif /* Utils_hpp */
