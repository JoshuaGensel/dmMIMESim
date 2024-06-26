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
    typedef unsigned long long id;
    typedef std::vector<id> ids; // species chunk ids

    long double nChoosek(const unsigned n, const unsigned k);

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
