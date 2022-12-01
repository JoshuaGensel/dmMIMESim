//
//  Utils.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 22.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Utils.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>

namespace utils
{

    long double nChoosek(unsigned n, unsigned k)
    {
        long double result;
        if (k > n)
            return 0;
        if (k == 0)
            return 1;
        if (k < n - k)
            result = k * std::beta(static_cast<long double>(k), static_cast<long double>(n - k + 1));
        else
            result = (n - k) * std::beta(static_cast<long double>(k + 1), static_cast<long double>(n - k));
        if (result == 0)
            throw std::overflow_error("Overflow error in nChoosek");
        result = 1 / result;

        return result;
    }

    std::string SampleIDStr(SampleID id)
    {
        std::string strings[8] = {"1", "2", "3", "4", "5", "6", "7", "8"};
        return strings[id];
    }

    void split_string(const std::string& s, char delim, std::vector<std::string>& elems)
    {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim))
        {
            elems.push_back(item);
        }
    }
}
