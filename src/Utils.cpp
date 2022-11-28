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

    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L)
    {
        // compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range
        // for each number of mutations
        std::vector<unsigned int> nMutRange(maxRange + 1);
        nMutRange[0] = 1;
        for (unsigned int i = 1; i <= maxRange; ++i)
        {
            nMutRange[i] = nMutRange[i - 1] + utils::nChoosek(L, i);
        }
        return nMutRange;
    }

    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q)
    {
        // compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range
        // for each number of mutations
        std::vector<unsigned int> nMutRange(maxRange + 1);
        nMutRange[0] = 1;
        for (unsigned int i = 1; i <= maxRange; ++i)
        {
            nMutRange[i] = nMutRange[i - 1] + (utils::nChoosek(L, i) * std::pow(q, i));
        }
        return nMutRange;
    }

    mutatedPositions specIdxToMutPos(const utils::id specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int>& nMutRange)
    {
        mutatedPositions mutPos{};
        // check if Id is valid
        if (specIdx <= nMutRange.back() && specIdx > 0)
        {
            // number of mutations in sequence
            unsigned int numMut = std::lower_bound(nMutRange.begin(), nMutRange.end(), specIdx) - nMutRange.begin();
            // id 1 = no mutations: don't start to compute anything
            if (numMut != 0)
            {
                unsigned int Lact = L;
                unsigned int numMutAct = numMut;
                // get the id within the range of number of mutations (substract the ids for the sequences with less
                // mutations)
                utils::id idAct = specIdx - nMutRange[numMut - 1];
                // determine the each mutations position seen from the mutations position before...
                unsigned int m = 0;
                do
                {
                    // for each possible positions within the length the actual mutation covers a range of ids
                    // depending on the residual mutations to follow

                    std::vector<long> cumSumRange(Lact - (numMutAct - 1));
                    auto n = Lact - 1;
                    for (unsigned int i = 0; i < cumSumRange.size(); ++i)
                    {
                        // the id range for each position when it is mutated (if pos 1 is mutated L-1 positions
                        // remain to
                        //  contain the mut-1 remaining mutations, each time there are q^nmut possibilies (it does
                        //  not change in every step)
                        cumSumRange[i] = utils::nChoosek(n, numMutAct - 1) * std::pow(numSymbols, numMut);
                        --n;
                    }

                    std::partial_sum(cumSumRange.begin(), cumSumRange.end(), cumSumRange.begin());

                    // find the id within the ranges and get the index (=position seen from the last mutated
                    // position)
                    int pos = std::lower_bound(cumSumRange.begin(), cumSumRange.end(), idAct) - cumSumRange.begin() + 1;
                    // if one mutation, the symbol can be simply find by modulo q, with 2 mutations it has to be
                    // divided by q first (since one symbol of the first mutation can be present with q symbols of
                    // the second one, with 3 mutation divided by q^2 etc
                    int mut = (int((idAct - 1) / std::pow(numSymbols, numMutAct - 1)) % numSymbols) + 1;
                    // put the position - symbol pair into the map. add the last cummulative position (=position
                    // seen from the beginning of the sequence)
                    unsigned int prePos = 0;
                    if (mutPos.begin() != mutPos.end())
                        prePos = mutPos.rbegin()->first;
                    mutPos.insert({pos + prePos, mut});
                    //  the residual length after the actual mutation
                    Lact = Lact - pos;
                    // the redsiudal number of mutations after the actual mutation
                    --numMutAct;
                    // the id within the residual length (-1 because the inidices start at 0, and another -1 because
                    // we substract the ids of the preceding mutation range changees for q symbols: the substracted
                    // ids have to be multiplied by q (since each position can have 3 values)
                    idAct = idAct - (pos == 1 ? 0 : cumSumRange[pos - 2]);
                    ++m;
                } while (m < numMut);
            }
        }
        else
        {
            std::cerr << "ERROR: id is too big for length of sequence and number of mutations" << std::endl;
        }
        return mutPos;
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
