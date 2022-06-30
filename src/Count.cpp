//
// Created by msmith on 15.01.20.
//

#include "Count.hpp"

#include <algorithm>
#include <fstream>

namespace count
{

    void counter_1::increment(const unsigned pos, const unsigned symbol)
    {
        ++counter_1::data[pos - 1][symbol - 1];
    }

    void counter_1::decrement(const unsigned pos, const unsigned symbol)
    {
        --counter_1::data[pos - 1][symbol - 1];
    }

    void counter_1::count(const unsigned pos, const unsigned symbol, const unsigned times)
    {
        counter_1::data[pos - 1][symbol - 1] += times;
    }

    const count_type counter_1::get(const unsigned pos, const unsigned symbol)
    {
        return counter_1::data[pos - 1][symbol - 1];
    }

    void counter_1::write_to_file(const std::string& out_file, const std::string& header)
    {
        std::ofstream outfile(out_file);

        if (outfile.good())
        {
            outfile << header;
            for (unsigned i = 0; i < data.size(); ++i)
            {
                // print position
                outfile << (i + 1);
                // print counts for all possible symbols
                std::for_each(data[i].cbegin(), data[i].cend(),
                              [&outfile](const auto& entry) { outfile << '\t' << entry; });
                outfile << '\n';
            }
        }
    }

    void counter_2::increment(const unsigned pos1, const unsigned pos2, const unsigned symbol1, const unsigned symbol2)
    {
        ++counter_2::data[counter_2::getPositionIndex(pos1, pos2)][counter_2::getSymbolIndex(symbol1, symbol2)];
    }

    void counter_2::decrement(const unsigned pos1, const unsigned pos2, const unsigned symbol1, const unsigned symbol2)
    {
        --counter_2::data[counter_2::getPositionIndex(pos1, pos2)][counter_2::getSymbolIndex(symbol1, symbol2)];
    }

    void counter_2::count(const unsigned int pos1, const unsigned int pos2, const unsigned int symbol1,
                          const unsigned int symbol2, const int times)
    {
        counter_2::data[counter_2::getPositionIndex(pos1, pos2)][counter_2::getSymbolIndex(symbol1, symbol2)] += times;
    }

    void counter_2::count(const unsigned pos, const unsigned symbol, const int times)
    {
        for (int i = 1; i <= L; ++i)
        {
            if (i < pos)
            {
                count(i, pos, 1, symbol, times);
            }
            else if (i > pos)
            {
                count(pos, i, symbol, 1, times);
            }
        }
    }

    unsigned counter_2::getPositionIndex(const unsigned pos1, const unsigned pos2)
    {
        return double(pos1 - 1) * (double(L) - double(pos1) / 2) + pos2 - pos1 - 1;
        // weg
        // return (pos1-1)*(L-(pos1-1))+pos2-pos1-1;
    }

    unsigned counter_2::getSymbolIndex(const unsigned symbol1, const unsigned symbol2)
    {
        return (symbol1 - 1) * q + symbol2 - 1;
    }

    void counter_2::write_to_file(const std::string& out_file, const std::string& header)
    {
        std::ofstream outfile(out_file);

        if (outfile.good())
        {
            outfile << header;

            unsigned line = 0;
            unsigned i = 1;
            unsigned j = 2;
            while (line < data.size())
            {
                outfile << i << '\t' << j;
                std::for_each(data[line].cbegin(), data[line].cend(),
                              [&outfile](const auto& entry) { outfile << '\t' << entry; });
                outfile << '\n';

                ++line;
                if (j == L)
                {
                    ++i;
                    j = i + 1;
                }
                else
                {
                    ++j;
                }
            }
        }
    }

}