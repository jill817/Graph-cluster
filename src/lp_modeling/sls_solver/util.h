#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#define __UNREACHABLE__
#define __ASSERT__
// #define __TRACE__

#ifdef __UNREACHABLE__
#define UNREACHABLE() std::cout << "UNREACHABLE, line: " << __LINE__ << std::endl
#endif

#ifdef __ASSERT__
#define ASSERT(T) if(!(T)) UNREACHABLE()
#define ASSERT_CODE(T, CODE) if(!(T)) { CODE } ((void) 0)
#endif

#ifdef __TRACE__
#define TRACE(CODE) { CODE } ((void) 0)
#else
#define TRACE(CODE) ((void) 0)
#endif

using bool_vector = std::vector<bool>;
using int_vector = std::vector<int>;
using double_vector = std::vector<double>;
using long_double_vector = std::vector<long double>;
using int_table = std::unordered_set<int>;
using int_table_vector = std::vector<int_table>;
using bool_pair = std::pair<bool, bool>;
using int_pair = std::pair<int, int>;
using long_double_pair = std::pair<long double, long double>;

struct pair_hasher {
    template<typename T, typename U>
    size_t operator()(std::pair<T, U> const & p) const {
        return std::hash<T>(p.first) ^ std::hash<U>(p.second);
    }
};

struct pair_comparator {
    template<typename T, typename U>
    bool operator()(std::pair<T, U> const & p1, std::pair<T, U> const & p2) const {
        return (p1.first == p2.first && p1.second == p2.second) 
            || (p1.first == p2.second && p1.second == p2.first);
    }
};

using int_pair_table = std::unordered_set<int_pair, pair_hasher, pair_comparator>;
using int_pair_vector = std::vector<int_pair>;