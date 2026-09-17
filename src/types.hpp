// Core data types shared across the pipeline.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "ankerl/unordered_dense.h"

namespace ar {

using GeneId = std::uint32_t;
using NodeId = std::uint32_t;

// Gene names interned to integers.
//
// The edge predicate asks whether two records share a V gene and a J gene. As
// strings that is repeated hashing and comparison of text; as sorted integer
// vectors it is a linear merge over a handful of elements. Names are kept so
// the usage table can be written back out.
class GeneTable {
public:
    GeneId intern(const std::string& name) {
        auto found = index_.find(name);
        if (found != index_.end()) return found->second;
        const auto id = static_cast<GeneId>(names_.size());
        names_.push_back(name);
        index_.emplace(name, id);
        return id;
    }

    const std::string& name(GeneId id) const { return names_[id]; }
    std::size_t size() const { return names_.size(); }

private:
    std::vector<std::string> names_;
    ankerl::unordered_dense::map<std::string, GeneId> index_;
};

// One sequencing read, normalised: allele suffixes stripped, gene calls split.
struct Read {
    std::string cdr3;
    std::string cdr3_aa;
    int cdr_length = 0;
    std::vector<GeneId> v_call;  // sorted, deduplicated
    std::vector<GeneId> j_call;  // sorted, deduplicated
};

struct Repertoire {
    std::vector<Read> reads;
    GeneTable genes;
};

// A distinct CDR3 sequence, carrying the union of gene calls across every read
// that reported it.
struct UniqueCdr3 {
    std::string cdr3;
    int cdr_length = 0;
    std::vector<GeneId> v_call;
    std::vector<GeneId> j_call;
};

}  // namespace ar
