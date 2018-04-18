#ifndef FAST_PTR_HASH_H
#define FAST_PTR_HASH_H

// std::unordered_map's default hashing function is slow.
template<typename Tval>
struct FashPtrHash {
  size_t operator()(const Tval* val) const {
    static const size_t shift = (size_t)log2(1 + sizeof(Tval));
    return (size_t)(val) >> shift;
  }
};

#endif