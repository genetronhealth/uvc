#ifndef Hashing_hpp_INLCUDED
#define Hashing_hpp_INLCUDED

#include "common.hpp"

template <class T>
inline
uvc1_hash_t
strnhash(const T *str, size_t n, const uvc1_hash_t base = 31UL) {
    uvc1_hash_t  ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * base + ((uvc1_hash_t)str[i]);
    }
    return ret;
}

template <class T>
inline
uvc1_hash_t
strnhash_rc(const T *str, size_t n, const uvc1_hash_t base = 31UL) {
    uvc1_hash_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * base + STATIC_REV_COMPLEMENT.data[((uvc1_hash_t)str[n-i-(size_t)1])];
    }
    return ret;
}

template<class T>
inline
uvc1_hash_t
strhash(const T *str, const uvc1_hash_t base = 31UL) {
    return strnhash(str, SIZE_MAX, base);
}

inline
uvc1_hash_t
hash2hash(uvc1_hash_t hash1, uvc1_hash_t hash2) {
    return hash1 * ((1UL<<(31UL)) - 1UL) + hash2;
}

#endif
