/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  2008-09-22 (0.1.0):

	* The initial version.

*/

#ifndef AC_KVEC_H
#define AC_KVEC_H

#include <stdlib.h>
#include "kalloc.h"

#define kv_roundup32(x)                                                                            \
  (--(x),                                                                                          \
   (x) |= (x) >> 1,                                                                                \
   (x) |= (x) >> 2,                                                                                \
   (x) |= (x) >> 4,                                                                                \
   (x) |= (x) >> 8,                                                                                \
   (x) |= (x) >> 16,                                                                               \
   ++(x))

#define kvec_t(type)                                                                               \
  struct                                                                                           \
  {                                                                                                \
    size_t n, m;                                                                                   \
    type* a;                                                                                       \
  }
#define kv_init(v) ((v).n = (v).m = 0, (v).a = 0)
#define kv_destroy(v) free((v).a)
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)

#define kv_resize(type, km, v, s)                                                                  \
  do {                                                                                             \
    if ((v).m < (s)) {                                                                             \
      (v).m = (s);                                                                                 \
      kv_roundup32((v).m);                                                                         \
      (v).a = (type*)krealloc((km), (v).a, sizeof(type) * (v).m);                                  \
    }                                                                                              \
  } while (0)

#define kv_copy(type, km, v1, v0)                                                                  \
  do {                                                                                             \
    if ((v1).m < (v0).n) kv_resize(type, (km), (v1), (v0).n);                                      \
    (v1).n = (v0).n;                                                                               \
    memcpy((v1).a, (v0).a, sizeof(type) * (v0).n);                                                 \
  } while (0)

#define kv_push(type, km, v, x)                                                                    \
  do {                                                                                             \
    if ((v).n == (v).m) {                                                                          \
      (v).m = (v).m ? (v).m << 1 : 2;                                                              \
      (v).a = (type*)krealloc((km), (v).a, sizeof(type) * (v).m);                                  \
    }                                                                                              \
    (v).a[(v).n++] = (x);                                                                          \
  } while (0)

#define kv_pushp(type, km, v, p)                                                                   \
  do {                                                                                             \
    if ((v).n == (v).m) {                                                                          \
      (v).m = (v).m ? (v).m << 1 : 2;                                                              \
      (v).a = (type*)krealloc((km), (v).a, sizeof(type) * (v).m);                                  \
    }                                                                                              \
    *(p) = &(v).a[(v).n++];                                                                        \
  } while (0)

#define kv_reverse(type, v, start)                                                                 \
  do {                                                                                             \
    if ((v).m > 0 && (v).n > (start)) {                                                            \
      size_t __i, __end = (v).n - (start);                                                         \
      type* __a = (v).a + (start);                                                                 \
      for (__i = 0; __i < __end >> 1; ++__i) {                                                     \
        type __t = __a[__end - 1 - __i];                                                           \
        __a[__end - 1 - __i] = __a[__i];                                                           \
        __a[__i] = __t;                                                                            \
      }                                                                                            \
    }                                                                                              \
  } while (0)

#endif
