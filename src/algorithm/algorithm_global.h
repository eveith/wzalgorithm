#ifndef ALGORITHM_GLOBAL_H
#define ALGORITHM_GLOBAL_H

#if defined(ALGORITHM_LIBRARY)
#  define ALGORTHMSHARED_EXPORT __attribute__((visibility("default")))
#else
#  define ALGORTHMSHARED_EXPORT __attribute__((visibility("default")))
#endif

#endif // ALGORITHM_GLOBAL_H
