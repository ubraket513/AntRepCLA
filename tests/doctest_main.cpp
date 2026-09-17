// doctest's implementation, compiled on its own.
//
// Kept apart from the tests so that editing a test recompiles only the test
// file. Pulling the implementation in alongside them cost ~6 s per rebuild.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
