#define CATCH_CONFIG_RUNNER
#include "TestsHeader.hpp"

int main(int argc, char const *argv[])
{
    int result = Catch::Session().run(argc, argv);

    return result;
}
