#include <gtest/gtest.h>
#include "math/math.hpp"

TEST(Math, Add) {
    EXPECT_EQ(mathx::add(2,2), 4);
}
TEST(Math, Mul) {
    EXPECT_EQ(mathx::mul(3,5), 15);
}
