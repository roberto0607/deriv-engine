cmake_minimum_required(VERSION 3.20)
project(deriv_engine CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.2
)
FetchContent_MakeAvailable(Catch2)

add_executable(deriv_engine src/main.cpp src/black_scholes.cpp src/greeks.cpp src/day_count.cpp src/tree_pricer.cpp src/monte_carlo.cpp src/tape.cpp src/hedging.cpp src/longstaff_schwartz.cpp src/pde_solver.cpp src/implied_vol.cpp)
target_include_directories(deriv_engine PUBLIC include)

enable_testing()
add_executable(tests tests/test_main.cpp src/black_scholes.cpp src/greeks.cpp src/day_count.cpp src/tree_pricer.cpp src/monte_carlo.cpp src/tape.cpp src/hedging.cpp src/longstaff_schwartz.cpp src/pde_solver.cpp src/implied_vol.cpp)
target_include_directories(tests PUBLIC include)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
add_test(NAME deriv_engine_tests COMMAND tests)