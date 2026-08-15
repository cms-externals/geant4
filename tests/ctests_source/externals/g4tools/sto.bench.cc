// sto.bench.cc

#include "tools/sto"

#include <benchmark/benchmark.h>

#include <cmath>
#include <iomanip>
#include <random>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//
// BENCHMARK FIXTURES
//---------------------------------------------------------------------------//
namespace G4Bench
{
/*!
 * \brief Selects n elements from a source vector, cycling through if needed.
 *
 * This function takes \p n elements from the \p from vector, cycling through
 * the input if \p n is greater than the size of \p from. The result is a vector
 * of size \p n, with elements selected in order from \p from, wrapping around as necessary.
 *
 * \param n Number of elements to take.
 * \param from Source vector to take elements from.
 * \return std::vector<std::string> Vector of n elements taken from \p from.
 */
std::vector<std::string> TakeFrom(size_t n, const std::vector<std::string> from)
{
  std::vector<std::string> result;
  result.reserve(n);
  for (size_t i = 0; i < n; ++i)
  {
    result.push_back(from[n % from.size()]);
  }
  return result;
}

/*!
 * \brief Join a vector of strings with a separator.
 *
 * Concatenates the elements of the input vector into a single string, inserting
 * the specified separator between each element. Useful for generating delimited
 * strings such as CSV rows or lists.
 *
 * \param input The vector of strings to join.
 * \param sep The separator string to insert between elements.
 * \return The joined string with separators.
 */
std::string Join(const std::vector<std::string>& input, const std::string& sep)
{
  std::string output;
  bool first = true;
  for (const auto& v : input)
  {
    if (!first)
    {
      output += sep;
    }
    output += v;
  }
  return output;
}

/*!
 * \brief Generate a set of varied string inputs for type T.
 *
 * This template function generates a vector of string representations suitable for parsing
 * into the type T. The generated strings are designed to provide a variety of cases for
 * benchmarking string-to-object conversion routines, and to prevent compiler optimizations
 * such as constant folding and branch prediction.
 *
 * Specialized versions exist for bool, int, unsigned, long, float, double, and long double.
 *
 * \tparam T The type for which to generate string inputs.
 * \param count The number of string inputs to generate.
 * \return Vector of string inputs for type T.
 */
template<typename T>
std::vector<std::string> GenerateStoInputs(size_t count);

template<>
std::vector<std::string> GenerateStoInputs<bool>(size_t count)
{
  std::vector<std::string> inputs = {"true", "false", "1",     "0",   "yes", "no", "on",
                                     "off",  "TRUE",  "FALSE", "Yes", "No",  "ON", "OFF"};
  return TakeFrom(count, inputs);
}

/*!
 * \brief Generate a set of varied string inputs for an integral type.
 *
 * This function generates a vector of string representations of random integral values
 * in the range [min, max], suitable for benchmarking string-to-integer conversion routines.
 *
 * \tparam T The integral type for which to generate string inputs.
 * \param count The number of string inputs to generate.
 * \param min The minimum value (inclusive) for generated numbers.
 * \param max The maximum value (inclusive) for generated numbers.
 * \return Vector of string inputs for type T.
 */
template<typename T>
std::vector<std::string> GenerateStoInputs_Integral(size_t count, T min, T max)
{
  static_assert(std::is_integral<T>::value, "GenerateIntegralTestInputs requires an integral type");
  std::vector<std::string> result;
  result.reserve(count);
  std::mt19937 gen(42);
  std::uniform_int_distribution<T> dist(min, max);
  for (size_t i = 0; i < count; ++i)
  {
    result.push_back(std::to_string(dist(gen)));
  }
  return result;
}

template<>
std::vector<std::string> GenerateStoInputs<int>(size_t count)
{
  return GenerateStoInputs_Integral(count, -100000, 100000);
}

template<>
std::vector<std::string> GenerateStoInputs<unsigned>(size_t count)
{
  return GenerateStoInputs_Integral(count, 0u, 200000u);
}

template<>
std::vector<std::string> GenerateStoInputs<long>(size_t count)
{
  return GenerateStoInputs_Integral(count, -1000000l, 1000000l);
}

/*!
 * \brief Generate a set of varied string inputs for float.
 *
 * This specialization generates a vector of string representations of random float values
 * suitable for benchmarking string-to-float conversion routines.
 *
 * The algorithm produces a realistic range of string representations by:
 *  - Randomizing the magnitude (from small decimals to large values) to vary the number of digits.
 *  - Randomly choosing between default, fixed-point, and scientific notation for each value.
 *  - Randomizing the sign and decimal precision to simulate real-world float input strings.
 *  - Ensuring the output covers both short and long strings, as well as different exponent formats.
 *
 * \param count The number of string inputs to generate.
 * \return Vector of string inputs for float.
 */
template<>
std::vector<std::string> GenerateStoInputs<float>(size_t count)
{
  std::vector<std::string> result;
  result.reserve(count);
  std::mt19937 gen(42);

  // Mix of different magnitudes to get varied string lengths
  std::uniform_int_distribution<int> magnitude_dist(0, 7);
  std::uniform_real_distribution<float> base_dist(0.1f, 9.9f);
  std::uniform_int_distribution<int> sign_dist(0, 1);
  std::uniform_int_distribution<int> format_dist(0, 2);

  for (size_t i = 0; i < count; ++i)
  {
    float base = base_dist(gen);
    int magnitude = magnitude_dist(gen);
    float value = base * std::pow(10.0f, magnitude - 3);  // Range: ~0.0001 to ~9900
    if (sign_dist(gen)) value = -value;

    // Vary the format: default, fixed precision, or scientific
    int format = format_dist(gen);
    if (format == 0)
    {
      // Default std::to_string
      result.push_back(std::to_string(value));
    }
    else if (format == 1)
    {
      // Fixed precision (short strings)
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(2) << value;
      result.push_back(oss.str());
    }
    else
    {
      // Scientific notation with varied exponents
      std::ostringstream oss;
      oss << std::scientific << std::setprecision(6) << value;
      result.push_back(oss.str());
    }
  }
  return result;
}

/*!
 * \brief Generate a set of varied string inputs for double.
 *
 * This specialization generates a vector of string representations of random double values
 * suitable for benchmarking string-to-double conversion routines.
 *
 * The algorithm is designed to provide a realistic and challenging set of string representations
 * by:
 *  - Covering the full range of double exponents (from -308 to +308) for scientific notation.
 *  - Generating values with many decimal places, large integers, and high-precision fixed-point
 * numbers.
 *  - Randomly selecting between fixed, scientific, and high-precision formats for each value.
 *  - Including both positive and negative values, and varying the number of digits and decimal
 * places.
 *  - Producing strings that mimic real-world double input, including edge cases for parsing
 * routines.
 *
 * \param count The number of string inputs to generate.
 * \return Vector of string inputs for double.
 */
template<>
std::vector<std::string> GenerateStoInputs<double>(size_t count)
{
  std::vector<std::string> result;
  result.reserve(count);
  std::mt19937 gen(42);

  // Much wider range for doubles to test longer strings
  std::uniform_int_distribution<int> magnitude_dist(-308, 308);  // Full double range
  std::uniform_real_distribution<double> base_dist(1.0, 9.999999999);
  std::uniform_int_distribution<int> sign_dist(0, 1);
  std::uniform_int_distribution<int> format_dist(0, 3);

  for (size_t i = 0; i < count; ++i)
  {
    int format = format_dist(gen);
    double value;

    if (format == 0)
    {
      // Small values with many decimal places (long strings)
      // e.g., "0.123456789012345"
      std::uniform_real_distribution<double> small_dist(0.0, 1.0);
      value = small_dist(gen);
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(15) << value;
      result.push_back(oss.str());
    }
    else if (format == 1)
    {
      // Large integers (long strings without decimals)
      // e.g., "123456789012345.000000"
      std::uniform_real_distribution<double> large_dist(1e10, 1e15);
      value = large_dist(gen);
      result.push_back(std::to_string(value));
    }
    else if (format == 2)
    {
      // Scientific notation with large exponents
      // e.g., "1.234567e+123" or "9.876543e-234"
      value = base_dist(gen);
      int exponent = magnitude_dist(gen);
      if (sign_dist(gen)) value = -value;
      std::ostringstream oss;
      oss << std::scientific << std::setprecision(12) << (value * std::pow(10.0, exponent));
      result.push_back(oss.str());
    }
    else
    {
      // Very high precision fixed format
      // e.g., "3.141592653589793238"
      std::uniform_real_distribution<double> med_dist(-1000.0, 1000.0);
      value = med_dist(gen);
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(18) << value;
      result.push_back(oss.str());
    }
  }
  return result;
}

template<>
std::vector<std::string> GenerateStoInputs<long double>(size_t count)
{
  return GenerateStoInputs<double>(count);
}

/*!
 * \brief Generate a set of invalid string inputs for testing.
 *
 * This function returns a vector of strings that are intentionally invalid for
 * numeric or boolean parsing. Useful for benchmarking error handling and robustness
 * of string-to-object conversion routines.
 *
 * \param count Number of invalid string inputs to generate.
 * \return Vector of invalid string inputs.
 */
std::vector<std::string> GenerateStoInputs_Invalid(size_t count)
{
  std::vector<std::string> inputs = {"not_a_number", "invalid", "abc123", "12.34.56",
                                     "hello",        "world",   "test",   "data"};
  return TakeFrom(count, inputs);
}

/*!
 * \brief Generate a set of special floating-point string inputs.
 *
 * This function returns a vector of strings representing special floating-point values
 * such as NaN, infinity, and common mathematical constants. Useful for testing the
 * handling of special cases in string-to-floating-point conversions.
 *
 * \param count Number of special floating-point string inputs to generate.
 * \return Vector of special floating-point string inputs.
 */
std::vector<std::string> GenerateStoInputs_SpecialFloat(size_t count)
{
  std::vector<std::string> inputs = {"nan",  "NaN",     "NAN",     "inf",     "Inf", "INF", "+inf",
                                     "-inf", "3.14159", "2.71828", "1.41421", "0.0", "-0.0"};
  return TakeFrom(count, inputs);
}

// Generate a vector of csv element strings (basically a column)
template<typename T>

/*!
 * \brief Generate a vector of delimited string rows (e.g., CSV column).
 *
 * For each row, generates a random number of elements (up to max), creates string
 * representations using GenerateStoInputs<T>, and joins them with the given separator.
 * This is useful for benchmarking parsing of delimited columns with varied element counts.
 *
 * \tparam T The type of elements to generate for each row.
 * \param rows Number of rows to generate.
 * \param max Maximum number of elements per row.
 * \param sep Separator to use between elements in each row.
 * \return Vector of delimited string rows.
 */
std::vector<std::string> GenerateColumn(size_t rows, size_t max, const std::string& sep)
{
  std::vector<std::string> output{rows, ""};
  std::mt19937 gen(42);
  std::uniform_int_distribution<size_t> elem_dist(0, max);

  for (auto& row : output)
  {
    auto data = GenerateStoInputs<T>(elem_dist(gen));
    row = Join(data, sep);
  }

  return output;
}

}  // namespace G4Bench
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// SNUMS BENCHMARKS
//---------------------------------------------------------------------------//
// Benchmark: multiple different invalid inputs
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_Invalid(benchmark::State& state)
{
  auto inputs = G4Bench::GenerateStoInputs_Invalid(1000);
  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_Invalid, bool);
BENCHMARK_TEMPLATE(BM_To_Invalid, int);
BENCHMARK_TEMPLATE(BM_To_Invalid, float);
BENCHMARK_TEMPLATE(BM_To_Invalid, double);

//---------------------------------------------------------------------------//
// Benchmark: partial parsing, e.g. "123abc" for integer
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_PartialParse(benchmark::State& state)
{
  auto inputs = G4Bench::GenerateStoInputs<T>(1000);
  for (auto& s : inputs)
  {
    s += "bleugh";  // append junk to trigger partial parse;
  }

  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_PartialParse, int);
BENCHMARK_TEMPLATE(BM_To_PartialParse, double);

//---------------------------------------------------------------------------//
// Benchmark: empty string
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_EmptyString(benchmark::State& state)
{
  T value;
  for (auto _ : state)
  {
    tools::to("", value);
    benchmark::DoNotOptimize(value);
  }
}
BENCHMARK_TEMPLATE(BM_To_EmptyString, int);
BENCHMARK_TEMPLATE(BM_To_EmptyString, double);

//---------------------------------------------------------------------------//
// Benchmark: numerical input from 2-10 digits
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_VaryingLength(benchmark::State& state)
{
  std::vector<std::string> inputs;
  std::mt19937 gen(42);

  if constexpr (std::is_integral<T>::value)
  {
    for (size_t i = 0; i < 1000; ++i)
    {
      int digits = 1 + (gen() % state.range(0));
      std::string s;
      for (int d = 0; d < digits; ++d)
      {
        s += '0' + (gen() % 10);
      }
      inputs.push_back(s);
    }
  }
  else if constexpr (std::is_floating_point<T>::value)
  {
    std::uniform_real_distribution<T> dist(-1000.0, 1000.0);
    for (size_t i = 0; i < 1000; ++i)
    {
      T val = dist(gen);
      int precision = 1 + (gen() % static_cast<int>(state.range(0)));
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(precision) << val;
      inputs.push_back(oss.str());
    }
  }

  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_VaryingLength, int)->Arg(2)->Arg(4)->Arg(8)->Arg(10);
BENCHMARK_TEMPLATE(BM_To_VaryingLength, double)->Arg(2)->Arg(4)->Arg(8)->Arg(10);

//---------------------------------------------------------------------------//
// Benchmark: basic varied input of given type
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_Varied(benchmark::State& state)
{
  auto inputs = G4Bench::GenerateStoInputs<T>(1000);
  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    benchmark::ClobberMemory();
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_Varied, bool);
BENCHMARK_TEMPLATE(BM_To_Varied, int);
BENCHMARK_TEMPLATE(BM_To_Varied, unsigned);
BENCHMARK_TEMPLATE(BM_To_Varied, long);
BENCHMARK_TEMPLATE(BM_To_Varied, float);
BENCHMARK_TEMPLATE(BM_To_Varied, double);
BENCHMARK_TEMPLATE(BM_To_Varied, long double);

//---------------------------------------------------------------------------//
// Benchmark: special FP values like "nan", "inf"
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_FPSpecialValues(benchmark::State& state)
{
  auto inputs = G4Bench::GenerateStoInputs_SpecialFloat(1000);
  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    benchmark::ClobberMemory();
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_FPSpecialValues, float);
BENCHMARK_TEMPLATE(BM_To_FPSpecialValues, double);

//---------------------------------------------------------------------------//
// Benchmark: Floating point input in 1.23E+12 format
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_FPScientific(benchmark::State& state)
{
  std::vector<std::string> inputs;
  std::mt19937 gen(42);
  std::uniform_real_distribution<T> mantissa(1.0, 9.999);
  std::uniform_int_distribution<int> exponent(-30, 30);
  for (size_t i = 0; i < 1000; ++i)
  {
    inputs.push_back(std::to_string(mantissa(gen)) + "e" + std::to_string(exponent(gen)));
  }
  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    benchmark::ClobberMemory();
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_FPScientific, double);
BENCHMARK_TEMPLATE(BM_To_FPScientific, double);

//---------------------------------------------------------------------------//
// Benchmark: realistic mixture of floating point inputs
//---------------------------------------------------------------------------//
template<typename T>
static void BM_To_FPRealisticMix(benchmark::State& state)
{
  // Mix of normal values (80%), special values (10%), and invalid (10%)
  auto normal = G4Bench::GenerateStoInputs<T>(800);
  auto special = G4Bench::GenerateStoInputs_SpecialFloat(100);
  auto invalid = G4Bench::GenerateStoInputs_Invalid(100);

  std::vector<std::string> inputs;
  inputs.insert(inputs.end(), normal.begin(), normal.end());
  inputs.insert(inputs.end(), special.begin(), special.end());
  inputs.insert(inputs.end(), invalid.begin(), invalid.end());

  // Shuffle to mix them up
  std::mt19937 gen(42);
  std::shuffle(inputs.begin(), inputs.end(), gen);

  T value;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(inputs[index % inputs.size()], value);
    benchmark::DoNotOptimize(value);
    benchmark::ClobberMemory();
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_To_FPRealisticMix, float);
BENCHMARK_TEMPLATE(BM_To_FPRealisticMix, double);

//---------------------------------------------------------------------------//
// SNUMS BENCHMARKS
//---------------------------------------------------------------------------//
// Benchmark: basic parse of cell of varying number of elements
//---------------------------------------------------------------------------//
template<typename T>
static void BM_Snums_VaryCellLength(benchmark::State& state)
{
  size_t vec_len = state.range(0);
  auto vec_input = G4Bench::GenerateStoInputs<T>(vec_len);
  std::string input = G4Bench::Join(vec_input, ";");

  std::vector<std::string> tmp;
  std::vector<T> values;
  for (auto _ : state)
  {
    tools::to(input, values, ";");
    benchmark::DoNotOptimize(values);
    benchmark::ClobberMemory();
  }
}
BENCHMARK_TEMPLATE(BM_Snums_VaryCellLength, int)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);
BENCHMARK_TEMPLATE(BM_Snums_VaryCellLength, float)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);
BENCHMARK_TEMPLATE(BM_Snums_VaryCellLength, double)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);

//---------------------------------------------------------------------------//
// Benchmark: 1000 cells each with range of vector lengths
//---------------------------------------------------------------------------//
template<typename T>
static void BM_Snums_RealisticCellMix(benchmark::State& state)
{
  size_t elems = state.range(0);
  auto column = G4Bench::GenerateColumn<T>(1000, elems, ";");

  std::vector<std::string> tmp;
  std::vector<T> values;
  size_t index = 0;
  for (auto _ : state)
  {
    tools::to(column[index % column.size()], values, ";");
    benchmark::DoNotOptimize(values);
    benchmark::ClobberMemory();
    ++index;
  }
}
BENCHMARK_TEMPLATE(BM_Snums_RealisticCellMix, int)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);
BENCHMARK_TEMPLATE(BM_Snums_RealisticCellMix, float)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);
BENCHMARK_TEMPLATE(BM_Snums_RealisticCellMix, double)->Arg(2)->Arg(10)->Arg(100)->Arg(1000);

// ============================================================================
// Main function
// ============================================================================
BENCHMARK_MAIN();
