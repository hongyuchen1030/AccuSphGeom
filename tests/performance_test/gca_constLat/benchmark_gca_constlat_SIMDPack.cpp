#include <Eigen/Dense>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_constlat_intersection.hpp"
#include "fp64_GCAconstLat.hh"

namespace {

using accusphgeom::constructions::accux_constlat;
using accusphgeom::constructions::try_gca_constlat_intersection;
using accusphgeom::numeric::Vec3;
using accusphgeom::performance_test::fp64_gca_constlat;
using accusphgeom::performance_test::fp64_try_gca_constlat_intersection;

constexpr int kMaxVecWidth = 8;
constexpr std::size_t kDefaultDataSize = 2000000;
constexpr std::size_t kDefaultNumTests = 100;
constexpr std::size_t kDefaultNumRepeats = 7;

constexpr char kRepeatsCsvPath[] =
    "tests/performance_test/gca_constLat/output/"
    "gca_constlat_SIMDPack_timing_repeats.csv";
constexpr char kSummaryCsvPath[] =
    "tests/performance_test/gca_constLat/output/"
    "gca_constlat_SIMDPack_timing_summary.csv";

struct IntersectionBuffers {
  Eigen::VectorXd pos_x;
  Eigen::VectorXd pos_y;
  Eigen::VectorXd pos_z;
  Eigen::VectorXd neg_x;
  Eigen::VectorXd neg_y;
  Eigen::VectorXd neg_z;

  explicit IntersectionBuffers(std::size_t size)
      : pos_x(static_cast<Eigen::Index>(size)),
        pos_y(static_cast<Eigen::Index>(size)),
        pos_z(static_cast<Eigen::Index>(size)),
        neg_x(static_cast<Eigen::Index>(size)),
        neg_y(static_cast<Eigen::Index>(size)),
        neg_z(static_cast<Eigen::Index>(size)) {}

  void setZero() {
    pos_x.setZero();
    pos_y.setZero();
    pos_z.setZero();
    neg_x.setZero();
    neg_y.setZero();
    neg_z.setZero();
  }
};

struct TryBuffers {
  Eigen::VectorXd point_x;
  Eigen::VectorXd point_y;
  Eigen::VectorXd point_z;
  Eigen::VectorXd status;

  explicit TryBuffers(std::size_t size)
      : point_x(static_cast<Eigen::Index>(size)),
        point_y(static_cast<Eigen::Index>(size)),
        point_z(static_cast<Eigen::Index>(size)),
        status(static_cast<Eigen::Index>(size)) {}

  void setZero() {
    point_x.setZero();
    point_y.setZero();
    point_z.setZero();
    status.setZero();
  }
};

struct SummaryStats {
  double min_time = 0.0;
  double median_time = 0.0;
  double mean_time = 0.0;
};

struct MethodSummaryRow {
  std::string method;
  int threads_num = 0;
  int vec_width = 0;
  SummaryStats stats{};
};

void create_output_directory(const std::string& directory) {
  struct stat info;
  if (stat(directory.c_str(), &info) != 0) {
    if (mkdir(directory.c_str(), 0777) == -1) {
      std::cerr << "Error creating directory: " << directory << '\n';
    }
  }
}

std::size_t parse_arg(const char* arg, const char* name) {
  const long long value = std::atoll(arg);
  if (value <= 0) {
    throw std::invalid_argument(std::string(name) + " must be positive");
  }
  return static_cast<std::size_t>(value);
}

template <typename T>
void do_not_optimize(const T& value) {
#if defined(__GNUC__) || defined(__clang__)
  asm volatile("" : : "g"(&value) : "memory");
#else
  (void)value;
#endif
}

void consume_output(const IntersectionBuffers& out) {
  do_not_optimize(out.pos_x.data());
  do_not_optimize(out.pos_y.data());
  do_not_optimize(out.pos_z.data());
  do_not_optimize(out.neg_x.data());
  do_not_optimize(out.neg_y.data());
  do_not_optimize(out.neg_z.data());
}

void consume_output(const TryBuffers& out) {
  do_not_optimize(out.point_x.data());
  do_not_optimize(out.point_y.data());
  do_not_optimize(out.point_z.data());
  do_not_optimize(out.status.data());
}

void validate_sizes(const Eigen::MatrixXd& ptsA,
                    const Eigen::MatrixXd& ptsB,
                    const Eigen::VectorXd& latitudes,
                    const IntersectionBuffers& out) {
  if ((ptsA.cols() != 3) || (ptsB.cols() != 3) ||
      (ptsB.rows() != ptsA.rows()) || (latitudes.rows() != ptsA.rows()) ||
      (out.pos_x.rows() != ptsA.rows()) || (out.pos_y.rows() != ptsA.rows()) ||
      (out.pos_z.rows() != ptsA.rows()) || (out.neg_x.rows() != ptsA.rows()) ||
      (out.neg_y.rows() != ptsA.rows()) || (out.neg_z.rows() != ptsA.rows())) {
    throw std::runtime_error("All passed arrays must have compatible sizes");
  }
}

void validate_sizes(const Eigen::MatrixXd& ptsA,
                    const Eigen::MatrixXd& ptsB,
                    const Eigen::VectorXd& latitudes,
                    const TryBuffers& out) {
  if ((ptsA.cols() != 3) || (ptsB.cols() != 3) ||
      (ptsB.rows() != ptsA.rows()) || (latitudes.rows() != ptsA.rows()) ||
      (out.point_x.rows() != ptsA.rows()) ||
      (out.point_y.rows() != ptsA.rows()) ||
      (out.point_z.rows() != ptsA.rows()) ||
      (out.status.rows() != ptsA.rows())) {
    throw std::runtime_error("All passed arrays must have compatible sizes");
  }
}

void store_scalar(
    const accusphgeom::constructions::GcaConstLatIntersections<double>& p,
    std::size_t i,
    IntersectionBuffers& out) {
  const auto ei = static_cast<Eigen::Index>(i);

  out.pos_x(ei) = p.point_pos[0];
  out.pos_y(ei) = p.point_pos[1];
  out.pos_z(ei) = p.point_pos[2];

  out.neg_x(ei) = p.point_neg[0];
  out.neg_y(ei) = p.point_neg[1];
  out.neg_z(ei) = p.point_neg[2];
}

template <int N>
void store_pack(
    const accusphgeom::constructions::GcaConstLatIntersections<EigenPack<N>>& p,
    std::size_t i,
    IntersectionBuffers& out) {
  const auto ei = static_cast<Eigen::Index>(i);

  out.pos_x.template segment<N>(ei).array() = p.point_pos[0];
  out.pos_y.template segment<N>(ei).array() = p.point_pos[1];
  out.pos_z.template segment<N>(ei).array() = p.point_pos[2];

  out.neg_x.template segment<N>(ei).array() = p.point_neg[0];
  out.neg_y.template segment<N>(ei).array() = p.point_neg[1];
  out.neg_z.template segment<N>(ei).array() = p.point_neg[2];
}

void store_scalar(
    const accusphgeom::constructions::GcaConstLatTryResult<double>& p,
    std::size_t i,
    TryBuffers& out) {
  const auto ei = static_cast<Eigen::Index>(i);

  out.point_x(ei) = p.point[0];
  out.point_y(ei) = p.point[1];
  out.point_z(ei) = p.point[2];
  out.status(ei) = p.status;
}

template <int N>
void store_pack(
    const accusphgeom::constructions::GcaConstLatTryResult<EigenPack<N>>& p,
    std::size_t i,
    TryBuffers& out) {
  const auto ei = static_cast<Eigen::Index>(i);

  out.point_x.template segment<N>(ei).array() = p.point[0];
  out.point_y.template segment<N>(ei).array() = p.point[1];
  out.point_z.template segment<N>(ei).array() = p.point[2];
  out.status.template segment<N>(ei).array() = p.status;
}

struct PureFpKernel {
  static constexpr const char* kMethod = "pure_fp";

  static auto eval(const Vec3<double>& a, const Vec3<double>& b, double z0) {
    return fp64_gca_constlat(a, b, z0);
  }

  template <int N>
  static auto eval(const Vec3<EigenPack<N>>& a,
                   const Vec3<EigenPack<N>>& b,
                   const EigenPack<N>& z0) {
    return fp64_gca_constlat<N>(a, b, z0);
  }
};

struct AccuxKernel {
  static constexpr const char* kMethod = "accux";

  static auto eval(const Vec3<double>& a, const Vec3<double>& b, double z0) {
    return accux_constlat(a, b, z0);
  }

  template <int N>
  static auto eval(const Vec3<EigenPack<N>>& a,
                   const Vec3<EigenPack<N>>& b,
                   const EigenPack<N>& z0) {
    return accux_constlat(a, b, z0);
  }
};

struct PureFpTryKernel {
  static constexpr const char* kMethod = "pure_fp_try";

  static auto eval(const Vec3<double>& a, const Vec3<double>& b, double z0) {
    return fp64_try_gca_constlat_intersection(a, b, z0);
  }

  template <int N>
  static auto eval(const Vec3<EigenPack<N>>& a,
                   const Vec3<EigenPack<N>>& b,
                   const EigenPack<N>& z0) {
    return fp64_try_gca_constlat_intersection(a, b, z0);
  }
};

struct AccuxTryKernel {
  static constexpr const char* kMethod = "accux_try";

  static auto eval(const Vec3<double>& a, const Vec3<double>& b, double z0) {
    return try_gca_constlat_intersection(a, b, z0);
  }

  template <int N>
  static auto eval(const Vec3<EigenPack<N>>& a,
                   const Vec3<EigenPack<N>>& b,
                   const EigenPack<N>& z0) {
    return try_gca_constlat_intersection(a, b, z0);
  }
};

template <typename Kernel, typename OutBuffers>
void compute_scalar(const Eigen::MatrixXd& ptsA,
                    const Eigen::MatrixXd& ptsB,
                    const Eigen::VectorXd& latitudes,
                    OutBuffers& out) {
  validate_sizes(ptsA, ptsB, latitudes, out);

  const std::size_t size = static_cast<std::size_t>(ptsA.rows());

#pragma omp parallel for schedule(static)
  for (long long ii = 0; ii < static_cast<long long>(size); ++ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    const auto ei = static_cast<Eigen::Index>(i);

    Vec3<double> a{};
    Vec3<double> b{};

    a[0] = ptsA(ei, 0);
    a[1] = ptsA(ei, 1);
    a[2] = ptsA(ei, 2);

    b[0] = ptsB(ei, 0);
    b[1] = ptsB(ei, 1);
    b[2] = ptsB(ei, 2);

    const double z0 = latitudes(ei);
    const auto p = Kernel::eval(a, b, z0);
    store_scalar(p, i, out);
  }
}

template <int N, typename Kernel, typename OutBuffers>
void compute_pack(const Eigen::MatrixXd& ptsA,
                  const Eigen::MatrixXd& ptsB,
                  const Eigen::VectorXd& latitudes,
                  OutBuffers& out) {
  validate_sizes(ptsA, ptsB, latitudes, out);

  const std::size_t size = static_cast<std::size_t>(ptsA.rows());
  assert(size % N == 0);

  const std::size_t num_packs = size / N;

#pragma omp parallel for schedule(static)
  for (long long pp = 0; pp < static_cast<long long>(num_packs); ++pp) {
    const std::size_t i = static_cast<std::size_t>(pp) * N;
    const auto ei = static_cast<Eigen::Index>(i);

    Vec3<EigenPack<N>> a{};
    Vec3<EigenPack<N>> b{};

    a[0] = ptsA.template block<N, 1>(ei, 0).array();
    a[1] = ptsA.template block<N, 1>(ei, 1).array();
    a[2] = ptsA.template block<N, 1>(ei, 2).array();

    b[0] = ptsB.template block<N, 1>(ei, 0).array();
    b[1] = ptsB.template block<N, 1>(ei, 1).array();
    b[2] = ptsB.template block<N, 1>(ei, 2).array();

    const EigenPack<N> z0 = latitudes.template segment<N>(ei).array();
    const auto p = Kernel::template eval<N>(a, b, z0);
    store_pack<N>(p, i, out);
  }
}

template <typename ComputeFn, typename OutBuffers>
double run_benchmark_trial(ComputeFn compute,
                           std::size_t num_tests,
                           const Eigen::MatrixXd& ptsA,
                           const Eigen::MatrixXd& ptsB,
                           const Eigen::VectorXd& latitudes,
                           OutBuffers& out,
                           int threads_num) {
  omp_set_num_threads(threads_num);

  out.setZero();
  compute(ptsA, ptsB, latitudes, out);
  consume_output(out);

  const auto start = std::chrono::high_resolution_clock::now();

  for (std::size_t t = 0; t < num_tests; ++t) {
    compute(ptsA, ptsB, latitudes, out);
    consume_output(out);
  }

  const std::chrono::duration<double> elapsed =
      std::chrono::high_resolution_clock::now() - start;

  return elapsed.count();
}

double seconds_per_point(double total_seconds,
                         std::size_t num_tests,
                         std::size_t data_size) {
  return total_seconds / static_cast<double>(num_tests * data_size);
}

SummaryStats summarize_times(const std::vector<double>& times) {
  if (times.empty()) {
    throw std::runtime_error("Cannot summarize an empty timing vector");
  }

  std::vector<double> sorted = times;
  std::sort(sorted.begin(), sorted.end());

  const double min_time = sorted.front();
  const double sum = std::accumulate(sorted.begin(), sorted.end(), 0.0);
  const double mean_time = sum / static_cast<double>(sorted.size());

  double median_time = 0.0;
  const std::size_t mid = sorted.size() / 2;
  if ((sorted.size() % 2) == 0) {
    median_time = 0.5 * (sorted[mid - 1] + sorted[mid]);
  } else {
    median_time = sorted[mid];
  }

  return SummaryStats{min_time, median_time, mean_time};
}

void write_repeat_row(std::ofstream& csv,
                      const char* method,
                      int threads_num,
                      int vec_width,
                      std::size_t repeat,
                      double seconds_per_point_value) {
  csv << method << "," << threads_num << "," << vec_width << "," << repeat
      << "," << seconds_per_point_value << "\n";
}

void write_summary_row(std::ofstream& csv,
                       const MethodSummaryRow& row) {
  csv << row.method << "," << row.threads_num << "," << row.vec_width << ","
      << row.stats.min_time << "," << row.stats.median_time << ","
      << row.stats.mean_time << "\n";
}

void print_timing(const char* method,
                  int threads_num,
                  int vec_width,
                  const SummaryStats& stats) {
  std::cout << "  " << std::setw(12) << method << " threads=" << std::setw(2)
            << threads_num << " width=" << vec_width << " median="
            << std::fixed << std::setprecision(3) << stats.median_time * 1e9
            << " ns/point min=" << stats.min_time * 1e9 << " ns/point\n";
}

template <typename Kernel, typename OutBuffers>
double run_width_trial(int vec_width,
                       std::size_t num_tests,
                       const Eigen::MatrixXd& ptsA,
                       const Eigen::MatrixXd& ptsB,
                       const Eigen::VectorXd& latitudes,
                       OutBuffers& out,
                       int threads_num) {
  switch (vec_width) {
    case 1:
      return seconds_per_point(
          run_benchmark_trial(compute_scalar<Kernel, OutBuffers>, num_tests,
                              ptsA, ptsB, latitudes, out, threads_num),
          num_tests, static_cast<std::size_t>(ptsA.rows()));
    case 2:
      return seconds_per_point(
          run_benchmark_trial(compute_pack<2, Kernel, OutBuffers>, num_tests,
                              ptsA, ptsB, latitudes, out, threads_num),
          num_tests, static_cast<std::size_t>(ptsA.rows()));
    case 4:
      return seconds_per_point(
          run_benchmark_trial(compute_pack<4, Kernel, OutBuffers>, num_tests,
                              ptsA, ptsB, latitudes, out, threads_num),
          num_tests, static_cast<std::size_t>(ptsA.rows()));
    case 8:
      return seconds_per_point(
          run_benchmark_trial(compute_pack<8, Kernel, OutBuffers>, num_tests,
                              ptsA, ptsB, latitudes, out, threads_num),
          num_tests, static_cast<std::size_t>(ptsA.rows()));
    default:
      throw std::invalid_argument("Unsupported vec_width");
  }
}

template <typename Kernel, typename OutBuffers>
double run_and_record_trial(int vec_width,
                            std::size_t repeat,
                            std::size_t num_tests,
                            const Eigen::MatrixXd& ptsA,
                            const Eigen::MatrixXd& ptsB,
                            const Eigen::VectorXd& latitudes,
                            OutBuffers& out,
                            int threads_num,
                            std::ofstream& repeats_csv) {
  const double time = run_width_trial<Kernel, OutBuffers>(
      vec_width, num_tests, ptsA, ptsB, latitudes, out, threads_num);
  write_repeat_row(repeats_csv, Kernel::kMethod, threads_num, vec_width, repeat,
                   time);
  return time;
}

template <typename PureKernel, typename AccuxKernelT, typename OutBuffers>
void run_kernel_pair_block(const char* block_title,
                           std::size_t data_size,
                           std::size_t num_tests,
                           std::size_t num_repeats,
                           const Eigen::MatrixXd& ptsA,
                           const Eigen::MatrixXd& ptsB,
                           const Eigen::VectorXd& latitudes,
                           const int* requested_threads,
                           int num_thread_entries,
                           const int* vec_widths,
                           int num_vec_width_entries,
                           int max_threads,
                           std::ofstream& repeats_csv,
                           std::ofstream& summary_csv) {
  std::cout << "\n" << block_title << "\n";

  for (int ti = 0; ti < num_thread_entries; ++ti) {
    const int threads_num = requested_threads[ti];
    if (threads_num > max_threads) {
      continue;
    }

    OutBuffers out(data_size);

    for (int wi = 0; wi < num_vec_width_entries; ++wi) {
      const int vec_width = vec_widths[wi];

      std::vector<double> pure_times;
      std::vector<double> accux_times;
      pure_times.reserve(num_repeats);
      accux_times.reserve(num_repeats);

      for (std::size_t repeat = 0; repeat < num_repeats; ++repeat) {
        if ((repeat % 2) == 0) {
          pure_times.push_back(run_and_record_trial<PureKernel>(
              vec_width, repeat, num_tests, ptsA, ptsB, latitudes, out,
              threads_num, repeats_csv));
          accux_times.push_back(run_and_record_trial<AccuxKernelT>(
              vec_width, repeat, num_tests, ptsA, ptsB, latitudes, out,
              threads_num, repeats_csv));
        } else {
          accux_times.push_back(run_and_record_trial<AccuxKernelT>(
              vec_width, repeat, num_tests, ptsA, ptsB, latitudes, out,
              threads_num, repeats_csv));
          pure_times.push_back(run_and_record_trial<PureKernel>(
              vec_width, repeat, num_tests, ptsA, ptsB, latitudes, out,
              threads_num, repeats_csv));
        }
      }

      const SummaryStats pure_stats = summarize_times(pure_times);
      const SummaryStats accux_stats = summarize_times(accux_times);

      const MethodSummaryRow pure_row{
          PureKernel::kMethod, threads_num, vec_width, pure_stats};
      const MethodSummaryRow accux_row{
          AccuxKernelT::kMethod, threads_num, vec_width, accux_stats};

      write_summary_row(summary_csv, pure_row);
      write_summary_row(summary_csv, accux_row);

      print_timing(PureKernel::kMethod, threads_num, vec_width, pure_stats);
      print_timing(AccuxKernelT::kMethod, threads_num, vec_width, accux_stats);

      const double ratio = accux_stats.median_time / pure_stats.median_time;
      std::cout << "  " << AccuxKernelT::kMethod << "/"
                << PureKernel::kMethod << " ratio threads=" << threads_num
                << " width=" << vec_width << " : " << std::fixed
                << std::setprecision(6) << ratio << "\n";
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
#ifndef NDEBUG
  std::cerr << "Warning: benchmark is not compiled with NDEBUG. "
               "Use Release mode for meaningful timings.\n";
#endif

  const std::size_t data_size =
      (argc > 1) ? parse_arg(argv[1], "data_size") : kDefaultDataSize;
  const std::size_t num_tests =
      (argc > 2) ? parse_arg(argv[2], "num_tests") : kDefaultNumTests;
  const std::size_t num_repeats =
      (argc > 3) ? parse_arg(argv[3], "num_repeats") : kDefaultNumRepeats;

  if ((data_size % kMaxVecWidth) != 0) {
    throw std::invalid_argument("data_size must be divisible by 8");
  }

  const int max_threads = omp_get_max_threads();
  const int requested_threads[] = {1, 2, 4, 8, 16};
  const int vec_widths[] = {1, 2, 4};

  std::srand(12345);

  Eigen::MatrixXd ptsA(data_size, 3);
  Eigen::MatrixXd ptsB(data_size, 3);
  Eigen::VectorXd latitudes(data_size);

  const double pointA_x_significand = 5028374390644146;
  const int pointA_x_exponent = -53;
  const double pointA_y_significand = -7472957205960756;
  const int pointA_y_exponent = -53;
  const double pointA_z_significand = 6254432438282003;
  const int pointA_z_exponent = -79;

  const double pointB_x_significand = 5167685454902838;
  const int pointB_x_exponent = -53;
  const double pointB_y_significand = -7377307466399399;
  const int pointB_y_exponent = -53;
  const double pointB_z_significand = 4525606513452550;
  const int pointB_z_exponent = -78;

  const double constZ_significand = 7998403280412384;
  const int constZ_exponent = -79;

  const Eigen::Vector3d pointA(
      pointA_x_significand * std::pow(2.0, pointA_x_exponent),
      pointA_y_significand * std::pow(2.0, pointA_y_exponent),
      pointA_z_significand * std::pow(2.0, pointA_z_exponent));

  const Eigen::Vector3d pointB(
      pointB_x_significand * std::pow(2.0, pointB_x_exponent),
      pointB_y_significand * std::pow(2.0, pointB_y_exponent),
      pointB_z_significand * std::pow(2.0, pointB_z_exponent));

  const double constZ = constZ_significand * std::pow(2.0, constZ_exponent);

  for (std::size_t i = 0; i < data_size; ++i) {
    ptsA.row(static_cast<Eigen::Index>(i)) = pointA;
    ptsB.row(static_cast<Eigen::Index>(i)) = pointB;
  }

  latitudes.setConstant(constZ);

  ptsA += 1e-12 * Eigen::MatrixXd::Random(ptsA.rows(), ptsA.cols());
  ptsB += 1e-12 * Eigen::MatrixXd::Random(ptsB.rows(), ptsB.cols());
  latitudes += 1e-12 * Eigen::VectorXd::Random(latitudes.size());

  create_output_directory("tests/performance_test/gca_constLat/output");

  std::ofstream repeats_csv(kRepeatsCsvPath);
  repeats_csv << std::scientific << std::setprecision(16);
  repeats_csv << "method,threadsNum,vec_width,repeat,time\n";

  std::ofstream summary_csv(kSummaryCsvPath);
  summary_csv << std::scientific << std::setprecision(16);
  summary_csv << "method,threadsNum,vec_width,min_time,median_time,mean_time\n";

  std::cout << "gca_constlat SIMD-pack timing benchmark\n";
  std::cout << "  data_size=" << data_size << ", num_tests=" << num_tests
            << ", num_repeats=" << num_repeats << "\n";
  std::cout << "  OpenMP max_threads=" << max_threads << "\n";

  run_kernel_pair_block<PureFpKernel, AccuxKernel, IntersectionBuffers>(
      "accux_constlat full-intersections timing benchmark", data_size,
      num_tests, num_repeats, ptsA, ptsB, latitudes, requested_threads,
      static_cast<int>(std::size(requested_threads)), vec_widths,
      static_cast<int>(std::size(vec_widths)), max_threads, repeats_csv,
      summary_csv);

  run_kernel_pair_block<PureFpTryKernel, AccuxTryKernel, TryBuffers>(
      "try_gca_constlat_intersection full-API timing benchmark", data_size,
      num_tests, num_repeats, ptsA, ptsB, latitudes, requested_threads,
      static_cast<int>(std::size(requested_threads)), vec_widths,
      static_cast<int>(std::size(vec_widths)), max_threads, repeats_csv,
      summary_csv);

  std::cout << "  repeats CSV: " << kRepeatsCsvPath << "\n";
  std::cout << "  summary CSV: " << kSummaryCsvPath << "\n";

  return EXIT_SUCCESS;
}
