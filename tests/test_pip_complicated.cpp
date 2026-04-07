#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

namespace {

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

void expect_equal(Location got,
                  Location expected,
                  const char* test_name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << test_name
              << ": expected " << static_cast<int>(expected)
              << ", got " << static_cast<int>(got) << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << test_name << '\n';
}

}  // namespace

int main() {
  // Visualization companion notebook for this test:
  //   tests/test_pip_complicated_visualization.nb
  //
  // That notebook visualizes the polygon and query points used here, together
  // with the query-to-perturbed-antipode arcs that correspond to the
  // implementation-generated rays used by point_in_polygon_sphere().
  // Q5 is the vertex case.
  const std::vector<std::array<double, 3>> poly = {
      {0.77114888623389370, -0.15726142646764130, 0.61692644537707060},
      {0.45249789144681710, -0.75061357063415830, 0.48148200985709080},
      {0.68946150885186746, -0.59933974587969335, 0.40673664307580021},
      {0.53398361424012150, -0.82144802877974800, 0.20021147753544170},
      {0.72547341102583852, -0.63064441484306173, 0.27563735581699919},
      {0.90662646752004000, -0.37288916572560260, 0.19743889808393390},
      {0.74736479846796566, -0.64967430761889954, 0.13917310096006544},
      {0.75468084319451650, -0.65603404827296060, -0.00872653549837396},
      {0.49138625363591330, -0.85368085756667700, -0.17253562867386300},
      {0.86555356123625300, -0.23932615843504300, -0.43993183849315200},
      {0.73819995144420940, -0.26096774566031860, -0.62205841157622660},
      {0.60166139617200880, -0.05234812405382043, -0.79703402578835670},
  };

  const std::array<double, 3> q1 = {
      0.75367527697268680, -0.65515992289232780, -0.05233595624294383};
  const std::array<double, 3> q2 = {
      0.92054211727315200, -0.38498585550407840, 0.06624274592780397};
  const std::array<double, 3> q3 = {
      0.53882393432914170, -0.82565565483991800, 0.16721694718218960};
  const std::array<double, 3> q4 = {
      0.63494819288856630, -0.65761549896072850, 0.40544130015845230};
  // Q5 is the vertex case: it is exactly P9.
  const std::array<double, 3> q5 = {
      0.49138625363591330, -0.85368085756667700, -0.17253562867386300};

  const Location loc_q1 = point_in_polygon_sphere(q1, poly);
  const Location loc_q2 = point_in_polygon_sphere(q2, poly);
  const Location loc_q3 = point_in_polygon_sphere(q3, poly);
  const Location loc_q4 = point_in_polygon_sphere(q4, poly);
  const Location loc_q5 = point_in_polygon_sphere(q5, poly);

  expect_equal(loc_q1, Location::Inside, "complicated polygon: Q1 inside");
  expect_equal(loc_q2, Location::Inside, "complicated polygon: Q2 inside");
  expect_equal(loc_q3, Location::Outside, "complicated polygon: Q3 outside");
  expect_equal(loc_q4, Location::Outside, "complicated polygon: Q4 outside");
  expect_equal(loc_q5, Location::OnVertex, "complicated polygon: Q5 on vertex");

  const std::vector<std::int64_t> vertex_ids = {
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
  };
  const std::int64_t q1_id = 200;
  const std::int64_t q2_id = 201;
  const std::int64_t q3_id = 202;
  const std::int64_t q4_id = 203;
  const std::int64_t q5_id = 204;
  const std::int64_t r1_id = 300;
  const std::int64_t r2_id = 301;
  const std::int64_t r3_id = 302;
  const std::int64_t r4_id = 303;
  const std::int64_t r5_id = 304;

  // Tier 1: full global robustness. Vertex IDs, query IDs, and the designated
  // outside points R1..R5 all participate in one explicit symbolic ordering.
  const std::array<double, 3> r1 = {
      -0.75367527657824340, 0.65515992254944344, 0.05233596621555329};
  const std::array<double, 3> r2 = {
      -0.92054211788294438, 0.38498585575910360, -0.06624273597168498};
  const std::array<double, 3> r3 = {
      -0.53882393523014660, 0.82565565622055415, -0.16721693746180466};
  const std::array<double, 3> r4 = {
      -0.63494819546290848, 0.65761550162697335, -0.40544129180227878};
  const std::array<double, 3> r5 = {
      -0.49138625278809683, 0.85368085609377320, 0.17253563837617750};
  const Location loc_q1_tier1 =
      point_in_polygon_sphere(q1, q1_id, r1, r1_id, poly, vertex_ids);
  const Location loc_q2_tier1 =
      point_in_polygon_sphere(q2, q2_id, r2, r2_id, poly, vertex_ids);
  const Location loc_q3_tier1 =
      point_in_polygon_sphere(q3, q3_id, r3, r3_id, poly, vertex_ids);
  const Location loc_q4_tier1 =
      point_in_polygon_sphere(q4, q4_id, r4, r4_id, poly, vertex_ids);
  const Location loc_q5_tier1 =
      point_in_polygon_sphere(q5, q5_id, r5, r5_id, poly, vertex_ids);
  expect_equal(loc_q1_tier1, Location::Inside,
               "complicated polygon tier1: Q1 inside");
  expect_equal(loc_q2_tier1, Location::Inside,
               "complicated polygon tier1: Q2 inside");
  expect_equal(loc_q3_tier1, Location::Outside,
               "complicated polygon tier1: Q3 outside");
  expect_equal(loc_q4_tier1, Location::Outside,
               "complicated polygon tier1: Q4 outside");
  expect_equal(loc_q5_tier1, Location::OnVertex,
               "complicated polygon tier1: Q5 on vertex");

  // Tier 2: semi-specified global robustness. Vertex IDs and query IDs are
  // explicit; the library infers R and assigns the internal R ID.
  const Location loc_q1_tier2 =
      point_in_polygon_sphere(q1, q1_id, poly, vertex_ids);
  const Location loc_q2_tier2 =
      point_in_polygon_sphere(q2, q2_id, poly, vertex_ids);
  const Location loc_q3_tier2 =
      point_in_polygon_sphere(q3, q3_id, poly, vertex_ids);
  const Location loc_q4_tier2 =
      point_in_polygon_sphere(q4, q4_id, poly, vertex_ids);
  const Location loc_q5_tier2 =
      point_in_polygon_sphere(q5, q5_id, poly, vertex_ids);
  expect_equal(loc_q1_tier2, Location::Inside,
               "complicated polygon tier2: Q1 inside");
  expect_equal(loc_q2_tier2, Location::Inside,
               "complicated polygon tier2: Q2 inside");
  expect_equal(loc_q3_tier2, Location::Outside,
               "complicated polygon tier2: Q3 outside");
  expect_equal(loc_q4_tier2, Location::Outside,
               "complicated polygon tier2: Q4 outside");
  expect_equal(loc_q5_tier2, Location::OnVertex,
               "complicated polygon tier2: Q5 on vertex");

  // Tier 3: local/internal robustness. Only the vertex IDs are explicit; the
  // library assigns internal IDs to both q and the inferred R.
  const Location loc_q1_tier3 = point_in_polygon_sphere(q1, poly, vertex_ids);
  const Location loc_q2_tier3 = point_in_polygon_sphere(q2, poly, vertex_ids);
  const Location loc_q3_tier3 = point_in_polygon_sphere(q3, poly, vertex_ids);
  const Location loc_q4_tier3 = point_in_polygon_sphere(q4, poly, vertex_ids);
  const Location loc_q5_tier3 = point_in_polygon_sphere(q5, poly, vertex_ids);
  expect_equal(loc_q1_tier3, Location::Inside,
               "complicated polygon tier3: Q1 inside");
  expect_equal(loc_q2_tier3, Location::Inside,
               "complicated polygon tier3: Q2 inside");
  expect_equal(loc_q3_tier3, Location::Outside,
               "complicated polygon tier3: Q3 outside");
  expect_equal(loc_q4_tier3, Location::Outside,
               "complicated polygon tier3: Q4 outside");
  expect_equal(loc_q5_tier3, Location::OnVertex,
               "complicated polygon tier3: Q5 on vertex");

  // Tier 4: no global IDs at all. This exercises the non-global-ID overload,
  // so SoS is unavailable and degeneracy handling must rely only on the
  // implementation's perturbed-antipode ray construction plus the half-open
  // crossing rule.
  const Location loc_q1_tier4 = point_in_polygon_sphere(q1, poly);
  const Location loc_q2_tier4 = point_in_polygon_sphere(q2, poly);
  const Location loc_q3_tier4 = point_in_polygon_sphere(q3, poly);
  const Location loc_q4_tier4 = point_in_polygon_sphere(q4, poly);
  const Location loc_q5_tier4 = point_in_polygon_sphere(q5, poly);
  expect_equal(loc_q1_tier4, Location::Inside,
               "complicated polygon tier4/no-vertex-id: Q1 inside");
  expect_equal(loc_q2_tier4, Location::Inside,
               "complicated polygon tier4/no-vertex-id: Q2 inside");
  expect_equal(loc_q3_tier4, Location::Outside,
               "complicated polygon tier4/no-vertex-id: Q3 outside");
  expect_equal(loc_q4_tier4, Location::Outside,
               "complicated polygon tier4/no-vertex-id: Q4 outside");
  expect_equal(loc_q5_tier4, Location::OnVertex,
               "complicated polygon tier4/no-vertex-id: Q5 on vertex");

  return EXIT_SUCCESS;
}
