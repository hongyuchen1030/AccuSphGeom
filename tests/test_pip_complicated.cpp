#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

namespace {

using spip::pip::Location;
using spip::pip::point_in_polygon_sphere;

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
  // Q5 is the vertex case, and Q6 is the edge case.
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

  // Q6 is the edge case from tests/test_pip_complicated_visualization.nb.
  // Use the notebook-derived long-double literals here so the great-circle
  // edge relation is preserved for the exact boundary classification.
  const std::vector<spip::V3_T<long double>> poly_q6 = {
      spip::V3_T<long double>(0.7711488862338939133931649L,
                              -0.1572614264676412392115167L,
                              0.6169264453770705474203694L),
      spip::V3_T<long double>(0.4524978914468171625447648L,
                              -0.7506135706341584205684991L,
                              0.4814820098570908587316732L),
      spip::V3_T<long double>(0.6894615088518674587796109L,
                              -0.5993397458796933456419420L,
                              0.4067366430758002077811529L),
      spip::V3_T<long double>(0.5339836142401214659925640L,
                              -0.8214480287797472635711628L,
                              0.2002114775354416563695869L),
      spip::V3_T<long double>(0.7254734110258385227200731L,
                              -0.6306444148430617328881034L,
                              0.2756373558169991856616615L),
      spip::V3_T<long double>(0.9066264675200399296428611L,
                              -0.3728891657256026464833856L,
                              0.1974388980839338736920788L),
      spip::V3_T<long double>(0.7473647984679656572099687L,
                              -0.6496743076188995401018013L,
                              0.1391731009600654440819083L),
      spip::V3_T<long double>(0.7546808431945165418323895L,
                              -0.6560340482729605666436731L,
                              -0.008726535498373935084119731L),
      spip::V3_T<long double>(0.4913862536359132189795005L,
                              -0.8536808575666768258180354L,
                              -0.1725356286738631755915772L),
      spip::V3_T<long double>(0.8655535612362530058770151L,
                              -0.2393261584350430547132437L,
                              -0.4399318384931519002251695L),
      spip::V3_T<long double>(0.7381999514442093898031194L,
                              -0.2609677456603186541559628L,
                              -0.6220584115762266047468948L),
      spip::V3_T<long double>(0.6016613961720088875036530L,
                              -0.05234812405382050469132365L,
                              -0.7970340257883566593348316L),
  };
  const spip::V3_T<long double> q6(0.7530936890789537183988137L,
                                   -0.6546543562494313614041203L,
                                   0.06540312923014306675410312L);
  const Location loc_q6 = point_in_polygon_sphere(q6, poly_q6);

  expect_equal(loc_q1, Location::Inside, "complicated polygon: Q1 inside");
  expect_equal(loc_q2, Location::Inside, "complicated polygon: Q2 inside");
  expect_equal(loc_q3, Location::Outside, "complicated polygon: Q3 outside");
  expect_equal(loc_q4, Location::Outside, "complicated polygon: Q4 outside");
  expect_equal(loc_q5, Location::OnVertex, "complicated polygon: Q5 on vertex");
  expect_equal(loc_q6, Location::OnEdge, "complicated polygon: Q6 on edge");

  const std::vector<std::int64_t> vertex_ids = {
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
  };
  const std::int64_t q1_id = 200;
  const std::int64_t q2_id = 201;
  const std::int64_t q3_id = 202;
  const std::int64_t q4_id = 203;
  const std::int64_t q5_id = 204;
  const std::int64_t q6_id = 205;
  const std::int64_t r1_id = 300;
  const std::int64_t r2_id = 301;
  const std::int64_t r3_id = 302;
  const std::int64_t r4_id = 303;
  const std::int64_t r5_id = 304;
  const std::int64_t r6_id = 305;

  // Tier 1: full global robustness. Vertex IDs, query IDs, and the designated
  // outside points R1..R6 all participate in one explicit symbolic ordering.
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
  const spip::V3_T<long double> r6(-0.7530936890789539095436567L,
                                   0.6546543562494314683064545L,
                                   -0.0654031292301420763557473L);

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
  const Location loc_q6_tier1 =
      point_in_polygon_sphere(q6, q6_id, r6, r6_id, poly_q6, vertex_ids);

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
  expect_equal(loc_q6_tier1, Location::OnEdge,
               "complicated polygon tier1: Q6 on edge");

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
  const Location loc_q6_tier2 =
      point_in_polygon_sphere(q6, q6_id, poly_q6, vertex_ids);

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
  expect_equal(loc_q6_tier2, Location::OnEdge,
               "complicated polygon tier2: Q6 on edge");

  // Tier 3: local/internal robustness. Only the vertex IDs are explicit; the
  // library assigns internal IDs to both q and the inferred R.
  const Location loc_q1_tier3 = point_in_polygon_sphere(q1, poly, vertex_ids);
  const Location loc_q2_tier3 = point_in_polygon_sphere(q2, poly, vertex_ids);
  const Location loc_q3_tier3 = point_in_polygon_sphere(q3, poly, vertex_ids);
  const Location loc_q4_tier3 = point_in_polygon_sphere(q4, poly, vertex_ids);
  const Location loc_q5_tier3 = point_in_polygon_sphere(q5, poly, vertex_ids);
  const Location loc_q6_tier3 =
      point_in_polygon_sphere(q6, poly_q6, vertex_ids);

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
  expect_equal(loc_q6_tier3, Location::OnEdge,
               "complicated polygon tier3: Q6 on edge");

  // EFT test on the same complicated polygon and the same query data.
  // Visualization companion notebook:
  //   tests/test_pip_complicated_visualization.nb
  const std::vector<spip::V3_T<long double>> poly_eft = {
      spip::V3_T<long double>(0.77114888623389370L, -0.15726142646764130L,
                              0.61692644537707060L),
      spip::V3_T<long double>(0.45249789144681710L, -0.75061357063415830L,
                              0.48148200985709080L),
      spip::V3_T<long double>(0.68946150885186746L, -0.59933974587969335L,
                              0.40673664307580021L),
      spip::V3_T<long double>(0.53398361424012150L, -0.82144802877974800L,
                              0.20021147753544170L),
      spip::V3_T<long double>(0.72547341102583852L, -0.63064441484306173L,
                              0.27563735581699919L),
      spip::V3_T<long double>(0.90662646752004000L, -0.37288916572560260L,
                              0.19743889808393390L),
      spip::V3_T<long double>(0.74736479846796566L, -0.64967430761889954L,
                              0.13917310096006544L),
      spip::V3_T<long double>(0.75468084319451650L, -0.65603404827296060L,
                              -0.00872653549837396L),
      spip::V3_T<long double>(0.49138625363591330L, -0.85368085756667700L,
                              -0.17253562867386300L),
      spip::V3_T<long double>(0.86555356123625300L, -0.23932615843504300L,
                              -0.43993183849315200L),
      spip::V3_T<long double>(0.73819995144420940L, -0.26096774566031860L,
                              -0.62205841157622660L),
      spip::V3_T<long double>(0.60166139617200880L, -0.05234812405382043L,
                              -0.79703402578835670L),
  };
  const spip::V3_T<long double> q1_eft(0.75367527697268680L,
                                       -0.65515992289232780L,
                                       -0.05233595624294383L);
  const spip::V3_T<long double> q2_eft(0.92054211727315200L,
                                       -0.38498585550407840L,
                                       0.06624274592780397L);
  const spip::V3_T<long double> q3_eft(0.53882393432914170L,
                                       -0.82565565483991800L,
                                       0.16721694718218960L);
  const spip::V3_T<long double> q4_eft(0.63494819288856630L,
                                       -0.65761549896072850L,
                                       0.40544130015845230L);
  const spip::V3_T<long double> q5_eft(0.49138625363591330L,
                                       -0.85368085756667700L,
                                       -0.17253562867386300L);
  const spip::V3_T<long double> r1_eft(-0.75367527657824340L,
                                       0.65515992254944344L,
                                       0.05233596621555329L);
  const spip::V3_T<long double> r2_eft(-0.92054211788294438L,
                                       0.38498585575910360L,
                                       -0.06624273597168498L);
  const spip::V3_T<long double> r3_eft(-0.53882393523014660L,
                                       0.82565565622055415L,
                                       -0.16721693746180466L);
  const spip::V3_T<long double> r4_eft(-0.63494819546290848L,
                                       0.65761550162697335L,
                                       -0.40544129180227878L);
  const spip::V3_T<long double> r5_eft(-0.49138625278809683L,
                                       0.85368085609377320L,
                                       0.17253563837617750L);

  // EFT no-global-ID path on the same data.
  const Location eft_loc_q1 = point_in_polygon_sphere(q1_eft, poly_eft);
  const Location eft_loc_q2 = point_in_polygon_sphere(q2_eft, poly_eft);
  const Location eft_loc_q3 = point_in_polygon_sphere(q3_eft, poly_eft);
  const Location eft_loc_q4 = point_in_polygon_sphere(q4_eft, poly_eft);
  const Location eft_loc_q5 = point_in_polygon_sphere(q5_eft, poly_eft);
  const Location eft_loc_q6 = point_in_polygon_sphere(q6, poly_q6);

  expect_equal(eft_loc_q1, Location::Inside,
               "complicated polygon eft/no-id: Q1 inside");
  expect_equal(eft_loc_q2, Location::Inside,
               "complicated polygon eft/no-id: Q2 inside");
  expect_equal(eft_loc_q3, Location::Outside,
               "complicated polygon eft/no-id: Q3 outside");
  expect_equal(eft_loc_q4, Location::Outside,
               "complicated polygon eft/no-id: Q4 outside");
  expect_equal(eft_loc_q5, Location::OnVertex,
               "complicated polygon eft/no-id: Q5 on vertex");
  expect_equal(eft_loc_q6, Location::OnEdge,
               "complicated polygon eft/no-id: Q6 on edge");

  // EFT Tier 1: full global robustness with explicit vertex IDs, query IDs,
  // and designated outside points R1..R6 from the same notebook dataset.
  const Location eft_loc_q1_tier1 =
      point_in_polygon_sphere(q1_eft, q1_id, r1_eft, r1_id, poly_eft, vertex_ids);
  const Location eft_loc_q2_tier1 =
      point_in_polygon_sphere(q2_eft, q2_id, r2_eft, r2_id, poly_eft, vertex_ids);
  const Location eft_loc_q3_tier1 =
      point_in_polygon_sphere(q3_eft, q3_id, r3_eft, r3_id, poly_eft, vertex_ids);
  const Location eft_loc_q4_tier1 =
      point_in_polygon_sphere(q4_eft, q4_id, r4_eft, r4_id, poly_eft, vertex_ids);
  const Location eft_loc_q5_tier1 =
      point_in_polygon_sphere(q5_eft, q5_id, r5_eft, r5_id, poly_eft, vertex_ids);
  const Location eft_loc_q6_tier1 =
      point_in_polygon_sphere(q6, q6_id, r6, r6_id, poly_q6, vertex_ids);

  expect_equal(eft_loc_q1_tier1, Location::Inside,
               "complicated polygon eft/tier1: Q1 inside");
  expect_equal(eft_loc_q2_tier1, Location::Inside,
               "complicated polygon eft/tier1: Q2 inside");
  expect_equal(eft_loc_q3_tier1, Location::Outside,
               "complicated polygon eft/tier1: Q3 outside");
  expect_equal(eft_loc_q4_tier1, Location::Outside,
               "complicated polygon eft/tier1: Q4 outside");
  expect_equal(eft_loc_q5_tier1, Location::OnVertex,
               "complicated polygon eft/tier1: Q5 on vertex");
  expect_equal(eft_loc_q6_tier1, Location::OnEdge,
               "complicated polygon eft/tier1: Q6 on edge");

  // EFT Tier 2: explicit vertex IDs and query IDs; the library infers R and
  // assigns the internal R ID.
  const Location eft_loc_q1_tier2 =
      point_in_polygon_sphere(q1_eft, q1_id, poly_eft, vertex_ids);
  const Location eft_loc_q2_tier2 =
      point_in_polygon_sphere(q2_eft, q2_id, poly_eft, vertex_ids);
  const Location eft_loc_q3_tier2 =
      point_in_polygon_sphere(q3_eft, q3_id, poly_eft, vertex_ids);
  const Location eft_loc_q4_tier2 =
      point_in_polygon_sphere(q4_eft, q4_id, poly_eft, vertex_ids);
  const Location eft_loc_q5_tier2 =
      point_in_polygon_sphere(q5_eft, q5_id, poly_eft, vertex_ids);
  const Location eft_loc_q6_tier2 =
      point_in_polygon_sphere(q6, q6_id, poly_q6, vertex_ids);

  expect_equal(eft_loc_q1_tier2, Location::Inside,
               "complicated polygon eft/tier2: Q1 inside");
  expect_equal(eft_loc_q2_tier2, Location::Inside,
               "complicated polygon eft/tier2: Q2 inside");
  expect_equal(eft_loc_q3_tier2, Location::Outside,
               "complicated polygon eft/tier2: Q3 outside");
  expect_equal(eft_loc_q4_tier2, Location::Outside,
               "complicated polygon eft/tier2: Q4 outside");
  expect_equal(eft_loc_q5_tier2, Location::OnVertex,
               "complicated polygon eft/tier2: Q5 on vertex");
  expect_equal(eft_loc_q6_tier2, Location::OnEdge,
               "complicated polygon eft/tier2: Q6 on edge");

  // EFT Tier 3: only vertex IDs are explicit; the library assigns internal
  // IDs to q and the inferred R.
  const Location eft_loc_q1_tier3 =
      point_in_polygon_sphere(q1_eft, poly_eft, vertex_ids);
  const Location eft_loc_q2_tier3 =
      point_in_polygon_sphere(q2_eft, poly_eft, vertex_ids);
  const Location eft_loc_q3_tier3 =
      point_in_polygon_sphere(q3_eft, poly_eft, vertex_ids);
  const Location eft_loc_q4_tier3 =
      point_in_polygon_sphere(q4_eft, poly_eft, vertex_ids);
  const Location eft_loc_q5_tier3 =
      point_in_polygon_sphere(q5_eft, poly_eft, vertex_ids);
  const Location eft_loc_q6_tier3 =
      point_in_polygon_sphere(q6, poly_q6, vertex_ids);

  expect_equal(eft_loc_q1_tier3, Location::Inside,
               "complicated polygon eft/tier3: Q1 inside");
  expect_equal(eft_loc_q2_tier3, Location::Inside,
               "complicated polygon eft/tier3: Q2 inside");
  expect_equal(eft_loc_q3_tier3, Location::Outside,
               "complicated polygon eft/tier3: Q3 outside");
  expect_equal(eft_loc_q4_tier3, Location::Outside,
               "complicated polygon eft/tier3: Q4 outside");
  expect_equal(eft_loc_q5_tier3, Location::OnVertex,
               "complicated polygon eft/tier3: Q5 on vertex");
  expect_equal(eft_loc_q6_tier3, Location::OnEdge,
               "complicated polygon eft/tier3: Q6 on edge");

  return EXIT_SUCCESS;
}
