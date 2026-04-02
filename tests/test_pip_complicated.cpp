#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

namespace {

void expect(Location got, Location expected, const char* name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << name << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << name << '\n';
}

}  // namespace

int main() {
  const std::vector<std::array<double, 3>> polygon = {
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
  const std::array<double, 3> q5 = {
      0.49138625363591330, -0.85368085756667700, -0.17253562867386300};

  expect(point_in_polygon_sphere(q1, polygon), Location::Inside, "Q1");
  expect(point_in_polygon_sphere(q2, polygon), Location::Inside, "Q2");
  expect(point_in_polygon_sphere(q3, polygon), Location::Outside, "Q3");
  expect(point_in_polygon_sphere(q4, polygon), Location::Outside, "Q4");
  expect(point_in_polygon_sphere(q5, polygon), Location::OnVertex, "Q5");

  const std::vector<std::int64_t> ids = {
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111};
  const std::array<double, 3> r1 = {
      -0.75367527657824340, 0.65515992254944344, 0.05233596621555329};

  expect(point_in_polygon_sphere(q1, polygon, ids), Location::Inside,
         "Q1 tier3");
  expect(point_in_polygon_sphere(q1, 200, polygon, ids), Location::Inside,
         "Q1 tier2");
  expect(point_in_polygon_sphere(q1, 200, r1, 300, polygon, ids),
         Location::Inside, "Q1 tier1");
  return EXIT_SUCCESS;
}
