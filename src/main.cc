#include "NetCDFData.hh"

typedef double real_t;
// "/media/zarnosch/2TB_SSD/FlowData/ensemble/0001/COMBINED_2011013100.nc"

static const int NC_ERR = 1;

int main(int argc, char* argv[]) {
  std::string netcdf_file_path;
  if (argc >= 2) {
    netcdf_file_path = argv[1];
  }

  std::unique_ptr<NetCDFData> data0 = std::make_unique<NetCDFData>();
  data0->setFilepath(netcdf_file_path);
  data0->loadDimensionsFromFile();
  data0->loadVariablesFromFile<4>();
  std::array<real_t, 4> position = {963930, 6.5, 20.0671, 38.0671};
  real_t lerped_value = data0->getLerpValue("U", position);
  std::cout << lerped_value << std::endl;
  return 0;
}