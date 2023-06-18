#include "util.hh"
using std::size_t;

std::vector<real_t> linspace(const real_t start, const real_t stop, const size_t num){
  assert(num > 1);
  std::vector<real_t> result = std::vector<real_t> (num);
  real_t delta = stop - start;
  real_t div = num-1;
  real_t step = delta / div;
  for (size_t i = 0; i < num; i++){
    result[i] = start + i * step;
  }
  return result;
}

std::vector<real_t> arange(const real_t start, const real_t stop, const real_t step){
  assert(step != 0);
  auto diff = stop-start;
  size_t num = static_cast<size_t>(diff/step)+1;
  std::vector<real_t> result = std::vector<real_t> (num);

  for (size_t i = 0; i < num; i++){
    result[i] = start + i * step;
  }
  return result;
}

size_t unravel_index(std::vector<size_t> indices, std::vector<size_t> shape){
  return _unravel_index(indices, shape, indices.size()-1);
}

size_t _unravel_index(std::vector<size_t> indices, std::vector<size_t> shape, size_t rek){
  if(rek <= 1){
    return indices[0];
  }
  else{
    size_t index = indices[rek-1];
    // std::cout << "index: " << index << std::endl;
    for(size_t i = 0; i < rek-1; i++){
      index = index * shape[i];
    }
    // std::cout << "index: " << index << std::endl;
    // indices.pop_back();
    // shape.pop_back();
    return index += _unravel_index(indices, shape, --rek);
  }
}
