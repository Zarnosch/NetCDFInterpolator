#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <netcdf>
#include <exception>
#include <memory>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include "util.hh"
#include "gcem.hpp"
#include <bitset>

using namespace netCDF;

typedef double real_t;

class NetCDFData {
  private:
    std::string _filepath;
    std::unique_ptr<NcFile> _ncFile;
    std::string _metadata;
    /// @brief Map of Dimension key and index
    std::unordered_map<std::string, size_t> _dim_keys;
    /// @brief Dimension data
    std::vector<std::vector<real_t>> _dimensions;
    /// @brief Ordering of the dimensions
    std::vector<bool> _dimensions_ordering;
    /// @brief Size of the dimensions
    std::vector<size_t> _dimensions_size;
    /// @brief Mean stride between two data points of dimensions
    std::vector<real_t> _dimensions_mean_stepsize;
    /// @brief minimum value of the dimensions
    std::vector<real_t> _dimensions_min;
    /// @brief maximum value of the dimensons
    std::vector<real_t> _dimensions_max;

    /// @brief map data index for a data key
    std::unordered_map<std::string, size_t> _var_keys;
    /// @brief currently loaded variable data
    std::shared_ptr<std::vector<std::vector<real_t>>> _var_data;
    /// @brief number of dimensions for variables
    std::vector<size_t> _var_nDims;
    /// @brief number of subvariables for variables
    std::vector<size_t> _var_nSubvar;
    /// @brief Dimension keys for variables
    std::vector<std::unique_ptr<std::string[]>> _var_dimensions;
    /// @brief Dimension sizes for variables
    std::vector<std::unique_ptr<size_t[]>> _var_shapes;
    /// @brief Dimension subspaces for variables
    std::vector<std::unique_ptr<size_t[]>> _var_subspaces;

    
    void _openFile(){
      try {
      _ncFile = std::unique_ptr<NcFile> (new NcFile(_filepath, netCDF::NcFile::read));
      } 
      catch (netCDF::exceptions::NcException &e) {
        std::cout << e.what() << std::endl;
      };
    };


    template<size_t D>
    std::array<size_t, D> _get_subspaces(std::array<std::string, D> dimensions){
      std::array<size_t, D> subspaces;
      subspaces[D-1] = 1;
      for (size_t i = D-1; i > 0; i--)
      {
        subspaces[i-1] = subspaces[i]*_dimensions_size[getDimIndex(dimensions[i])];
      }
      return subspaces;
    }

    size_t* _getSubspaces(std::string var_key){return _var_subspaces[_var_keys[var_key]].get();}; 
    size_t* _getSubspaces(size_t var_index){return _var_subspaces[var_index].get();};
    size_t* _getShape(std::string var_key) {return _var_shapes[_var_keys[var_key]].get();};

    void _init(){
      _dim_keys = std::unordered_map<std::string, size_t> ();
      _var_data = std::make_shared<std::vector<std::vector<real_t>>>();
      _var_keys = std::unordered_map<std::string, size_t> ();
      _var_nDims = std::vector<size_t> ();
      _var_dimensions = std::vector<std::unique_ptr<std::string[]>>();
      _var_shapes = std::vector<std::unique_ptr<size_t[]>>();
      _var_subspaces = std::vector<std::unique_ptr<size_t[]>>();
    }

  public:
    NetCDFData(){_init();};
    ~NetCDFData(){};
    void setFilepath(std::string filepath)
    {
      _filepath = filepath;
      _openFile();
    };

    size_t const& getDimIndex(std::string dim_key){return _dim_keys[dim_key];}

    void loadDimensionsFromFile(){
      auto dims = _ncFile->getDims();
      for (const auto& pair : dims)
      {
        std::string key = pair.first;
        NcDim dim = pair.second;
        std::vector<real_t> dim_values(dim.getSize());
        if(_ncFile->getVar(key).getId() >= 0)
        {
          _ncFile->getVar(key).getVar(dim_values.data());
          addDimension(key, dim_values);
        }
      }
    };

    /// @brief Load all Variables from File with a given number of Dimensions
    /// @tparam D number of Dimensions
    template<size_t D>
    void loadVariablesFromFile(){
      auto vars = _ncFile->getVars();
      for (const auto& pair : vars)
      {
        std::string key = pair.first;
        NcVar var = pair.second;
        auto dims = var.getDims();
        // only load variables with correct number of dimensions
        if (dims.size() != D)
          continue;

        std::array<std::string, D> var_dims;
        for (size_t i = 0; i < D; i++)
        {
          auto dim = dims[i];
          if(_dim_keys.count(dim.getName()) == 0){
            break;
          }
          else{
            var_dims[i] = dim.getName();
          }
        }
        // calculate the required space to load the data
        size_t varSize = _get_subspaces(var_dims)[0]*_dimensions_size[getDimIndex(var_dims[0])];
        std::vector<real_t> var_values(varSize);
        if(_ncFile->getVar(key).getId() >= 0)
        {
          _ncFile->getVar(key).getVar(var_values.data());
          addVar(key, var_values, var_dims);
        }
      }
    }

    template<size_t P>
    inline std::array<size_t, P> ravel_index(const size_t* subspaces, const size_t &index){
      std::array<size_t, P> indices;
      size_t remainder = index;
      for(size_t rek = 0; rek < P; rek++){
        indices[rek] = remainder/subspaces[rek];
        remainder = remainder%subspaces[rek];
      }
      return indices;
    }

    template<size_t P>
    inline size_t unravel_index(const size_t* subspaces, const std::array<size_t, P> &indices){
      size_t index = 0;
      for(size_t rek = 0; rek != P; rek++){
        index += (indices[rek] * subspaces[rek]);
      }
      return index;
    }

    size_t addDimension(std::string key, std::vector<real_t> &data){
      size_t dimension_index = _dimensions.size();
      _dim_keys[key] = dimension_index;
      _dimensions_ordering.push_back(std::is_sorted(data.begin(), data.end()));
      _dimensions.push_back(data);
      _dimensions_size.push_back(data.size());
      if(_dimensions_ordering[dimension_index]){
        _dimensions_min.push_back(data[0]);
        _dimensions_max.push_back(data[data.size()-1]);
      }
      else{
        _dimensions_min.push_back(data[data.size()-1]);
        _dimensions_max.push_back(data[0]);
      }
      _dimensions_mean_stepsize.push_back((_dimensions_max[dimension_index]-_dimensions_min[dimension_index])/(_dimensions_size[dimension_index]-1));
      // std::cout << "Added Dimension " << key << " at index: " << dimension_index << " and size " << data.size() << std::endl;
      return dimension_index;
    };

    template<size_t D>
    bool addVar(std::string key, std::vector<real_t> &data, std::array<std::string, D> var_dimensions){
      // std::cout << "add data with key: " << key.c_str() << " at position: " << (*_griddata).size() << std::endl;
      size_t sub_vars = 1;

      // check if datasize fits the dimensions
      size_t dimensions_size = 1;
      for (size_t i = 0; i < D; i++)
      {
        dimensions_size*=_dimensions[getDimIndex(var_dimensions[i])].size();
      }

      size_t var_index = (*_var_data).size();
      // check if key is already presend and only update if this is the case
      auto key_pos = _var_keys.find(key);
      if(key_pos != _var_keys.end()){
        var_index = _var_keys[key];
        _var_keys[key] = var_index;
        (*_var_data)[var_index] = data;
        _var_nDims[var_index] = D;
        _var_nSubvar[var_index] = sub_vars;
        _var_dimensions[var_index] = std::make_unique<std::string[]>(D);
        _var_shapes[var_index] = std::make_unique<size_t[]>(D);
        _var_subspaces[var_index] = std::make_unique<size_t[]>(D);
      }
      else{ // add data at the end if the key is not present already
        _var_keys[key] = var_index;
        (*_var_data).push_back(data);
        _var_nDims.push_back(D);
        _var_nSubvar.push_back(sub_vars);
        _var_dimensions.push_back(std::make_unique<std::string[]>(D));
        _var_shapes.push_back(std::make_unique<size_t[]>(D));
        _var_subspaces.push_back(std::make_unique<size_t[]>(D));
      }

      for (size_t i = 0; i < D; i++)
      {
        _var_dimensions[var_index][i] = var_dimensions[i];
        _var_shapes[var_index][i] = _dimensions[getDimIndex(var_dimensions[i])].size();
      }

      _var_subspaces[var_index][D-1] = 1;
      for (size_t i = D-1; i > 0; i--)
      {
        _var_subspaces[var_index][i-1] = _var_subspaces[var_index][i]*_var_shapes[var_index][i];
      }
    
      return true;
    }

    template<size_t D>
    real_t getLerpValue(const std::string &key, std::array<real_t, D> position){

      constexpr size_t pow2D = gcem::pow(2, D); // calculate  the number of required samples at compiletime
      auto shape = _getShape(key); // get the length of each dimension for later
      auto subspaces = _getSubspaces(key); // get the subspaces of the dimensions for later

      // calculate interpolation factors for every dimension first
      std::array<real_t, D> interp_factors;
      std::array<size_t, D> position_indices;
      for (size_t i = 0; i < D; i++)
      {
        size_t dim_i = getDimIndex(_var_dimensions[_var_keys[key]][i]); // get the internal index of this dimension
        size_t upper_bound = own_upper_bound(dim_i, position[i]); // binary search depending on dimension ordering
        auto start_index = upper_bound -1;
        position_indices[i] = start_index;
        interp_factors[i] = (position[i]-_dimensions[dim_i][start_index])/
          (_dimensions[dim_i][start_index+1] -_dimensions[dim_i][start_index]);
        if (shape[i] == 1)
        {
          position_indices[i] = 0;
          interp_factors[i] = 1;
        }
      }

      std::array<real_t, pow2D> values; // the original sample points from which we interpolate
      std::array<std::array<size_t, D>, pow2D> changed_indices;
      std::bitset<D> bits; // Bitmask to shift the indices to acess the 1D bulk of data
      // get the hyperslab with all needed values
      for (size_t i = 0; i < pow2D; i++)
      {
        bits = std::bitset<D>(i);
        for (size_t j = 0; j < D; j++)
        {
          if (shape[j] == 1){
            changed_indices[i][j] = position_indices[j];
          }
          else{
            changed_indices[i][j] = position_indices[j] + bits[j];
          }
        }
        values[i] = (*_var_data)[_var_keys[key]][unravel_index<D>(subspaces, changed_indices[i])];
      }
      for (size_t j = 0; j < D; j++)
      {
        for (size_t i = 0; i < std::pow(2, D-j); i+=2)
        {
          values[i/2] = std::lerp(values[i], values[i+1], interp_factors[j]);
        }
      }
      return values[0];
    }

    inline size_t own_upper_bound(const size_t dim_i, const real_t value) const{
      size_t index = 0;
      if(_dimensions_ordering[dim_i]){
        index = std::distance(_dimensions[dim_i].begin(), std::upper_bound(_dimensions[dim_i].begin(), _dimensions[dim_i].end(), value));
      }
      else{
        index = std::distance(_dimensions[dim_i].begin(), std::upper_bound(_dimensions[dim_i].begin(), _dimensions[dim_i].end(), value, std::greater<double>()));
      }
      if(index == _dimensions[dim_i].size()){
        index--;
      }
      return index;
    }

    /*
     * @brief brief Get the Lerp Values object
     * @tparam D 
     * @param keys 
     * @param position 
     * @return std::vector<real_t> 
     * 
     * calculates the interpolation of multiple fields with the same dimensions at the identical position
     * 
     */
    template<size_t D, size_t V>
    std::array<real_t, V> getLerpValues(const std::array<size_t, V> data_key_indices, std::array<real_t, D> position){
      constexpr size_t pow2D = gcem::pow(2, D);
      constexpr size_t pow2D2 = gcem::pow(2, D)/2;

      // calculate interpolation factors for each dimension first
      std::array<real_t, D> interp_factors;
      std::array<size_t, D> position_indices;
      std::array<std::array<real_t, pow2D2>, V> values;
      std::array<std::array<size_t, D>, pow2D> changed_indices;
      std::array<real_t, V> result;
      std::bitset<D> bits;
      // #pragma omp parallel for num_threads(D)
      for (size_t i = 0; i < D; i++)
      {
        size_t dim_i = getDimIndex(_var_dimensions[data_key_indices[0]][i]);
        size_t upper_bound = own_upper_bound(dim_i, position[i]);
        size_t start_index = upper_bound -1;
        position_indices[i] = start_index;
        interp_factors[i] = (position[i]-_dimensions[dim_i][start_index])/(_dimensions[dim_i][start_index+1] -_dimensions[dim_i][start_index]);

        if (_var_shapes[data_key_indices[0]][i] == 1)
        {
          position_indices[i] = 0;
          interp_factors[i] = 1;
        }
      }

      auto subspaces = getSubspaces(data_key_indices[0]);
      // get the hyperslab with all needed values
      for (size_t i = 0; i < pow2D; i++)
      {
        bits = std::bitset<D>(i);
        for (size_t j = 0; j < D; j++)
        {
          if (_var_shapes[data_key_indices[0]][j] == 1){
            changed_indices[i][j] = position_indices[j];
          }
          else{
            changed_indices[i][j] = position_indices[j] + bits[j];
          }         
        }
      }


      for (size_t di = 0; di < V; di++)
      {
        for (size_t i = 0; i < pow2D; i+=2)
        {
          real_t value_0 = (*_var_data)[data_key_indices[di]][unravel_index<D>(subspaces, changed_indices[i])];
          real_t value_1 = (*_var_data)[data_key_indices[di]][unravel_index<D>(subspaces, changed_indices[i+1])];
          values[di][i/2] = std::lerp(value_0, value_1, interp_factors[0]);
        }
        for (size_t j = 1; j < D; j++)
        {
          for (size_t i = 0; i < gcem::pow(2, D-j); i+=2)
          {
            values[di][i/2] = std::lerp(values[di][i], values[di][i+1], interp_factors[j]);
          }
        }
        result[di] = values[di][0];
      }
      
      return result;
    }

    template<size_t D>
    std::shared_ptr<std::vector<real_t>> getHyperslab(const std::string key, const std::array<size_t, D> position, const std::array<size_t, D> stride){
      size_t newsize = 1;
      for (size_t i = 0; i < D; i++)
      {
        if(stride[i] > 0)
          newsize*=stride[i];
      }
      
      auto result =std::make_shared<std::vector<real_t>>(newsize);
      std::array<size_t, D> position_indices;
      std::array<size_t, D> final_position_indices;
      size_t counter = 0;

      for (size_t i = 0; i < D; i++)
      {
        position_indices[i] = position[i];
        final_position_indices[i] = position[i] + stride[i]-1;
      }

      size_t idx = 0;
      auto subspace = _getSubspaces(key);
      while (position_indices != final_position_indices)
      {
        (*result)[idx++] = (*_var_data)[_var_keys[key]][unravel_index(subspace, position_indices)];
        int dim_i = D-1;
        bool increased = false;
        while (!increased && dim_i >= 0)
        {
          if(position_indices[dim_i] < final_position_indices[dim_i]){
            position_indices[dim_i] +=1;
            increased = true;
          }
          else{
            position_indices[dim_i] = position[dim_i];
            dim_i--;
          }
        }
      }
      (*result)[idx] = (*_var_data)[_var_keys[key]][unravel_index(subspace, final_position_indices)];
      return result;
    }
};