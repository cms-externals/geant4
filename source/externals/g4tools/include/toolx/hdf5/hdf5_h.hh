// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef toolx_hdf5_h
#define toolx_hdf5_h

#include <hdf5.h>

// HDF5 1.8+ only: use H5Dcreate2, H5Gcreate2, H5Acreate2, etc. with correct argument lists
#define toolx_H5Dopen(loc, name) H5Dopen2((loc), (name), H5P_DEFAULT)

#define toolx_H5Dcreate(loc, name, type, space, plist) \
  H5Dcreate2((loc), (name), (type), (space), H5P_DEFAULT, (plist), H5P_DEFAULT)

#define toolx_H5Acreate(obj, name, type, space, plist) \
  H5Acreate2((obj), (name), (type), (space), H5P_DEFAULT, (plist))

#define toolx_H5Tarray_create(base_type_id, ndims, dims) \
  H5Tarray_create2((base_type_id), (ndims), (dims))

#define toolx_H5Tget_array_dims(type_id, dims) H5Tget_array_dims((type_id), (dims))

#define toolx_H5Gcreate(loc, name, plist) \
  H5Gcreate2((loc), (name), H5P_DEFAULT, (plist), H5P_DEFAULT)

#define toolx_H5Gopen(loc, name) H5Gopen2((loc), (name), H5P_DEFAULT)

#define toolx_H5Aiterate(loc_id, idx, op, op_data) \
  H5Aiterate2((loc_id), H5_INDEX_NAME, H5_ITER_INC, (hsize_t*)(idx), (op), (op_data))

#define toolx_H5free_memory H5free_memory

// Additional wrappers for consistency and future-proofing
#define toolx_H5Giterate(loc_id, name, idx, op, op_data) \
  H5Giterate((loc_id), (name), (idx), (op), (op_data))

#define toolx_H5Gget_objinfo(loc_id, name, follow_link, statbuf) \
  H5Gget_objinfo((loc_id), (name), (follow_link), (statbuf))

#define toolx_H5Aopen_name(obj_id, attr_name) \
  H5Aopen_by_name((obj_id), ".", (attr_name), H5P_DEFAULT, H5P_DEFAULT)

#define toolx_H5Aread(attr_id, mem_type_id, buf) H5Aread((attr_id), (mem_type_id), (buf))

#define toolx_H5Aclose(attr_id) H5Aclose((attr_id))

#define toolx_H5Aget_type(attr_id) H5Aget_type((attr_id))

#define toolx_H5Tclose(type_id) H5Tclose((type_id))

#endif
