// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/** \file stereo_blend.cc

The purpose of this tool is to blend the boundaries of extra-large
stereo_corr tiles so that no seams are visible.  This is only required
when running the SGM algorithm on large data sets using parallel_stereo.
*/

#include <asp/Tools/stereo.h>
#include <vw/Stereo/DisparityMap.h>
#include <asp/Sessions/ResourceLoader.h>
//#include <xercesc/util/PlatformUtils.hpp>
#include <boost/filesystem.hpp>

using namespace vw;
using namespace vw::stereo;
using namespace asp;
using namespace std;

//namespace vw {
//  template<> struct PixelFormatID<PixelMask<Vector<float, 5> > >   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_6_CHANNEL; };
//}

typedef DiskImageView<PixelMask<Vector2f> > DiskImageType;
typedef ImageView    <PixelMask<Vector2f> > DispImageType;
typedef ImageView    <double              > WeightsType;

const size_t NUM_NEIGHBORS = 8;
enum Position {TL = 0, T = 1, TR = 2,
                L = 3, M = 8,  R = 4,
               BL = 5, B = 6, BR = 7};
// M is the central non-buffered area.

/// Debugging aid
std::string position_string(int p) {
  switch(p) {
    case TL: return "TL";
    case T:  return "T";
    case TR: return "TR";
    case L:  return "L";
    case R:  return "R";
    case BL: return "BL";
    case B:  return "B";
    case BR: return "BR";
    default: return "M";
  };
}

/// Returns the opposite position (what it is in the neighbor)
Position get_opposed_position(Position p) {
  switch(p) {
    case TL: return BR;
    case T:  return B;
    case TR: return BL;
    case L:  return R;
    case R:  return L;
    case BL: return TR;
    case B:  return T;
    case BR: return TL;
    default: return M;
  };
}

/// Given "out-2048_0_1487_2048" return "2048_0_1487_2048"
std::string extract_process_folder_bbox_string(std::string const& s) {
  size_t num_start = s.rfind("-");
  if (num_start == std::string::npos)
    vw_throw( ArgumentErr() << "Error parsing folder string: " << s );
  return s.substr(num_start+1);
}

/// Constructs a BBox2i from a parallel_stereo formatted folder.
BBox2i bbox_from_folder(std::string const& s) {
  std::string cropped = extract_process_folder_bbox_string(s);
  int x, y, width, height;
  sscanf(cropped.c_str(), "%d_%d_%d_%d", &x, &y, &width, &height);
  return BBox2i(x, y, width, height);
}


/// Returns one of eight possible ROI locations for the given tile.
/// - If get_buffer is set, fetch the ROI from the buffer region.
///   Otherwise get it from the non-buffer region at that location.
/// - Some tiles do not have all buffers available, if one of these
///   is requested the function will return false.
/// - Set buffers_stripped if you want the output ROI in reference to
///   an image with the buffers removed.  This will always fail if combined
///   with get_buufer==true.
/// - Generally you would get the non-buffer region for the main tile,
///   and the opposed buffer region for the neighboring tile.
bool get_roi_from_tile(std::string const& tile_path, Position pos,
                              int buffer_size, bool get_buffer,
                              BBox2i &output_roi,
                              bool buffers_stripped=false) {

  // Get the ROI of this tile and of the entire processing job
  BBox2i tile_roi = bbox_from_folder(tile_path);
  BBox2i full_roi = stereo_settings().left_image_crop_win;
  
  // Determine if this tile sits on the border of the tiles.
  bool left_edge  = (tile_roi.min().x() == full_roi.min().x());
  bool top_edge   = (tile_roi.min().y() == full_roi.min().y());
  bool right_edge = (tile_roi.max().x() == full_roi.max().x());
  bool bot_edge   = (tile_roi.max().y() == full_roi.max().y());
  if (buffers_stripped) {
    left_edge = top_edge = right_edge = bot_edge = true;
  }
  
  // Get the size of the buffers on each edge
  int left_offset  = (left_edge ) ? 0 : buffer_size;
  int right_offset = (right_edge) ? 0 : buffer_size;
  int top_offset   = (top_edge  ) ? 0 : buffer_size;
  int bot_offset   = (bot_edge  ) ? 0 : buffer_size;
  
  // Adjust the size of the output ROI according to the available buffer area
  int roi_width  = tile_roi.width()  - (left_offset + right_offset);
  int roi_height = tile_roi.height() - (top_offset  + bot_offset  );
  
  // The three sizes of bboxes that will be used.
  Vector2i corner_size         (buffer_size, buffer_size);
  Vector2i horizontal_edge_size(roi_width,   buffer_size);
  Vector2i vertical_edge_size  (buffer_size, roi_height);
  Vector2i dummy(0,0);

  switch(pos) {
  case TL: 
    if (get_buffer) {
      if (left_edge || top_edge)
        return false;
      output_roi = BBox2i(Vector2i(0, 0), dummy);
    } else
      output_roi = BBox2i(Vector2i(left_offset, top_offset), dummy);
    output_roi.set_size(corner_size);
    break;
  case T:
    if (get_buffer) {
      if (top_edge)
        return false;
      output_roi = BBox2i(Vector2i(left_offset, 0), dummy);
    } else
      output_roi = BBox2i(Vector2i(left_offset, top_offset), dummy);
    output_roi.set_size(horizontal_edge_size);
    break;             
  case TR:
    if (get_buffer) {
      if (right_edge || top_edge)
          return false;
      output_roi = BBox2i(Vector2i(tile_roi.width()-buffer_size, 0), dummy);
    } else
      output_roi = BBox2i(Vector2i(tile_roi.width()-buffer_size-right_offset, top_offset), dummy);
    output_roi.set_size(corner_size);
    break;
  case L:
    if (get_buffer) {
      if (left_edge)
        return false;
      output_roi = BBox2i(Vector2i(0, top_offset), dummy);
    } else
      output_roi = BBox2i(Vector2i(left_offset, top_offset), dummy);
    output_roi.set_size(vertical_edge_size);
    break;
  case R:
    if (get_buffer) {
      if (right_edge)
        return false;
      output_roi = BBox2i(Vector2i(tile_roi.width()-buffer_size, 0), dummy);
    } else
      output_roi = BBox2i(Vector2i(tile_roi.width()-buffer_size-right_offset, top_offset), dummy);
    output_roi.set_size(vertical_edge_size);
    break;
  case BL:
    if (get_buffer) {
      if (left_edge || bot_edge)
        return false;
      output_roi = BBox2i(Vector2i(0, tile_roi.height()-buffer_size), dummy);
    } else
      output_roi = BBox2i(Vector2i(left_offset, tile_roi.height()-buffer_size-bot_offset), dummy);
    output_roi.set_size(corner_size);
    break;
  case B:
    if (get_buffer) {
      if (bot_edge)
        return false;
      output_roi = BBox2i(Vector2i(left_offset, tile_roi.height()-buffer_size), dummy);
    } else
      output_roi = BBox2i(Vector2i(left_offset, tile_roi.height()-buffer_size-bot_offset), dummy);
    output_roi.set_size(horizontal_edge_size);
    break;
  case BR:
    if (get_buffer) {
      if (right_edge || bot_edge) 
        return false;
      output_roi = BBox2i(Vector2i(tile_roi.width()-buffer_size, tile_roi.height()-buffer_size), dummy);
    } else
      output_roi = BBox2i(Vector2i(tile_roi.width() -buffer_size-right_offset, 
                                 tile_roi.height()-buffer_size-bot_offset), dummy);
    output_roi.set_size(corner_size);
    break;
  default: // M (central area, everything except for the buffers)
    if (get_buffer)
      return false; // Central area is never a buffer
    output_roi = BBox2i(Vector2i(left_offset, top_offset), dummy);
    output_roi.set_size(Vector2i(roi_width, roi_height));
    break;
  };
  return true;
}

/*
/// Compute the bounding boxes for blending with one neigboring tile
/// - Note that the edge ROIs include the corner ROIs
void compute_rois(BBox2i const& input_bbox, ///< BBox in the input (center) tile.
                  Position pos, 
                  Vector2i const& tile_size, int buffer_size, 
                  BBox2i& tile_roi, BBox2i& input_roi) {

  // No tile should ever be smaller than buffer_size in width and height
  
  // The three sizes of bboxes that will be used.
  Vector2i corner_size         (buffer_size,        buffer_size);
  Vector2i horizontal_edge_size(input_bbox.width(), buffer_size);
  Vector2i vertical_edge_size  (buffer_size,        input_bbox.height());
  Vector2i dummy(0,0);

  switch(pos) {
  case TL: 
    tile_roi  = BBox2i(Vector2i(tile_size[0]-buffer_size, tile_size[1]-buffer_size), dummy);
    input_roi = BBox2i(input_bbox.min(), dummy);
    tile_roi.set_size (corner_size);
    input_roi.set_size(corner_size);
    break;
  case T:
    tile_roi  = BBox2i(Vector2i(buffer_size, tile_size[1]-buffer_size), dummy);
    input_roi = BBox2i(input_bbox.min(), dummy);
    tile_roi.set_size (horizontal_edge_size);
    input_roi.set_size(horizontal_edge_size);
    break;             
  case TR:
    tile_roi  = BBox2i(Vector2i(0, tile_size[1]-buffer_size),  dummy);
    input_roi = BBox2i(input_bbox.min() + Vector2i(input_bbox.width()-buffer_size, 0), dummy);
    tile_roi.set_size (corner_size);
    input_roi.set_size(corner_size);
    break;
  case L:
    tile_roi  = BBox2i(Vector2i(tile_size[0]-buffer_size, buffer_size), dummy);
    input_roi = BBox2i(input_bbox.min(), dummy);
    tile_roi.set_size (vertical_edge_size);
    input_roi.set_size(vertical_edge_size);
    break;
  case R:
    tile_roi  = BBox2i(Vector2i(0, buffer_size), dummy);
    input_roi = BBox2i(input_bbox.min() + Vector2i(input_bbox.width()-buffer_size, 0), dummy);
    tile_roi.set_size (vertical_edge_size);
    input_roi.set_size(vertical_edge_size);
    break;
  case BL:
    tile_roi  = BBox2i(Vector2i(tile_size[0]-buffer_size, 0), dummy);
    input_roi = BBox2i(input_bbox.min() + Vector2i(0, input_bbox.height()-buffer_size), dummy);
    tile_roi.set_size (corner_size);
    input_roi.set_size(corner_size);
    break;
  case B:
    tile_roi  = BBox2i(Vector2i(buffer_size, 0), dummy);
    input_roi = BBox2i(input_bbox.min() + Vector2i(0, input_bbox.height()-buffer_size), dummy);
    tile_roi.set_size (horizontal_edge_size);
    input_roi.set_size(horizontal_edge_size);
    break;
  default: // BR
    tile_roi  = BBox2i(Vector2i(0, 0), dummy);
    input_roi = BBox2i(input_bbox.max() - Vector2i(buffer_size, buffer_size), dummy);
    tile_roi.set_size (corner_size);
    input_roi.set_size(corner_size);
    break;
  };
}
*/
/// Load the desired portion of a disparity tile and associated image weights.
bool load_image_and_weights(std::string const& file_path, BBox2i const& roi,
                            DispImageType & image, WeightsType & weights) {
  // Verify image exists
  if (file_path == "")
    return false;
    
  // Load the image from disk
  DispImageType full_image = DiskImageType(file_path);
  
  // Compute the desired weights
  centerline_weights(full_image, weights, roi);
  
  // Extract the desired portion of the full image
  image = crop(full_image, roi);
  return true;
}

/// Perform the blending of 
void blend_tile_region(DispImageType      & main_image, WeightsType      & main_weights,     BBox2i const& main_roi,
                       DispImageType const& neighbor,   WeightsType const& neighbor_weights, BBox2i const& neighbor_roi) {

  // Multiply-accumulate the image values, then accumulate the weights.
  crop(main_image,   main_roi) += neighbor * neighbor_weights;
  crop(main_weights, main_roi) += neighbor_weights;
}

struct BlendOptions {
  std::string main_path;
  BBox2i      main_roi;
  std::string tile_paths[NUM_NEIGHBORS];
  BBox2i      rois      [NUM_NEIGHBORS]; // TODO: Keep?
  int sgm_collar_size;
  
  // Helper functions to indicate which directions tiles are present in
  bool has_upper_tiles() const { 
    return ((tile_paths[TL] != "") || (tile_paths[T] != "") || (tile_paths[TR] != ""));
  };
  bool has_lower_tiles() const { 
    return ((tile_paths[BL] != "") || (tile_paths[B] != "") || (tile_paths[BR] != ""));
  };
  bool has_left_tiles() const { 
    return ((tile_paths[TL] != "") || (tile_paths[L] != "") || (tile_paths[BL] != ""));
  };
  bool has_right_tiles() const { 
    return ((tile_paths[TR] != "") || (tile_paths[R] != "") || (tile_paths[BR] != ""));
  };
};





/// Blend the borders of an input disparity tile using the adjacent disparity tiles.
DispImageType
tile_blend( DispImageType const& input_image,
            BlendOptions const& opt) {

  // The amount of padding applied to each tile.
  int buff_size = opt.sgm_collar_size;

  const bool GET_BUFFER = true;
  const bool NOT_BUFFER = false;

  // Retrieve the output bounding box in the input image
  BBox2i output_bbox;
  get_roi_from_tile(opt.main_path, M, buff_size, NOT_BUFFER, output_bbox);

  // The bbox of the input image and the bbox that will be written as
  //  output (with the padding removed)
  BBox2i input_bbox  = bounding_box(input_image);  
  std::cout << "Input bbox = " << input_bbox << std::endl;
  std::cout << "Output bbox = " << output_bbox << std::endl;
  
  // Allocate the output image as a copy of the input image.
  // - This sets the non-blended portion of the image.
  DispImageType output_image = crop(input_image, output_bbox);

  // Compute weights for the main tile
  WeightsType main_weights;
  centerline_weights(input_image, main_weights, output_bbox);

  write_image("main_weights.tif", main_weights);

  // Load the neighboring eight tiles
  //Vector2i      tile_sizes [NUM_NEIGHBORS];
  BBox2i        tile_rois  [NUM_NEIGHBORS]; // ROIs in the neighbors
  BBox2i        input_rois [NUM_NEIGHBORS]; // ROIs in the main image
  DispImageType images     [NUM_NEIGHBORS]; // Contains cropped regions
  WeightsType   weights    [NUM_NEIGHBORS]; // Contains cropped regions

  // Figure out the ROIs, load images, and initialize output blend values.
  for (size_t i=0; i<NUM_NEIGHBORS; ++i) {
    if (opt.tile_paths[i] == "")
      continue;
    //tile_sizes[i] = asp::file_image_size(opt.tile_paths[i]);
    //compute_rois(output_bbox, Position(i), tile_sizes[i], buff_size, tile_rois[i], input_rois[i]);
    
    // Get the ROI from the cropped input image, then from the neighboring tile
    get_roi_from_tile(opt.main_path, Position(i), buff_size, NOT_BUFFER, input_rois[i], true);
    get_roi_from_tile(opt.tile_paths[i], get_opposed_position(Position(i)), 
                      buff_size, GET_BUFFER, tile_rois[i]);
    
    std::cout << "For tile " << position_string(i) << ", tile roi = " 
              << tile_rois[i] << ", input_roi = " << input_rois[i] << std::endl;
    //std::cout << opt.tile_paths[i] << std::endl;
    
    load_image_and_weights(opt.tile_paths[i], tile_rois[i], images[i], weights[i]);
  }

  std::cout << "Premultiply...\n";

//vw_throw( ArgumentErr() << "DEBUG!" );

  std::cout << "BBox output = " << bounding_box(output_image) << std::endl;
  std::cout << "BBox weights = " << bounding_box(main_weights) << std::endl;
  
  // Multiply the main image values by blend weights
  output_image *= main_weights;
  
  std::cout << "Performing blending...\n";

  // Blend in the neighbors one section at a time.
  for (size_t i=0; i<NUM_NEIGHBORS; ++i) {
    if (opt.tile_paths[i] == "")
      continue;
    std::cout << "Blending tile " << position_string(i) << std::endl;
    blend_tile_region(output_image, main_weights, input_rois[i],
                      images[i],    weights[i],   tile_rois[i]);
  }

  std::cout << "Postmultiply...\n";
  
  // Normalize the main image values to account for the applied weighting.
  output_image /= main_weights;

  return output_image;
}


void fill_blend_options(ASPGlobalOptions const& opt, BlendOptions & blend_options) {

  blend_options.main_path = opt.out_prefix + "-Dnosym.tif";

  //std::cout << "left_image_crop_win  = " << stereo_settings().left_image_crop_win << std::endl;
  //std::cout << "right_image_crop_win = " << stereo_settings().right_image_crop_win << std::endl;
  //std::cout << "trans_crop_win = " << stereo_settings().trans_crop_win << std::endl;
  
  boost::filesystem::path parallel_stereo_folder(opt.out_prefix);
  
  std::cout << "parallel_stereo_folder = " << parallel_stereo_folder.string() << std::endl;
  
  parallel_stereo_folder = parallel_stereo_folder.parent_path().parent_path();
  
  std::cout << "parallel_stereo_folder = " << parallel_stereo_folder.string() << std::endl;
  
  // Get a list of all the folders in the parallel_stereo output directory
  std::vector<std::string> folder_list;
  for (boost::filesystem::directory_iterator iter(parallel_stereo_folder); 
       iter!=boost::filesystem::directory_iterator(); ++iter) {
    //std::cout << itr->path();;
    if (boost::filesystem::is_directory(iter->status())) {
      //std::cout << " is a folder";
      folder_list.push_back(iter->path().filename().string());
    }
    //std::cout << '\n';
  }
  
  // Get the main tile bbox
  boost::filesystem::path mpath(blend_options.main_path);
  std::string mbb = mpath.parent_path().filename().string();
  std::cout << "mbb = " << mbb << std::endl;
  blend_options.main_roi = bbox_from_folder(mbb);
  BBox2i main_bbox = blend_options.main_roi;
  std::cout << "mbb = " << main_bbox << std::endl;
   
  // Figure out where each folder goes
  for (size_t i=0; i<folder_list.size(); ++i) {
    BBox2i bbox = bbox_from_folder(folder_list[i]);
    std::cout << folder_list[i] << " ---> " << bbox << "\n";
    
    std::string bbox_string = extract_process_folder_bbox_string(folder_list[i]);
    
    const std::string abs_path = parallel_stereo_folder.string() + 
        "/" + folder_list[i] + "/" + bbox_string + "-Dnosym.tif";
    
    if (bbox == main_bbox)
      continue;
    if (bbox.min().x() < main_bbox.min().x()) { // Tiles to the left
      if (bbox.min().y() < main_bbox.min().y()) { // Top left
        blend_options.tile_paths[TL] = abs_path;
        blend_options.rois      [TL] = bbox;
        continue;
      }
      if (bbox.min().y() > main_bbox.min().y()) { // Bot left
        blend_options.tile_paths[BL] = abs_path;
        blend_options.rois      [BL] = bbox;
        continue;
      }
      // Left
      blend_options.tile_paths[L] = abs_path;
      blend_options.rois      [L] = bbox;
      continue;
    }
    if (bbox.min().x() > main_bbox.min().x()) { // Tiles to the right
      if (bbox.min().y() < main_bbox.min().y()) { // Top right
        blend_options.tile_paths[TR] = abs_path;
        blend_options.rois      [TR] = bbox;
        continue;
      }
      if (bbox.min().y() > main_bbox.min().y()) { // Bot right
        blend_options.tile_paths[BR] = abs_path;
        blend_options.rois      [BR] = bbox;
        continue;
      }
      // Right
      blend_options.tile_paths[R] = abs_path;
      blend_options.rois      [R] = bbox;
      continue;
    }
    // Top and bottom tiles
    if (bbox.min().y() < main_bbox.min().y()) { // Top
      blend_options.tile_paths[T] = abs_path;
      blend_options.rois      [T] = bbox;
      continue;
    }
    if (bbox.min().y() > main_bbox.min().y()) { // Bottom
      blend_options.tile_paths[B] = abs_path;
      blend_options.rois      [B] = bbox;
      continue;
    }
    vw_throw( ArgumentErr() << "Unrecognized folder location: " << folder_list[i] );
  }
  
  //vw_throw( ArgumentErr() << "DEBUG!" );
  
  blend_options.sgm_collar_size = stereo_settings().sgm_collar_size;
}

void stereo_blending( ASPGlobalOptions const& opt ) {

  BlendOptions blend_options;
  fill_blend_options(opt, blend_options);

  // This tool is only intended to run as part of parallel_stereo, which
  //  renames the normal -D.tif file to -Dnosym.tif.

  // Since this tool is only for follow-up processing of SGM results in
  //  parallel_stereo, it can be safely assumed that the input images are small enough
  //  to load entirely into memory.
  DispImageType integer_disp;

  try {

    // Verify that the input correlation file is float, indicating SGM processing.
    // - No need to run the blend operation on integer files!
    boost::scoped_ptr<SrcImageResource> rsrc(DiskImageResource::open(blend_options.main_path));
    ChannelTypeEnum disp_data_type = rsrc->channel_type();
    if (disp_data_type == VW_CHANNEL_INT32)
      vw_throw( ArgumentErr() << "Error: stereo_blend should only be called after SGM correlation." );
    integer_disp = DiskImageType(blend_options.main_path);
    
  } catch (IOErr const& e) {
    vw_throw( ArgumentErr() << "\nUnable to start at blending stage -- could not read input files.\n" 
                            << e.what() << "\nExiting.\n\n" );
  }
  
  cartography::GeoReference left_georef;
  bool   has_left_georef = read_georeference(left_georef,  opt.out_prefix + "-L.tif");
  bool   has_nodata      = false;
  double nodata          = -32768.0;

  DispImageType output = tile_blend(integer_disp, blend_options);

  string rd_file = opt.out_prefix + "-RD.tif";
  vw_out() << "Writing: " << rd_file << "\n";
  vw::cartography::block_write_gdal_image(rd_file, output,
                                          has_left_georef, left_georef,
                                          has_nodata, nodata, opt,
                                          TerminalProgressCallback("asp", "\t--> Blending :") );
}

int main(int argc, char* argv[]) {

  //try {
    //xercesc::XMLPlatformUtils::Initialize();

    vw_out() << "\n[ " << current_posix_time_string()
             << " ] : Stage 2 --> BLENDING \n";

    stereo_register_sessions();

    bool verbose = false;
    vector<ASPGlobalOptions> opt_vec;
    string output_prefix;
    asp::parse_multiview(argc, argv, SubpixelDescription(),
                         verbose, output_prefix, opt_vec);
    ASPGlobalOptions opt = opt_vec[0];

    // Subpixel refinement uses smaller tiles.
    //---------------------------------------------------------
    int ts = ASPGlobalOptions::rfne_tile_size();
    opt.raster_tile_size = Vector2i(ts, ts);

    // Internal Processes
    //---------------------------------------------------------
    stereo_blending( opt );

    vw_out() << "\n[ " << current_posix_time_string()
             << " ] : BLENDING FINISHED \n";

    //xercesc::XMLPlatformUtils::Terminate();
  //} ASP_STANDARD_CATCHES;

  return 0;
}
