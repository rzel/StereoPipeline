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


/// \file StereoSession.cc
///
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Math/Vector.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

#include <asp/asp_config.h>
#include <asp/Core/StereoSettings.h>
#include <asp/Core/Common.h>
#include <asp/Core/AffineEpipolar.h>
#include <asp/Sessions/DG/StereoSessionDG.h>
#include <asp/Sessions/DGMapRPC/StereoSessionDGMapRPC.h>
#include <asp/Sessions/ISIS/StereoSessionIsis.h>
#include <asp/Sessions/ISISMapISIS/StereoSessionISISMapISIS.h>
#include <asp/Sessions/NadirPinhole/StereoSessionNadirPinhole.h>
#include <asp/Sessions/Pinhole/StereoSessionPinhole.h>
#include <asp/Sessions/RPC/StereoSessionRPC.h>
#include <asp/Sessions/StereoSession.h>

#include <map>
#include <utility>
#include <string>
#include <ostream>
#include <limits>

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

using namespace vw;

// This creates an anonymous namespace where the lookup table for
// stereo sessions lives.
namespace {
  typedef std::map<std::string,asp::StereoSession::construct_func> ConstructMapType;
  ConstructMapType *stereo_session_construct_map = 0;
}

// Allows FileIO to correctly read/write these pixel types
namespace vw {
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
}

namespace asp {

  // Pass over all the string variables we use
  void StereoSession::initialize( BaseOptions const& options,
                                  std::string const& left_image_file,
                                  std::string const& right_image_file,
                                  std::string const& left_camera_file,
                                  std::string const& right_camera_file,
                                  std::string const& out_prefix,
                                  std::string const& input_dem,
                                  std::string const& extra_argument1,
                                  std::string const& extra_argument2,
                                  std::string const& extra_argument3) {
    m_options = options;
    m_left_image_file = left_image_file;
    m_right_image_file = right_image_file;
    m_left_camera_file = left_camera_file;
    m_right_camera_file = right_camera_file;
    m_out_prefix = out_prefix;
    m_input_dem = input_dem;
    m_extra_argument1 = extra_argument1;
    m_extra_argument2 = extra_argument2;
    m_extra_argument3 = extra_argument3;
  }

  static void register_default_session_types() {
    static bool already = false;
    if ( already ) return;
    already = true;
    StereoSession::register_session_type( "dg", &StereoSessionDG::construct );
    StereoSession::register_session_type( "dgmaprpc", &StereoSessionDGMapRPC::construct );
    StereoSession::register_session_type( "nadirpinhole", &StereoSessionNadirPinhole::construct );
    StereoSession::register_session_type( "pinhole", &StereoSessionPinhole::construct );
    StereoSession::register_session_type( "rpc", &StereoSessionRPC::construct );
#if defined(ASP_HAVE_PKG_ISISIO) && ASP_HAVE_PKG_ISISIO == 1
    asp::StereoSession::register_session_type( "isis", &asp::StereoSessionIsis::construct);
    asp::StereoSession::register_session_type( "isismapisis", &asp::StereoSessionISISMapISIS::construct);
#endif
  }

  StereoSession* StereoSession::create( std::string & session_type, // in-out variable
                                        BaseOptions const& options,
                                        std::string const& left_image_file,
                                        std::string const& right_image_file,
                                        std::string const& left_camera_file,
                                        std::string const& right_camera_file,
                                        std::string const& out_prefix,
                                        std::string const& input_dem,
                                        std::string const& extra_argument1,
                                        std::string const& extra_argument2,
                                        std::string const& extra_argument3) {
    register_default_session_types();

    // Known user session types are:
    // DG, RPC, ISIS, Pinhole, NadirPinhole
    //
    // Hidden sessions are:
    // DGMapRPC, ISISMapISIS, Blank (Guessing)

    // Try to guess the session if not provided
    std::string actual_session_type = session_type;
    boost::to_lower(actual_session_type);
    if ( actual_session_type.empty() ) {
      
      if ( asp::has_cam_extension( left_camera_file ) ||
           asp::has_cam_extension( right_camera_file ) ) {
        actual_session_type = "pinhole";
      }
      if ( boost::iends_with(boost::to_lower_copy(left_image_file), ".cub")   ||
           boost::iends_with(boost::to_lower_copy(right_image_file), ".cub")  ||
           boost::iends_with(boost::to_lower_copy(left_camera_file), ".cub")  ||
           boost::iends_with(boost::to_lower_copy(right_camera_file), ".cub")
           ) {
        actual_session_type = "isis";
      }
      if (boost::iends_with(boost::to_lower_copy(left_camera_file), ".xml") ||
          boost::iends_with(boost::to_lower_copy(right_camera_file), ".xml") ) {
        actual_session_type = "dg";
      }
    }

    // If a DEM is provided, need to use a map-projected session
    if ( !input_dem.empty() && actual_session_type == "dg" ) {
      // User says DG .. but also gives a DEM.
      actual_session_type = "dgmaprpc";
      VW_OUT(DebugMessage,"asp")
        << "Changing session type to: dgmaprpc" << std::endl;
    }
    if ( !input_dem.empty() && actual_session_type == "isis" ) {
      // User says ISIS .. but also gives a DEM.
      actual_session_type = "isismapisis";
      VW_OUT(DebugMessage,"asp")
        << "Changing session type to: isismapisis" << std::endl;
    }
    
    try {
      if ( actual_session_type.empty() ) {
        // RPC can be in the main file or it can be in the camera file
        // DG sessions are always RPC sessions because they contain that
        // as an extra camera model. Thus this RPC check must happen
        // last.
        StereoSessionRPC session;
        boost::shared_ptr<camera::CameraModel>
          left_model = session.camera_model( left_image_file, left_camera_file ),
          right_model = session.camera_model( right_image_file, right_camera_file );
        actual_session_type = "rpc";
      }
    } catch ( vw::NotFoundErr const& e ) {
      // If it throws, it wasn't RPC
    } catch ( ... ) {
      // It didn't even have XML!
    }

    // We should know the session type by now.
    VW_ASSERT( !actual_session_type.empty(),
               ArgumentErr() << "Could not determine stereo session type. "
               << "Please set it explicitly using the -t switch.\n"
               << "Options include: [pinhole isis dg rpc].\n" );
    VW_OUT(DebugMessage,"asp") << "Guessed session type to be " << actual_session_type << std::endl;
    
    if( stereo_session_construct_map ) {
      ConstructMapType::const_iterator i =
        stereo_session_construct_map->find( actual_session_type );
      if( i != stereo_session_construct_map->end() ) {
        StereoSession* session_new = i->second();
        session_new->initialize( options,
                                 left_image_file, right_image_file,
                                 left_camera_file, right_camera_file,
                                 out_prefix, input_dem, extra_argument1,
                                 extra_argument2, extra_argument3 );
        session_type = session_new->name(); // we count on this in the caller
        return session_new;
      }
    }
    vw_throw( vw::NoImplErr() << "Unsuppported stereo session type: "
              << session_type );
    return 0; // never reached
  }

  void StereoSession::register_session_type( std::string const& id,
                                             StereoSession::construct_func func) {
    if( ! stereo_session_construct_map )
      stereo_session_construct_map = new ConstructMapType();
    stereo_session_construct_map->insert( std::make_pair( id, func ) );
  }

  // Base class implementation of processing steps. All of these don't
  // perform anything, they're just place holders.
  void StereoSession::camera_models(boost::shared_ptr<vw::camera::CameraModel> &cam1,
                                    boost::shared_ptr<vw::camera::CameraModel> &cam2) {
    cam1 = camera_model(m_left_image_file, m_left_camera_file);
    cam2 = camera_model(m_right_image_file, m_right_camera_file);
  }
  
  bool StereoSession::ip_matching( std::string const& match_filename,
                                   double left_nodata_value,
                                   double right_nodata_value ){
    vw_throw( NoImplErr() << "ip_matching not implemented for session: " << name()
              <<"\n." );
  }
  
  // The default hook assumes that the images are tifs, applies alignment
  // and scales them.
  void StereoSession::pre_preprocessing_hook(std::string const& left_input_file,
                                             std::string const& right_input_file,
                                             std::string &left_output_file,
                                             std::string &right_output_file) {
    
    boost::shared_ptr<DiskImageResource>
      left_rsrc( DiskImageResource::open(m_left_image_file) ),
      right_rsrc( DiskImageResource::open(m_right_image_file) );

    float left_nodata_value, right_nodata_value;
    get_nodata_values(left_rsrc, right_rsrc, left_nodata_value, right_nodata_value);

    // Load the unmodified images
    DiskImageView<float> left_disk_image( left_rsrc ), right_disk_image( right_rsrc );

    // Filenames of normalized images
    left_output_file = m_out_prefix + "-L.tif";
    right_output_file = m_out_prefix + "-R.tif";

    // If these files already exist, don't bother writting them again.
    bool rebuild = false;
    try {
      vw_log().console_log().rule_set().add_rule(-1,"fileio");
      DiskImageView<float> test_left(left_output_file);
      DiskImageView<float> test_right(right_output_file);
      vw_settings().reload_config();
    } catch (vw::IOErr const& e) {
      vw_settings().reload_config();
      rebuild = true;
    } catch (vw::ArgumentErr const& e ) {
      // Throws on a corrupted file.
      vw_settings().reload_config();
      rebuild = true;
    }

    if (!rebuild) {
      vw_out() << "\t--> Using cached L and R files.\n";
      return;
    }

    ImageViewRef< PixelMask<float> > left_masked_image
      = create_mask_less_or_equal(left_disk_image, left_nodata_value);
    ImageViewRef< PixelMask<float> > right_masked_image
      = create_mask_less_or_equal(right_disk_image, right_nodata_value);

    Vector4f left_stats  = gather_stats( left_masked_image,  "left" );
    Vector4f right_stats = gather_stats( right_masked_image, "right" );

    ImageViewRef< PixelMask<float> > Limg, Rimg;
    std::string lcase_file = boost::to_lower_copy(m_left_camera_file);

    if ( stereo_settings().alignment_method == "homography" ||
         stereo_settings().alignment_method == "affineepipolar" ) {
      std::string match_filename
        = ip::match_filename(m_out_prefix, left_input_file, right_input_file);

      if (!fs::exists(match_filename)) {
        // This is calling an internal virtualized method.
        bool inlier =
          ip_matching( match_filename,
                       left_nodata_value,
                       right_nodata_value );
        if ( !inlier ) {
          fs::remove( match_filename );
          vw_throw( IOErr() << "Unable to match left and right images." );
        }
      } else {
        vw_out() << "\t--> Using cached match file: " << match_filename << "\n";
      }

      std::vector<ip::InterestPoint> left_ip, right_ip;
      ip::read_binary_match_file( match_filename, left_ip, right_ip  );
      Matrix<double> align_left_matrix = math::identity_matrix<3>(),
        align_right_matrix = math::identity_matrix<3>();

      Vector2i left_size = file_image_size( left_input_file ),
        right_size = file_image_size( right_input_file );

      if ( stereo_settings().alignment_method == "homography" ) {
        left_size =
          homography_rectification( left_size, right_size, left_ip, right_ip,
                                    align_left_matrix, align_right_matrix );
        vw_out() << "\t--> Aligning right image to left using matrices:\n"
                 << "\t      " << align_left_matrix << "\n"
                 << "\t      " << align_right_matrix << "\n";
      } else {
        left_size =
          affine_epipolar_rectification( left_size, right_size,
                                         left_ip, right_ip,
                                         align_left_matrix,
                                         align_right_matrix );
        vw_out() << "\t--> Aligning left and right images using affine matrices:\n"
                 << "\t      " << submatrix(align_left_matrix,0,0,2,3) << "\n"
                 << "\t      " << submatrix(align_right_matrix,0,0,2,3) << "\n";
      }
      write_matrix( m_out_prefix + "-align-L.exr", align_left_matrix );
      write_matrix( m_out_prefix + "-align-R.exr", align_right_matrix );

      // Applying alignment transform
      Limg = transform(left_masked_image,
                       HomographyTransform(align_left_matrix),
                       left_size.x(), left_size.y() );
      Rimg = transform(right_masked_image,
                       HomographyTransform(align_right_matrix),
                       left_size.x(), left_size.y() );
    } else if ( stereo_settings().alignment_method == "epipolar" ) {
      vw_throw( NoImplErr() << "Session " << name() << " does not support epipolar rectification" );
    } else {
      // Do nothing just provide the original files.
      Limg = left_masked_image;
      Rimg = right_masked_image;
    }

    // Apply our normalization options.
    if ( stereo_settings().force_use_entire_range > 0 ) {
      if ( stereo_settings().individually_normalize > 0 ) {
        vw_out() << "\t--> Individually normalize images to their respective min max\n";
        Limg = normalize( Limg, left_stats[0], left_stats[1], 0.0, 1.0 );
        Rimg = normalize( Rimg, right_stats[0], right_stats[1], 0.0, 1.0 );
      } else {
        float low = std::min(left_stats[0], right_stats[0]);
        float hi  = std::max(left_stats[1], right_stats[1]);
        vw_out() << "\t--> Normalizing globally to: [" << low << " " << hi << "]\n";
        Limg = normalize( Limg, low, hi, 0.0, 1.0 );
        Rimg = normalize( Rimg, low, hi, 0.0, 1.0 );
      }
    } else {
      if ( stereo_settings().individually_normalize > 0 ) {
        vw_out() << "\t--> Individually normalize images to their respective 4 std dev window\n";
        Limg = normalize( Limg, left_stats[2] - 2*left_stats[3],
                          left_stats[2] + 2*left_stats[3], 0.0, 1.0 );
        Rimg = normalize( Rimg, right_stats[2] - 2*right_stats[3],
                          right_stats[2] + 2*right_stats[3], 0.0, 1.0 );
      } else {
        float low = std::min(left_stats[2] - 2*left_stats[3],
                             right_stats[2] - 2*right_stats[3]);
        float hi  = std::max(left_stats[2] + 2*left_stats[3],
                             right_stats[2] + 2*right_stats[3]);
        vw_out() << "\t--> Normalizing globally to: [" << low << " " << hi << "]\n";
        Limg = normalize( Limg, low, hi, 0.0, 1.0 );
        Rimg = normalize( Rimg, low, hi, 0.0, 1.0 );
      }
    }


    // The output no-data value must be < 0 as we scale the images to [0, 1].
    float output_nodata = -32768.0;

    // Enforce no predictor in compression, it works badly with L.tif and R.tif.
    asp::BaseOptions options = m_options;
    options.gdal_options["PREDICTOR"] = "1";

    vw_out() << "\t--> Writing pre-aligned images.\n";
    block_write_gdal_image( left_output_file, apply_mask(Limg, output_nodata),
                            output_nodata, options,
                            TerminalProgressCallback("asp","\t  L:  ") );

    if ( stereo_settings().alignment_method == "none" )
      block_write_gdal_image( right_output_file, apply_mask(Rimg, output_nodata),
                              output_nodata, options,
                              TerminalProgressCallback("asp","\t  R:  ") );
    else
      block_write_gdal_image( right_output_file,
                              apply_mask(crop(edge_extend(Rimg, ConstantEdgeExtension()), bounding_box(Limg)), output_nodata),
                              output_nodata, options,
                              TerminalProgressCallback("asp","\t  R:  ") );
  }


  // Other processing hooks. The default is to do nothing.

  void StereoSession::post_preprocessing_hook(std::string const& input_file1,
                                              std::string const& input_file2,
                                              std::string &output_file1,
                                              std::string &output_file2) {
    output_file1 = input_file1;
    output_file2 = input_file2;
  }

  void StereoSession::pre_correlation_hook(std::string const& input_file1,
                                           std::string const& input_file2,
                                           std::string &output_file1,
                                           std::string &output_file2) {
    output_file1 = input_file1;
    output_file2 = input_file2;
  }

  void StereoSession::post_correlation_hook(std::string const& input_file,
                                            std::string & output_file) {
    output_file = input_file;
  }

  void StereoSession::pre_filtering_hook(std::string const& input_file,
                                         std::string & output_file) {
    output_file = input_file;
  }

  void StereoSession::post_filtering_hook(std::string const& input_file,
                                          std::string & output_file) {
    output_file = input_file;
  }

  ImageViewRef<PixelMask<Vector2f> >
  StereoSession::pre_pointcloud_hook(std::string const& input_file) {
    return DiskImageView<PixelMask<Vector2f> >( input_file );
  }

  void StereoSession::post_pointcloud_hook(std::string const& input_file,
                                           std::string & output_file) {
    output_file = input_file;
  }

  void StereoSession::get_nodata_values(boost::shared_ptr<vw::DiskImageResource> left_rsrc,
                                        boost::shared_ptr<vw::DiskImageResource> right_rsrc,
                                        float & left_nodata_value,
                                        float & right_nodata_value){

    // The no-data value read from options overrides the value present
    // in the image files.
    left_nodata_value = std::numeric_limits<float>::quiet_NaN();
    right_nodata_value = std::numeric_limits<float>::quiet_NaN();
    if ( left_rsrc->has_nodata_read()  ) left_nodata_value  = left_rsrc->nodata_read();
    if ( right_rsrc->has_nodata_read() ) right_nodata_value = right_rsrc->nodata_read();
    float opt_nodata = stereo_settings().nodata_value;
    if (!std::isnan(opt_nodata)){

      if ( opt_nodata < left_nodata_value )
        vw_out(WarningMessage) << "It appears that the user-supplied no-data value is less than the no-data value of left image. This may not be what was intended.\n";
      if ( opt_nodata < right_nodata_value )
        vw_out(WarningMessage) << "It appears that the user-supplied no-data value is less than the no-data value of right image. This may not be what was intended.\n";

      left_nodata_value  = opt_nodata;
      right_nodata_value = opt_nodata;
    }

    return;
  }

}
