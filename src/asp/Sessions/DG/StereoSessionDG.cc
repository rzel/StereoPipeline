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


/// \file StereoSessionDG.cc
///
#include <vw/Image/ImageMath.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/Transform.h>
#include <vw/Camera/Extrinsics.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Math/Matrix.h>
#include <vw/Cartography/Datum.h>

#include <asp/Core/StereoSettings.h>
#include <asp/Core/InterestPointMatching.h>
#include <asp/Core/AffineEpipolar.h>
#include <asp/Sessions/DG/LinescanDGModel.h>
#include <asp/Sessions/DG/StereoSessionDG.h>
#include <asp/Sessions/RPC/RPCModel.h>
#include <asp/Sessions/DG/XML.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <xercesc/util/PlatformUtils.hpp>

using namespace vw;
using namespace asp;
using namespace xercesc;

namespace pt = boost::posix_time;
namespace fs = boost::filesystem;

// Allows FileIO to correctly read/write these pixel types
namespace vw {
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<Vector2f>  { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };
}

// Helper class for converting to floating point seconds based on a
// given reference.
class SecondsFrom {
  pt::ptime m_reference;
public:
  SecondsFrom( pt::ptime const& time ) : m_reference(time) {}

  double operator()( pt::ptime const& time ) const {
    return double( (time - m_reference).total_microseconds() ) / 1e6;
  }
};


namespace asp {

  pt::ptime parse_time(std::string str){
    try{
      return pt::time_from_string(str);
    }catch(...){
      vw_throw(ArgumentErr() << "Failed to parse time from string: "
               << str << "\n");
    }
    return pt::time_from_string(str); // never reached
  }


  // These are initializers and closers for Xercesc since we use it to
  // read our RPC models.
  StereoSessionDG::StereoSessionDG() {
    XMLPlatformUtils::Initialize();
  }

  StereoSessionDG::~StereoSessionDG() {
    XMLPlatformUtils::Terminate();
  }

  // Provide our camera model
  boost::shared_ptr<camera::CameraModel>
  StereoSessionDG::camera_model( std::string const& /*image_file*/,
                                 std::string const& camera_file ) {

    GeometricXML geo;
    AttitudeXML att;
    EphemerisXML eph;
    ImageXML img;
    RPCXML rpc;
    read_xml( camera_file, geo, att, eph, img, rpc );

    // Convert measurements in millimeters to pixels.
    geo.principal_distance /= geo.detector_pixel_pitch;
    geo.detector_origin /= geo.detector_pixel_pitch;

    bool correct_velocity_aberration = !stereo_settings().disable_correct_velocity_aberration;

    // Convert all time measurements to something that boost::date_time can read.
    boost::replace_all( eph.start_time, "T", " " );
    boost::replace_all( img.tlc_start_time, "T", " " );
    boost::replace_all( img.first_line_start_time, "T", " " );
    boost::replace_all( att.start_time, "T", " " );

    // Convert UTC time measurements to line measurements. Ephemeris
    // start time will be our reference frame to calculate seconds
    // against.
    SecondsFrom convert( parse_time( eph.start_time ) );

    // I'm going make the assumption that EPH and ATT are sampled at the
    // same rate and time.
    VW_ASSERT( eph.position_vec.size() == att.quat_vec.size(),
               MathErr() << "Ephemeris and Attitude don't have the same number of samples." );
    VW_ASSERT( eph.start_time == att.start_time && eph.time_interval == att.time_interval,
               MathErr() << "Ephemeris and Attitude don't seem to sample with the same t0 or dt." );

    // Convert ephemeris to be position of camera. Change attitude to be
    // be the rotation from camera frame to world frame. We also add an
    // additional rotation to the camera frame so X is the horizontal
    // direction to the picture and +Y points down the image (in the
    // direction of flight).
    Quat sensor_coordinate = math::euler_xyz_to_quaternion(Vector3(0,0,geo.detector_rotation * M_PI/180.0 - M_PI/2));
    for ( size_t i = 0; i < eph.position_vec.size(); i++ ) {
      eph.position_vec[i] += att.quat_vec[i].rotate( geo.perspective_center );
      att.quat_vec[i] = att.quat_vec[i] * geo.camera_attitude * sensor_coordinate;
    }

    // Load up the time interpolation class. If the TLCList only has
    // one entry ... then we have to manually drop in the slope and
    // offset.
    if ( img.tlc_vec.size() == 1 ) {
      double direction = 1;
      if ( boost::to_lower_copy( img.scan_direction ) !=
           "forward" ) {
        direction = -1;
      }
      img.tlc_vec.push_back( std::make_pair(img.tlc_vec.front().first +
                                            img.avg_line_rate, direction) );
    }

    // Build the TLCTimeInterpolation object and do a quick sanity check.
    camera::TLCTimeInterpolation tlc_time_interpolation( img.tlc_vec,
                                                         convert( parse_time( img.tlc_start_time ) ) );
    VW_ASSERT( fabs( convert( parse_time( img.first_line_start_time ) ) -
                     tlc_time_interpolation( 0 ) ) < fabs( 1.0 / (10.0 * img.avg_line_rate ) ),
               MathErr() << "First Line Time and output from TLC lookup table do not agree of the ephemeris time for the first line of the image." );

    typedef LinescanDGModel<camera::PiecewiseAPositionInterpolation, camera::LinearPiecewisePositionInterpolation, camera::SLERPPoseInterpolation, camera::TLCTimeInterpolation> camera_type;
    typedef boost::shared_ptr<camera::CameraModel> result_type;
    double et0 = convert( parse_time( eph.start_time ) );
    double at0 = convert( parse_time( att.start_time ) );
    double edt = eph.time_interval;
    double adt = att.time_interval;
    return result_type(new camera_type(camera::PiecewiseAPositionInterpolation(eph.position_vec, eph.velocity_vec,
                                                                               et0, edt ),
                                       camera::LinearPiecewisePositionInterpolation(eph.velocity_vec, et0, edt),
                                       camera::SLERPPoseInterpolation(att.quat_vec, at0, adt),
                                       tlc_time_interpolation, img.image_size,
                                       subvector(inverse(sensor_coordinate).rotate(Vector3(geo.detector_origin[0],
                                                                                           geo.detector_origin[1],
                                                                                           0)), 0, 2),
                                       geo.principal_distance, correct_velocity_aberration)
                      );
  }

  bool StereoSessionDG::ip_matching( std::string const& match_filename,
                                     double left_nodata_value,
                                     double right_nodata_value ) {
    // Load the unmodified images
    DiskImageView<float> left_disk_image( m_left_image_file ),
      right_disk_image( m_right_image_file );
    
    boost::shared_ptr<camera::CameraModel> left_cam, right_cam;
    camera_models( left_cam, right_cam );

    return
      ip_matching_w_alignment( left_cam.get(), right_cam.get(),
                               left_disk_image, right_disk_image,
                               cartography::Datum("WGS84"), match_filename,
                               left_nodata_value,
                               right_nodata_value);
  }
  
  StereoSessionDG::left_tx_type
  StereoSessionDG::tx_left() const {
    Matrix<double> tx = math::identity_matrix<3>();
    if ( stereo_settings().alignment_method == "homography" ||
         stereo_settings().alignment_method == "affineepipolar" ) {
      read_matrix( tx, m_out_prefix + "-align-L.exr" );
    }
    return left_tx_type( tx );
  }

  StereoSessionDG::right_tx_type
  StereoSessionDG::tx_right() const {
    Matrix<double> tx = math::identity_matrix<3>();
    if ( stereo_settings().alignment_method == "homography" ||
         stereo_settings().alignment_method == "affineepipolar" ) {
      read_matrix( tx, m_out_prefix + "-align-R.exr" );
    }
    return right_tx_type( tx );
  }

}
