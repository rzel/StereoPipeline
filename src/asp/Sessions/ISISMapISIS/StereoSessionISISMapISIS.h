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


/// \file StereoSessionISISMapISIS.h
///
/// This a session that support ISIS Mapproject ISIS images. It is built
/// entirely so that left and right TX are objects and not
/// TransformRefs.

#ifndef __STEREO_SESSION_ISISMAPISIS_H__
#define __STEREO_SESSION_ISISMAPISIS_H__

#include <asp/Sessions/ISIS/StereoSessionIsis.h>
#include <vw/Cartography/Map2CamTrans.h>
#include <vw/Image/Transform.h>

namespace vw{ namespace camera{
  class CameraModel;
}}
  
namespace asp {
  class StereoSessionISISMapISIS : public StereoSessionIsis {
  public:
    StereoSessionISISMapISIS(){};
    virtual ~StereoSessionISISMapISIS(){};

    // Initializer verifies that the input is map projected
    virtual void initialize(BaseOptions const& options,
                            std::string const& left_image_file,
                            std::string const& right_image_file,
                            std::string const& left_camera_file,
                            std::string const& right_camera_file,
                            std::string const& out_prefix,
                            std::string const& input_dem,
                            std::string const& extra_argument1,
                            std::string const& extra_argument2,
                            std::string const& extra_argument3);

    virtual std::string name() const { return "isismapisis"; }

    virtual void pre_preprocessing_hook(std::string const& input_file1,
                                        std::string const& input_file2,
                                        std::string &output_file1,
                                        std::string &output_file2);
    
    // Allows specialization of how matches are captured.
    virtual bool ip_matching( std::string const& match_filename,
                              double left_nodata_value,
                              double right_nodata_value );

    // For reversing the arithmetic applied in preprocessing plus the
    // map projection.
    typedef CompositionTransformPassBBox<vw::cartography::Map2CamTrans,vw::HomographyTransform> left_tx_type;
    typedef CompositionTransformPassBBox<vw::cartography::Map2CamTrans,vw::HomographyTransform> right_tx_type;
    typedef vw::stereo::StereoModel stereo_model_type;
    left_tx_type tx_left() const;
    right_tx_type tx_right() const;

    static StereoSession* construct() { return new StereoSessionISISMapISIS; }

    boost::shared_ptr<vw::camera::CameraModel> m_left_model, m_right_model;
    
  };

}

#endif//__STEREO_SESSION_ISISMAPISIS_H__
