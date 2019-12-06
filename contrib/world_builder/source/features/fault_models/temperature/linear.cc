/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <world_builder/utilities.h>
#include <world_builder/assert.h>
#include <world_builder/nan.h>
#include <world_builder/parameters.h>

#include <world_builder/types/double.h>
#include <world_builder/types/string.h>
#include <world_builder/types/object.h>
#include <world_builder/features/fault_models/temperature/linear.h>


namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace FaultModels
    {
      namespace Temperature
      {
        Linear::Linear(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          top_temperature(NaN::DSNAN),
          bottom_temperature(NaN::DSNAN),
          operation("")
        {
          this->world = world_;
          this->name = "linear";
        }

        Linear::~Linear()
        { }

        void
        Linear::declare_entries(Parameters &prm, const std::string &)
        {
          // Add max depth to the required parameters.
          prm.declare_entry("", Types::Object({"max distance fault center"}), "Temperature model object");

          prm.declare_entry("min distance fault center", Types::Double(0),
                            "The minimum distance to the center of the fault. This determines where the linear temperature starts.");

          prm.declare_entry("max distance fault center", Types::Double(std::numeric_limits<double>::max()),
                            "The minimum distance to the center of the fault. This determines where the linear temperature end.");

          prm.declare_entry("top temperature", Types::Double(293.15),
                            "The temperature at the top in degree Kelvin of this feature."
                            "If the value is below zero, the an adiabatic temperature is used.");

          prm.declare_entry("bottom temperature", Types::Double(-1),
                            "The temperature at the top in degree Kelvin of this feature. "
                            "If the value is below zero, an adiabatic temperature is used.");

          prm.declare_entry("operation", Types::String("replace"),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value (add).");
        }

        void
        Linear::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance fault center");
          max_depth = prm.get<double>("max distance fault center");
          WBAssert(max_depth >= min_depth, "max depth needs to be larger or equal to min depth.");
          top_temperature = prm.get<double>("top temperature");
          bottom_temperature = prm.get<double>("bottom temperature");
        }


        double
        Linear::get_temperature(const Point<3> &,
                                const double /*depth*/,
                                const double gravity_norm,
                                double temperature,
                                const double /*feature_min_depth*/,
                                const double /*feature_max_depth*/,
                                const std::map<std::string,double> &distance_from_planes) const
        {

          if (std::fabs(distance_from_planes.at("distanceFromPlane")) <= max_depth && std::fabs(distance_from_planes.at("distanceFromPlane")) >= min_depth)
            {
              const double min_depth_local = min_depth;
              const double max_depth_local = max_depth;

              double top_temperature_local = top_temperature;
              if (top_temperature_local < 0)
                {
                  top_temperature_local =  this->world->potential_mantle_temperature *
                                           std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                     this->world->specific_heat) * min_depth_local);
                }

              double bottom_temperature_local = bottom_temperature;
              if (bottom_temperature_local < 0)
                {
                  bottom_temperature_local =  this->world->potential_mantle_temperature *
                                              std::exp(((this->world->thermal_expansion_coefficient * gravity_norm) /
                                                        this->world->specific_heat) * max_depth_local);
                }

              return top_temperature +
                     (std::fabs(distance_from_planes.at("distanceFromPlane")) - min_depth_local) * ((bottom_temperature_local - top_temperature_local) / (max_depth_local - min_depth_local));

            }

          return temperature;
        }

        WB_REGISTER_FEATURE_FAULT_TEMPERATURE_MODEL(Linear, linear)
      }
    }
  }
}
