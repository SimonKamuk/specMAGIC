#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <netcdf.h>

#include "navigation.hpp"
#include "tables.hpp"
#include "types.hpp"

// Fields extracted from the zarr archives by py_utils/extract_zarr_inputs.py, already on
// the output grid. They replace the climatologies (and the albedo) pixel by pixel, wherever
// they have a value. Every field is optional: a missing file, a missing variable or a
// missing pixel falls back to the climatology, and the run as a whole carries on.

namespace ZarrInputs {

    // Stands in for NaN once the fields are read: -ffast-math lets the compiler assume
    // that NaN never occurs, so the per-pixel tests must not rely on it
    constexpr float MISSING = -9999.f;

    struct Field {
        std::vector<float> values;
        std::string file;

        bool present() const { return !values.empty(); }
    };

    // How the inputs of the day pixels were sourced, summed over threads after the loop
    struct Tally {
        uint64_t primary = 0;
        uint64_t backup = 0;
        uint64_t climatology = 0;
        uint64_t clipped = 0;

        Tally& operator+=(const Tally& o) {
            primary += o.primary; backup += o.backup;
            climatology += o.climatology; clipped += o.clipped;
            return *this;
        }
    };

    struct Usage {
        Tally aod, o3, h2o, albedo;

        Usage& operator+=(const Usage& o) {
            aod += o.aod; o3 += o.o3; h2o += o.h2o; albedo += o.albedo;
            return *this;
        }
    };

    struct Inputs {

        Field aod, o3, h2o, h2o_cams, albedo_vis, albedo_nir;
        int nlon = 0;

        // Concentration axes of the LUTs. Values from the archives are clipped to these,
        // because unlike the smooth monthly climatologies they do leave them
        MAGIC_REAL aod_min = 0, aod_max = 0;
        MAGIC_REAL h2o_min = 0, h2o_max = 0;
        MAGIC_REAL o3_min = 0, o3_max = 0;

        bool any() const {
            return aod.present() || o3.present() || h2o.present() || h2o_cams.present()
                || albedo_vis.present() || albedo_nir.present();
        }

        bool hasAlbedo() const { return albedo_vis.present() && albedo_nir.present(); }

        // Reads every file named in the list (one path per line)
        void load(const std::string& list_path, const Geography& geo) {

            nlon = geo.nlon;

            std::ifstream in(list_path);
            if (!in) {
                std::printf("[WARNING]: Could not open zarr input list %s. Using climatologies. \n",
                    list_path.c_str());
                return;
            }

            std::string path;
            while (std::getline(in, path)) {
                if (!path.empty()) readFile(path, geo);
            }

            if (!any()) std::printf("[WARNING]: No zarr inputs could be read. Using climatologies. \n");
        }

        void setLimits(const Tables::RTM& rtm, const Tables::Correction::Absorber& water,
            const Tables::Correction::Absorber& ozone) {

            aod_min = rtm.aod[0];       aod_max = rtm.aod[rtm.aoddim - 1];
            h2o_min = water.concen[0];  h2o_max = water.concen[water.nrows - 1];
            o3_min = ozone.concen[0];   o3_max = ozone.concen[ozone.nrows - 1];
        }

        // Overrides the climatology values of this pixel wherever an archive has one
        void apply(PixelClimate& clim, const Area& a, bool count_albedo, Usage& u) const {

            override(clim.aod, aod, nullptr, aod_min, aod_max, u.aod, a);
            override(clim.o3, o3, nullptr, o3_min, o3_max, u.o3, a);

            // DINI first, CAMS outside the DINI domain
            override(clim.h2o, h2o, &h2o_cams, h2o_min, h2o_max, u.h2o, a);

            if (count_albedo) {
                if (at(albedo_vis, a) != MISSING && at(albedo_nir, a) != MISSING) u.albedo.primary++;
                else u.albedo.climatology++;
            }
        }

        // LSA SAF visible covers 0.4-0.7 um and near-infrared 0.7-4.0 um. Returns -1 where
        // there is no value, which sends the caller to the land-use albedo.
        MAGIC_REAL getKatoAlbedo(int zero_based_kato_band, const Area& a) const {

            // Band 16 (684.2-704.4 nm) is the last that is mostly visible
            return albedo(zero_based_kato_band <= 15 ? albedo_vis : albedo_nir, a);
        }

        MAGIC_REAL getSatelliteAlbedo(double wavelength_nm, const Area& a) const {
            return albedo(wavelength_nm < 700. ? albedo_vis : albedo_nir, a);
        }

        // For the global attributes of the output files
        std::vector<std::pair<std::string, std::string>> describe(const Usage& u, bool lsasaf) const {

            std::vector<std::pair<std::string, std::string>> attrs;

            const std::pair<const char*, const Field*> fields[] = {
                {"zarr_aod_file", &aod}, {"zarr_o3_file", &o3}, {"zarr_h2o_file", &h2o},
                {"zarr_h2o_backup_file", &h2o_cams},
                {"zarr_albedo_vis_file", &albedo_vis}, {"zarr_albedo_nir_file", &albedo_nir}};

            for (const auto& f : fields) {
                if (f.second->present()) attrs.emplace_back(f.first, f.second->file);
            }

            attrs.emplace_back("zarr_aod_usage", summary(u.aod));
            attrs.emplace_back("zarr_o3_usage", summary(u.o3));
            attrs.emplace_back("zarr_h2o_usage", summary(u.h2o));
            if (lsasaf) attrs.emplace_back("zarr_albedo_usage", summary(u.albedo));

            return attrs;
        }

        void printUsage(const Usage& u, bool lsasaf) const {

            report("aod", u.aod);
            report("o3", u.o3);
            report("h2o", u.h2o);
            if (lsasaf) report("albedo", u.albedo);
        }

    private:

        float at(const Field& f, const Area& a) const {
            if (!f.present()) return MISSING;
            return f.values[static_cast<size_t>(a.nlat) * nlon + a.nlon];
        }

        // Takes the primary field's value, else the backup's, clipped to the LUT axis.
        // Leaves the climatology value in place when neither has one.
        void override(MAGIC_REAL& target, const Field& primary, const Field* backup,
            MAGIC_REAL lo, MAGIC_REAL hi, Tally& t, const Area& a) const {

            float v = at(primary, a);
            if (v != MISSING) {
                t.primary++;
            } else if (backup && (v = at(*backup, a)) != MISSING) {
                t.backup++;
            } else {
                t.climatology++;
                return;
            }

            MAGIC_REAL x = static_cast<MAGIC_REAL>(v);
            if (x < lo || x > hi) t.clipped++;
            target = std::min(std::max(x, lo), hi);
        }

        MAGIC_REAL albedo(const Field& f, const Area& a) const {
            float v = at(f, a);
            if (v == MISSING || v < 0.f || v > 1.f) return (MAGIC_REAL)-1;
            return static_cast<MAGIC_REAL>(v);
        }

        static bool isNan(float v) {
            uint32_t bits;
            std::memcpy(&bits, &v, sizeof bits);
            return (bits & 0x7f800000u) == 0x7f800000u && (bits & 0x007fffffu) != 0;
        }

        static std::string summary(const Tally& t) {
            uint64_t total = t.primary + t.backup + t.climatology;
            if (total == 0) return "no day pixels";

            char buf[160];
            std::snprintf(buf, sizeof buf,
                "%.1f%% zarr, %.1f%% zarr backup, %.1f%% climatology, %.1f%% clipped to LUT",
                100. * t.primary / total, 100. * t.backup / total,
                100. * t.climatology / total, 100. * t.clipped / total);
            return buf;
        }

        static void report(const char* name, const Tally& t) {
            uint64_t total = t.primary + t.backup + t.climatology;
            if (total == 0) return;

            const char* level = (t.climatology > 0 || t.clipped > 0) ? "[WARNING]: " : "";
            std::printf("%sZarr input %s: %s \n", level, name, summary(t).c_str());
        }

        bool sameAxis(const std::vector<double>& file_axis, const MAGIC_EXACT* geo_axis) const {
            for (size_t i = 0; i < file_axis.size(); i++) {
                if (std::fabs(file_axis[i] - geo_axis[i]) > 1e-3) return false;
            }
            return true;
        }

        void readFile(const std::string& path, const Geography& geo) {

            int ncid;
            if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
                std::printf("[WARNING]: Could not open zarr input %s. Skipping it. \n", path.c_str());
                return;
            }

            // The file must be on exactly the grid we compute on
            int dimid;
            size_t file_nlat = 0, file_nlon = 0;
            bool ok = nc_inq_dimid(ncid, "lat", &dimid) == NC_NOERR
                && nc_inq_dimlen(ncid, dimid, &file_nlat) == NC_NOERR
                && nc_inq_dimid(ncid, "lon", &dimid) == NC_NOERR
                && nc_inq_dimlen(ncid, dimid, &file_nlon) == NC_NOERR
                && file_nlat == geo.nlat && file_nlon == geo.nlon;

            if (ok) {
                std::vector<double> file_lat(file_nlat), file_lon(file_nlon);
                int varid;
                ok = nc_inq_varid(ncid, "lat", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lat.data()) == NC_NOERR
                    && nc_inq_varid(ncid, "lon", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lon.data()) == NC_NOERR
                    && sameAxis(file_lat, geo.lat) && sameAxis(file_lon, geo.lon);
            }

            if (!ok) {
                std::printf("[WARNING]: Zarr input %s is not on the output grid. Skipping it. \n",
                    path.c_str());
                nc_close(ncid);
                return;
            }

            const std::pair<const char*, Field*> fields[] = {
                {"aod", &aod}, {"o3", &o3}, {"h2o", &h2o}, {"h2o_cams", &h2o_cams},
                {"albedo_vis", &albedo_vis}, {"albedo_nir", &albedo_nir}};

            for (const auto& f : fields) {

                int varid;
                if (nc_inq_varid(ncid, f.first, &varid) != NC_NOERR) continue;

                if (f.second->present()) {
                    std::printf("[WARNING]: Zarr input %s found again in %s. Keeping the one from %s. \n",
                        f.first, path.c_str(), f.second->file.c_str());
                    continue;
                }

                std::vector<float> values(file_nlat * file_nlon);
                if (nc_get_var_float(ncid, varid, values.data()) != NC_NOERR) {
                    std::printf("[WARNING]: Could not read %s from %s. Skipping it. \n",
                        f.first, path.c_str());
                    continue;
                }

                for (float& v : values) {
                    if (isNan(v)) v = MISSING;
                }

                f.second->values = std::move(values);
                f.second->file = path;
            }

            nc_close(ncid);
        }
    };

}
