/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../../src/inout/hdf5_time_series_exporter_mpi.hpp"

namespace fs = std::filesystem;

namespace
{
constexpr double     kPi               = 3.141592653589793238462643383279502884;
constexpr const char* kMeshRelativePath = "../tests_data/test_time_export.poly";
constexpr int        kNumSteps         = 6;
constexpr cusfloat   kDt               = static_cast<cusfloat>(0.15);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int exit_code = 0;

    try
    {
        const fs::path default_root = fs::path(__FILE__).parent_path() / "../../aux_data/tests_output/sample_hdf5_time_mpi";
        const std::string output_root = (argc > 1)
            ? std::string(argv[1])
            : fs::weakly_canonical(default_root).string();

        if ( !fs::exists(output_root) )
        {
            fs::create_directories(output_root);
        }

        const fs::path _mesh_path   = fs::path(__FILE__).parent_path( ) / kMeshRelativePath;
        const std::string mesh_path = fs::weakly_canonical(_mesh_path).string( );
        if (!fs::exists(mesh_path))
        {
            throw std::runtime_error("Mesh file not found at " + mesh_path);
        }

        cusfloat cog[3] = {0.0F, 0.0F, 0.0F};
        Mesh mesh(
                        mesh_path,
                        "surface",
                        cog,
                        false,
                        PanelTypeE::DIFFRAC,
                        false
                    );

        MpiConfig mpi_config(rank, size, 0, MPI_COMM_WORLD);
        HDF5TimeSeriesExporterMPI exporter(output_root, &mesh, &mpi_config);
 
        exporter.add_field("Pressure", 1);
        exporter.add_field("Velocity", 3);

        int node_start = 0;
        int node_end   = 0;
        mpi_config.get_1d_bounds(mesh.nodes_np, node_start, node_end);
        const std::size_t local_np = static_cast<std::size_t>(node_end - node_start);

        cut::CusTensor<cusfloat> pressure(local_np);
        cut::CusTensor<cusfloat> velocity({local_np, static_cast<std::size_t>(3)});

        const cusfloat period = static_cast<cusfloat>(2.0);
        const cusfloat omega  = static_cast<cusfloat>(2.0 * kPi) / period;
        const cusfloat wave_k = (omega * omega) / static_cast<cusfloat>(9.81);

        for (int step = 0; step < kNumSteps; ++step)
        {
            const cusfloat time = static_cast<cusfloat>(step) * kDt;
            for (std::size_t local = 0; local < local_np; ++local)
            {
                const std::size_t gindex = static_cast<std::size_t>(node_start) + local;
                const cusfloat phase_x = wave_k * mesh.x[gindex] - omega * time;
                const cusfloat phase_y = wave_k * mesh.y[gindex] - omega * time;
                const cusfloat cos_x   = std::cos(phase_x);
                const cusfloat cos_y   = std::cos(phase_y);
                const cusfloat sin_x   = std::sin(phase_x);
                const cusfloat sin_y   = std::sin(phase_y);

                pressure[local]    = cos_x * cos_y;
                velocity(local, 0) = -wave_k * sin_x * cos_y;
                velocity(local, 1) = -wave_k * cos_x * sin_y;
                velocity(local, 2) = 0.0F;
            }

            exporter.append_time(time);
            exporter.append_step("Pressure", &pressure);
            exporter.append_step("Velocity", &velocity);
        }

        exporter.write_xdmf();
        if (rank == 0)
        {
            std::cout << "MPI HDF5/XDMF files written to " << output_root << '\n';
        }
    }
    catch (const std::exception& ex)
    {
        if (rank == 0)
        {
            std::cerr << "MPI exporter example failed: " << ex.what() << '\n';
        }
        exit_code = 1;
    }

    MPI_Finalize();
    return exit_code;
}
