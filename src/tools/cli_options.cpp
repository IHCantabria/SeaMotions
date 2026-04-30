/*
 * Copyright (c) 2025 Sergio Fernandez Ruano / IHCantabria
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

#include "cli_options.hpp"

#include <ostream>

CliOptions parse_cli_options(int argc, char* argv[])
{
    CliOptions options;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i] ? argv[i] : "";

        if (arg == "-h" || arg == "--help")
        {
            options.show_help = true;
            return options;
        }

        if (arg == "--mem-report" || arg == "--memory-report" || arg == "--ram-report")
        {
            options.out_memory_report = true;
            continue;
        }

        if (arg == "-c" || arg == "--case")
        {
            if (i + 1 >= argc)
            {
                options.error = "Missing value for --case";
                return options;
            }
            options.case_path = argv[++i];
            continue;
        }

        if (!arg.empty() && arg[0] == '-')
        {
            options.error = "Unknown option: " + arg;
            return options;
        }

        if (options.case_path.empty())
        {
            options.case_path = arg;
        }
        else
        {
            options.error = "Unexpected extra argument: " + arg;
            return options;
        }
    }

    return options;
}

void print_cli_usage(std::ostream& out, const std::string& program_name)
{
    out << "Usage: " << program_name << " [options] <case_folder>\n";
    out << "Options:\n";
    out << "  -c, --case <path>     Case folder path\n";
    out << "  --mem-report          Write RAM usage CSV reports\n";
    out << "  -h, --help            Show this help\n";
}
