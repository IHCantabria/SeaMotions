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

#pragma once

#include <iosfwd>
#include <string>

struct CliOptions
{
    std::string case_path;
    bool        out_memory_report   = false;
    bool        show_help           = false;
    std::string error;
};

CliOptions parse_cli_options(int argc, char* argv[]);
void print_cli_usage(std::ostream& out, const std::string& program_name);
